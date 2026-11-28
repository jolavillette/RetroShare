/*******************************************************************************
 * idle/idle.cpp                                                               *
 *                                                                             *
 * Copyright (C) 2003  Justin Karneges <retroshare.project@gmail.com>          *
 *                                                                             *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU Affero General Public License as              *
 * published by the Free Software Foundation, either version 3 of the          *
 * License, or (at your option) any later version.                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful,             *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                *
 * GNU Affero General Public License for more details.                         *
 *                                                                             *
 * You should have received a copy of the GNU Affero General Public License    *
 * along with this program. If not, see <https://www.gnu.org/licenses/>.       *
 *                                                                             *
 *******************************************************************************/

#include "idle.h"

//linux
#ifdef HAVE_XSS

#include <qapplication.h>
#if QT_VERSION >= QT_VERSION_CHECK (6, 0, 0)
#include <QGuiApplication>
#else
#include <QX11Info>
#endif

#if QT_VERSION < QT_VERSION_CHECK(5,0,0)
#include <QDesktopWidget>
#endif

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/extensions/scrnsaver.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

static XErrorHandler old_handler = 0;
extern "C" int xerrhandler(Display* dpy, XErrorEvent* err)
{
	if(err->error_code == BadDrawable)
		return 0;

	return (*old_handler)(dpy, err);
}

// Display of the Qt application, or null when Qt does not run on X11.
static Display *qtDisplay()
{
#if QT_VERSION >= QT_VERSION_CHECK (6, 0, 0)
	QNativeInterface::QX11Application *nativeInterface = qGuiApp->nativeInterface<QNativeInterface::QX11Application>();
	return nativeInterface ? nativeInterface->display() : nullptr;
#elif QT_VERSION >= QT_VERSION_CHECK(5,0,0)
	return QX11Info::isPlatformX11() ? QX11Info::display() : nullptr;
#else
	return QApplication::desktop()->screen()->x11Info().display();
#endif
}

// XScreenSaverQueryInfo() is a round trip to the X server: the caller waits for
// the reply. Asked every second from the GUI thread, it froze the whole window
// for as long as the server took to answer (seconds on a busy remote display).
// The question is therefore asked by a worker thread, on a connection of its
// own, and the GUI thread only reads the last answer.
class IdlePlatform::Private
{
public:
	Private() : display(nullptr), ss_info(nullptr), stop(false), idleSeconds(0) {}

	void run()
	{
		for(;;) {
			int idle = 0;
			if(XScreenSaverQueryInfo(display, DefaultRootWindow(display), ss_info))
				idle = ss_info->idle / 1000;
			idleSeconds.store(idle);

			std::unique_lock<std::mutex> lock(mtx);
			if(cv.wait_for(lock, std::chrono::seconds(1), [this]{ return stop; }))
				return;
		}
	}

	Display *display;          // own connection, used by the worker only once it runs
	XScreenSaverInfo *ss_info; // idem

	std::thread worker;
	std::mutex mtx;            // protects stop
	std::condition_variable cv;
	bool stop;

	std::atomic<int> idleSeconds;
};

IdlePlatform::IdlePlatform()
{
	d = new Private;
}

IdlePlatform::~IdlePlatform()
{
	if(d->worker.joinable()) {
		{
			std::lock_guard<std::mutex> lock(d->mtx);
			d->stop = true;
		}
		d->cv.notify_one();
		d->worker.join();
	}
	if(d->ss_info)
		XFree(d->ss_info);
	if(d->display)
		XCloseDisplay(d->display);
	if(old_handler) {
		XSetErrorHandler(old_handler);
		old_handler = 0;
	}
	delete d;
}

bool IdlePlatform::init()
{
	if(d->ss_info)
		return true;

	Display *appDisplay = qtDisplay();
	if(!appDisplay)
		return false;

	// The connection is opened, and the extension registered on it, here on the
	// GUI thread: the worker then only sends its query on it.
	d->display = XOpenDisplay(DisplayString(appDisplay));
	if(!d->display)
		return false;

	old_handler = XSetErrorHandler(xerrhandler);

	int event_base, error_base;
	if(!XScreenSaverQueryExtension(d->display, &event_base, &error_base))
		return false;

	d->ss_info = XScreenSaverAllocInfo();
	if(!d->ss_info)
		return false;

	d->worker = std::thread(&Private::run, d);
	return true;
}

int IdlePlatform::secondsIdle()
{
	return d->idleSeconds.load();
}

#else // windows
#ifdef WINDOWS_SYS

#include <QLibrary>
#include <windows.h>

#ifndef tagLASTINPUTINFO
typedef struct __tagLASTINPUTINFO {
	UINT cbSize;
	DWORD dwTime;
 } __LASTINPUTINFO, *__PLASTINPUTINFO;
#endif

class IdlePlatform::Private
{
public:
	Private()
	{
		GetLastInputInfo = NULL;
		lib = 0;
	}

	BOOL (__stdcall * GetLastInputInfo)(__PLASTINPUTINFO);
	DWORD (__stdcall * IdleUIGetLastInputTime)(void);
	QLibrary *lib;
};

IdlePlatform::IdlePlatform()
{
	d = new Private;
}

IdlePlatform::~IdlePlatform()
{
	delete d->lib;
	delete d;
}

bool IdlePlatform::init()
{
	if(d->lib)
		return true;

#if QT_VERSION >= QT_VERSION_CHECK(5,0,0)
	QFunctionPointer p;
#else
	void *p;
#endif

	// try to find the built-in Windows 2000 function
	d->lib = new QLibrary("user32");
	if(d->lib->load() && (p = d->lib->resolve("GetLastInputInfo"))) {
		d->GetLastInputInfo = (BOOL (__stdcall *)(__PLASTINPUTINFO))p;
		return true;
	} else {
		delete d->lib;
		d->lib = 0;
	}

	// fall back on idleui
	d->lib = new QLibrary("idleui");
	if(d->lib->load() && (p = d->lib->resolve("IdleUIGetLastInputTime"))) {
		d->IdleUIGetLastInputTime = (DWORD (__stdcall *)(void))p;
		return true;
	} else {
		delete d->lib;
		d->lib = 0;
	}

	return false;
}

int IdlePlatform::secondsIdle()
{
	int i;
	if(d->GetLastInputInfo) {
		__LASTINPUTINFO li;
		li.cbSize = sizeof(__LASTINPUTINFO);
		bool ok = d->GetLastInputInfo(&li);
		if(!ok)
			return 0;
		i = li.dwTime;
	} else if (d->IdleUIGetLastInputTime) {
		i = d->IdleUIGetLastInputTime();
	} else
		return 0;

	return (GetTickCount() - i) / 1000;
}

#else //mac
#ifdef MAC_IDLE

#include <Carbon/Carbon.h>


// Why does Apple have to make this so complicated?
static OSStatus LoadFrameworkBundle(CFStringRef framework, CFBundleRef *bundlePtr) {
	OSStatus  err;
	FSRef   frameworksFolderRef;
	CFURLRef baseURL;
	CFURLRef bundleURL;

	if ( bundlePtr == nil ) return( -1 );

	*bundlePtr = nil;

	baseURL = nil;
	bundleURL = nil;

	err = FSFindFolder(kOnAppropriateDisk, kFrameworksFolderType, true, &frameworksFolderRef);
	if (err == noErr) {
		baseURL = CFURLCreateFromFSRef(kCFAllocatorSystemDefault, &frameworksFolderRef);
		if (baseURL == nil) {
			err = coreFoundationUnknownErr;
		}
	}
	if (err == noErr) {
		bundleURL = CFURLCreateCopyAppendingPathComponent(kCFAllocatorSystemDefault, baseURL, framework, false);
		if (bundleURL == nil) {
			err = coreFoundationUnknownErr;
		}
	}
	if (err == noErr) {
		*bundlePtr = CFBundleCreate(kCFAllocatorSystemDefault, bundleURL);
		if (*bundlePtr == nil) {
			err = coreFoundationUnknownErr;
		}
	}
	if (err == noErr) {
		if ( ! CFBundleLoadExecutable( *bundlePtr ) ) {
			err = coreFoundationUnknownErr;
		}
	}

	// Clean up.
	if (err != noErr && *bundlePtr != nil) {
		CFRelease(*bundlePtr);
		*bundlePtr = nil;
	}
	if (bundleURL != nil) {
		CFRelease(bundleURL);
	}
	if (baseURL != nil) {
		CFRelease(baseURL);
	}

	return err;
}


class IdlePlatform::Private {
public:
	EventLoopTimerRef mTimerRef;
	int mSecondsIdle;

	Private() : mTimerRef(0), mSecondsIdle(0) {}

	static pascal void IdleTimerAction(EventLoopTimerRef, EventLoopIdleTimerMessage inState, void* inUserData);

};


pascal void IdlePlatform::Private::IdleTimerAction(EventLoopTimerRef, EventLoopIdleTimerMessage inState, void* inUserData) {
	switch (inState) {
		case kEventLoopIdleTimerStarted:
		case kEventLoopIdleTimerStopped:
    		// Get invoked with this constant at the start of the idle period,
			// or whenever user activity cancels the idle.
		   ((IdlePlatform::Private*)inUserData)->mSecondsIdle = 0;
			break;
		case kEventLoopIdleTimerIdling:
			// Called every time the timer fires (i.e. every second).
		   ++((IdlePlatform::Private*)inUserData)->mSecondsIdle;
			break;
	}
}


IdlePlatform::IdlePlatform() {
	d = new Private();
}

IdlePlatform::~IdlePlatform() {
	RemoveEventLoopTimer(d->mTimerRef);
	delete d;
}


// Typedef for the function we're getting back from CFBundleGetFunctionPointerForName.
typedef OSStatus (*InstallEventLoopIdleTimerPtr)(EventLoopRef inEventLoop,
																 EventTimerInterval   inFireDelay,
																 EventTimerInterval   inInterval,
																 EventLoopIdleTimerUPP    inTimerProc,
																 void *               inTimerData,
																 EventLoopTimerRef *  outTimer);


bool IdlePlatform::init() {
	// May already be init'ed.
	if (d->mTimerRef) {
		return true;
	}

	// According to the docs, InstallEventLoopIdleTimer is new in 10.2.
	// According to the headers, it has been around since 10.0.
	// One of them is lying.  We'll play it safe and weak-link the function.

	// Load the "Carbon.framework" bundle.
	CFBundleRef carbonBundle;
	if (LoadFrameworkBundle( CFSTR("Carbon.framework"), &carbonBundle ) != noErr) {
		return false;
	}

	// Load the Mach-O function pointers for the routine we will be using.
	InstallEventLoopIdleTimerPtr myInstallEventLoopIdleTimer = (InstallEventLoopIdleTimerPtr)CFBundleGetFunctionPointerForName(carbonBundle, CFSTR("InstallEventLoopIdleTimer"));
	if (myInstallEventLoopIdleTimer == 0) {
		return false;
	}

	EventLoopIdleTimerUPP timerUPP = NewEventLoopIdleTimerUPP(Private::IdleTimerAction);
	if ((*myInstallEventLoopIdleTimer)(GetMainEventLoop(), kEventDurationSecond, kEventDurationSecond, timerUPP, 0, &d->mTimerRef)) {
		return true;
	}

	return false;
}


int IdlePlatform::secondsIdle() {
	return d->mSecondsIdle;
}

#else

IdlePlatform::IdlePlatform() {}
IdlePlatform::~IdlePlatform() {}
bool IdlePlatform::init() { return false; }
int IdlePlatform::secondsIdle() { return 0; }

#endif
#endif
#endif
