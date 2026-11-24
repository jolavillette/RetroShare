/*******************************************************************************
 * plugins/VOIP/gui/VOIPChatWidgetHolder.cpp                                   *
 *                                                                             *
 * Copyright (C) 2012 by Retroshare Team <retroshare.project@gmail.com>        *
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

//C++
#include <time.h>
//Qt
#include <QIcon>
#include <QLayout>
#include <QPropertyAnimation>
#include <QToolButton>
//VOIP
#include <gui/audiodevicehelper.h>
#include "interface/rsVOIP.h"
#include "VOIPChatWidgetHolder.h"
#include "VideoProcessor.h"
#include "QVideoDevice.h"
//retroshare GUI
#include "gui/SoundManager.h"
#include "util/HandleRichText.h"
#include "gui/common/StatusDefs.h"
#include "gui/chat/ChatWidget.h"
//libretroshare
#include <retroshare/rsstatus.h>
#include <retroshare/rspeers.h>
#include <retroshare/rschats.h>
#include <retroshare/rsidentity.h>

#include <QGuiApplication>
#include <QKeyEvent>
#include <QScreen>

#define CALL_START ":/images/call-start.png"
#define CALL_STOP  ":/images/call-stop.png"
#define CALL_HOLD  ":/images/call-hold.png"

// A microphone that has delivered nothing this long after the call started is
// reported.
static const int AUDIO_INPUT_CHECK_DELAY_MS = 5000 ;

// An unanswered call stops ringing after this long.
static const int CALL_RING_TIMEOUT_MS = 60000 ;


VOIPChatWidgetHolder::VOIPChatWidgetHolder(ChatWidget *chatWidget, VOIPNotify *notify)
  : QObject(), ChatWidgetHolder(chatWidget), mVOIPNotify(notify)
{
	int S = QFontMetricsF(chatWidget->font()).height() ;
	QSize iconSize = QSize(3*S,3*S);
	QSize buttonSize = QSize(iconSize + QSize(3,3));

	QIcon iconaudioListenToggleButton ;
	iconaudioListenToggleButton.addPixmap(QPixmap(":/images/speaker_mute.png")) ;
	iconaudioListenToggleButton.addPixmap(QPixmap(":/images/speaker.png"),QIcon::Normal,QIcon::On) ;
	iconaudioListenToggleButton.addPixmap(QPixmap(":/images/speaker.png"),QIcon::Disabled,QIcon::On) ;
	iconaudioListenToggleButton.addPixmap(QPixmap(":/images/speaker.png"),QIcon::Active,QIcon::On) ;
	iconaudioListenToggleButton.addPixmap(QPixmap(":/images/speaker.png"),QIcon::Selected,QIcon::On) ;

	audioListenToggleButton = new QToolButton ;
	audioListenToggleButton->setIcon(iconaudioListenToggleButton) ;
	audioListenToggleButton->setIconSize(iconSize) ;
	audioListenToggleButton->setMinimumSize(buttonSize) ;
	audioListenToggleButton->setMaximumSize(buttonSize) ;
	audioListenToggleButton->setCheckable(true);
	audioListenToggleButton->setAutoRaise(true) ;
	audioListenToggleButton->setText(QString()) ;
	audioListenToggleButton->setToolTip(tr("Mute"));
	audioListenToggleButton->setEnabled(false);

	QIcon iconaudioCaptureToggleButton ;
	iconaudioCaptureToggleButton.addPixmap(QPixmap(":/images/phone.png")) ;
	iconaudioCaptureToggleButton.addPixmap(QPixmap(":/images/microphone_mute.png"),QIcon::Normal,QIcon::On) ;
	iconaudioCaptureToggleButton.addPixmap(QPixmap(":/images/microphone_mute.png"),QIcon::Disabled,QIcon::On) ;
	iconaudioCaptureToggleButton.addPixmap(QPixmap(":/images/microphone_mute.png"),QIcon::Active,QIcon::On) ;
	iconaudioCaptureToggleButton.addPixmap(QPixmap(":/images/microphone_mute.png"),QIcon::Selected,QIcon::On) ;

	audioCaptureToggleButton = new QToolButton ;
	audioCaptureToggleButton->setIcon(iconaudioCaptureToggleButton) ;
	audioCaptureToggleButton->setIconSize(iconSize) ;
	audioCaptureToggleButton->setMinimumSize(buttonSize) ;
	audioCaptureToggleButton->setMaximumSize(buttonSize) ;
	audioCaptureToggleButton->setCheckable(true) ;
	audioCaptureToggleButton->setAutoRaise(true) ;
	audioCaptureToggleButton->setText(QString()) ;
	audioCaptureToggleButton->setToolTip(tr("Start Call"));

	QIcon iconvideoCaptureToggleButton ;
	iconvideoCaptureToggleButton.addPixmap(QPixmap(":/images/filmcam.png")) ;
	iconvideoCaptureToggleButton.addPixmap(QPixmap(":/images/filmcam-off.png"),QIcon::Normal,QIcon::On) ;
	iconvideoCaptureToggleButton.addPixmap(QPixmap(":/images/filmcam-off.png"),QIcon::Disabled,QIcon::On) ;
	iconvideoCaptureToggleButton.addPixmap(QPixmap(":/images/filmcam-off.png"),QIcon::Active,QIcon::On) ;
	iconvideoCaptureToggleButton.addPixmap(QPixmap(":/images/filmcam-off.png"),QIcon::Selected,QIcon::On) ;

	videoCaptureToggleButton = new QToolButton ;
	videoCaptureToggleButton->setIcon(iconvideoCaptureToggleButton) ;
	videoCaptureToggleButton->setIconSize(iconSize) ;
	videoCaptureToggleButton->setMinimumSize(buttonSize) ;
	videoCaptureToggleButton->setMaximumSize(buttonSize) ;
	videoCaptureToggleButton->setCheckable(true) ;
	videoCaptureToggleButton->setAutoRaise(true) ;
	videoCaptureToggleButton->setText(QString()) ;
	videoCaptureToggleButton->setToolTip(tr("Start Video Call"));

	hangupButton = new QToolButton ;
	hangupButton->setIcon(QIcon(":/images/phone_hangup.png")) ;
	hangupButton->setIconSize(iconSize) ;
	hangupButton->setMinimumSize(buttonSize) ;
	hangupButton->setMaximumSize(buttonSize) ;
	hangupButton->setCheckable(false) ;
	hangupButton->setAutoRaise(true) ;
	hangupButton->setText(QString()) ;
	hangupButton->setToolTip(tr("Hangup Call"));
	hangupButton->hide();

	QIcon iconhideChatTextToggleButton ;
	iconhideChatTextToggleButton.addPixmap(QPixmap(":/images/chat-bubble.png")) ;
	iconhideChatTextToggleButton.addPixmap(QPixmap(":/images/chat-bubble.png"),QIcon::Normal,QIcon::On) ;
	iconhideChatTextToggleButton.addPixmap(QPixmap(":/images/chat-bubble.png"),QIcon::Disabled,QIcon::On) ;
	iconhideChatTextToggleButton.addPixmap(QPixmap(":/images/chat-bubble.png"),QIcon::Active,QIcon::On) ;
	iconhideChatTextToggleButton.addPixmap(QPixmap(":/images/chat-bubble.png"),QIcon::Selected,QIcon::On) ;

	hideChatTextToggleButton = new QToolButton ;
	hideChatTextToggleButton->setIcon(iconhideChatTextToggleButton) ;
	hideChatTextToggleButton->setIconSize(iconSize) ;
	hideChatTextToggleButton->setMinimumSize(buttonSize) ;
	hideChatTextToggleButton->setMaximumSize(buttonSize) ;
	hideChatTextToggleButton->setCheckable(true) ;
	hideChatTextToggleButton->setAutoRaise(true) ;
	hideChatTextToggleButton->setText(QString()) ;
	hideChatTextToggleButton->setToolTip(tr("Hide Chat Text"));
	hideChatTextToggleButton->setEnabled(false) ;

	QIcon iconfullscreenToggleButton ;
	iconfullscreenToggleButton.addPixmap(QPixmap(":/images/fullscreen_arrows.png")) ;
	iconfullscreenToggleButton.addPixmap(QPixmap(":/images/fullscreen.png"),QIcon::Normal,QIcon::On) ;
	iconfullscreenToggleButton.addPixmap(QPixmap(":/images/fullscreen.png"),QIcon::Disabled,QIcon::On) ;
	iconfullscreenToggleButton.addPixmap(QPixmap(":/images/fullscreen.png"),QIcon::Active,QIcon::On) ;
	iconfullscreenToggleButton.addPixmap(QPixmap(":/images/fullscreen.png"),QIcon::Selected,QIcon::On) ;

	fullscreenToggleButton = new QToolButton ;
	fullscreenToggleButton->setIcon(iconfullscreenToggleButton) ;
	fullscreenToggleButton->setIconSize(iconSize) ;
	fullscreenToggleButton->setMinimumSize(buttonSize) ;
	fullscreenToggleButton->setMaximumSize(buttonSize) ;
	fullscreenToggleButton->setCheckable(true) ;
	fullscreenToggleButton->setAutoRaise(true) ;
	fullscreenToggleButton->setText(QString()) ;
	fullscreenToggleButton->setToolTip(tr("Fullscreen mode"));
	fullscreenToggleButton->setEnabled(false) ;

	QIcon iconMicMute ;
	iconMicMute.addPixmap(QPixmap(":/images/microphone.png")) ;
	iconMicMute.addPixmap(QPixmap(":/images/microphone_mute.png"),QIcon::Normal,QIcon::On) ;
	iconMicMute.addPixmap(QPixmap(":/images/microphone_mute.png"),QIcon::Active,QIcon::On) ;

	QIcon iconCamera ;
	iconCamera.addPixmap(QPixmap(":/images/filmcam.png")) ;
	iconCamera.addPixmap(QPixmap(":/images/filmcam-off.png"),QIcon::Normal,QIcon::On) ;
	iconCamera.addPixmap(QPixmap(":/images/filmcam-off.png"),QIcon::Active,QIcon::On) ;

	auto makeCallBarButton = [&](QWidget *parent, const QIcon& icon, const QString& tip) {
		QToolButton *b = new QToolButton(parent) ;
		b->setIcon(icon) ;
		b->setIconSize(iconSize) ;
		b->setMinimumSize(buttonSize) ;
		b->setMaximumSize(buttonSize) ;
		b->setCheckable(true) ;
		b->setAutoRaise(true) ;
		b->setText(QString()) ;
		b->setToolTip(tip) ;
		return b ;
	};
	micMuteButton = makeCallBarButton(NULL, iconMicMute, tr("Mute the microphone")) ;
	cameraToggleButton = makeCallBarButton(NULL, iconCamera, tr("Turn the camera off")) ;
	answerNoCameraButton = makeCallBarButton(NULL, QIcon(":/images/phone.png"), tr("Answer without camera")) ;
	answerNoCameraButton->setCheckable(false) ;
	audioStartedWithVideo = false ;
	startWithCameraOff = false ;

	connect(answerNoCameraButton, SIGNAL(clicked()), this , SLOT(answerWithoutCamera()));

	connect(micMuteButton, SIGNAL(clicked()), this , SLOT(toggleMicMute()));
	connect(cameraToggleButton, SIGNAL(clicked()), this , SLOT(toggleCamera()));
	connect(audioListenToggleButton, SIGNAL(clicked()), this , SLOT(toggleAudioListen()));
	connect(audioCaptureToggleButton, SIGNAL(clicked()), this , SLOT(audioCallButtonClicked()));
	connect(videoCaptureToggleButton, SIGNAL(clicked()), this , SLOT(videoCallButtonClicked()));
	connect(hangupButton, SIGNAL(clicked()), this , SLOT(hangupOrDecline()));
	connect(hideChatTextToggleButton, SIGNAL(clicked()), this , SLOT(toggleHideChatText()));
	connect(fullscreenToggleButton, SIGNAL(clicked()), this , SLOT(toggleFullScreen()));

	// At rest: audio call, video call. Ringing in: green answer buttons, red
	// decline. In a call: microphone, camera, full screen, hang-up (see
	// updateCallButtons()). A gap keeps them apart from the chat's own buttons.
	QWidget *titleBarGap = new QWidget ;
	titleBarGap->setFixedWidth(S) ;
	mChatWidget->addTitleBarWidget(titleBarGap) ;
	mChatWidget->addTitleBarWidget(audioCaptureToggleButton) ;
	mChatWidget->addTitleBarWidget(videoCaptureToggleButton) ;
	mChatWidget->addTitleBarWidget(answerNoCameraButton) ;
	mChatWidget->addTitleBarWidget(micMuteButton) ;
	mChatWidget->addTitleBarWidget(audioListenToggleButton) ;
	mChatWidget->addTitleBarWidget(cameraToggleButton) ;
	mChatWidget->addTitleBarWidget(hideChatTextToggleButton) ;
	mChatWidget->addTitleBarWidget(fullscreenToggleButton) ;
	mChatWidget->addTitleBarWidget(hangupButton) ;

	outputAudioProcessor = NULL ;
	outputAudioDevice = NULL ;
	inputAudioProcessor = NULL ;
	inputAudioDevice = NULL ;

	inputVideoDevice = new QVideoInputDevice(mChatWidget) ; // not started yet ;-)
	connect(inputVideoDevice, SIGNAL(cameraErrorString(QString)), this, SLOT(showCameraProblem(QString))) ;

	audioInputBytesAtStart = 0 ;
	audioInputProblemShown = false ;
	audioOutputProblemShown = false ;
	timerAudioInputCheck = new QTimer(this) ;
	timerAudioInputCheck->setInterval(AUDIO_INPUT_CHECK_DELAY_MS) ;
	timerAudioInputCheck->setSingleShot(true) ;
	connect(timerAudioInputCheck, SIGNAL(timeout()), this, SLOT(audioInputCheckTimeOut())) ;
	videoProcessor = new VideoProcessor ;

	// Make a widget with two video devices, one for echo, and one for the talking peer.
	videoWidget = new QWidget(mChatWidget) ;
	videoWidget->setLayout(new QVBoxLayout()) ;
	videoWidget->layout()->addWidget(outputVideoDevice = new QVideoOutputDevice(videoWidget)) ;
	videoWidget->layout()->addWidget(echoVideoDevice = new QVideoOutputDevice(videoWidget)) ;
	videoWidget->hide();

	connect(inputVideoDevice, SIGNAL(networkPacketReady()), this, SLOT(sendVideoData()));

	echoVideoDevice->setMinimumSize(320,240) ;//4/3
	outputVideoDevice->setMinimumSize(320,240) ;//4/3

	echoVideoDevice->showFrameOff();
	outputVideoDevice->showFrameOff();
	
	echoVideoDevice->setStyleSheet("border: 4px solid #CCCCCC; border-radius: 4px;");
	outputVideoDevice->setStyleSheet("border: 4px solid #CCCCCC; border-radius: 4px;");

	/// FULLSCREEN ///
	fullScreenFrame = new QFrame();

	outputVideoDeviceFS = new QVideoOutputDevice(fullScreenFrame);
	outputVideoDeviceFS->setGeometry(QRect(QPoint(0,0),fullScreenFrame->geometry().size()));
	outputVideoDeviceFS->showFrameOff();

	echoVideoDeviceFS = new QVideoOutputDevice(fullScreenFrame);
	echoVideoDeviceFS->setGeometry(QRect(QPoint(fullScreenFrame->width(), fullScreenFrame->height()) - QPoint(320,240), QSize(320,240)));
	echoVideoDeviceFS->showFrameOff();

	toolBarFS = new QFrame(fullScreenFrame);
	QHBoxLayout *toolBarFSLayout = new QHBoxLayout(toolBarFS);

	audioListenToggleButtonFS = new QToolButton(fullScreenFrame) ;
	audioListenToggleButtonFS->setIcon(iconaudioListenToggleButton) ;
	audioListenToggleButtonFS->setIconSize(iconSize) ;
	audioListenToggleButtonFS->setMinimumSize(buttonSize) ;
	audioListenToggleButtonFS->setMaximumSize(buttonSize) ;
	audioListenToggleButtonFS->setCheckable(true);
	audioListenToggleButtonFS->setAutoRaise(true) ;
	audioListenToggleButtonFS->setText(QString()) ;
	audioListenToggleButtonFS->setToolTip(tr("Mute"));

	audioCaptureToggleButtonFS = new QToolButton(fullScreenFrame) ;
	audioCaptureToggleButtonFS->setIcon(iconaudioCaptureToggleButton) ;
	audioCaptureToggleButtonFS->setIconSize(iconSize) ;
	audioCaptureToggleButtonFS->setMinimumSize(buttonSize) ;
	audioCaptureToggleButtonFS->setMaximumSize(buttonSize) ;
	audioCaptureToggleButtonFS->setCheckable(true) ;
	audioCaptureToggleButtonFS->setAutoRaise(true) ;
	audioCaptureToggleButtonFS->setText(QString()) ;
	audioCaptureToggleButtonFS->setToolTip(tr("Start Call"));

	videoCaptureToggleButtonFS = new QToolButton(fullScreenFrame) ;
	videoCaptureToggleButtonFS->setIcon(iconvideoCaptureToggleButton) ;
	videoCaptureToggleButtonFS->setIconSize(iconSize) ;
	videoCaptureToggleButtonFS->setMinimumSize(buttonSize) ;
	videoCaptureToggleButtonFS->setMaximumSize(buttonSize) ;
	videoCaptureToggleButtonFS->setCheckable(true) ;
	videoCaptureToggleButtonFS->setAutoRaise(true) ;
	videoCaptureToggleButtonFS->setText(QString()) ;
	videoCaptureToggleButtonFS->setToolTip(tr("Start Video Call"));

	hangupButtonFS = new QToolButton(fullScreenFrame) ;
	hangupButtonFS->setIcon(QIcon(":/images/phone_hangup.png")) ;
	hangupButtonFS->setIconSize(iconSize) ;
	hangupButtonFS->setMinimumSize(buttonSize) ;
	hangupButtonFS->setMaximumSize(buttonSize) ;
	hangupButtonFS->setCheckable(false) ;
	hangupButtonFS->setAutoRaise(true) ;
	hangupButtonFS->setText(QString()) ;
	hangupButtonFS->setToolTip(tr("Hangup Call"));
	hangupButtonFS->hide();

	fullscreenToggleButtonFS = new QToolButton(fullScreenFrame);
	fullscreenToggleButtonFS->setIcon(iconfullscreenToggleButton);
	fullscreenToggleButtonFS->setIconSize(iconSize);
	fullscreenToggleButtonFS->setMinimumSize(buttonSize);
	fullscreenToggleButtonFS->setMaximumSize(buttonSize);
	fullscreenToggleButtonFS->setCheckable(true);
	fullscreenToggleButtonFS->setAutoRaise(true);
	fullscreenToggleButtonFS->setText(QString());
	fullscreenToggleButtonFS->setToolTip(tr("Fullscreen mode"));
	fullscreenToggleButtonFS->setEnabled(false);

	micMuteButtonFS = makeCallBarButton(fullScreenFrame, iconMicMute, tr("Mute the microphone")) ;
	cameraToggleButtonFS = makeCallBarButton(fullScreenFrame, iconCamera, tr("Turn the camera off")) ;
	connect(micMuteButtonFS, SIGNAL(clicked()), this , SLOT(toggleMicMuteFS()));
	connect(cameraToggleButtonFS, SIGNAL(clicked()), this , SLOT(toggleCameraFS()));

	connect(audioListenToggleButtonFS, SIGNAL(clicked()), this , SLOT(toggleAudioListenFS()));
	connect(audioCaptureToggleButtonFS, SIGNAL(clicked()), this , SLOT(toggleAudioCaptureFS()));
	connect(videoCaptureToggleButtonFS, SIGNAL(clicked()), this , SLOT(toggleVideoCaptureFS()));
	connect(hangupButtonFS, SIGNAL(clicked()), this , SLOT(hangupOrDecline()));
	connect(fullscreenToggleButtonFS, SIGNAL(clicked()), this , SLOT(toggleFullScreenFS()));

	toolBarFSLayout->setDirection(QBoxLayout::LeftToRight);
	toolBarFSLayout->setSpacing(2);
	toolBarFSLayout->setSizeConstraint(QLayout::SetMinAndMaxSize);
	toolBarFSLayout->addWidget(audioCaptureToggleButtonFS);
	toolBarFSLayout->addWidget(videoCaptureToggleButtonFS);
	toolBarFSLayout->addWidget(micMuteButtonFS);
	toolBarFSLayout->addWidget(audioListenToggleButtonFS);
	toolBarFSLayout->addWidget(cameraToggleButtonFS);
	toolBarFSLayout->addWidget(fullscreenToggleButtonFS);
	toolBarFSLayout->addWidget(hangupButtonFS);
	toolBarFS->setLayout(toolBarFSLayout);

	fullScreenFrame->setParent(0);
#if defined(Q_OS_MACOS)
	// No native full screen on macOS: it puts the window in a Space of its own,
	// and hiding the window without leaving that Space first left a black
	// screen. A frameless window over the whole screen is shown and hidden
	// instead (toggleFullScreen()).
	fullScreenFrame->setWindowFlags(Qt::Window | Qt::FramelessWindowHint | Qt::WindowStaysOnTopHint);
#else
	fullScreenFrame->setWindowFlags(Qt::WindowStaysOnTopHint);
	fullScreenFrame->setWindowState(Qt::WindowFullScreen);
#endif
	fullScreenFrame->setFocusPolicy(Qt::StrongFocus);
	fullScreenFrame->hide();
	fullScreenFrame->installEventFilter(this);

	mChatWidget->addChatHorizontalWidget(videoWidget) ;

	inputVideoDevice->setEchoVideoTarget(echoVideoDevice) ;
	inputVideoDevice->setVideoProcessor(videoProcessor) ;
	videoProcessor->setDisplayTarget(outputVideoDevice) ;

	//Ring
	pbAudioRing = new QProgressBar();
	pbAudioRing->setOrientation(Qt::Horizontal);
	pbAudioRing->setRange(0, 99);
	pbAudioRing->setTextVisible(false);
	pbAudioRing->setHidden(true);
	pbVideoRing = new QProgressBar();
	pbVideoRing->setOrientation(Qt::Horizontal);
	pbVideoRing->setRange(0, 99);
	pbVideoRing->setTextVisible(false);
	pbVideoRing->setHidden(true);
	mChatWidget->addChatBarWidget(pbAudioRing);
	mChatWidget->addChatBarWidget(pbVideoRing);

	sendAudioRingTime = -1;
	sendVideoRingTime = -1;
	recAudioRingTime = -1;
	recVideoRingTime = -1;
	updateCallButtons();

	timerAudioRing = new QTimer(this);
	timerAudioRing->setInterval(300);
	timerAudioRing->setSingleShot(true);
	connect(timerAudioRing, SIGNAL(timeout()), this, SLOT(timerAudioRingTimeOut()));
	timerVideoRing = new QTimer(this);
	timerVideoRing->setInterval(300);
	timerVideoRing->setSingleShot(true);
	connect(timerVideoRing, SIGNAL(timeout()), this, SLOT(timerVideoRingTimeOut()));

	lastTimePlayOccurs = time(NULL);
}

VOIPChatWidgetHolder::~VOIPChatWidgetHolder()
{
	hangupCall();

	if(inputAudioDevice != NULL)
		inputAudioDevice->stop() ;

	delete inputVideoDevice ;
	delete videoProcessor ;
	deleteButtonMap();

	// stop and delete timers
	timerAudioRing->stop();
	delete(timerAudioRing);
	timerVideoRing->stop();
	delete(timerVideoRing);
}

void VOIPChatWidgetHolder::deleteButtonMap(int flags)
{
	button_map::iterator it = buttonMapTakeCall.begin();
	while (it != buttonMapTakeCall.end()) {
		if (((it.key().left(1) == "a") && (flags & RS_VOIP_FLAGS_AUDIO_DATA))
		    || ((it.key().left(1) == "v") && (flags & RS_VOIP_FLAGS_VIDEO_DATA)) ) {
			QPair<RSButtonOnText*,RSButtonOnText*> pair = it.value();
			delete pair.second;
			delete pair.first;
			if (flags & RS_VOIP_FLAGS_AUDIO_DATA) recAudioRingTime = -1;
			if (flags & RS_VOIP_FLAGS_VIDEO_DATA) recVideoRingTime = -1;
			it = buttonMapTakeCall.erase(it);
		} else {
			++it;
		}
	}
}

void VOIPChatWidgetHolder::addNewAudioButtonMap(const RsPeerId &peer_id)
{
	if (!mChatWidget)
		return;

	// Answered or declined with the buttons at the top of the window, green and
	// red while it rings. The chat only keeps a line of it.
	const bool already_ringing = recAudioRingTime >= 0;
	recAudioRingTime = 0;
	timerAudioRingTimeOut();

	if (!already_ringing)
		addVoipStatusMessage(tr("%1 is calling you (audio call). Answer or decline with the buttons at the top of the window.").arg(peerName(peer_id)));

	updateCallButtons();
}

void VOIPChatWidgetHolder::addNewVideoButtonMap(const RsPeerId &peer_id)
{
	if (!mChatWidget)
		return;

	// Answered or declined with the buttons at the top of the window, green and
	// red while it rings. The chat only keeps a line of it.
	const bool already_ringing = recVideoRingTime >= 0;
	recVideoRingTime = 0;
	timerVideoRingTimeOut();

	if (!already_ringing)
		addVoipStatusMessage(tr("%1 is calling you (video call). Answer or decline with the buttons at the top of the window.").arg(peerName(peer_id)));

	updateCallButtons();
}

bool VOIPChatWidgetHolder::eventFilter(QObject *obj, QEvent *event)
{
	if (obj == fullScreenFrame) {
		if (event->type() == QEvent::Close || event->type() == QEvent::MouseButtonDblClick) {
			showNormalView();
		}
		// Escape leaves full screen, as everywhere else.
		if (event->type() == QEvent::KeyPress && static_cast<QKeyEvent*>(event)->key() == Qt::Key_Escape) {
			showNormalView();
			return true;
		}
		if (event->type() == QEvent::Resize) {
			replaceFullscreenWidget();
		}

	}
	// pass the event on to the parent class
	return QObject::eventFilter(obj, event);
}

void VOIPChatWidgetHolder::hangupCall()
{
	hangupCallAudio();
	hangupCallVideo();
	hangupButton->hide();
	hangupButtonFS->hide();
	deleteButtonMap();
}

void VOIPChatWidgetHolder::hangupCallAudio()
{
	bool atLeastOneChecked = false;
	if (audioCaptureToggleButton->isChecked()) {
		audioCaptureToggleButton->setChecked(false);
		toggleAudioCapture();
		atLeastOneChecked = true;
	}
	if (!atLeastOneChecked) {
		//Decline button ,Friend hang up or chat close
		if (recAudioRingTime != -1) {
			rsVOIP->sendVoipHangUpCall(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_AUDIO_DATA);
			deleteButtonMap(RS_VOIP_FLAGS_AUDIO_DATA);
		}
		sendAudioRingTime = -1;
		recAudioRingTime = -1;
	}
	updateCallButtons();
}

void VOIPChatWidgetHolder::hangupCallVideo()
{
	bool atLeastOneChecked = false;
	if (videoCaptureToggleButton->isChecked()) {
		videoCaptureToggleButton->setChecked(false);
		toggleVideoCapture();
		atLeastOneChecked = true;
	}
	if (fullscreenToggleButton->isChecked()) {
		fullscreenToggleButton->setChecked(false);
		toggleFullScreen();
		atLeastOneChecked = true;
	}
	if (hideChatTextToggleButton->isChecked()) {
		hideChatTextToggleButton->setChecked(false);
		toggleHideChatText();
		atLeastOneChecked = true;
	}
	if (!atLeastOneChecked) {
		//Decline button ,Friend hang up or chat close
		if (recVideoRingTime != -1) {
			rsVOIP->sendVoipHangUpCall(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_VIDEO_DATA);
			deleteButtonMap(RS_VOIP_FLAGS_VIDEO_DATA);
		}
		sendVideoRingTime = -1;
		recVideoRingTime = -1;
	}
	updateCallButtons();
}

void VOIPChatWidgetHolder::toggleAudioListenFS()
{
	audioListenToggleButton->setChecked(audioListenToggleButtonFS->isChecked());
	toggleAudioListen();
}

void VOIPChatWidgetHolder::toggleAudioListen()
{
    // Checked: the speakers play the peer.
    if (outputAudioDevice)
        outputAudioDevice->setVolume(audioListenToggleButton->isChecked() ? 1.0 : 0.0);

    if (audioListenToggleButton->isChecked()) {
        audioListenToggleButton->setToolTip(tr("Mute the speakers"));
    } else {
        audioListenToggleButton->setToolTip(tr("Unmute the speakers"));
        //audioListenToggleButton->setChecked(false);
        /*if (outputAudioDevice) {
            outputAudioDevice->stop();
        }*/
    }
    audioListenToggleButtonFS->setChecked(audioListenToggleButton->isChecked());
    audioListenToggleButtonFS->setToolTip(audioListenToggleButton->toolTip());
}

void VOIPChatWidgetHolder::startAudioCapture()
{
	addVoipStatusMessage(tr("You answered the audio call."));
	recAudioRingTime = -2;
	audioCaptureToggleButton->setChecked(true);
	toggleAudioCapture();
}

void VOIPChatWidgetHolder::toggleAudioCaptureFS()
{
	audioCaptureToggleButton->setChecked(audioCaptureToggleButtonFS->isChecked());
	toggleAudioCapture();
}

void VOIPChatWidgetHolder::toggleAudioCapture()
{
	if (audioCaptureToggleButton->isChecked()) {
		if (recAudioRingTime == -1) {
			if (sendAudioRingTime == -1) {
				sendAudioRingTime = 0;
				timerAudioRingTimeOut();
				rsVOIP->sendVoipRinging(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_AUDIO_DATA);
				return; //Start Audio when accept received
			}
		}
		if (recAudioRingTime != -1)
			rsVOIP->sendVoipAcceptCall(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_AUDIO_DATA);
		recAudioRingTime = -1; //Stop ringing

		//activate buttons
		audioListenToggleButton->setEnabled(true);
		audioListenToggleButton->setChecked(true);
		audioListenToggleButtonFS->setEnabled(true);
		audioListenToggleButtonFS->setChecked(true);
		audioCaptureToggleButton->setToolTip(tr("Hold Call"));
		hangupButton->show();
		hangupButtonFS->show();

		//activate audio input
		if (!inputAudioProcessor) {
			inputAudioProcessor = new QtSpeex::SpeexInputProcessor();
			if (outputAudioProcessor) {
				connect(outputAudioProcessor, SIGNAL(playingFrame(QByteArray*)), inputAudioProcessor, SLOT(addEchoFrame(QByteArray*)));
			}
			inputAudioProcessor->open(QIODevice::WriteOnly | QIODevice::Unbuffered);
		}
		if (!inputAudioDevice) {
			inputAudioDevice = AudioDeviceHelper::getPreferedInputDevice();
		}
		connect(inputAudioProcessor, SIGNAL(networkPacketReady()), this, SLOT(sendAudioData()));
		inputAudioDevice->start(inputAudioProcessor);
		checkAudioDevicesAtCallStart();


		deleteButtonMap(RS_VOIP_FLAGS_AUDIO_DATA);
	} else {
		//desactivate buttons
		audioListenToggleButton->setEnabled(false);
		audioListenToggleButton->setChecked(false);
		audioListenToggleButtonFS->setEnabled(false);
		audioListenToggleButtonFS->setChecked(false);
		audioCaptureToggleButton->setToolTip(tr("Resume Call"));
		if (!videoCaptureToggleButton->isChecked()) {
			hangupButton->hide();
			hangupButtonFS->hide();
		}

		if (recAudioRingTime <= -1){
			//desactivate audio input
			disconnect(inputAudioProcessor, SIGNAL(networkPacketReady()), this, SLOT(sendAudioData()));
			if (inputAudioDevice) {
				inputAudioDevice->stop();
			}


			rsVOIP->sendVoipHangUpCall(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_AUDIO_DATA);
		}

		sendAudioRingTime = -1;
		recAudioRingTime = -1;

	}
	audioCaptureToggleButtonFS->setChecked(audioCaptureToggleButton->isChecked());
	audioCaptureToggleButtonFS->setToolTip(audioCaptureToggleButton->toolTip());
	updateCallButtons();
}

void VOIPChatWidgetHolder::startVideoCapture()
{
	addVoipStatusMessage(startWithCameraOff ? tr("You answered the video call, with your camera off.")
	                                        : tr("You answered the video call."));
	recVideoRingTime = -2;
	videoCaptureToggleButton->setChecked(true);
	toggleVideoCapture();
}

void VOIPChatWidgetHolder::toggleVideoCaptureFS()
{
	videoCaptureToggleButton->setChecked(videoCaptureToggleButtonFS->isChecked());
	toggleVideoCapture();
}

void VOIPChatWidgetHolder::toggleVideoCapture()
{
	if (videoCaptureToggleButton->isChecked()) 
	{
		if (recVideoRingTime == -1) {
			if (sendVideoRingTime == -1) {
				sendVideoRingTime = 0;
				timerVideoRingTimeOut();
				rsVOIP->sendVoipRinging(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_VIDEO_DATA);
				return; //Start Video when accept received
			}
		}
		if (recVideoRingTime != -1)
			rsVOIP->sendVoipAcceptCall(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_VIDEO_DATA);
		recVideoRingTime = -1; //Stop ringing

		//activate buttons
		hideChatTextToggleButton->setEnabled(true);
		fullscreenToggleButton->setEnabled(true);
		fullscreenToggleButtonFS->setEnabled(true);
		videoCaptureToggleButton->setToolTip(tr("Shut camera off"));
		hangupButton->show();
    hangupButtonFS->show();

		//activate video input
		videoWidget->show();
		videoProcessor->setMaximumBandwidth(rsVOIP->getVoipVideoMaxBandwidth()) ;	// user setting (wizard / settings page)
		if (!startWithCameraOff)
			inputVideoDevice->start(QVideoInputDevice::preferredCamera()) ;
		videoCallStarted();


		deleteButtonMap(RS_VOIP_FLAGS_VIDEO_DATA);
	} else {
		//desactivate buttons
		hideChatTextToggleButton->setEnabled(false);
		hideChatTextToggleButton->setChecked(false);
		toggleHideChatText();
		fullscreenToggleButton->setEnabled(false);
		fullscreenToggleButton->setChecked(false);
		fullscreenToggleButtonFS->setEnabled(false);
		fullscreenToggleButtonFS->setChecked(false);
		toggleFullScreen();
		videoCaptureToggleButton->setToolTip(tr("Activate camera"));
		if (!audioCaptureToggleButton->isChecked()) {
			hangupButton->hide();
			hangupButtonFS->hide();
		}

		if (recVideoRingTime<=-1){
			//desactivate video input
			inputVideoDevice->stop() ;
			outputVideoDevice->showFrameOff();
			videoWidget->hide();
			videoCallEnded();


			rsVOIP->sendVoipHangUpCall(mChatWidget->getChatId().toPeerId(), RS_VOIP_FLAGS_VIDEO_DATA);
		}

		sendVideoRingTime = -1;
		recVideoRingTime = -1;

	}
	videoCaptureToggleButtonFS->setChecked(videoCaptureToggleButton->isChecked());
	videoCaptureToggleButtonFS->setToolTip(videoCaptureToggleButton->toolTip());
	updateCallButtons();
}

void VOIPChatWidgetHolder::addVideoData(const RsPeerId &/*peer_id*/, QByteArray* array)
{
	// Only the peer's ring invites us. Video outside a call is what the peer
	// sent before it saw our hang-up (a tunnel holds many such frames): drop it.
	if (!videoCaptureToggleButton->isChecked())
		return;

	sendVideoRingTime = -2;//Receive Video so Accepted

	RsVOIPDataChunk chunk ;
	chunk.type = RsVOIPDataChunk::RS_VOIP_DATA_TYPE_VIDEO ;
	chunk.size = array->size() ;
	chunk.data = array->data() ;

	videoProcessor->receiveEncodedData(chunk) ;

}

void VOIPChatWidgetHolder::toggleHideChatText()
{
	QBoxLayout *layout = static_cast<QBoxLayout*>(videoWidget->layout());

	if (hideChatTextToggleButton->isChecked()) {
		mChatWidget->hideChatText(true);
		if (layout) layout->setDirection(QBoxLayout::LeftToRight);
		hideChatTextToggleButton->setToolTip(tr("Show Chat Text"));
	} else {
		mChatWidget->hideChatText(false);
		if (layout) layout->setDirection(QBoxLayout::TopToBottom);
		hideChatTextToggleButton->setToolTip(tr("Hide Chat Text"));
		fullscreenToggleButton->setChecked(false);
		toggleFullScreen();
	}
}

void VOIPChatWidgetHolder::toggleFullScreenFS()
{
	fullscreenToggleButton->setChecked(fullscreenToggleButtonFS->isChecked());
	toggleFullScreen();
}

void VOIPChatWidgetHolder::toggleFullScreen()
{
	if (fullscreenToggleButton->isChecked()) {
		fullscreenToggleButton->setToolTip(tr("Return to normal view."));
		inputVideoDevice->setEchoVideoTarget(echoVideoDeviceFS) ;
		videoProcessor->setDisplayTarget(outputVideoDeviceFS) ;
#if defined(Q_OS_MACOS)
		// Cover the screen the chat window is on.
		QScreen *screen = QGuiApplication::screenAt(mChatWidget->mapToGlobal(mChatWidget->rect().center()));
		if (screen == NULL)
			screen = QGuiApplication::primaryScreen();
		fullScreenFrame->setGeometry(screen->geometry());
		fullScreenFrame->show();
		fullScreenFrame->raise();
		fullScreenFrame->activateWindow();
#else
		fullScreenFrame->show();
#endif
	} else {
		mChatWidget->hideChatText(false);
		fullscreenToggleButton->setToolTip(tr("Fullscreen mode"));
		inputVideoDevice->setEchoVideoTarget(echoVideoDevice) ;
		videoProcessor->setDisplayTarget(outputVideoDevice) ;
		fullScreenFrame->hide();
	}
	fullscreenToggleButtonFS->setChecked(fullscreenToggleButton->isChecked());
	fullscreenToggleButtonFS->setToolTip(fullscreenToggleButton->toolTip());
}

void VOIPChatWidgetHolder::replaceFullscreenWidget()
{
	if (QSize(toolBarFS->geometry().size() - fullScreenFrame->geometry().size()).isValid()){
		QRect fsRect = fullScreenFrame->geometry();
		fsRect.setSize(toolBarFS->geometry().size());
		fullScreenFrame->setGeometry(fsRect);
	}

	outputVideoDeviceFS->setGeometry(QRect(QPoint(0,0),fullScreenFrame->geometry().size()));
	echoVideoDeviceFS->setGeometry(QRect(QPoint(fullScreenFrame->width(), fullScreenFrame->height()) - QPoint(320,240), QSize(320,240)));
	QRect toolBarFSGeo = QRect( (fullScreenFrame->width() - toolBarFS->geometry().width()) / 2
	                            , fullScreenFrame->height() - (toolBarFS->geometry().height() * 2)
	                            , toolBarFS->geometry().width(), toolBarFS->geometry().height());
	toolBarFS->setGeometry(toolBarFSGeo);

	if (!videoCaptureToggleButton->isChecked()) {
		outputVideoDeviceFS->showFrameOff();
		echoVideoDeviceFS->showFrameOff();
	}
}

void VOIPChatWidgetHolder::showNormalView()
{
	hideChatTextToggleButton->setChecked(false);
	toggleHideChatText();
	fullscreenToggleButton->setChecked(false);
	fullscreenToggleButtonFS->setChecked(false);
	toggleFullScreen();
}

void VOIPChatWidgetHolder::botMouseEnterTake()
{
	RSButtonOnText *source = qobject_cast<RSButtonOnText *>(QObject::sender());
	if (source){
		source->setStyleSheet(QString("border: 1px solid #333333;")
		                      .append("font-size: 12pt; color: white;")
		                      .append("min-width: 128px; min-height: 24px;")
		                      .append("border-radius: 6px;")
		                      .append("padding: 3px;")
		                      .append("background-color: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 0.67, "
	                                  "stop: 0 #444444, stop: 1 #222222);")

		                      );
		//source->setDown(true);
	}
}

void VOIPChatWidgetHolder::botMouseLeaveTake()
{
	RSButtonOnText *source = qobject_cast<RSButtonOnText *>(QObject::sender());
	if (source){
		source->setStyleSheet(QString("border: 1px solid #116a06;")
		                      .append("font-size: 12pt; color: white;")
		                      .append("min-width: 128px; min-height: 24px;")
		                      .append("border-radius: 6px;")
		                      .append("padding: 3px;")		                      
		                      .append("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 rgba(153, 240, 143, 255)," 
   	                               "stop:0.494318 rgba(59, 201, 40, 255), stop:0.505682 rgba(46, 172, 29, 255), stop:1 rgba(30, 116, 20, 255));")
		                      );
		//source->setDown(false);
	}
}

void VOIPChatWidgetHolder::botMouseEnterDecline()
{
	RSButtonOnText *source = qobject_cast<RSButtonOnText *>(QObject::sender());
	if (source){
		source->setStyleSheet(QString("border: 1px solid #333333;")
		                      .append("font-size: 12pt; color: white;")
		                      .append("min-width: 128px; min-height: 24px;")
		                      .append("border-radius: 6px;")
		                      .append("padding: 3px;")
		                      .append("background-color:  qlineargradient(x1: 0, y1: 0, x2: 0, y2: 0.67, "
	                                  "stop: 0 #444444, stop: 1 #222222);")
		                      );
		//source->setDown(true);
	}
}

void VOIPChatWidgetHolder::botMouseLeaveDecline()
{
	RSButtonOnText *source = qobject_cast<RSButtonOnText *>(QObject::sender());
	if (source){
		source->setStyleSheet(QString("border: 1px solid #6a1106;")
		                      .append("font-size: 12pt; color: white;")
		                      .append("min-width: 128px; min-height: 24px;")
		                      .append("border-radius: 6px;")			                       
		                      .append("padding: 3px;")
	                          .append("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:0, y2:1, stop:0 rgba(240, 154, 143, 255), "
	                                  "stop:0.494318 rgba(201, 57, 40, 255), stop:0.505682 rgba(172, 45, 29, 255), stop:1 rgba(116, 30, 20, 255));")
		                      );
		//source->setDown(false);
	}
}

void VOIPChatWidgetHolder::setAcceptedBandwidth(uint32_t bytes_per_sec)
{
	// The peer reports the video bandwidth it actually received. That is just what we
	// sent (smoothed, and starting from zero), not what the connection can carry, so
	// using it as the encoder limit only ever lowered the quality. The encoder uses
	// the bandwidth chosen by the user instead.
	std::cerr << "VOIP: peer reports " << bytes_per_sec << " B/s of incoming video (target: "
	          << rsVOIP->getVoipVideoMaxBandwidth() << " B/s)" << std::endl;
}

void VOIPChatWidgetHolder::addVoipStatusMessage(const QString& text)
{
	if (mChatWidget)
		mChatWidget->addChatMsg(true, tr("VoIP Status"), QDateTime::currentDateTime(), QDateTime::currentDateTime(), text, ChatWidget::MSGTYPE_SYSTEM);
}

QString VOIPChatWidgetHolder::peerName(const RsPeerId& peer_id)
{
	// In a distant chat the correspondent is an identity and peer_id a tunnel id
	// that rsPeers does not know. Its name is looked up now: the window title is
	// set when the window opens, to the identity id if the name was not loaded.
	if (mChatWidget && mChatWidget->getChatId().isDistantChatId())
	{
		DistantChatPeerInfo info;
		RsIdentityDetails details;

		if (rsChats->getDistantChatStatus(mChatWidget->getChatId().toDistantChatId(), info)
		        && rsIdentity->getIdDetails(info.to_id, details) && !details.mNickname.empty())
			return QString::fromUtf8(details.mNickname.c_str());

		return mChatWidget->getTitle();
	}

	QString name;
	if (!peer_id.isNull()) name = QString::fromUtf8(rsPeers->getPeerName(peer_id).c_str());
	if (name.isEmpty() && mChatWidget) name = mChatWidget->getTitle();
	if (name.isEmpty() && !peer_id.isNull()) name = QString::fromStdString(peer_id.toStdString());
	return name;
}

void VOIPChatWidgetHolder::showCameraProblem(const QString& msg)
{
	addVoipStatusMessage(tr("Cannot use the camera: %1").arg(msg));
}

void VOIPChatWidgetHolder::showAudioInputProblem(const QString& text)
{
	if (audioInputProblemShown)
		return;

	audioInputProblemShown = true;
	addVoipStatusMessage(text);
}

void VOIPChatWidgetHolder::checkAudioDevicesAtCallStart()
{
	audioInputProblemShown = false;
	audioOutputProblemShown = false;	// the speakers open with the first sound received
	micMuteButton->setChecked(false);
	micMuteButtonFS->setChecked(false);
	micMuteButton->setToolTip(tr("Mute the microphone"));
	micMuteButtonFS->setToolTip(micMuteButton->toolTip());

	if (inputAudioDevice->state() == QAudio::StoppedState && inputAudioDevice->error() != QAudio::NoError)
	{
		showAudioInputProblem(tr("Cannot open the microphone. Check the audio input device in the VOIP settings."));
		return;
	}
	audioInputBytesAtStart = inputAudioProcessor->bytesReceived();
	timerAudioInputCheck->start();
}

void VOIPChatWidgetHolder::audioInputCheckTimeOut()
{
	// A microphone held by another application, or to which the system refuses
	// access, can open without any error and still deliver nothing.
	if (inputAudioDevice == NULL || inputAudioProcessor == NULL || !audioCaptureToggleButton->isChecked()
	        || micMuteButton->isChecked())
		return;

	if (inputAudioProcessor->bytesReceived() == audioInputBytesAtStart)
		showAudioInputProblem(tr("The microphone sends no sound. Another application may be using it, or the system may refuse access to it."));
}

void VOIPChatWidgetHolder::updateCallButtons()
{
	// Ringing out: only the hang-up button, which cancels the call. Once the
	// call is up: microphone and camera on/off, full screen for a video call,
	// hang-up -- what every chat application offers. The call buttons hold the
	// call state and stay checked while hidden; the video call button stays
	// offered during an audio-only call. The speakers and chat-text buttons are
	// left out.
	const bool ringingOut = sendAudioRingTime >= 0 || sendVideoRingTime >= 0;
	const bool audioOn = audioCaptureToggleButton->isChecked() && sendAudioRingTime < 0;
	const bool videoOn = videoCaptureToggleButton->isChecked() && sendVideoRingTime < 0;
	const bool inCall = audioOn || videoOn;
	// Ringing in: the call button of that kind answers (green), the red one
	// declines. A video call can also be answered without our camera.
	const bool audioRingingIn = !audioCaptureToggleButton->isChecked() && recAudioRingTime >= 0;
	const bool videoRingingIn = !videoCaptureToggleButton->isChecked() && recVideoRingTime >= 0;
	const bool ringingIn = audioRingingIn || videoRingingIn;

	const QString green("QToolButton { background-color: #2fa84f; border-radius: 6px; }");
	audioCaptureToggleButton->setStyleSheet(audioRingingIn ? green : QString());
	videoCaptureToggleButton->setStyleSheet(videoRingingIn ? green : QString());
	answerNoCameraButton->setStyleSheet(green);

	audioCaptureToggleButton->setVisible(!inCall && !ringingOut && !videoRingingIn);
	videoCaptureToggleButton->setVisible(!videoOn && !ringingOut && !audioRingingIn);
	answerNoCameraButton->setVisible(videoRingingIn);
	micMuteButton->setVisible(audioOn);
	audioListenToggleButton->setVisible(false);
	cameraToggleButton->setVisible(videoOn);
	hideChatTextToggleButton->setVisible(false);
	fullscreenToggleButton->setVisible(videoOn);
	hangupButton->setVisible(inCall || ringingOut || ringingIn);
	hangupButton->setToolTip((ringingIn && !inCall) ? tr("Decline") : tr("Hangup Call"));

	audioCaptureToggleButtonFS->setVisible(!inCall && !ringingOut);
	videoCaptureToggleButtonFS->setVisible(!videoOn && !ringingOut);
	micMuteButtonFS->setVisible(audioOn);
	audioListenToggleButtonFS->setVisible(false);
	cameraToggleButtonFS->setVisible(videoOn);
	fullscreenToggleButtonFS->setVisible(videoOn);
	hangupButtonFS->setVisible(inCall || ringingOut);
}

void VOIPChatWidgetHolder::videoCallStarted()
{
	// Checked camera button: our camera is off.
	cameraToggleButton->setChecked(startWithCameraOff);
	cameraToggleButtonFS->setChecked(startWithCameraOff);
	cameraToggleButton->setToolTip(startWithCameraOff ? tr("Turn the camera on") : tr("Turn the camera off"));
	cameraToggleButtonFS->setToolTip(cameraToggleButton->toolTip());
	startWithCameraOff = false;

	if (audioCaptureToggleButton->isChecked())
		return;	// already in an audio call: the video was added to it

	// Both ends start their sound when the video call starts; the answer to the
	// video ring stands for the audio one.
	audioStartedWithVideo = true;
	audioCaptureToggleButton->setChecked(true);
	sendAudioRingTime = -2;
	toggleAudioCapture();
}

void VOIPChatWidgetHolder::videoCallEnded()
{
	if (!audioStartedWithVideo)
		return;

	audioStartedWithVideo = false;
	if (audioCaptureToggleButton->isChecked())
	{
		audioCaptureToggleButton->setChecked(false);
		toggleAudioCapture();
	}
}

void VOIPChatWidgetHolder::toggleMicMute()
{
	// The call goes on: the microphone stops delivering, so nothing is sent.
	if (inputAudioDevice)
	{
		if (micMuteButton->isChecked())
			inputAudioDevice->suspend();
		else
			inputAudioDevice->resume();
	}
	micMuteButton->setToolTip(micMuteButton->isChecked() ? tr("Unmute the microphone") : tr("Mute the microphone"));
	micMuteButtonFS->setChecked(micMuteButton->isChecked());
	micMuteButtonFS->setToolTip(micMuteButton->toolTip());
}

void VOIPChatWidgetHolder::toggleMicMuteFS()
{
	micMuteButton->setChecked(micMuteButtonFS->isChecked());
	toggleMicMute();
}

void VOIPChatWidgetHolder::toggleCamera()
{
	// The call goes on: the peer's picture still shows, ours stops.
	if (cameraToggleButton->isChecked())
		inputVideoDevice->stop();
	else
		inputVideoDevice->start(QVideoInputDevice::preferredCamera());

	cameraToggleButton->setToolTip(cameraToggleButton->isChecked() ? tr("Turn the camera on") : tr("Turn the camera off"));
	cameraToggleButtonFS->setChecked(cameraToggleButton->isChecked());
	cameraToggleButtonFS->setToolTip(cameraToggleButton->toolTip());
}

void VOIPChatWidgetHolder::declineAudioCall()
{
	addVoipStatusMessage(tr("You declined the audio call."));
	hangupCallAudio();
}

void VOIPChatWidgetHolder::declineVideoCall()
{
	addVoipStatusMessage(tr("You declined the video call."));
	hangupCallVideo();
}

void VOIPChatWidgetHolder::hangupOrDecline()
{
	const bool inCall = (audioCaptureToggleButton->isChecked() && sendAudioRingTime < 0)
	                 || (videoCaptureToggleButton->isChecked() && sendVideoRingTime < 0);
	const bool audioRingingIn = !audioCaptureToggleButton->isChecked() && recAudioRingTime >= 0;
	const bool videoRingingIn = !videoCaptureToggleButton->isChecked() && recVideoRingTime >= 0;

	if (!inCall && (audioRingingIn || videoRingingIn))
	{
		if (audioRingingIn) declineAudioCall();
		if (videoRingingIn) declineVideoCall();
		return;
	}
	addVoipStatusMessage(inCall ? tr("You ended the call.") : tr("You cancelled the call."));
	hangupCall();
}

void VOIPChatWidgetHolder::audioCallButtonClicked()
{
	// Clicking has already checked the button.
	if (audioCaptureToggleButton->isChecked() && recAudioRingTime >= 0)
		startAudioCapture();	// it rang in: answer
	else
		toggleAudioCapture();
}

void VOIPChatWidgetHolder::videoCallButtonClicked()
{
	if (videoCaptureToggleButton->isChecked() && recVideoRingTime >= 0)
		startVideoCapture();	// it rang in: answer
	else
		toggleVideoCapture();
}

void VOIPChatWidgetHolder::answerWithoutCamera()
{
	startWithCameraOff = true;
	startVideoCapture();	// answers the video call; videoCallStarted() leaves the camera off
}

void VOIPChatWidgetHolder::toggleCameraFS()
{
	cameraToggleButton->setChecked(cameraToggleButtonFS->isChecked());
	toggleCamera();
}

void VOIPChatWidgetHolder::addAudioData(const RsPeerId &peer_id, QByteArray* array)
{
    // Only the peer's ring invites us. Sound outside a call is what the peer
    // sent before it saw our hang-up: drop it.
    if (!audioCaptureToggleButton->isChecked())
        return;

    sendAudioRingTime = -2;//Receive Audio so Accepted

    if (!outputAudioDevice) {
        outputAudioDevice = AudioDeviceHelper::getPreferedOutputDevice();
    }

    if (!outputAudioProcessor) {
        //start output audio device
        outputAudioProcessor = new QtSpeex::SpeexOutputProcessor();
        if (inputAudioProcessor) {
            connect(outputAudioProcessor, SIGNAL(playingFrame(QByteArray*)), inputAudioProcessor, SLOT(addEchoFrame(QByteArray*)));
        }
        outputAudioProcessor->open(QIODevice::ReadOnly | QIODevice::Unbuffered);
        outputAudioDevice->start(outputAudioProcessor);
        outputAudioDevice->setVolume(audioListenToggleButton->isChecked() ? 1.0 : 0.0);
    }

    if (outputAudioDevice && outputAudioDevice->error() != QAudio::NoError) {
        std::cerr << "Restarting output device. Error before reset " << outputAudioDevice->error() << " buffer size : " << outputAudioDevice->bufferSize() << std::endl;
        outputAudioDevice->stop();
        outputAudioDevice->reset();
        if (outputAudioDevice->error() == QAudio::UnderrunError)
            outputAudioDevice->setBufferSize(20);
        outputAudioDevice->start(outputAudioProcessor);
    }

    if (outputAudioDevice && !audioOutputProblemShown
            && outputAudioDevice->state() == QAudio::StoppedState && outputAudioDevice->error() != QAudio::NoError)
    {
        audioOutputProblemShown = true;
        addVoipStatusMessage(tr("Cannot open the speakers or headphones. Check the audio output device in the VOIP settings."));
    }
    outputAudioProcessor->putNetworkPacket(QString::fromStdString(peer_id.toStdString()), *array);

    //check the input device for errors
    if (inputAudioDevice && inputAudioDevice->error() != QAudio::NoError) {
        std::cerr << "Restarting input device. Error before reset " << inputAudioDevice->error() << std::endl;
        inputAudioDevice->stop();
        inputAudioDevice->reset();
        inputAudioDevice->start(inputAudioProcessor);
    }
}

void VOIPChatWidgetHolder::sendVideoData()
{
	RsVOIPDataChunk chunk ;

	// Switch to H.264 once the peer has said it decodes it. Until then (older
	// plugin, or its answer not received yet) video is sent as MPEG-4.
	if(videoProcessor && (rsVOIP->getVoipPeerCapabilities(mChatWidget->getChatId().toPeerId()) & RS_VOIP_CAPS_H264_VIDEO))
		videoProcessor->setPeerDecodesH264(true) ;

	while(inputVideoDevice && inputVideoDevice->getNextEncodedPacket(chunk))
		rsVOIP->sendVoipData(mChatWidget->getChatId().toPeerId(),chunk) ;
}

void VOIPChatWidgetHolder::sendAudioData()
{
    while(inputAudioProcessor && inputAudioProcessor->hasPendingPackets()) {
        QByteArray qbarray = inputAudioProcessor->getNetworkPacket();
        RsVOIPDataChunk chunk;
        chunk.size = qbarray.size();
        chunk.data = (void*)qbarray.constData();
        chunk.type = RsVOIPDataChunk::RS_VOIP_DATA_TYPE_AUDIO ;
        rsVOIP->sendVoipData(mChatWidget->getChatId().toPeerId(),chunk);
    }
}

void VOIPChatWidgetHolder::updateStatus(RsStatusValue status)
{
	bool enabled = (status != RsStatusValue::RS_STATUS_OFFLINE);

	audioListenToggleButton->setEnabled(audioCaptureToggleButton->isChecked() && enabled);
	audioListenToggleButtonFS->setEnabled(audioCaptureToggleButton->isChecked() && enabled);
	audioCaptureToggleButton->setEnabled(enabled);
	audioCaptureToggleButtonFS->setEnabled(enabled);
	videoCaptureToggleButton->setEnabled(enabled);
	videoCaptureToggleButtonFS->setEnabled(enabled);
	hideChatTextToggleButton->setEnabled(videoCaptureToggleButton->isChecked() && enabled);
	fullscreenToggleButton->setEnabled(videoCaptureToggleButton->isChecked() && enabled);
	fullscreenToggleButtonFS->setEnabled(videoCaptureToggleButton->isChecked() && enabled);
	hangupButton->setEnabled(enabled);
	hangupButtonFS->setEnabled(enabled);
}

void VOIPChatWidgetHolder::ReceivedInvitation(const RsPeerId &peer_id, int flags)
{
	switch(flags){
		case RS_VOIP_FLAGS_AUDIO_DATA: {
			if (audioCaptureToggleButton->isChecked()) {
				if (recAudioRingTime != -1)
					toggleAudioCapture();
			} else {
				addNewAudioButtonMap(peer_id);
			}
		}
		break;
		case RS_VOIP_FLAGS_VIDEO_DATA: {
			if (videoCaptureToggleButton->isChecked()) {
				if (recVideoRingTime != -1)
					toggleVideoCapture();
			} else {
				addNewVideoButtonMap(peer_id);
			}
		}
		break;
		default:
			std::cerr << "VOIPChatWidgetHolder::ReceivedInvitation(): Received unknown flags item # " << flags << ": not handled yet ! Sorry" << std::endl;
		break ;
	}
}

void VOIPChatWidgetHolder::ReceivedVoipHangUp(const RsPeerId &peer_id, int flags)
{
	switch(flags){
		case RS_VOIP_FLAGS_AUDIO_DATA | RS_VOIP_FLAGS_VIDEO_DATA: {
			if (mChatWidget) {
				if (!videoCaptureToggleButton->isChecked() && !audioCaptureToggleButton->isChecked()
				        && (recAudioRingTime != -1 || recVideoRingTime != -1))
					addVoipStatusMessage(tr("Missed call from %1.").arg(peerName(peer_id)));
				if (videoCaptureToggleButton->isChecked() || audioCaptureToggleButton->isChecked())
					addVoipStatusMessage((sendAudioRingTime >= 0 || sendVideoRingTime >= 0)
					                     ? tr("%1 declined your call.").arg(peerName(peer_id))
					                     : tr("%1 ended the call.").arg(peerName(peer_id)));
				hangupCall();
			}
		}
		break;
		case RS_VOIP_FLAGS_AUDIO_DATA: {
			if (mChatWidget) {
				if (!audioCaptureToggleButton->isChecked() && recAudioRingTime != -1)
					addVoipStatusMessage(tr("Missed call from %1.").arg(peerName(peer_id)));
				// The sound of a video call ends with it: the video hang-up says so.
				if (audioCaptureToggleButton->isChecked() && !audioStartedWithVideo)
					addVoipStatusMessage((sendAudioRingTime >= 0)
					                     ? tr("%1 declined your audio call.").arg(peerName(peer_id))
					                     : tr("%1 ended the call.").arg(peerName(peer_id)));
				hangupCallAudio();
			}
		}
		break;
		case RS_VOIP_FLAGS_VIDEO_DATA: {
			if (mChatWidget) {
				if (!videoCaptureToggleButton->isChecked() && recVideoRingTime != -1)
					addVoipStatusMessage(tr("Missed call from %1.").arg(peerName(peer_id)));
				if (videoCaptureToggleButton->isChecked())
					addVoipStatusMessage((sendVideoRingTime >= 0)
					                     ? tr("%1 declined your video call.").arg(peerName(peer_id))
					                     : tr("%1 ended the call.").arg(peerName(peer_id)));
				hangupCallVideo();
			}
		}
		break;
		default:
			std::cerr << "VOIPChatWidgetHolder::ReceivedVoipHangUp(): Received unknown flags item # " << flags << ": not handled yet ! Sorry" << std::endl;
		break ;
	}
	//deleteButtonMap();
}

void VOIPChatWidgetHolder::ReceivedVoipAccept(const RsPeerId &peer_id, int flags)
{
	switch(flags){
		case RS_VOIP_FLAGS_AUDIO_DATA: {
			if (mChatWidget) {
				sendAudioRingTime = -2;
				QString name = peerName(peer_id);
				mChatWidget->addChatMsg(true, tr("VoIP Status"), QDateTime::currentDateTime(), QDateTime::currentDateTime()
				                        , tr("%1 answered your audio call.").arg(name), ChatWidget::MSGTYPE_SYSTEM);
				if (audioCaptureToggleButton->isChecked())
					toggleAudioCapture();
			}
		}
		break;
		case RS_VOIP_FLAGS_VIDEO_DATA: {
			if (mChatWidget) {
				sendVideoRingTime = -2;
				QString name = peerName(peer_id);
				mChatWidget->addChatMsg(true, tr("VoIP Status"), QDateTime::currentDateTime(), QDateTime::currentDateTime()
				                        , tr("%1 answered your video call.").arg(name), ChatWidget::MSGTYPE_SYSTEM);
				if (videoCaptureToggleButton->isChecked())
					toggleVideoCapture();
			}
		}
		break;
		default:
			std::cerr << "VOIPChatWidgetHolder::ReceivedVoipHangUp(): Received unknown flags item # " << flags << ": not handled yet ! Sorry" << std::endl;
		break ;
	}
}

void VOIPChatWidgetHolder::timerAudioRingTimeOut()
{
	//Sending or receiving (-2 connected, -1 reseted, >=0 in progress)
	if (sendAudioRingTime >= 0) {
		//Sending
		if (pbAudioRing->isHidden())	// first tick of this call
		{
			ringClock.start();
			addVoipStatusMessage(tr("Audio call to %1, waiting for an answer...").arg(peerName()));
			updateCallButtons();	// the hang-up button cancels the call
		}
		if (ringClock.isValid() && ringClock.elapsed() > CALL_RING_TIMEOUT_MS)
		{
			addVoipStatusMessage(tr("%1 did not answer.").arg(peerName()));
			hangupCallAudio();	// the peer reads it as a missed call
			return;
		}
		++sendAudioRingTime;
		if (sendAudioRingTime == 100) sendAudioRingTime = 0;
		pbAudioRing->setValue(sendAudioRingTime);
		pbAudioRing->setToolTip(tr("Waiting for your friend to respond to your audio call."));
		pbAudioRing->setVisible(true);

		if (time(NULL) > lastTimePlayOccurs) {
			soundManager->play(VOIP_SOUND_OUTGOING_AUDIO_CALL);
			lastTimePlayOccurs = time(NULL) + 1;
		}

		timerAudioRing->start();
	} else if(recAudioRingTime >= 0) {
		//Receiving
		if (pbAudioRing->isHidden())	// first tick of this call (the counter wraps every 30 s)
			ringClock.start();
		else if (ringClock.isValid() && ringClock.elapsed() > CALL_RING_TIMEOUT_MS + 10000)
		{
			// The caller gives up first; it may be gone without a word.
			addVoipStatusMessage(tr("Missed call from %1.").arg(peerName()));
			hangupCallAudio();
			return;
		}
		++recAudioRingTime;
		if (recAudioRingTime == 100) recAudioRingTime = 0;
		pbAudioRing->setValue(recAudioRingTime);
		pbAudioRing->setToolTip(tr("Your friend is calling you for audio. Respond."));
		pbAudioRing->setVisible(true);

		//launch an animation. Don't launch it if already animating
		if (!audioCaptureToggleButton->graphicsEffect()
		    || (audioCaptureToggleButton->graphicsEffect()->inherits("QGraphicsOpacityEffect")
		        && ((QGraphicsOpacityEffect*)audioCaptureToggleButton->graphicsEffect())->opacity() == 1)
		    ) {
			QGraphicsOpacityEffect *effect = new QGraphicsOpacityEffect(audioListenToggleButton);
			audioCaptureToggleButton->setGraphicsEffect(effect);
			QPropertyAnimation *anim = new QPropertyAnimation(effect, "opacity", effect);
			anim->setStartValue(1);
			anim->setKeyValueAt(0.5,0);
			anim->setEndValue(1);
			anim->setDuration(timerAudioRing->interval());
			anim->start();
		}
		audioCaptureToggleButton->setToolTip(tr("Answer"));

		if (time(NULL) > lastTimePlayOccurs) {
			soundManager->play(VOIP_SOUND_INCOMING_AUDIO_CALL);
			lastTimePlayOccurs = time(NULL) + 1;
		}

		if (mVOIPNotify) mVOIPNotify->notifyReceivedVoipAudioCall(mChatWidget->getChatId().toPeerId());

		timerAudioRing->start();
	} else {
		//Nothing to do, reset stat
		pbAudioRing->setHidden(true);
		updateCallButtons();
		pbAudioRing->setValue(0);
		pbAudioRing->setToolTip("");
		audioCaptureToggleButton->setGraphicsEffect(0);
	}
}

void VOIPChatWidgetHolder::timerVideoRingTimeOut()
{
	//Sending or receiving (-2 connected, -1 reseted, >=0 in progress)
	if (sendVideoRingTime >= 0) {
		//Sending
		if (pbVideoRing->isHidden())	// first tick of this call
		{
			ringClock.start();
			addVoipStatusMessage(tr("Video call to %1, waiting for an answer...").arg(peerName()));
			updateCallButtons();	// the hang-up button cancels the call
		}
		if (ringClock.isValid() && ringClock.elapsed() > CALL_RING_TIMEOUT_MS)
		{
			addVoipStatusMessage(tr("%1 did not answer.").arg(peerName()));
			hangupCallVideo();	// the peer reads it as a missed call
			return;
		}
		++sendVideoRingTime;
		if (sendVideoRingTime == 100) sendVideoRingTime = 0;
		pbVideoRing->setValue(sendVideoRingTime);
		pbVideoRing->setToolTip(tr("Waiting for your friend to respond to your video call."));
		pbVideoRing->setVisible(true);

		if (time(NULL) > lastTimePlayOccurs) {
			soundManager->play(VOIP_SOUND_OUTGOING_VIDEO_CALL);
			lastTimePlayOccurs = time(NULL) + 1;
		}

		timerVideoRing->start();
	} else if(recVideoRingTime >= 0) {
		//Receiving
		if (pbVideoRing->isHidden())	// first tick of this call (the counter wraps every 30 s)
			ringClock.start();
		else if (ringClock.isValid() && ringClock.elapsed() > CALL_RING_TIMEOUT_MS + 10000)
		{
			// The caller gives up first; it may be gone without a word.
			addVoipStatusMessage(tr("Missed call from %1.").arg(peerName()));
			hangupCallVideo();
			return;
		}
		++recVideoRingTime;
		if (recVideoRingTime == 100) recVideoRingTime = 0;
		pbVideoRing->setValue(recVideoRingTime);
		pbVideoRing->setToolTip(tr("Your friend is calling you for video. Respond."));
		pbVideoRing->setVisible(true);

		//launch an animation. Don't launch it if already animating
		if (!videoCaptureToggleButton->graphicsEffect()
		    || (videoCaptureToggleButton->graphicsEffect()->inherits("QGraphicsOpacityEffect")
		        && ((QGraphicsOpacityEffect*)videoCaptureToggleButton->graphicsEffect())->opacity() == 1)
		    ) {
			QGraphicsOpacityEffect *effect = new QGraphicsOpacityEffect(audioListenToggleButton);
			videoCaptureToggleButton->setGraphicsEffect(effect);
			QPropertyAnimation *anim = new QPropertyAnimation(effect, "opacity", effect);
			anim->setStartValue(1);
			anim->setKeyValueAt(0.5,0);
			anim->setEndValue(1);
			anim->setDuration(timerVideoRing->interval());
			anim->start();
		}
		videoCaptureToggleButton->setToolTip(tr("Answer"));

		if (time(NULL) > lastTimePlayOccurs) {
			soundManager->play(VOIP_SOUND_INCOMING_VIDEO_CALL);
			lastTimePlayOccurs = time(NULL) + 1;
		}

		if (mVOIPNotify) mVOIPNotify->notifyReceivedVoipVideoCall(mChatWidget->getChatId().toPeerId());

		timerVideoRing->start();
	} else {
		//Nothing to do, reset stat
		pbVideoRing->setHidden(true);
		updateCallButtons();
		pbVideoRing->setValue(0);
		pbVideoRing->setToolTip("");
		videoCaptureToggleButton->setGraphicsEffect(0);
	}
}
