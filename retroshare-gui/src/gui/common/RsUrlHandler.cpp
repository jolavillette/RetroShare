/*******************************************************************************
 * retroshare-gui/src/gui/common/RSUrlHandler.h                                *
 *                                                                             *
 * Copyright (C) 2011 Cyril Soler     <retroshare.project.com>                 *
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

#include <stdexcept>
#include <QDesktopServices>
#include <QUrl>
#include "RsCollectionDialog.h"
#include "RsUrlHandler.h"

#ifdef Q_OS_WIN
#include <string>
#include <thread>
#include <QDir>
#include "util/rsdebug.h"
#include <windows.h>
#include <shellapi.h>
#include <objbase.h>
#endif

bool RsUrlHandler::openUrl(const QUrl& url)
{
	if(url.scheme() == QString("file") && url.toLocalFile().endsWith("."+RsCollection::ExtensionString))
        return RsCollectionDialog::openExistingCollection(url.toLocalFile());

#ifdef Q_OS_WIN
	// Opening a local file goes through ShellExecute, which returns only once
	// Windows is done with it: a UAC or SmartScreen prompt for an executable,
	// an "Open with" chooser, a virus scan. On the GUI thread, that froze the
	// RetroShare window for as long as the prompt stayed up. Launch from a
	// short-lived thread instead. The outcome is not known when this returns,
	// so a failure can only be logged.
	if(url.isLocalFile())
	{
		const std::wstring path = QDir::toNativeSeparators(url.toLocalFile()).toStdWString();

		std::thread([path]()
		{
			// The shell may hand the file over to COM-based extensions.
			const HRESULT com = CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);

			SHELLEXECUTEINFOW info = {};
			info.cbSize = sizeof(info);
			info.fMask = SEE_MASK_NOASYNC;	// this thread ends right after the call
			info.lpFile = path.c_str();
			info.nShow = SW_SHOWNORMAL;

			if(!ShellExecuteExW(&info))
			{
				const DWORD err = GetLastError();
				if(err != ERROR_CANCELLED)	// not when the user declined a prompt
					RsWarn() << "Could not open " << QString::fromStdWString(path).toStdString()
					         << ": Windows error " << err;
			}

			if(SUCCEEDED(com))
				CoUninitialize();
		}).detach();

		return true;
	}
#endif

	return QDesktopServices::openUrl(url) ;
}
