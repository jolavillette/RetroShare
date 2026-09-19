/*******************************************************************************
 * util/retroshare-service_path.cpp                                            *
 *                                                                             *
 * Copyright (C) 2026  Retroshare Team       <retroshare.project@gmail.com>    *
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

#ifdef WINDOWS_SYS
#include <windows.h>
#include <vector>
#include "util/rsdebug.h"
#endif

#include "retroshare-service_path.h"

/** Get file path of the main executable */
std::string getMainExecutablePath(const char *argv0)
{
	std::string mainExecutablePath;

#ifdef WINDOWS_SYS
	// Use a dynamic buffer to handle paths longer than MAX_PATH safely
	std::vector<wchar_t> filename(MAX_PATH);
	while (true) {
		// Clear error state before the call to avoid stale error codes causing infinite loops
		SetLastError(ERROR_SUCCESS);

		// Passing NULL retrieves the path of the current executable (.exe)
		DWORD result = GetModuleFileNameW(NULL, filename.data(), static_cast<DWORD>(filename.size()));

		// Check for function failure
		if (result == 0) {
			RsFatal() << "Failed to calculate string size. Error code: " << GetLastError() << std::endl;
			return std::string();
		}

		// If the buffer was too small or truncated, double its size and try again
		if (result == filename.size() || GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
			filename.resize(filename.size() * 2);
			continue;
		}

		// Success: path was fully copied
		filename.resize(result); // Resize vector to the actual path length
		break;
	}

	// Calculate the required buffer size for the UTF-8 string
	int sizeNeeded = WideCharToMultiByte(CP_UTF8, 0, filename.data(), -1, NULL, 0, NULL, NULL);
	if (sizeNeeded == 0) {
		RsFatal() << "Failed to calculate string size. Error code: " << GetLastError() << std::endl;
		return std::string();
	}

	// Instantiate a std::string with the matching size
	mainExecutablePath.resize(sizeNeeded - 1, 0); // -1 excludes the null-terminator

	// Convert UTF-16 to UTF-8
	int result = WideCharToMultiByte(CP_UTF8, 0, filename.data(), -1, &mainExecutablePath[0], sizeNeeded, NULL, NULL);
	if (result == 0) {
		RsFatal() << "Failed to convert string to UTF-8. Error code: " << GetLastError() << std::endl;
		return std::string();
	}

#else
	if (argv0 != nullptr) {
		mainExecutablePath = std::string(argv0);
	}
#endif

	return mainExecutablePath;
}
