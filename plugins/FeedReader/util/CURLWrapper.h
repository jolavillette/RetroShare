/*******************************************************************************
 * plugins/FeedReader/util/CURLWrapper.h                                       *
 *                                                                             *
 * Copyright (C) 2012 by Thunder <retroshare.project@gmail.com>                *
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

#ifndef CURLWRAPPER
#define CURLWRAPPER

#include <functional>
#include <string>
#include <vector>
#include <curl/curl.h>

class CURLWrapper
{
public:
	/* abortCheck is polled while a transfer runs (libcurl calls the progress
	 * callback about once a second even when no byte arrives); returning true
	 * aborts the transfer with CURLE_ABORTED_BY_CALLBACK. */
	CURLWrapper(const std::string &proxy, std::function<bool()> abortCheck = std::function<bool()>());
	~CURLWrapper();

	/* libcurl keeps a pointer to this object for the progress callback */
	CURLWrapper(const CURLWrapper &) = delete;
	CURLWrapper &operator=(const CURLWrapper &) = delete;

	CURLcode downloadText(const std::string &link, std::string &data);
	CURLcode downloadBinary(const std::string &link, std::vector<unsigned char> &data);

	long responseCode() { return longInfo(CURLINFO_RESPONSE_CODE); }
	std::string contentType() { return stringInfo(CURLINFO_CONTENT_TYPE); }
	std::string effectiveUrl() { return stringInfo(CURLINFO_EFFECTIVE_URL); }

protected:
	long longInfo(CURLINFO info);
	std::string stringInfo(CURLINFO info);

private:
	static int xferInfoCallback(void *clientp, curl_off_t dltotal, curl_off_t dlnow, curl_off_t ultotal, curl_off_t ulnow);

	CURL *mCurl;
	std::function<bool()> mAbortCheck;
};

#endif 
