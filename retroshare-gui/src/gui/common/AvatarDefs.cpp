/*******************************************************************************
 * retroshare-gui/src/gui/common/AvatarDefs.cpp                                *
 *                                                                             *
 * Copyright (C) 2012, Robert Fernie <retroshare.project@gmail.com>            *
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

#include <QPixmap>
#include <QMutex>
#include <QMutexLocker>

#include <ctime>
#include <map>

#include <retroshare/rschats.h>
#include <retroshare/rsevents.h>
#include <retroshare/rspeers.h>
#include <retroshare/rsidentity.h>
#include <gui/gxs/GxsIdDetails.h>

#include "AvatarDefs.h"
#include "gui/common/FilesDefs.h"

void AvatarDefs::getOwnAvatar(QPixmap &avatar, const QString& defaultImage)
{
	unsigned char *data = NULL;
	int size = 0;

	/* get avatar */
    rsChats->getOwnNodeAvatarData(data, size);

	if (size == 0) {
        avatar = FilesDefs::getPixmapFromQtResourcePath(defaultImage);
		return;
	}

	/* load image */
	GxsIdDetails::loadPixmapFromData(data, size, avatar,GxsIdDetails::ORIGINAL) ;

	free(data);
}
// The avatar of a node is asked by the item views for every row they lay out or paint (the friend list asks it
// from its decoration role, which the size hint of each row calls). Each request copies the picture out of the
// chat service and hashes it to find the decoded image: fine once, seconds of cpu for a whole list. So the
// result is kept per node, and dropped when the chat service says the avatar of that node changed. "No avatar"
// is only kept for a short time, because asking is also what makes the chat service request it from the peer.

namespace
{
    struct CachedNodeAvatar
    {
        bool    stale ;		// set by the events thread; the entry is refreshed by the next request
        bool    has_avatar ;
        QPixmap pixmap ;
        time_t  ts ;
    };

    const time_t NO_AVATAR_RETRY_DELAY = 60 ; // seconds

    QMutex                             node_avatars_mutex ;
    std::map<RsPeerId,CachedNodeAvatar> node_avatars ;
    bool                               node_avatars_handler_registered = false ;

    void registerNodeAvatarsHandler()
    {
        if(node_avatars_handler_registered || !rsEvents)
            return ;

        node_avatars_handler_registered = true ;

        RsEventsHandlerId_t handler_id = 0 ;

        // Runs in the events thread: it only marks the entry, it never touches the pixmap.

        rsEvents->registerEventsHandler( [](std::shared_ptr<const RsEvent> event)
        {
            const RsFriendListEvent *e = dynamic_cast<const RsFriendListEvent*>(event.get());

            if(!e || (e->mEventCode != RsFriendListEventCode::NODE_AVATAR_CHANGED && e->mEventCode != RsFriendListEventCode::OWN_AVATAR_CHANGED))
                return ;

            QMutexLocker lock(&node_avatars_mutex) ;

            auto it = node_avatars.find(e->mSslId) ;

            if(it != node_avatars.end())
                it->second.stale = true ;
        }, handler_id, RsEventType::FRIEND_LIST );
    }
}

bool AvatarDefs::getAvatarFromSslId(const RsPeerId& sslId, QPixmap &avatar, const QString& defaultImage)
{
    registerNodeAvatarsHandler() ;

    const time_t now = time(NULL) ;
    bool known = false ;
    bool has_avatar = false ;

    {
        QMutexLocker lock(&node_avatars_mutex) ;

        auto it = node_avatars.find(sslId) ;

        if(it != node_avatars.end() && !it->second.stale && (it->second.has_avatar || it->second.ts + NO_AVATAR_RETRY_DELAY > now))
        {
            known = true ;
            has_avatar = it->second.has_avatar ;

            if(has_avatar)
                avatar = it->second.pixmap ;
        }
    }

    if(!known)
    {
        unsigned char *data = NULL;
        int size = 0;

        /* get avatar */
        rsChats->getAvatarData(RsPeerId(sslId), data, size);

        QPixmap pixmap ;

        if(size != 0)
        {
            /* load image */
            has_avatar = GxsIdDetails::loadPixmapFromData(data, size, pixmap, GxsIdDetails::LARGE) ;
            free(data);
        }

        if(has_avatar)
            avatar = pixmap ;

        QMutexLocker lock(&node_avatars_mutex) ;

        CachedNodeAvatar& c(node_avatars[sslId]) ;

        c.stale = false ;
        c.has_avatar = has_avatar ;
        c.pixmap = pixmap ;
        c.ts = now ;
    }

    if(!has_avatar)
    {
        if (!defaultImage.isEmpty()) {
            avatar = GxsIdDetails::makeDefaultGroupIconFromString(QString::fromStdString(sslId.toStdString()), ":icons/person.png", GxsIdDetails::LARGE);
        }
        return false;
    }

    return true;
}
bool AvatarDefs::getAvatarFromGxsId(const RsGxsId& gxsId, QPixmap &avatar, const QString& defaultImage)
{
    //int size = 0;

    /* get avatar */
    RsIdentityDetails details ;

    if(!rsIdentity->getIdDetails(gxsId, details))
    {
        avatar = FilesDefs::getPixmapFromQtResourcePath(defaultImage);
        return false;
    }

    /* load image */

        if(details.mAvatar.mSize == 0 || !GxsIdDetails::loadPixmapFromData(details.mAvatar.mData, details.mAvatar.mSize, avatar,GxsIdDetails::LARGE))
            avatar = GxsIdDetails::makeDefaultIcon(gxsId,GxsIdDetails::LARGE);

        return true;
}

bool AvatarDefs::getAvatarFromGpgId(const RsPgpId& gpgId, QPixmap &avatar, const QString& defaultImage)
{
	unsigned char *data = NULL;
	int size = 0;

    if (gpgId == rsPeers->getGPGOwnId()) {
		/* Its me */
        rsChats->getOwnNodeAvatarData(data,size);
	} else {
		/* get the first available avatar of one of the ssl ids */
        std::list<RsPeerId> sslIds;
        if (rsPeers->getAssociatedSSLIds(gpgId, sslIds)) {
            std::list<RsPeerId>::iterator sslId;
			for (sslId = sslIds.begin(); sslId != sslIds.end(); ++sslId) {
                rsChats->getAvatarData(*sslId, data, size);
				if (size) {
					break;
				}
			}
		}
	}

    if (size == 0) {
        if (!defaultImage.isEmpty()) {
            avatar = GxsIdDetails::makeDefaultGroupIconFromString(QString::fromStdString(gpgId.toStdString()), ":icons/person.png", GxsIdDetails::LARGE);
        }
        return false;
    }

	/* load image */
	GxsIdDetails::loadPixmapFromData(data, size, avatar, GxsIdDetails::LARGE);

	free(data);

    return true;
}
