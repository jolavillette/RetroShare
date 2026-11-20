/*******************************************************************************
 * gui/statusbar/hashingstatus.cpp                                             *
 *                                                                             *
 * Copyright (c) 2008 Retroshare Team <retroshare.project@gmail.com>           *
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

#include <algorithm>

#include <QFileInfo>
#include <QLayout>
#include <QLabel>
#include <QMovie>
#include <QToolButton>

#include "retroshare/rsfiles.h"
#include "hashingstatus.h"
#include "gui/common/ElidedLabel.h"
#include "util/qtthreadsutils.h"
#include "util/misc.h"
#include "gui/common/FilesDefs.h"

//#define DEBUG_HASHING_STATUS 1

HashingStatus::HashingStatus(QWidget *parent)
 : QWidget(parent)
{
    QHBoxLayout *hbox = new QHBoxLayout(this);
    hbox->setContentsMargins(0, 0, 0, 0);
    hbox->setSpacing(6);
        
    movie = new QMovie(":/images/loader/indicator-16.gif");
    movie->setSpeed(80); // 2x speed
    hashloader = new QLabel(this);
    hashloader->setMovie(movie);
    hbox->addWidget(hashloader);

    movie->jumpToNextFrame(); // to calculate the real width
    statusHashing = new ElidedLabel(this);
    hbox->addWidget(statusHashing);

    _compactMode = false;
    mPaused = false;

    setLayout(hbox);

    hashloader->hide();
    statusHashing->hide();

	mEventHandlerId=0;
    rsEvents->registerEventsHandler( [this](std::shared_ptr<const RsEvent> event) { RsQThreadUtils::postToObject( [this,event]() { handleEvent_main_thread(event); }) ;}, mEventHandlerId, RsEventType::SHARED_DIRECTORIES );
}

void HashingStatus::handleEvent_main_thread(std::shared_ptr<const RsEvent> event)
{
    if(event->mType != RsEventType::SHARED_DIRECTORIES)
        return;

    const RsSharedDirectoriesEvent *fe = dynamic_cast<const RsSharedDirectoriesEvent*>(event.get());

	if(!fe)
        return;

 	switch (fe->mEventCode)
    {
    default:
        break;

    case RsSharedDirectoriesEventCode::HASHING_PROCESS_RESUMED:
        mPaused = false;
        statusHashing->setText(mLastText);						// fallthrough
#ifdef DEBUG_HASHING_STATUS
        std::cerr << "HashStatusEventHandler: received event " << (void*)fe->mEventCode << ": " ;
        std::cerr << "HASHING RESUMED" << std::endl;
#endif

    case RsSharedDirectoriesEventCode::HASHING_PROCESS_STARTED:
    {
#ifdef DEBUG_HASHING_STATUS
        std::cerr << "HashStatusEventHandler: received event " << (void*)fe->mEventCode << ": " ;
        std::cerr << "HASHING STARTED" << std::endl;
#endif
        hashloader->show() ;
        hashloader->setMovie(movie) ;
        movie->start() ;

        // A run starts with the first file queued by the scan of the shared directories, and nothing is hashed
        // until that scan is over. Not on a resume, which comes through here with its text already set.
        if(statusHashing->text().isEmpty())
            statusHashing->setText(tr("Scanning shared directories..."));

        statusHashing->setVisible(!_compactMode) ;
    }
    break;

    case RsSharedDirectoriesEventCode::HASHING_PROCESS_PAUSED:
    {
#ifdef DEBUG_HASHING_STATUS
        std::cerr << "HashStatusEventHandler: received event " << (void*)fe->mEventCode << ": " ;
        std::cerr << "HASHING PAUSED" << std::endl;
#endif
        movie->stop() ;
        hashloader->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":/images/resume.png")) ;

        mPaused = true;
        mLastText = statusHashing->text();
        statusHashing->setText(QObject::tr("[Hashing is paused]"));
        setToolTip(QObject::tr("Click to resume the hashing process"));
    }
    break;

    case RsSharedDirectoriesEventCode::HASHING_PROCESS_FINISHED:
    {
#ifdef DEBUG_HASHING_STATUS
        std::cerr << "HashStatusEventHandler: received event " << (void*)fe->mEventCode << ": " ;
        std::cerr << "HASHING FINISHED" << std::endl;
#endif
        movie->stop() ;
        mPaused = false;
        statusHashing->setText(QString());
        statusHashing->hide() ;
        setToolTip(QString());
        hashloader->hide() ;
    }
    break;

    case RsSharedDirectoriesEventCode::HASHING_FILE:
    {
        // All sizes in bytes. The totals cover the whole hashing run, the files already done included.

        const uint64_t total_size = fe->mTotalSizeToHash;
        const uint64_t done_size = std::min(fe->mTotalHashedSize, total_size);
        const uint64_t file_number = fe->mHashCounter+1;

        const QString path = QString::fromUtf8(fe->mFilePath.c_str());
        const QString separator = QString(" ") + QChar(0x00B7) + " ";
        QStringList parts;

        // Last event of a run: every file is done, and the hashing thread waits a few seconds for more before it
        // stops. Showing the last file at 100 % with the animation would read as a file still being hashed.
        if(fe->mHashCounter >= fe->mTotalFilesToHash)
        {
            parts << (fe->mTotalFilesToHash == 1 ? tr("1 file") : tr("%1 files").arg(fe->mTotalFilesToHash));
            parts << misc::friendlyUnit(total_size);

            if(fe->mTotalIsProvisional)
                parts << tr("scanning...");

            const QString text = tr("Hashing done") + separator + parts.join(separator);

            if(mPaused)
            {
                mLastText = text;
                break;
            }

            if(!fe->mTotalIsProvisional)
            {
                movie->stop() ;
                hashloader->hide() ;
            }
            statusHashing->setText(text);
            setToolTip(text);
            break;
        }

        if(fe->mTotalIsProvisional)
        {
            // The shared directories are still being scanned, so the totals still grow: a ratio to them
            // would go backwards. Only what is done is shown.

            parts << QString::number(file_number);
            parts << tr("%1 done").arg(misc::friendlyUnit(fe->mTotalHashedSize));
        }
        else
        {
            parts << QString::number(std::min(file_number, fe->mTotalFilesToHash)) + "/" + QString::number(fe->mTotalFilesToHash);
            parts << misc::friendlyUnit(done_size) + " / " + misc::friendlyUnit(total_size);

            if(total_size > 0)
                parts << QString::number(int(done_size * 100.0 / total_size)) + " %";
        }

        if(fe->mHashingSpeedBytes > 0)
        {
            parts << misc::friendlyUnit(fe->mHashingSpeedBytes) + "/s";

            if(!fe->mTotalIsProvisional && total_size > done_size)
                parts << tr("%1 left").arg(misc::userFriendlyDuration((total_size - done_size) / fe->mHashingSpeedBytes));
        }

        if(fe->mTotalIsProvisional)
            parts << tr("scanning...");

        const QString progress = parts.join(separator);

        // The label cuts what does not fit at its right end. The name comes first so that it is the last thing to
        // go, and a long one is shortened in its middle so that it does not push the figures out.
        QString name = QFileInfo(path).fileName();

        if(name.length() > 40)
            name = name.left(20) + QChar(0x2026) + name.right(19);

        const QString msg = name + separator + progress;
        const QString text = tr("Hashing file") + " " + msg;

        // The file being hashed when the pause is asked is hashed to its end, and keeps sending its progress.
        if(mPaused)
        {
            mLastText = text;
            break;
        }

        // The animation was stopped if a previous batch of this run was shown as done.
        if(movie->state() != QMovie::Running)
        {
            hashloader->show() ;
            hashloader->setMovie(movie) ;
            movie->start() ;
        }

        statusHashing->setText(text);
        setToolTip(progress + "\n" + path + "\n" + QObject::tr("Click to pause the hashing process"));
#ifdef DEBUG_HASHING_STATUS
        std::cerr << "HashStatusEventHandler: received event " << (void*)fe->mEventCode << ": " ;
        std::cerr << "HASHING FILE " << msg.toStdString() << std::endl;
#endif
    }
    break;

	case RsSharedDirectoriesEventCode::SAVING_FILE_INDEX:
    {
#ifdef DEBUG_HASHING_STATUS
        std::cerr << "HashStatusEventHandler: received event " << (void*)fe->mEventCode << ": " ;
        std::cerr << "SAVING FILE INDEX" << std::endl;
#endif
        statusHashing->setText(tr("Saving file index..."));
    }
    break;

    };
}

HashingStatus::~HashingStatus()
{
    rsEvents->unregisterEventsHandler(mEventHandlerId);
    delete(movie);
}

void HashingStatus::mousePressEvent(QMouseEvent *)
{
	rsFiles->togglePauseHashingProcess() ;
}

