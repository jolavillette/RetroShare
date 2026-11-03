/*******************************************************************************
 * gui/statusbar/dhtstatus.cpp                                                 *
 *                                                                             *
 * Copyright (c) 2009 Retroshare Team <retroshare.project@gmail.com>           *
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

#include "dhtstatus.h"

#include <QLayout>
#include <QLabel>
#include <QIcon>
#include <QPixmap>

#include "retroshare/rsconfig.h"
#include "retroshare/rsdht.h"
#include "retroshare/rspeers.h"
#include "gui/common/FilesDefs.h"

#include "util/misc.h"

#include <iomanip>

DHTStatus::DHTStatus(QWidget *parent)
 : QWidget(parent)
{
    QHBoxLayout *hbox = new QHBoxLayout();
    hbox->setContentsMargins(0, 0, 0, 0);
    hbox->setSpacing(6);
       
    statusDHT = new QLabel("<strong>" + tr("DHT") + ":</strong>", this );
	 statusDHT->setToolTip(tr("<p>Retroshare uses Bittorrent's DHT as a proxy for connexions. It does not \"store\" your IP in the DHT. \
				                        Instead the DHT is used by your trusted nodes to reach you while processing standard DHT requests. \
												The status bullet will turn green as soon as Retroshare gets a DHT response from one of your trusted nodes.</p>")) ;
    hbox->addWidget(statusDHT);
    
    dhtstatusLabel = new QLabel( this );
    dhtstatusLabel->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":/icons/bullet_grey_128.png"));
    hbox->addWidget(dhtstatusLabel);
    
    spaceLabel = new QLabel( "|", this );
    spaceLabel->setVisible(false);
    hbox->addWidget(spaceLabel);

    dhtnetworkLabel = new QLabel( this );
    dhtnetworkLabel->setVisible(false);
    int S = QFontMetricsF(dhtnetworkLabel->font()).height();
    dhtnetworkLabel->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":/images/dht32.png").scaledToHeight(S,Qt::SmoothTransformation));
    hbox->addWidget(dhtnetworkLabel);

    dhtnetworksizeLabel = new QLabel( "0 (0) ",this );
    dhtnetworksizeLabel->setVisible(false);
    dhtnetworksizeLabel->setTextInteractionFlags(Qt::TextSelectableByMouse);
    hbox->addWidget(dhtnetworksizeLabel);

    hbox->addSpacing(2);

    _compactMode = false;

    setLayout( hbox );
}

/* "203–282": what we saw – what we estimate; one number when they agree */
static QString sizeRange(uint32_t lo, uint32_t hi)
{
	if (lo == 0 && hi == 0)
		return QString(QChar(0x2013));
	if (hi <= lo)
		return misc::userFriendlyUnit(lo, 1);
	return misc::userFriendlyUnit(lo, 1) + QString(QChar(0x2013)) + misc::userFriendlyUnit(hi, 1);
}

void DHTStatus::setSizeEstimate(const QString &prefix)
{
	RsDhtSizeEstimate est;
	if (rsDht)
		rsDht->getDhtSizeEstimate(est);

	dhtnetworksizeLabel->setText(QString("%1 (%2)").arg(sizeRange(est.rsMin, est.rsMax),
	                                                    sizeRange(est.allMin, est.allMax)));

	QString tip = prefix + "<p><b>" + tr("RetroShare nodes in the DHT") + "</b><br>";
	if (est.rsMin == 0)
		tip += tr("None seen yet.");
	else
		tip += tr("%1 seen in our DHT routing table over the last 24 hours; about %2 "
		          "when the nodes we have not crossed yet are estimated (%3 hours "
		          "observed over the last 7 days).")
		       .arg(est.rsMin).arg(est.rsMax).arg(est.rsObservedSecs / 3600);
	tip += "</p><p><b>" + tr("Whole DHT") + "</b><br>";
	if (est.allMin == 0)
		tip += tr("Not measured yet: this needs a few random lookups (about half an hour).");
	else
		tip += tr("%1 nodes answer our random lookups; about %2 when the nodes we "
		          "hear of but that never answer are counted (most DHT nodes do not "
		          "answer an unknown node). From the last %3 lookups.")
		       .arg(misc::userFriendlyUnit(est.allMin, 1), misc::userFriendlyUnit(est.allMax, 1))
		       .arg(est.allSamples);
	tip += "</p>";
	dhtnetworksizeLabel->setToolTip(tip);
}

void DHTStatus::getDHTStatus()
{
	statusDHT->setVisible(!_compactMode);
	QString text = _compactMode?statusDHT->text():"";

	/* now the extra bit .... switch on check boxes */
	RsConfigNetStatus config;
	rsConfig->getConfigNetStatus(config);

    int S = QFontMetricsF(dhtstatusLabel->font()).height();
	if (!(config.DHTActive))
	{
		// GRAY.
        dhtstatusLabel->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":/icons/bullet_grey_128.png").scaledToHeight(S,Qt::SmoothTransformation));
		dhtstatusLabel->setToolTip( text + tr("DHT Off"));

		spaceLabel->setVisible(false);
		dhtnetworkLabel->setVisible(false);
		dhtnetworksizeLabel->setVisible(false);
		
		dhtnetworksizeLabel->setText("");
		dhtnetworksizeLabel->setToolTip("");
	}
	else
	{
		if (config.netDhtOk)
		{
#define MIN_RS_NET_SIZE		10
			// YELLOW or GREEN.
			if (config.netDhtRsNetSize < MIN_RS_NET_SIZE)
			{
                dhtstatusLabel->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":/icons/bullet_yellow_128.png").scaledToHeight(S,Qt::SmoothTransformation));
                dhtstatusLabel->setToolTip( text + tr("DHT Searching for RetroShare Peers"));

				spaceLabel->setVisible(true);
				dhtnetworkLabel->setVisible(true);
				dhtnetworksizeLabel->setVisible(true);

				setSizeEstimate(text);
			}
			else
			{
                dhtstatusLabel->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":/icons/bullet_green_128.png").scaledToHeight(S,Qt::SmoothTransformation));
				dhtstatusLabel->setToolTip( text + tr("DHT Good"));

				spaceLabel->setVisible(true);
				dhtnetworkLabel->setVisible(true);
				dhtnetworksizeLabel->setVisible(true);

				setSizeEstimate(text);
			}
		}
		else
		{
			// RED - some issue.
                dhtstatusLabel->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":/icons/bullet_red_128.png").scaledToHeight(S,Qt::SmoothTransformation));
            dhtstatusLabel->setToolTip( text + tr("No peer found in DHT"));

			spaceLabel->setVisible(false);
			dhtnetworkLabel->setVisible(false);
			dhtnetworksizeLabel->setVisible(false);
			
			dhtnetworksizeLabel->setText("");
			dhtnetworksizeLabel->setToolTip("");
		}
	}
}
