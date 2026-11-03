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

/* "180–260", one number when both ends agree */
static QString sizeRange(uint32_t lo, uint32_t hi)
{
	/* userFriendlyUnit() leaves "199 " when there is no unit */
	if (hi <= lo)
		return misc::userFriendlyUnit(lo, 1).trimmed();
	return misc::userFriendlyUnit(lo, 1).trimmed() + QChar(0x2013) + misc::userFriendlyUnit(hi, 1).trimmed();
}

/* label value: the range, flagged when it rests on too few lookups yet */
static QString sizeValue(bool available, uint32_t lo, uint32_t hi, bool provisional)
{
	if (!available)
		return DHTStatus::tr("n/a");
	if (hi == 0)
		return DHTStatus::tr("measuring");
	QString v = sizeRange(lo, hi);
	if (provisional)
		v += " " + DHTStatus::tr("provisional");
	return v;
}

static QString sizeDetails(uint32_t lo, uint32_t hi, bool provisional, uint32_t lookups,
                           uint32_t reach, uint32_t reachErr, uint32_t heard, uint32_t heardErr)
{
	if (hi == 0)
		return DHTStatus::tr("Not measured yet: the first random lookups take a few minutes.");
	QString t = DHTStatus::tr("Range %1, from %2 random lookups%3.")
	        .arg(sizeRange(lo, hi)).arg(lookups)
	        .arg(provisional ? DHTStatus::tr(" (provisional: fewer than 12 lookups, it will still move)") : QString());
	t += "<br>" + DHTStatus::tr("Nodes that answer us: %1 (±%2%). With the silent ones at least two nodes report: %3 (±%4%).")
	        .arg(misc::userFriendlyUnit(reach, 1)).arg(reachErr)
	        .arg(misc::userFriendlyUnit(heard, 1)).arg(heardErr);
	return t;
}

void DHTStatus::setSizeEstimate(const QString &prefix)
{
	RsDhtSizeEstimate est;
	if (rsDht)
		rsDht->getDhtSizeEstimate(est);

	dhtnetworksizeLabel->setText(tr("RS: %1 · all: %2")
	        .arg(sizeValue(est.lookupsAvailable, est.rsMin, est.rsMax, est.rsProvisional),
	             sizeValue(est.lookupsAvailable, est.allMin, est.allMax, est.allProvisional)));

	QString tip = prefix + "<p><b>" + tr("RetroShare nodes online now") + "</b><br>";
	if (!est.lookupsAvailable)
		tip += tr("This build cannot measure it (its BitDHT lacks network size lookups).");
	else
		tip += sizeDetails(est.rsMin, est.rsMax, est.rsProvisional, est.rsLookups,
		                   est.rsReachable, est.rsReachableErr, est.rsHeard, est.rsHeardErr);
	tip += "<br>" + tr("Checks: %1 seen in our DHT routing table over the last 24 h")
	        .arg(est.rsSeen24h);
	if (est.rsDistinct24h)
		tip += tr("; about %1 distinct ones online during the last 24 h (routing table crossed with the lookups)")
		        .arg(est.rsDistinct24h);
	if (est.rsFriendCalibrated)
		tip += tr("; about %1 online now, calibrated on our %2 friends present in the DHT (a low estimate: our table favours our friends)")
		        .arg(est.rsFriendCalibrated).arg(est.rsFriendsInDht);
	tip += ".</p><p><b>" + tr("Whole DHT") + "</b><br>";
	if (!est.lookupsAvailable)
		tip += tr("This build cannot measure it (its BitDHT lacks network size lookups).");
	else
		tip += sizeDetails(est.allMin, est.allMax, est.allProvisional, est.allLookups,
		                   est.allReachable, est.allReachableErr, est.allHeard, est.allHeardErr)
		       + "<br>" + tr("Most DHT nodes do not answer a node they do not know: the truth lies towards the high end.");
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
