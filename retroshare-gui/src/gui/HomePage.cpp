/*******************************************************************************
 * gui/HomePage.cpp                                                            *
 *                                                                             *
 * Copyright (C) 2016 Defnax          <retroshare.project@gmail.com>           *
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

#include "HomePage.h"
#include "ui_HomePage.h"

#include "retroshare/rsinit.h"

#include "util/qtthreadsutils.h"
#include "util/misc.h"

#include "gui/common/QrCodeDialog.h"
#include "gui/common/FilesDefs.h"
#include "gui/RetroShareLink.h"
#include "gui/msgs/MessageComposer.h"
#include "gui/connect/ConnectFriendWizard.h"
#include "gui/connect/ConfCertDialog.h"
#include <gui/QuickStartWizard.h>
#include "settings/rsharesettings.h"

#if QT_VERSION >= QT_VERSION_CHECK(5,0,0)
#include <QUrlQuery>
#endif

#include <iostream>
#include <string>
#include <QTime>
#include <QMenu>
#include <QPushButton>
#include <QDesktopServices>
#include <QMessageBox>
#include <QFileDialog>
#include <QImageReader>
#include <QMimeData>
#include <QClipboard>
#include <QTextStream>
#include <QTextCodec>

HomePage::HomePage(QWidget *parent) :
    MainPage(parent),
    ui(new Ui::HomePage)
{
    ui->setupUi(this);

    connect(ui->pasteButton, SIGNAL(clicked()), this, SLOT(pasteFriendId()));
    connect(ui->openFileButton, SIGNAL(clicked()), this, SLOT(openFriendFile()));
    connect(ui->copyIDButton, SIGNAL(clicked()), this, SLOT(copyId()));
    connect(ui->saveButton, SIGNAL(clicked()), this, SLOT(saveCert()));
    connect(ui->emailButton, SIGNAL(clicked()), this, SLOT(runEmailClient()));

    // The two format buttons exclude each other: one signal is enough.
    connect(ui->fullCertButton, SIGNAL(toggled(bool)), this, SLOT(updateOwnCert()));

    connect(ui->qrCodeButton, &QPushButton::clicked, this, [this]()
    {
        // The code carries what the page shows: the short ID, or the whole
        // certificate when the full format is chosen.
        QString invite, description;
        getOwnCert(invite, description);
        QrCodeDialog::show(this, tr("Retroshare ID"), invite);
    });

    // A hidden node has no address to choose: the options button stays away.
    ui->optionsButton->setVisible(!RsAccounts::isHiddenNode());

    QMenu *menu = new QMenu(this);

    if(!RsAccounts::isHiddenNode())
    {
        mIncludeLocIPact = new QAction(QIcon(), tr("Include current local IP"),this);
        connect(mIncludeLocIPact, SIGNAL(triggered()), this, SLOT(updateOwnCert()));
        mIncludeLocIPact->setCheckable(true);
        mIncludeLocIPact->setChecked(true);
        menu->addAction(mIncludeLocIPact);

        mIncludeExtIPact = new QAction(QIcon(), tr("Include current external IP"),this);
        connect(mIncludeExtIPact, SIGNAL(triggered()), this, SLOT(updateOwnCert()));
        mIncludeExtIPact->setCheckable(true);
        mIncludeExtIPact->setChecked(true);
        menu->addAction(mIncludeExtIPact);

        mIncludeDNSact = new QAction(QIcon(), tr("Include my DNS"),this);
        connect(mIncludeDNSact, SIGNAL(triggered()), this, SLOT(updateOwnCert()));
        mIncludeDNSact->setCheckable(true);
        mIncludeDNSact->setChecked(true);
        menu->addAction(mIncludeDNSact);

        mIncludeIPHistoryact = new QAction(QIcon(), tr("Include all IPs history"),this);
        connect(mIncludeIPHistoryact, SIGNAL(triggered()), this, SLOT(updateOwnCert()));
        mIncludeIPHistoryact->setCheckable(true);
        mIncludeIPHistoryact->setChecked(false);
        menu->addAction(mIncludeIPHistoryact);
    }

    ui->optionsButton->setMenu(menu);

    connect(ui->openwebhelp,SIGNAL(clicked()), this,SLOT(openWebHelp())) ;

	int H = misc::getFontSizeFactor("HelpButton").height();
	QString help_str = tr(
	    "<h1><img width=\"%1\" src=\":/icons/help_64.png\">&nbsp;&nbsp;Welcome to Retroshare!</h1>"
	    "<p>You need to <b>make friends</b>! After you create a network of friends or join an existing network,"
	    "   you'll be able to exchange files, chat, talk in forums, etc. </p>"
	    "<div align=\"center\"><IMG width=\"%2\" height=\"%3\" src=\":/images/network_map.png\" style=\"display: block; margin-left: auto; margin-right: auto; \"/></div>"
	    "<p>To do so, copy your Retroshare ID on this page and send it to friends, and add your friends' Retroshare ID.</p>"
	    "<p>Another option is to search the internet for \"Retroshare chat servers\" (independently administrated). These servers allow you to exchange"
	    "   Retroshare ID with a dedicated Retroshare node, through which"
	    "   you will be able to anonymously meet other people.</p>"
	                     ).arg(QString::number(2*H), QString::number(width()*0.5), QString::number(width()*0.5*(337.0/800.0)));//<img> needs height and width defined.
	registerHelpButton(ui->helpButton,help_str,"HomePage") ;

	// register a event handler to catch IP updates

    mEventHandlerId = 0;
    rsEvents->registerEventsHandler( [this](std::shared_ptr<const RsEvent> event) { handleEvent(event); }, mEventHandlerId, RsEventType::NETWORK );

    updateOwnCert();

    updateHomeLogo();
}

void HomePage::handleEvent(std::shared_ptr<const RsEvent> e)
{
    if(e->mType != RsEventType::NETWORK)
        return;

    const RsNetworkEvent *ne = dynamic_cast<const RsNetworkEvent*>(e.get());

    if(!ne)
        return;

    // in any case we update the IPs

    switch(ne->mNetworkEventCode)
    {
    case RsNetworkEventCode::LOCAL_IP_UPDATED:  // [fallthrough]
    case RsNetworkEventCode::EXTERNAL_IP_UPDATED:  // [fallthrough]
    case RsNetworkEventCode::DNS_UPDATED:  // [fallthrough]
                RsQThreadUtils::postToObject( [=]()
                {
                    updateOwnCert();
                },this);
        break;
    default:
        break;
    }
}

#ifdef DEAD_CODE
void HomePage::certContextMenu(QPoint /*point*/)
{
    QMenu menu(this) ;

    QAction *CopyAction = new QAction(QIcon(),tr("Copy your Cert to Clipboard"), this);
    connect(CopyAction, SIGNAL(triggered()), this, SLOT(copyCert()));

    QAction *SaveAction = new QAction(QIcon(),tr("Save your Cert into a File"), this);
    connect(SaveAction, SIGNAL(triggered()), this, SLOT(saveCert()));

    menu.addAction(CopyAction);
    menu.addAction(SaveAction);

    QAction *shortFormatAct = new QAction(QIcon(), tr("Use new (short) certificate format"),this);
    connect(shortFormatAct, SIGNAL(triggered()), this, SLOT(toggleUseShortFormat()));
    shortFormatAct->setCheckable(true);
    shortFormatAct->setChecked(mUseShortFormat);

    menu.addAction(shortFormatAct);

    if(!RsAccounts::isHiddenNode())
    {
        QAction *includeIPsAct = new QAction(QIcon(), tr("Include all your known IPs"),this);
        connect(includeIPsAct, SIGNAL(triggered()), this, SLOT(toggleIncludeAllIPs()));
        includeIPsAct->setCheckable(true);
        includeIPsAct->setChecked(mIncludeAllIPs);

        menu.addAction(includeIPsAct);
    }

    menu.exec(QCursor::pos());
}
#endif

HomePage::~HomePage()
{
    rsEvents->unregisterEventsHandler(mEventHandlerId);
    delete ui;
}

RetroshareInviteFlags HomePage::currentInviteFlags() const
{
    RetroshareInviteFlags invite_flags = RetroshareInviteFlags::NOTHING;

    if(!RsAccounts::isHiddenNode())
    {
        if(mIncludeLocIPact->isChecked())
            invite_flags |= RetroshareInviteFlags::CURRENT_LOCAL_IP;

        if(mIncludeExtIPact->isChecked())
            invite_flags |= RetroshareInviteFlags::CURRENT_EXTERNAL_IP;

        if(mIncludeDNSact->isChecked())
            invite_flags |= RetroshareInviteFlags::DNS;

        if(mIncludeIPHistoryact->isChecked())
            invite_flags |= RetroshareInviteFlags::FULL_IP_HISTORY;
    }

    return invite_flags;
}
void HomePage::getOwnCert(QString& invite,QString& description) const
{
    RsPeerDetails detail;

    if (!rsPeers->getPeerDetails(rsPeers->getOwnId(), detail))
    {
        std::cerr << "(EE) Cannot retrieve information about own certificate. That is a real problem!!" << std::endl;
        return ;
    }
    auto invite_flags = currentInviteFlags();

    invite_flags |= RetroshareInviteFlags::SLICE_TO_80_CHARS;

    if(!ui->fullCertButton->isChecked())
    {
        std::string short_invite;
        rsPeers->getShortInvite(short_invite,rsPeers->getOwnId(),invite_flags | RetroshareInviteFlags::RADIX_FORMAT);
        invite = QString::fromStdString(short_invite);
    }
    else
        invite = QString::fromStdString(rsPeers->GetRetroshareInvite(detail.id,invite_flags));

    description = ConfCertDialog::getCertificateDescription(detail,false,!ui->fullCertButton->isChecked(),invite_flags);
}

void HomePage::updateOwnCert()
{
    QString certificate, description;

    getOwnCert(certificate,description);

    if(!ui->fullCertButton->isChecked())	// in this case we have to split the cert for a better display
    {
        QString S;
        QString txt;

        for(int i=0;i<certificate.size();)
            if(S.length() < 64)
                S += certificate[i++];
            else
            {
                txt += S + "\n";
                S.clear();
            }

        txt += S;
        certificate = txt;	// the "\n" is here to make some space
    }

    ui->retroshareid->setText("\n"+certificate+"\n");
    ui->retroshareid->setToolTip(description);
}

static void sendMail(const QString &sAddress, const QString &sSubject, const QString &sBody)
{
    QUrl url = QUrl("mailto:" + sAddress);

#if QT_VERSION >= QT_VERSION_CHECK(5,0,0)
    QUrlQuery urlQuery;
#else
    QUrl &urlQuery(url);
#endif

    urlQuery.addQueryItem("subject", sSubject);
#ifdef Q_OS_WIN
    /* search and replace the end of lines with: "%0D%0A" */
    urlQuery.addQueryItem("body", QString(sBody).replace("\n", "%0D%0A"));
#else
    urlQuery.addQueryItem("body", sBody);
#endif

#if QT_VERSION >= QT_VERSION_CHECK(5,0,0)
    url.setQuery(urlQuery);
#endif

    std::cerr << "MAIL STRING:" << (std::string)url.toEncoded().constData() << std::endl;

    /* pass the url directly to QDesktopServices::openUrl */
    QDesktopServices::openUrl (url);
}

void HomePage::runEmailClient()
{
    sendMail("", tr("RetroShare Invite"), ui->retroshareid->text());
}

void HomePage::copyCert()
{
    QClipboard *clipboard = QApplication::clipboard();
    clipboard->setText(ui->retroshareid->text());
    QMessageBox::information(this, "RetroShare", tr("Your Retroshare certificate is copied to Clipboard, paste and send it to your friend via email or some other way"));
}

void HomePage::copyId()
{
    QClipboard *clipboard = QApplication::clipboard();
    clipboard->setText(ui->retroshareid->text());
    QMessageBox::information(this, "RetroShare", tr("Your Retroshare ID is copied to Clipboard, paste and send it to your friend via email or some other way"));
}

void HomePage::saveCert()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save as..."), "", tr("RetroShare Certificate (*.rsc );;All Files (*)"));
    if (fileName.isEmpty())
        return;

    QFile file(fileName);
    if (!file.open(QFile::WriteOnly))
        return;

    //Todo: move save to file to p3Peers::SaveCertificateToFile

    QTextStream ts(&file);
#if QT_VERSION >= QT_VERSION_CHECK (6, 0, 0)
    ts.setEncoding(QStringConverter::Utf8);
#else
    ts.setCodec(QTextCodec::codecForName("UTF-8"));
#endif
    ts << ui->retroshareid->text();
}

void HomePage::pasteFriendId()
{
    const QClipboard *clipboard = QApplication::clipboard();
    const QMimeData *data = clipboard->mimeData();

    if(data && data->hasImage())
    {
        QString text;
        if(QrCodeDialog::decode(clipboard->image(), text))
            addFriendFromText(text);
        else
            QMessageBox::information(this, tr("Add a friend"), tr("The picture in the clipboard holds no readable QR code."));
    }
    else if(!clipboard->text().trimmed().isEmpty())
        addFriendFromText(clipboard->text());
    else
        QMessageBox::information(this, tr("Add a friend"), tr("The clipboard holds no text and no picture. Copy the Retroshare ID of your friend first."));
}

void HomePage::openFriendFile()
{
    QString fileName;
    if(!misc::getOpenFileName(this, RshareSettings::LASTDIR_CERT, tr("Open a Retroshare ID"),
                              tr("Retroshare ID (*.rsc *.txt *.png *.jpg *.jpeg *.bmp *.gif);;All Files (*)"), fileName)
       || fileName.isEmpty())
        return;

    // The content decides, not the name: a picture is searched for a QR code,
    // anything else is read as the text of an ID or certificate.
    QImageReader reader(fileName);
    reader.setDecideFormatFromContent(true);

    if(reader.canRead())
    {
        QString text;
        if(QrCodeDialog::decode(reader.read(), text))
            addFriendFromText(text);
        else
            QMessageBox::information(this, tr("Add a friend"), tr("This picture holds no readable QR code."));
        return;
    }

    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly))
    {
        QMessageBox::warning(this, tr("Add a friend"), tr("Cannot read this file."));
        return;
    }

    // Far above the size of any certificate.
    addFriendFromText(QString::fromUtf8(file.read(1024*1024)));
}

void HomePage::addFriendFromText(const QString& text)
{
    QString certificate = text.trimmed();

    RetroShareLink link(certificate);
    if(link.valid() && link.type() == RetroShareLink::TYPE_CERTIFICATE)
        certificate = link.radix();

    // The wizard opens on the details of the friend, or on what is wrong
    // with the text.
    ConnectFriendWizard connwiz (this);

    connwiz.setCertificate(certificate, false);
    connwiz.exec ();
}

// void HomePage::loadCert()
// {
//     ConnectFriendWizard connwiz (this);
//
//     connwiz.setStartId(ConnectFriendWizard::Page_Cert);
//     connwiz.exec ();
// }

void HomePage::openWebHelp()
{
    QDesktopServices::openUrl(QUrl(QString("https://retrosharedocs.readthedocs.io/en/latest/")));
}

void HomePage::showEvent(QShowEvent *event)
{
	if (!event->spontaneous()) {
		updateHomeLogo();
	}
}

void HomePage::updateHomeLogo()
{
	if (Settings->getSheetName() == ":Standard_Dark")
		ui->label->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":images/logo/logo_web_nobackground_black.png"));
	else
		ui->label->setPixmap(FilesDefs::getPixmapFromQtResourcePath(":images/logo/logo_web_nobackground.png"));
}
