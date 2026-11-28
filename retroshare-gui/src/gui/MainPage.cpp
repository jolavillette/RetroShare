/*******************************************************************************
 * gui/MainPage.cpp                                                            *
 *                                                                             *
 * Copyright (c) 2006 Retroshare Team  <retroshare.project@gmail.com>          *
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

#include <retroshare-gui/mainpage.h>

#include "common/FloatingHelpBrowser.h"
#include "gui/help/browser/helpbrowser.h"
#include "gui/settings/rsharesettings.h"
#include "util/misc.h"

#include <QTextDocument>
#include <QToolButton>
#include <QTimer>

MainPage::MainPage(QWidget *parent , Qt::WindowFlags flags ) : QWidget(parent, flags)
{
	mHelpBrowser = NULL ;
	mIcon = QIcon();
	mName = "";
	mHelp = "";
    mUserNotify = nullptr;
}

UserNotify *MainPage::getUserNotify()
{
    if(!mUserNotify)
        mUserNotify = createUserNotify(this);

    return mUserNotify;
}

void MainPage::registerHelpButton(QToolButton *button, const QString& help_html_text, const QString &code_name)
{
	mHelpCodeName = code_name ;

	// The short help text is shown in a bubble the first time the page is
	// displayed (see showEvent()).
	if (mHelpBrowser == nullptr)
		mHelpBrowser = new FloatingHelpBrowser(this, nullptr) ;

	mHelpBrowser->setHelpText(help_html_text) ;

	// The help button opens the manual on the topic of the page. A page that
	// has no topic in the manual shows its bubble instead.
	int H = misc::getFontSizeFactor("HelpButton").height();
	button->setIconSize(QSize(H, H)) ;//Square Icon
	button->setCheckable(false) ;
	button->setToolTip(tr("Open the manual on this page (F1)")) ;

	button->disconnect(this) ;
	connect(button, &QToolButton::clicked, this, [this]()
	{
		const QString topic = helpTopic() ;

		if (!topic.isEmpty())
			HelpBrowser::showWindow(topic) ;
		else if (!mHelpBrowser->document()->isEmpty())
			mHelpBrowser->showHelp(true) ;
	});
}

void MainPage::setHelpTopic(const QString& topic)
{
	HelpBrowser::setWidgetTopic(this, topic) ;
}

QString MainPage::helpTopic() const
{
	return HelpBrowser::widgetTopic(this) ;
}

void MainPage::showEvent(QShowEvent *s)
{
    if(!Settings->getPageAlreadyDisplayed(mHelpCodeName) && mHelpBrowser!=NULL && !mHelpBrowser->document()->isEmpty())
    {
        // I use a timer to make sure that the GUI is able to process that.
        FloatingHelpBrowser *helpBrowser = mHelpBrowser ;
        QTimer::singleShot(1000, helpBrowser, [helpBrowser]() { helpBrowser->showHelp(true) ; });

        Settings->setPageAlreadyDisplayed(mHelpCodeName,true);
    }

    QWidget::showEvent(s);
}
