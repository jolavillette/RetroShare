/*******************************************************************************
 * retroshare-gui/src/gui/common/FloatingHelpBrowser.cpp                       *
 *                                                                             *
 * Copyright (C) 2013, Retroshare Team <retroshare.project@gmail.com>          *
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

#include <QGraphicsDropShadowEffect>
#include <QAbstractButton>
#include <QShowEvent>
#include <QMouseEvent>

#include <iostream>

#include "FloatingHelpBrowser.h"
#include "gui/help/browser/helpbrowser.h"

#define HELP_LINK_SCHEME "help:"

FloatingHelpBrowser::FloatingHelpBrowser(QWidget *parent, QAbstractButton *button)
	:QTextBrowser(parent), mButton(button)
{
	QGraphicsDropShadowEffect *effect = new QGraphicsDropShadowEffect(this);
	effect->setBlurRadius(30.0);
	setGraphicsEffect(effect);

	setSizePolicy(QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum));
	setOpenLinks(false);
	hide();

	connect(this, SIGNAL(textChanged()), this, SLOT(textChanged()));

	if (mButton) {
		connect(mButton, SIGNAL(toggled(bool)), this, SLOT(showHelp(bool)));
	}
}

void FloatingHelpBrowser::showEvent(QShowEvent */*event*/)
{
	// The help topic of the page may have been set after the help text.
	updateText();

	if (mButton) {
		mButton->setChecked(true);
	}
}

void FloatingHelpBrowser::hideEvent(QHideEvent */*event*/)
{
	if (mButton) {
		mButton->setChecked(false);
	}
}

void FloatingHelpBrowser::mousePressEvent(QMouseEvent *event)
{
	const QString link = anchorAt(event->pos());

	hide();

	if (link.startsWith(HELP_LINK_SCHEME)) {
		HelpBrowser::showWindow(link.mid(QString(HELP_LINK_SCHEME).length()));
	}
}

void FloatingHelpBrowser::textChanged()
{
	if (mButton) {
		mButton->setVisible(!document()->isEmpty());
	}
}

void FloatingHelpBrowser::showHelp(bool state)
{
	QWidget *p = parentWidget();
	if (!p) {
		return;
	}

	resize(p->size() * 0.75);
	move(p->width() / 2 - width() / 2, p->height() / 2 - height() / 2);
	update();

	std::cerr << "Toggling help to " << state << std::endl;

	setVisible(state);
}

void FloatingHelpBrowser::setHelpText(const QString &helpText)
{
	mHelpText = helpText;
	updateText();
}

// Shows the help text, followed by a link to the manual when the page that
// owns this browser has a help topic.
void FloatingHelpBrowser::updateText()
{
	if (mHelpText.isEmpty()) {
		clear();
		return;
	}

	QString html = "<div align=\"justify\">" + mHelpText + "</div>";

	const QString topic = HelpBrowser::widgetTopic(parentWidget());
	if (!topic.isEmpty()) {
		html += "<p align=\"right\"><a href=\"" HELP_LINK_SCHEME + topic + "\">" + tr("Read more in the manual") + "</a></p>";
	}

	setHtml(html);
}
