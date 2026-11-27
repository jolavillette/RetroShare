/*******************************************************************************
 * retroshare-gui/src/gui/help/browser/helptextbrowser.h                       *
 *                                                                             *
 * Copyright (c) 2008, defnax          <retroshare.project@gmail.com>          *
 * Copyright (c) 2008, Matt Edman, Justin Hipple                               *
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

/*
** \file helptextbrowser.h
** \brief Displays a help page stored in the application resources
*/

#ifndef _HELPTEXTBROWSER_H
#define _HELPTEXTBROWSER_H

#include <QTextBrowser>
#include <QVariant>


class HelpTextBrowser : public QTextBrowser
{
  Q_OBJECT

public:
  /** Default constructor. */
  HelpTextBrowser(QWidget *parent = 0);
  /** Loads a resource into the browser. Help pages are read from
   * ":/help/<language>/", English being used for the pages that have no
   * translation. */
  QVariant loadResource(int type, const QUrl &name);

  /** Returns the resource path of the given help file, in the language of
   * the GUI when the file exists there, in English otherwise. */
  static QString resourcePath(const QString &fileName);

#if QT_VERSION >= QT_VERSION_CHECK(6,0,0)
protected:
  /** Called when the displayed document is changed. An external link is
   * opened in the Web browser after the user agreed, a help page is loaded. */
  virtual void doSetSource(const QUrl &url, QTextDocument::ResourceType type = QTextDocument::UnknownResource);
#else
public:
  /** Called when the displayed document is changed. An external link is
   * opened in the Web browser after the user agreed, a help page is loaded. */
  virtual void setSource(const QUrl &url);
#endif

private:
  /** Returns true when <b>url</b> points outside of the help pages. */
  static bool isExternalLink(const QUrl &url);
  /** Asks the user and opens <b>url</b> in the default Web browser. */
  void openExternalLink(const QUrl &url);
};

#endif

