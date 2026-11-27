/*******************************************************************************
 * retroshare-gui/src/gui/help/browser/helptextbrowser.cpp                     *
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
** \file helptextbrowser.cpp
** \brief Displays a help page stored in the application resources
*/

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QDesktopServices>
#include "gui/common/vmessagebox.h"
#include "gui/common/rshtml.h"
#include <rshare.h>

#include "helptextbrowser.h"

/* The help pages are written in Markdown. Qt reads Markdown since 5.14; with
 * an older Qt, or a Qt built without the Markdown reader, the pages are shown
 * as plain text. */
#if QT_VERSION >= QT_VERSION_CHECK(5,14,0)
#if QT_CONFIG(textmarkdownreader)
#define HELP_MARKDOWN_READER
#endif
#endif

static bool isMarkdownFile(const QUrl &url)
{
  return url.path().endsWith(".md", Qt::CaseInsensitive);
}

/** Default constructor. */
HelpTextBrowser::HelpTextBrowser(QWidget *parent)
  : QTextBrowser(parent)
{
  /* Links are followed here, so that every one of them goes through the
   * external link check, whatever the Qt version. */
  setOpenExternalLinks(false);
  setOpenLinks(false);
  connect(this, &QTextBrowser::anchorClicked, this, [this](const QUrl &url) { setSource(url); });
}

/** Returns the resource path of the given help file, in the language of the
 * GUI when the file exists there, in English otherwise. */
QString
HelpTextBrowser::resourcePath(const QString &fileName)
{
  QString path = ":/help/" + RsApplication::language() + "/" + fileName;
  if (!QFile::exists(path))
    path = ":/help/en/" + fileName;
  return path;
}

/** Loads a resource into the browser. Help pages are read as UTF-8 from the
 * application resources. */
QVariant
HelpTextBrowser::loadResource(int type, const QUrl &name)
{
  if (type == QTextDocument::HtmlResource
#ifdef HELP_MARKDOWN_READER
      || type == QTextDocument::MarkdownResource
#endif
     ) {
    QFile file(resourcePath(QFileInfo(name.path()).fileName()));
    if (!file.open(QIODevice::ReadOnly)) {
      return tr("Error opening help file:")+" "+ name.path();
    }
    QString text = QString::fromUtf8(file.readAll());

#ifndef HELP_MARKDOWN_READER
    if (isMarkdownFile(name)) {
      /* No Markdown reader: the source is readable as it is. */
      return "<pre style=\"white-space: pre-wrap;\">" + text.toHtmlEscaped() + "</pre>";
    }
#endif
    return text;
  }
  /* Everything else, just let QTextBrowser take care of it. */
  return QTextBrowser::loadResource(type, name);
}

/** Returns true when <b>url</b> points outside of the help pages. */
bool
HelpTextBrowser::isExternalLink(const QUrl &url)
{
  return url.scheme() != "qrc" && !url.isRelative();
}

/** Asks the user and opens <b>url</b> in the default Web browser. */
void
HelpTextBrowser::openExternalLink(const QUrl &url)
{
  int ret = VMessageBox::question(this,
              tr("Opening External Link"),
              p(tr("RetroShare can open the link you selected in your default "
                   "Web browser. If your browser is not currently "
                   "configured to use Tor then the request will not be "
                   "anonymous.")) +
              p(tr("Do you want Retroshare to open the link in your Web "
                   "browser?")),
              VMessageBox::Yes|VMessageBox::Default,
              VMessageBox::Cancel|VMessageBox::Cancel);

  if (ret == VMessageBox::Cancel)
    return;

  bool ok = QDesktopServices::openUrl(url);
  if (!ok) {
    VMessageBox::information(this,
      tr("Unable to Open Link"),
      tr("RetroShare was unable to open the selected link in your Web browser. "
         "You can still copy the URL and paste it into your browser."),
      VMessageBox::Ok);
  }
}

#if QT_VERSION >= QT_VERSION_CHECK(6,0,0)
void
HelpTextBrowser::doSetSource(const QUrl &url, QTextDocument::ResourceType /*type*/)
{
  if (isExternalLink(url)) {
    openExternalLink(url);
    return;
  }
#ifdef HELP_MARKDOWN_READER
  if (isMarkdownFile(url)) {
    QTextBrowser::doSetSource(url, QTextDocument::MarkdownResource);
    return;
  }
#endif
  QTextBrowser::doSetSource(url, QTextDocument::HtmlResource);
}
#else
void
HelpTextBrowser::setSource(const QUrl &url)
{
  if (isExternalLink(url)) {
    openExternalLink(url);
    return;
  }
#if QT_VERSION >= QT_VERSION_CHECK(5,14,0)
#ifdef HELP_MARKDOWN_READER
  if (isMarkdownFile(url)) {
    QTextBrowser::setSource(url, QTextDocument::MarkdownResource);
    return;
  }
#endif
  QTextBrowser::setSource(url, QTextDocument::HtmlResource);
#else
  QTextBrowser::setSource(url);
#endif
}
#endif

