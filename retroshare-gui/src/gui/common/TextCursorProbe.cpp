/*******************************************************************************
 * gui/common/TextCursorProbe.cpp                                              *
 *                                                                             *
 * Copyright (C) 2026 RetroShare Team <retroshare.project@gmail.com>           *
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
#include "TextCursorProbe.h"

#include <QDateTime>
#include <QEvent>
#include <QHash>
#include <QMetaEnum>
#include <QObject>
#include <QSet>
#include <QTextCursor>
#include <QTextDocument>
#include <QTextEdit>
#include <QTimer>

#include <private/qtextcursor_p.h>
#include <private/qtextdocument_p.h>

#include <cstdio>
#include <cstdlib>

#ifdef _WIN32
#	include <windows.h>
#endif

#include "util/stacktrace.h"

namespace {

/* Read two private members without patching Qt: an explicit template
 * instantiation may name a private member ([temp.explicit]/12), and the friend
 * it defines hands the member pointer out. */
template<typename Tag, typename Tag::type M>
struct Rob { friend typename Tag::type get(Tag) { return M; } };

struct DocCursors {
	using type = QSet<QTextCursorPrivate *> QTextDocumentPrivate::*;
	friend type get(DocCursors);
};
template struct Rob<DocCursors, &QTextDocumentPrivate::cursors>;

struct CursorD {
	using type = QSharedDataPointer<QTextCursorPrivate> QTextCursor::*;
	friend type get(CursorD);
};
template struct Rob<CursorD, &QTextCursor::d>;

struct State
{
	const QTextCursorPrivate *cursorPriv = nullptr;
	const QTextDocumentPrivate *docPriv = nullptr;
	bool ok = true;
	QString lastGood = "attach";
	QString name;
	QObject *filter = nullptr;
};

QHash<QTextEdit *, State> gStates;

QString now() { return QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss.zzz"); }

void report(QTextEdit *edit, const State &s, const QString &msg)
{
	const QString name = edit->objectName().isEmpty() ? s.name : edit->objectName();
	fprintf(stderr, "%s TEXTCURSOR-PROBE %s edit=%p window=\"%s\" %s\n",
	        now().toUtf8().constData(), name.toUtf8().constData(),
	        static_cast<void *>(edit),
	        edit->window()->windowTitle().toUtf8().constData(),
	        msg.toUtf8().constData());
	fflush(stderr);
}

QString ptr(const void *p) { return QString("0x%1").arg(quintptr(p), 0, 16); }

class Filter : public QObject
{
public:
	Filter(QTextEdit *edit) : QObject(edit), mEdit(edit) {}
	bool eventFilter(QObject *, QEvent *e) override
	{
		const char *n = QMetaEnum::fromType<QEvent::Type>().valueToKey(e->type());
		TextCursorProbe::check(mEdit, n ? n : "event");
		return false;
	}
private:
	QTextEdit *mEdit;
};

} // namespace

void TextCursorProbe::attach(QTextEdit *edit, const char *name)
{
	State &s = gStates[edit];
	s.name = name;
	s.filter = new Filter(edit);
	edit->installEventFilter(s.filter);

	QObject::connect(edit, &QTextEdit::textChanged, s.filter, [edit]() { check(edit, "textChanged"); });
	QObject::connect(edit, &QTextEdit::cursorPositionChanged, s.filter, [edit]() { check(edit, "cursorPositionChanged"); });
	QObject::connect(edit, &QTextEdit::selectionChanged, s.filter, [edit]() { check(edit, "selectionChanged"); });
	QObject::connect(edit, &QTextEdit::currentCharFormatChanged, s.filter, [edit]() { check(edit, "currentCharFormatChanged"); });

	auto *timer = new QTimer(s.filter);
	QObject::connect(timer, &QTimer::timeout, s.filter, [edit]() { check(edit, "timer"); });
	timer->start(1000);

	report(edit, s, "attached");
	check(edit, "attach");
}

void TextCursorProbe::check(QTextEdit *edit, const char *where)
{
	auto it = gStates.find(edit);
	if (it == gStates.end())
		return;
	State &s = *it;

	/* textCursor() returns a copy that shares the control's private: read it
	 * through constData() so that nothing detaches. */
	const QTextCursor c = edit->textCursor();
	const QTextCursorPrivate *cp = (c.*get(CursorD())).constData();
	const QTextDocumentPrivate *dp = QTextDocumentPrivate::get(edit->document());

	const bool sameDoc = cp && cp->priv == dp;
	const bool registered = cp && cp->priv
	        && (cp->priv->*get(DocCursors())).contains(const_cast<QTextCursorPrivate *>(cp));
	/* The copy c itself holds no extra private: it shares cp. */
	const bool ok = cp && sameDoc && registered;

	if (cp != s.cursorPriv || dp != s.docPriv) {
		if (!s.ok || !ok)
			report(edit, s, QString("cursor private %1 -> %2, document private %3 -> %4, at %5")
			       .arg(ptr(s.cursorPriv), ptr(cp), ptr(s.docPriv), ptr(dp), where));
		s.cursorPriv = cp;
		s.docPriv = dp;
	}

	if (ok) {
		if (!s.ok)
			report(edit, s, QString("RECOVERED at %1").arg(where));
		s.ok = true;
		s.lastGood = where;
		return;
	}
	if (!s.ok)
		return;

	s.ok = false;
	report(edit, s, QString("LOST at %1, last good at %2: cursor private %3 cursor->priv %4 "
	                        "document private %5 sameDoc=%6 registered=%7 setSize=%8 "
	                        "inContentsChange=%9 position=%10 chars=%11")
	       .arg(where, s.lastGood, ptr(cp), ptr(cp ? cp->priv : nullptr), ptr(dp))
	       .arg(sameDoc).arg(registered)
	       .arg(cp && cp->priv ? (cp->priv->*get(DocCursors())).size() : -1)
	       .arg(dp ? dp->inContentsChange : false)
	       .arg(c.position()).arg(edit->document()->characterCount()));
	print_stacktrace(true, stderr);

#ifdef _WIN32
	if (IsDebuggerPresent() && getenv("RS_TEXTCURSOR_PROBE_BREAK"))
		DebugBreak();
#endif
}

void TextCursorProbe::detach(QTextEdit *edit)
{
	auto it = gStates.find(edit);
	if (it == gStates.end())
		return;
	check(edit, "destructor");
	report(edit, *it, it->ok ? "destroyed, cursor registered" : "destroyed with an UNREGISTERED cursor: the crash follows");
	edit->removeEventFilter(it->filter);
	gStates.erase(it);
}
