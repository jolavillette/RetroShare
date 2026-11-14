/*******************************************************************************
 * gui/common/TextCursorProbe.h                                                *
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
#pragma once

class QTextEdit;

/* Diagnostics only (RS_TEXTCURSOR_PROBE). Watches the cursor a QTextEdit keeps
 * internally and reports the moment it stops being registered in its
 * document's cursor set: Qt then never resets it when the document dies, and
 * the edit crashes in ~QTextCursorPrivate when it is destroyed. Every report
 * is one "TEXTCURSOR-PROBE" line on stderr. */
namespace TextCursorProbe
{
	void attach(QTextEdit *edit, const char *name);
	void check(QTextEdit *edit, const char *where);
	void detach(QTextEdit *edit);
}
