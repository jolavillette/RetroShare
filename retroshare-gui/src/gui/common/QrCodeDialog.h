/*******************************************************************************
 * gui/common/QrCodeDialog.h                                                   *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
 * GNU Affero General Public License for more details.                         *
 *                                                                             *
 * You should have received a copy of the GNU Affero General Public License    *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.      *
 *                                                                             *
 *******************************************************************************/

#pragma once

#include <QString>

class QImage;
class QWidget;

namespace QrCodeDialog
{
/** Shows @p text as a QR code in a dialog, with a button to save it as an
 * image. Line breaks and spaces of @p text are left out of the code: a
 * certificate is cut in lines for display only. Tells the user when the text
 * is too long for a QR code. */
void show(QWidget *parent, const QString& title, const QString& text);

/** Reads the first QR code found in @p image into @p text.
 * @return false when the image holds no readable QR code. */
bool decode(const QImage& image, QString& text);
}
