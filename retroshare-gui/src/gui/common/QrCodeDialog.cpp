/*******************************************************************************
 * gui/common/QrCodeDialog.cpp                                                 *
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

#include "QrCodeDialog.h"

#include <QCoreApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QImage>
#include <QLabel>
#include <QMessageBox>
#include <QPixmap>
#include <QPushButton>
#include <QRegularExpression>
#include <QVBoxLayout>

#include <stdexcept>

#include "util/qrcodegen/qrcodegen.hpp"

namespace
{
QString tr(const char *text) { return QCoreApplication::translate("QrCodeDialog", text); }

/* Black modules on white, with the 4 modules of quiet zone a reader needs. */
QImage codeImage(const qrcodegen::QrCode& code, int targetSize)
{
	const int border = 4;
	const int modules = code.getSize() + 2 * border;
	const int scale = qMax(2, targetSize / modules);

	QImage image(modules * scale, modules * scale, QImage::Format_RGB32);
	image.fill(Qt::white);

	for (int y = 0; y < code.getSize(); ++y)
		for (int x = 0; x < code.getSize(); ++x)
			if (code.getModule(x, y))
				for (int dy = 0; dy < scale; ++dy)
					for (int dx = 0; dx < scale; ++dx)
						image.setPixel( (x + border) * scale + dx,
						                (y + border) * scale + dy, qRgb(0, 0, 0) );
	return image;
}
}

void QrCodeDialog::show(QWidget *parent, const QString& title, const QString& text)
{
	QString data = text;
	data.remove(QRegularExpression("\\s"));

	if (data.isEmpty())
		return;

	QImage image;
	try
	{
		image = codeImage( qrcodegen::QrCode::encodeText(
		                       data.toUtf8().constData(),
		                       qrcodegen::QrCode::Ecc::LOW ), 560 );
	}
	catch (const std::length_error&)
	{
		QMessageBox::information( parent, title,
		        tr("This is too long for a QR code. Use the short format, or copy the text.") );
		return;
	}

	QDialog dialog(parent);
	dialog.setWindowTitle(title);

	QLabel *label = new QLabel(&dialog);
	label->setPixmap(QPixmap::fromImage(image));
	label->setAlignment(Qt::AlignCenter);

	QDialogButtonBox *buttons = new QDialogButtonBox(QDialogButtonBox::Close, &dialog);
	QPushButton *save = buttons->addButton(tr("Save as image"), QDialogButtonBox::ActionRole);

	QObject::connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
	QObject::connect(save, &QPushButton::clicked, &dialog, [&dialog, image]()
	{
		QString file = QFileDialog::getSaveFileName( &dialog, tr("Save as image"),
		                                             "retroshare-id.png", tr("Images (*.png)") );
		if (!file.isEmpty() && !image.save(file, "PNG"))
			QMessageBox::warning(&dialog, tr("Save as image"), tr("Cannot write this file."));
	});

	QVBoxLayout *layout = new QVBoxLayout(&dialog);
	layout->addWidget(label);
	layout->addWidget(buttons);

	dialog.exec();
}
