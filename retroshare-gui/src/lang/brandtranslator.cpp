/*******************************************************************************
 * lang/brandtranslator.cpp                                                    *
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

#include <QRegularExpression>

#include "brandtranslator.h"

/* Matches the product name as a standalone word only. Case-sensitive on
 * purpose: lowercase "retroshare" belongs to identifiers that must not change
 * (retroshare:// links, libretroshare, retroshare.cc). The lookbehind skips
 * URL and path parts (github.com/RetroShare/RetroShare), the lookahead skips
 * identifiers (RetroShareLink) and file names (RetroShare.conf). */
static const QRegularExpression &upstreamName()
{
	static const QRegularExpression re(
	        QStringLiteral("(?<![\\w/.])(?:RetroShare|Retroshare|(RETROSHARE))(?![\\w]|\\.conf)"));
	return re;
}

BrandTranslator::BrandTranslator(const QList<QTranslator*> &wrapped, QObject *parent)
    : QTranslator(parent), mWrapped(wrapped)
{
}

QString BrandTranslator::translate(const char *context, const char *sourceText,
                                   const char *disambiguation, int n) const
{
	for (QTranslator *translator : mWrapped) {
		QString result = translator->translate(context, sourceText, disambiguation, n);
		if (!result.isNull())
			return rebrand(result);
	}

	/* No translation: answer with the rebranded source text, or null when
	 * there is nothing to rebrand so the next translators get their turn. */
	const QString source = QString::fromUtf8(sourceText);
	const QString rebranded = rebrand(source);
	return rebranded == source ? QString() : rebranded;
}

QString BrandTranslator::rebrand(QString text)
{
	if (!text.contains(QLatin1String("etro"), Qt::CaseInsensitive))
		return text;

	static const QString brand = QStringLiteral(BRAND_NAME);
	static const QString brandUpper = brand.toUpper();

	QString out;
	int last = 0;
	QRegularExpressionMatchIterator it = upstreamName().globalMatch(text);
	while (it.hasNext()) {
		const QRegularExpressionMatch m = it.next();
		out += text.mid(last, m.capturedStart() - last);
		out += m.capturedLength(1) ? brandUpper : brand;
		last = m.capturedEnd();
	}
	if (last == 0)
		return text;

	out += text.mid(last);
	return out;
}
