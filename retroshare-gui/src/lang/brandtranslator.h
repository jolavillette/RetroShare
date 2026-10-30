/*******************************************************************************
 * lang/brandtranslator.h                                                      *
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

#ifndef _BRANDTRANSLATOR_H
#define _BRANDTRANSLATOR_H

#include <QTranslator>
#include <QList>

#define BRAND_NAME "NeoShare"

/**
 * Rebrands every translated UI string: it asks the wrapped translators (the
 * RetroShare one and the plugin ones) for the translation, falls back to the
 * source text, and replaces the upstream product name with BRAND_NAME.
 *
 * Qt queries the most recently installed translator first and stops at the
 * first non-null answer, so this one must be installed last. It answers null
 * for strings that do not mention the product and that no wrapped translator
 * knows, so Qt's own translator still handles its widgets.
 *
 * The wrapped list is fixed at construction: LanguageSupport uninstalls the
 * BrandTranslator before deleting any wrapped translator and installs a new
 * one afterwards (QCoreApplication::removeTranslator waits for in-flight
 * translate() calls).
 */
class BrandTranslator : public QTranslator
{
public:
	/** @param wrapped translators to consult, most recently installed first */
	explicit BrandTranslator(const QList<QTranslator*> &wrapped, QObject *parent = nullptr);

	QString translate(const char *context, const char *sourceText,
	                  const char *disambiguation = nullptr, int n = -1) const override;
	bool isEmpty() const override { return false; }

	/** Replaces the upstream product name with BRAND_NAME in @p text. */
	static QString rebrand(QString text);

private:
	QList<QTranslator*> mWrapped;
};

#endif
