/*******************************************************************************
 * retroshare-gui/src/gui/common/RsEdlideLabelItemDelegate.cpp                 *
 *                                                                             *
 * Copyright (C) 2010, Retroshare Team <retroshare.project@gmail.com>          *
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

#include "RSElidedItemDelegate.h"

#include "gui/common/ElidedLabel.h"
#include "util/rsdebug.h"

#include <QAbstractItemView>
#include <QApplication>
#include <QHash>
#include <QPainter>
#include <QTextDocument>
#include <QTextLayout>
#include <QTimer>
#include <QToolTip>

#include <QtMath>

#include <cmath>
#include <chrono>

//#define DEBUG_EID_PAINT 1

/* To test it you can make an empty.qss file with:
QTreeView::item, QTreeWidget::item, QListWidget::item{
	color: #AB0000;
	background-color: #00DC00;
}
QTreeView::item:selected, QTreeWidget::item:selected, QListWidget::item:selected{
	color: #00CD00;
	background-color: #0000BA;
}
QTreeView::item:hover, QTreeWidget::item:hover, QListWidget::item:hover{
	color: #0000EF;
	background-color: #FEDCBA;
}
QTreeView::item:selected:hover, QTreeWidget::item:selected:hover, QListWidget::item:selected:hover{
	color: #ABCDEF;
	background-color: #FE0000;
}

ForumsDialog, GxsForumThreadWidget
{
	qproperty-textColorRead: darkgray;
	qproperty-textColorUnread: white;
	qproperty-textColorUnreadChildren: red;
	qproperty-textColorNotSubscribed: white;
	qproperty-textColorMissing: darkred;
}
*/

RSElidedItemDelegate::RSElidedItemDelegate(QObject *parent)
  : RSStyledItemDelegate(parent)
  , mOnlyPlainText(false), mPaintRoundedRect(true)
{
}

QSize RSElidedItemDelegate::sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const
{
	QStyleOptionViewItem ownOption (option);
	initStyleOption(&ownOption, index);

	const QWidget* widget = option.widget;
	QStyle* ownStyle =  widget ? widget->style() : QApplication::style();

	//Only need "…" for text
	ownOption.text = "…";
	QRect checkRect = ownStyle->subElementRect(QStyle::SE_ItemViewItemCheckIndicator, &ownOption, widget);
	QRect iconRect = ownStyle->subElementRect(QStyle::SE_ItemViewItemDecoration, &ownOption, widget);
	QRect textRect = ownStyle->subElementRect(QStyle::SE_ItemViewItemText, &ownOption, widget);

	QSize contSize = ownStyle->sizeFromContents( QStyle::CT_ItemViewItem,&ownOption
	                                            ,QSize( checkRect.width()+iconRect.width()+textRect.width()
	                                                  , qMax(checkRect.height(),qMax(iconRect.height(),textRect.height()))
	                                                  ), widget ) ;

	contSize += QSize( 2*spacing().width()
	                 , qMax(checkRect.height(),iconRect.height()) > textRect.height()
	                   ? 0 : 2*spacing().height() );

	return contSize;
}

// WCAG relative luminance and contrast ratio (1:1 to 21:1).
static qreal relativeLuminance(const QColor& c)
{
	auto lin = [](qreal v) { return v <= 0.03928 ? v / 12.92 : std::pow((v + 0.055) / 1.055, 2.4); };
	return 0.2126 * lin(c.redF()) + 0.7152 * lin(c.greenF()) + 0.0722 * lin(c.blueF());
}

static qreal contrastRatio(const QColor& a, const QColor& b)
{
	const qreal la = relativeLuminance(a), lb = relativeLuminance(b);
	return (qMax(la, lb) + 0.05) / (qMin(la, lb) + 0.05);
}

// fg painted over an opaque bg, using fg's alpha.
static QColor blendOver(const QColor& fg, const QColor& bg)
{
	const qreal a = fg.alphaF();
	return QColor::fromRgbF(fg.redF()   * a + bg.redF()   * (1 - a),
	                        fg.greenF() * a + bg.greenF() * (1 - a),
	                        fg.blueF()  * a + bg.blueF()  * (1 - a));
}

inline QColor getImagePixelColor(QImage img, int x, int y)
{
#if QT_VERSION >= QT_VERSION_CHECK(5,6,0)
#ifdef DEBUG_EID_PAINT
	//RsDbg(" RSEID: Found Color ", img.pixelColor(x,y).name(QColor::HexArgb).toStdString(), " at ", x, ",", y);
#endif
	return img.pixelColor(x,y);
#else
	return img.pixel(x,y);
#endif
}

void RSElidedItemDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
	if(!index.isValid())
	{
		RS_ERR(" attempt to draw an invalid index.");
		return ;
	}
	painter->save();
	QStyleOptionViewItem ownOption (option);
	initStyleOption(&ownOption, index);
	//Prefer use icon from option
	if (!option.icon.isNull())
		ownOption.icon = option.icon;

#ifdef DEBUG_EID_PAINT
	RS_DBG("\n RSEID: Enter for item with text:", ownOption.text.toStdString());
#endif

	const QWidget* widget = option.widget;
	QStyle* ownStyle =  widget ? widget->style() : QApplication::style();

	if (!mOnlyPlainText)
	{
		QTextDocument td;
		td.setHtml(ownOption.text);
		ownOption.text = td.toPlainText();
	}
	// Get Font as option.font is not accurate. Merge it into the font of the view instead
	// of replacing it: a font stored in the model usually only sets a couple of attributes
	// (bold, italic) and expects the rest (family, size) to be inherited from the view.
	// QFont::resolve() keeps exactly what the model really set and takes the rest from the
	// view, which is what QStyledItemDelegate::initStyleOption() does as well.
	if (index.data(Qt::FontRole).type() == QVariant::Font) {
		QFont font = index.data(Qt::FontRole).value<QFont>().resolve(ownOption.font);
		ownOption.font = font;
		ownOption.fontMetrics = QFontMetrics(font);
#ifdef DEBUG_EID_PAINT
		RsDbg(" RSEID: Found font in model:", font.family().toStdString(), " size:", font.pointSize());
#endif
	}
	// Get Text color from model if one exists
	QColor textColor;
	if (index.data(Qt::ForegroundRole).isValid()) {
		//textColor = QColor(index.data(Qt::TextColorRole).toString());//Needs to pass from string else loose RBG format.
		textColor = index.data(Qt::ForegroundRole).value<QColor>();
#ifdef DEBUG_EID_PAINT
		RsDbg(" RSEID: Found text color in model:", textColor.name().toStdString());
#endif
	}
	// Get Brush from model if one exists
	QBrush bgBrush;
	bgBrush.setColor(QColor());// To get color().spec()==QColor::Invalid)
	const bool bgFromModel = index.data(Qt::BackgroundRole).isValid();
	if (bgFromModel) {
		bgBrush = index.data(Qt::BackgroundRole).value<QBrush>();
#ifdef DEBUG_EID_PAINT
		RsDbg(" RSEID: Found bg brush in model:", bgBrush.color().name().toStdString());
#endif
	}

	// If we get text and bg color from model, no need to retrieve it from base
	if ( (bgBrush.color().spec()==QColor::Invalid) || (textColor.spec()!=QColor::Invalid) )
	{
#ifdef DEBUG_EID_PAINT
		RsDbg( " RSEID:"
		     , ((bgBrush.color().spec()==QColor::Invalid) ? " Brush not defined" : "")
		     , ((textColor.spec()==QColor::Invalid) ? " Text Color not defined" : "")
		     , " so get it from base image.");
#endif
		// QPalette is not updated by StyleSheet all occurs in internal class. (QRenderRule)
		// https://code.woboq.org/qt5/qtbase/src/widgets/styles/qstylesheetstyle.cpp.html#4138
		//   void QStyleSheetStyle::drawControl(ControlElement ce, const QStyleOption *opt, QPainter *p, const QWidget *w) const
		//     case CE_ItemViewItem:
		// So we have to print it in Image to get colors by pixel
		QSize moSize=sizeHint(option,index);
		if (moSize.width() <= 20)
			moSize.setWidth(20);
#ifdef DEBUG_EID_PAINT
		RsDbg(" RSEID: for item size = ", moSize.width(), "x", moSize.height());
#endif

		QImage moImg(moSize,QImage::Format_ARGB32_Premultiplied);
		QPainter moPnt;
		moPnt.begin(&moImg);
		moPnt.setCompositionMode (QPainter::CompositionMode_Source);
		moPnt.fillRect(moImg.rect(), Qt::transparent);
		moPnt.setCompositionMode (QPainter::CompositionMode_SourceOver);

		QStyleOptionViewItem moOption (option);
		// Define option to get only what we want
		{
			moOption.rect = QRect(QPoint(0,0),moSize);
			moOption.state = ownOption.state;
			moOption.text = " ████████████████";//Add a blank char to get BackGround Color at top left
			// Remove unwanted info. Yes it can draw without that all public data ...
			moOption.backgroundBrush = QBrush();
			moOption.checkState = Qt::Unchecked;
			moOption.decorationAlignment = Qt::AlignLeft;
			moOption.decorationPosition = QStyleOptionViewItem::Left;
			moOption.decorationSize = QSize();
			moOption.displayAlignment = Qt::AlignLeft | Qt::AlignTop;
			moOption.features=QStyleOptionViewItem::ViewItemFeatures();
			moOption.font = QFont();
			moOption.icon = QIcon();
			moOption.index = QModelIndex();
			moOption.locale = QLocale();
			moOption.showDecorationSelected = false;
			moOption.textElideMode = Qt::ElideNone;
			moOption.viewItemPosition = QStyleOptionViewItem::Middle;
			//moOption.widget = nullptr; //Needed.

			moOption.direction = Qt::LayoutDirectionAuto;
			moOption.fontMetrics = QFontMetrics(QFont());
			// Keep the view's palette: it is what the real item is painted with. The
			// application palette can differ, e.g. on macOS item views get their own
			// (white selected text on the accent blue) while the application one has
			// black HighlightedText, which turned the selected row's text dark on blue.
			moOption.palette = option.palette;
			moOption.styleObject = nullptr;
		}
		//QStyledItemDelegate::paint(&moPnt, moOption, QModelIndex(index));//This update option now.
		ownStyle->drawControl(QStyle::CE_ItemViewItem, &moOption, &moPnt, widget);

		//// But these lines doesn't works.
		{
			//QStyleOptionViewItem moOptionsState;
			//moOptionsState.initFrom(moOption.widget);
			//moOptionsState.rect = QRect(QPoint(0,0),moSize);
			//moOptionsState.state = QStyle::State_MouseOver | QStyle::State_Enabled | QStyle::State_Sibling;
			//moOptionsState.text = "████████";
			//moOptionsState.widget = option.widget;
			//QStyledItemDelegate::paint(&moPnt, moOptionsState, QModelIndex());
		}

		moPnt.end();
#ifdef DEBUG_EID_PAINT
		// To save what it paint in application path
		moImg.save("image.png");
#endif
		// Get Color in this rect.
		{
			QColor moBGColor=getImagePixelColor(moImg,1,1); // BackGround may be paint.
			// The text color is the most frequent one that is not the background: the
			// probe text is made of full blocks, so their solid interior outnumbers the
			// anti-aliased glyph edges. Picking the first differing pixel could land on
			// such an edge and return a semi-transparent color, fading the item text.
			QColor moColor;
			QHash<QRgb, int> moCount;
			int moBest = 0;
			const QRgb moBGRgba = moBGColor.rgba();
			int moWidth = moImg.size().width(), moHeight = moImg.size().height();
			for (int x = 0; x < moWidth; x++)
				for (int y = 0; y < moHeight; y++)
				{
					const QRgb c = getImagePixelColor(moImg,x,y).rgba();
					if (c == moBGRgba)
						continue;
					const int n = ++moCount[c];
					if (n > moBest)
					{
						moBest = n;
						moColor = QColor::fromRgba(c);
					}
				}

			// If not found color is same as BackGround.
			if (moColor.spec() == QColor::Invalid)
				moColor = moBGColor;

			if (bgBrush.color().spec()==QColor::Invalid)
			{
				bgBrush = QBrush(moBGColor);
#ifdef DEBUG_EID_PAINT
				RsDbg(" RSEID: bg brush setted to ", moBGColor.name(QColor::HexArgb).toStdString());
#endif
			}
			if (textColor.spec()==QColor::Invalid)
			{
				textColor = moColor;
#ifdef DEBUG_EID_PAINT
				RsDbg(" RSEID: text color setted to ", moColor.name(QColor::HexArgb).toStdString());
#endif
			}
		}
	}

	// Last guard against an unreadable item. The probe above is a heuristic and can
	// get it wrong: on macOS (Qt 6.12, Fusion, Retina display, built-in default sheet)
	// it returned a shade of the selection band itself as the selected text color, so
	// the selected row's text vanished in the People and Preferences lists. Work out
	// the background the text really lands on and, when the chosen text contrasts
	// less than 3:1 with it, fall back to the palette's own text color for that state
	// (HighlightedText when selected, Text otherwise), then to black or white,
	// whichever reads better.
	{
		const bool selected = ownOption.state & QStyle::State_Selected;
		const QPalette::ColorGroup cg = !(ownOption.state & QStyle::State_Enabled) ? QPalette::Disabled
		                              : (ownOption.state & QStyle::State_Active) ? QPalette::Normal
		                                                                         : QPalette::Inactive;
		const QColor base = option.palette.color(cg, QPalette::Base);
		const QColor brushColor = bgBrush.color().isValid() ? bgBrush.color() : QColor(Qt::transparent);
		QColor seenBg;
		if (brushColor.alpha() == 255)
			seenBg = brushColor;
		else if (selected && !bgFromModel)
			// Not opaque: PE_PanelItemViewItem below paints the native selection band.
			seenBg = blendOver(option.palette.color(cg, QPalette::Highlight), base);
		else
			seenBg = blendOver(brushColor, base);

		if (contrastRatio(blendOver(textColor, seenBg), seenBg) < 3.0)
		{
			const QColor paletteText = option.palette.color(cg, selected ? QPalette::HighlightedText : QPalette::Text);
			if (contrastRatio(blendOver(paletteText, seenBg), seenBg) >= 3.0)
				textColor = paletteText;
			else
				textColor = contrastRatio(Qt::black, seenBg) >= contrastRatio(Qt::white, seenBg)
				            ? QColor(Qt::black) : QColor(Qt::white);
		}
	}

	painter->setPen(textColor);
	painter->setBrush(bgBrush);
	ownOption.backgroundBrush = bgBrush;

	// Code from: https://code.woboq.org/qt5/qtbase/src/widgets/styles/qcommonstyle.cpp.html#2271
	QRect checkRect = ownStyle->subElementRect(QStyle::SE_ItemViewItemCheckIndicator, &ownOption, widget);
	QRect iconRect = ownStyle->subElementRect(QStyle::SE_ItemViewItemDecoration, &ownOption, widget);
	QRect textRect = ownStyle->subElementRect(QStyle::SE_ItemViewItemText, &ownOption, widget);

	// Draw the background
	if (bgBrush.color().alpha() == 0)
		// No BackGround Color found, use default delegate to draw it.
		ownStyle->proxy()->drawPrimitive(QStyle::PE_PanelItemViewItem, &option, painter, widget);// This prefer draw StyleSheet bg than item one.
	else
		painter->fillRect(ownOption.rect,bgBrush);

#ifdef DEBUG_EID_PAINT
	{
		QStyleOptionViewItem tstOption = option;
		// Reduce rect to get this item bg color external and base internal
		tstOption.rect.adjust(2,2,-2,-2);
		// To draw with base for debug purpose
		RSStyledItemDelegate::paint(painter, tstOption, index);
	}
#endif

	// draw the check mark
	if (ownOption.features & QStyleOptionViewItem::HasCheckIndicator) {
		QStyleOptionViewItem cmOption(*&ownOption);
		cmOption.rect = checkRect;
		cmOption.state = cmOption.state & ~QStyle::State_HasFocus;
		switch (ownOption.checkState) {
			case Qt::Unchecked:
				cmOption.state |= QStyle::State_Off;
			break;
			case Qt::PartiallyChecked:
				cmOption.state |= QStyle::State_NoChange;
			break;
			case Qt::Checked:
				cmOption.state |= QStyle::State_On;
			break;
		}
		ownStyle->drawPrimitive(QStyle::PE_IndicatorItemViewItemCheck, &cmOption, painter, widget);
	}
	// draw the icon
	{
		if (!ownOption.icon.isNull())
		{
			QString status;
			if (index.data(Qt::StatusTipRole).canConvert(QMetaType::QString))
				status = index.data(Qt::StatusTipRole).toString();

			// Draw item Icon
			{
				QIcon::Mode mode = QIcon::Normal;
				if (!(ownOption.state & QStyle::State_Enabled))
					mode = QIcon::Disabled;
				else if (ownOption.state & QStyle::State_Selected)
					mode = QIcon::Selected;
				QIcon::State state = (ownOption.state & QStyle::State_Open) ? QIcon::On : QIcon::Off;
				ownOption.icon.paint(painter, iconRect, ownOption.decorationAlignment, mode, state);
			}
			// Then overlay with waiting
			if (status.toLower() == "waiting")
			{
				const QAbstractItemView* aiv = dynamic_cast<const QAbstractItemView*>(option.widget);
				if (aiv)
					QTimer::singleShot(200, aiv->viewport(), SLOT(update()));

				QSize waitSize=iconRect.size();
				qreal diag = qMin(waitSize.height(),waitSize.height())*std::sqrt(2);
				auto now = std::chrono::system_clock::now().time_since_epoch();
				auto s = std::chrono::duration_cast<std::chrono::seconds>(now).count();
				auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - std::chrono::duration_cast<std::chrono::seconds>(now)).count();
				int duration = 3;// Time (s) to make a revolution.
				auto time = (s%duration)*1000 + ms;
				qreal angle = 360.0*(time/(duration*1000.0));
				qreal add = 120*(time/(duration*1000.0))*std::abs(sin(qDegreesToRadians(angle/2)));
				painter->setPen(QPen(QBrush(ownOption.palette.color(QPalette::Normal, QPalette::WindowText)),diag/10,Qt::DotLine,Qt::RoundCap));
				painter->drawEllipse( iconRect.x()+iconRect.width() /2 + (diag/4)*cos(qDegreesToRadians(angle      ))
				                      , iconRect.y()+iconRect.height()/2 + (diag/4)*sin(qDegreesToRadians(angle      )), 1, 1);
				painter->setPen(QPen(QBrush(ownOption.palette.color(QPalette::Normal, QPalette::Midlight)),diag/10,Qt::DotLine,Qt::RoundCap));
				painter->drawEllipse( iconRect.x()+iconRect.width() /2 + (diag/4)*cos(qDegreesToRadians(angle-  add))
				                      , iconRect.y()+iconRect.height()/2 + (diag/4)*sin(qDegreesToRadians(angle-  add)), 1, 1);
				painter->setPen(QPen(QBrush(ownOption.palette.color(QPalette::Normal, QPalette::Window)),diag/10,Qt::DotLine,Qt::RoundCap));
				painter->drawEllipse( iconRect.x()+iconRect.width() /2 + (diag/4)*cos(qDegreesToRadians(angle-2*add))
				                      , iconRect.y()+iconRect.height()/2 + (diag/4)*sin(qDegreesToRadians(angle-2*add)), 1, 1);
			}
		}
	}
	// draw the text
	if (!ownOption.text.isEmpty()) {
#ifdef DEBUG_EID_PAINT
		// To draw text near base one.
		ownOption.text = ownOption.text.prepend("__");
#endif

		QTextOption::WrapMode wm = (ownOption.features & QStyleOptionViewItem::WrapText) ? QTextOption::WordWrap : QTextOption::NoWrap;
		const int textHMargin = ownStyle->pixelMetric(QStyle::PM_FocusFrameHMargin, nullptr, widget) + 1;
		const int textVMargin = ownStyle->pixelMetric(QStyle::PM_FocusFrameVMargin, nullptr, widget) + 1;
		textRect = textRect.adjusted(textHMargin, textVMargin, -textHMargin, -textVMargin); // remove width padding

		ElidedLabel::paintElidedLine(painter,ownOption.text,textRect,ownOption.font,ownOption.displayAlignment,wm,mPaintRoundedRect);
	}
	painter->restore();
#ifdef DEBUG_EID_PAINT
	RsDbg(" RSEID: Finished");
#endif
}

bool RSElidedItemDelegate::editorEvent(QEvent *event, QAbstractItemModel *model, const QStyleOptionViewItem &option, const QModelIndex &index)
{
	if (event->type() == QEvent::MouseButtonPress) {
		QMouseEvent *ev = static_cast<QMouseEvent *>(event);
		if (ev) {
			if (ev->buttons()==Qt::LeftButton) {
#ifdef DEBUG_EID_PAINT
				QVariant var = index.data();
				Q_UNUSED(var);
#endif
				if (index.data().type() == QVariant::String) {
					QString text = index.data().toString();
					if (!text.isEmpty()) {

						QStyleOptionViewItem ownOption (option);
						initStyleOption(&ownOption, index);

						const QWidget* widget = option.widget;
						QStyle* ownStyle =  widget ? widget->style() : QApplication::style();

						if (!mOnlyPlainText)
						{
							QTextDocument td;
							td.setHtml(ownOption.text);
							ownOption.text = td.toPlainText();
						}
						//Get Font as option.font is not accurate, merged into the one of the view (see paint())
						if (index.data(Qt::FontRole).type() == QVariant::Font) {
							QFont font = index.data(Qt::FontRole).value<QFont>().resolve(ownOption.font);
							ownOption.font = font;
							ownOption.fontMetrics = QFontMetrics(font);
						}
						QRect textRect = ownStyle->subElementRect(QStyle::SE_ItemViewItemText, &ownOption, widget);
						const int textHMargin = ownStyle->pixelMetric(QStyle::PM_FocusFrameHMargin, nullptr, widget) + 1;
						const int textVMargin = ownStyle->pixelMetric(QStyle::PM_FocusFrameVMargin, nullptr, widget) + 1;
						textRect = textRect.adjusted(textHMargin, textVMargin, -textHMargin, -textVMargin); // remove width padding

						QTextLayout textLayout(text, ownOption.font);
						QTextOption to = textLayout.textOption();
						to.setWrapMode((ownOption.features & QStyleOptionViewItem::WrapText) ? QTextOption::WordWrap : QTextOption::NoWrap);

						//Update RSElidedItemDelegate as only one delegate for all items
						QRect rectElision;
						bool elided = ElidedLabel::paintElidedLine(nullptr,text,textRect,ownOption.font,ownOption.displayAlignment,to.wrapMode(),true,&rectElision);
						if (elided && (rectElision.contains(ev->pos()))){
							QToolTip::showText(ev->globalPos(),QString("<FONT>") + text + QString("</FONT>"));
							return true; // eat event
						}
					}
				}
			}
		}
	}
	return RSStyledItemDelegate::editorEvent(event, model, option, index);
}
