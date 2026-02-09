/*******************************************************************************
 * gui/elastic/elnode.cpp                                                      *
 *                                                                             *
 * Copyright (c) 2012, RetroShare Team <retroshare.project@gmail.com>          *
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

// This code is inspired from http://doc.qt.io/qt-5/qtwidgets-graphicsview-elasticnodes-node-cpp.html

#include "gui/common/FilesDefs.h"
#include "gui/settings/rsharesettings.h"
#include <math.h>

#include <QApplication>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QAction>
#include <QMenu>
#include <QStyleOption>
#include <iostream>

#include <gui/connect/PGPKeyDialog.h>

#include <retroshare/rspeers.h>
#include "edge.h"
#include "elnode.h"
#include "graphwidget.h"
#include "util/RsQtVersion.h"

#define IMAGE_AUTHED         ":/images/accepted16.png"
#define IMAGE_DENIED         ":/images/denied16.png"
#define IMAGE_TRUSTED        ":/images/rs-2.png"
#define IMAGE_MAKEFRIEND     ":/images/user/add_user16.png"

// Node *Node::_selected_node = NULL ; // Supprimé pour isolation

Node::Node(const std::string& node_string,GraphWidget::NodeType type,GraphWidget::AuthType auth,GraphWidget *graphWidget,const RsPeerId& ssl_id,const RsPgpId& gpg_id)
    : graph(graphWidget),_desc_string(node_string),_type(type),_auth(auth),_ssl_id(ssl_id),_gpg_id(gpg_id)
{
#ifdef DEBUG_ELASTIC
	std::cerr << "Created node type " << type << ", string=" << node_string << std::endl ;
#endif
    setFlag(ItemIsMovable);
#if QT_VERSION >= 0x040600
    setFlag(ItemSendsGeometryChanges);
#endif
    setCacheMode(DeviceCoordinateCache);
    setZValue(1);
	 mDeterminedBB = false ;
	 mBBWidth = 0 ;
	 mNodeDrawSize = 20;

	 _speedx=_speedy=0;
	 _steps=0;

	if(_type == GraphWidget::ELASTIC_NODE_TYPE_OWN)
		_auth = GraphWidget::ELASTIC_NODE_AUTH_FULL ;
}

Node::~Node()
{
    foreach (Edge *edge, edgeList) {
        if (edge->sourceNode() == this) edge->setSourceNode(NULL);
        if (edge->destNode() == this) edge->setDestNode(NULL);
    }
}

const float Node::MASS_FACTOR = 50.0f ;
const float Node::FRICTION_FACTOR = 20.0f ; 
const float Node::REPULSION_FACTOR = 10.0f ;
const float Node::NODE_DISTANCE = 100.0f ;

void Node::addEdge(Edge *edge)
{
    if (edge && !edgeList.contains(edge)) {
        edgeList << edge;
        edge->adjust();
    }
}

void Node::removeEdge(Edge *edge)
{
    edgeList.removeAll(edge);
}

const QList<Edge *>& Node::edges() const
{
    return edgeList;
}

//static double interpolate(const double *map,int W,int H,float x,float y)
//{
//	if(x>W-2) x=W-2 ;
//	if(y>H-2) y=H-2 ;
//	if(x<0  ) x=0   ;
//	if(y<0  ) y=0   ;
//
//	int i=(int)floor(x) ;
//	int j=(int)floor(y) ;
//	double di = x-i ;
//	double dj = y-j ;
//
//	return (1-di)*( (1-dj)*map[2*(i+W*j)] + dj*map[2*(i+W*(j+1))])
//				+di *( (1-dj)*map[2*(i+1+W*j)] + dj*map[2*(i+1+W*(j+1))]) ;
//}

void Node::calculateForces(const double *map, int /*width*/, int /*height*/, int W, int H, float x, float y, float friction_factor)
{
	if (!scene() || scene()->mouseGrabberItem() == this)
	{
		newPos = pos();
		return;
	}

	// Sum up all forces pushing this item away
	qreal xforce = 0;
	qreal yforce = 0;

	float dei = 0.0f;
	float dej = 0.0f;

	const int KS = 5;
	float *e = graph->repulsionKernel();

	// 1. Force de répulsion FFT (basée sur la grille FFT)
	for (int i = -KS; i <= KS; ++i)
		for (int j = -KS; j <= KS; ++j)
		{
			int X = std::min(W - 1, std::max(0, (int)rint(x)));
			int Y = std::min(H - 1, std::max(0, (int)rint(y)));

			float val = map[2 * ((i + X + W) % W + W * ((j + Y + H) % H))];

			dei += i * e[i + KS + (2 * KS + 1) * (j + KS)] * val;
			dej += j * e[i + KS + (2 * KS + 1) * (j + KS)] * val;
		}

	xforce = REPULSION_FACTOR * dei / 25.0;
	yforce = REPULSION_FACTOR * dej / 25.0;

	// Sécurité anti-explosion (NaN/Inf)
	if (!std::isfinite(xforce)) xforce = 0;
	if (!std::isfinite(yforce)) yforce = 0;

	// 2. Force d'attraction des liens (Ressorts)
	double weight = (n_edges() + 1);

	foreach (Edge *edge, edgeList)
	{
		QPointF edge_pos;
		double w2;
		if (edge->sourceNode() == this)
		{
			edge_pos = mapFromItem(edge->destNode(), 0, 0);
			w2 = sqrtf(std::min(n_edges(), edge->destNode()->n_edges()));
		}
		else
		{
			edge_pos = mapFromItem(edge->sourceNode(), 0, 0);
			w2 = sqrtf(std::min(n_edges(), edge->sourceNode()->n_edges()));
		}

		float dist = sqrtf(edge_pos.x() * edge_pos.x() + edge_pos.y() * edge_pos.y());
		float val = dist - graph->edgeLength() * w2;

		xforce += 0.10 * edge_pos.x() * val / weight;
		yforce += 0.10 * edge_pos.y() * val / weight;
	}

	// 3. Friction (Amortissement)
	// On utilise le friction_factor qui augmente avec le temps pour stabiliser le graphe
	xforce -= FRICTION_FACTOR * friction_factor * _speedx;
	yforce -= FRICTION_FACTOR * friction_factor * _speedy;

	// 4. Physique du cadre et Gravité
	QPointF current_pos = pos();
	QRectF scene_rect = scene()->sceneRect();
	QPointF center = scene_rect.center();

	// Gravité vers le centre de la vue (pas (0,0) qui est en haut à gauche)
	float k_gravity = 0.05f;
	xforce += (center.x() - current_pos.x()) * k_gravity;
	yforce += (center.y() - current_pos.y()) * k_gravity;

	// Murs de sécurité basés sur le sceneRect
	float margin = 50.0f;
	float k_wall = 50.0f;

	if (current_pos.x() < scene_rect.left() + margin)
		xforce += k_wall * (scene_rect.left() + margin - current_pos.x());
	if (current_pos.x() > scene_rect.right() - margin)
		xforce -= k_wall * (current_pos.x() - (scene_rect.right() - margin));
	if (current_pos.y() < scene_rect.top() + margin)
		yforce += k_wall * (scene_rect.top() + margin - current_pos.y());
	if (current_pos.y() > scene_rect.bottom() - margin)
		yforce -= k_wall * (current_pos.y() - (scene_rect.bottom() - margin));

	// Application de l'accélération avec garde (Military Grade : checking all components)
	if (std::isfinite(xforce) && std::isfinite(yforce) && std::isfinite(_speedx) && std::isfinite(_speedy))
	{
		_speedx += xforce / MASS_FACTOR;
		_speedy += yforce / MASS_FACTOR;
	}
	else
	{
		_speedx = 0 ;
		_speedy = 0 ;
	}

	// Limitation de la vitesse (Aérospatiale) : Clamping strict
	const float max_speed = 30.0f;
	_speedx = std::max(-max_speed, std::min(max_speed, (float)_speedx));
	_speedy = std::max(-max_speed, std::min(max_speed, (float)_speedy));

	newPos = current_pos + QPointF(_speedx, _speedy);
}

bool Node::progress()
{
	if(_type == GraphWidget::ELASTIC_NODE_TYPE_OWN)
		return false;

	float f = std::max(fabs(newPos.x() - pos().x()), fabs(newPos.y() - pos().y())) ;

    setPos(newPos);

	// Garde ultime : Si la position devient invalide, reset au centre
	if(!std::isfinite(newPos.x()) || !std::isfinite(newPos.y()))
	{
		setPos(scene() ? scene()->sceneRect().center() : QPointF(0,0));
		_speedx = _speedy = 0 ;
	}

    return f > 0.5;
}

QRectF Node::boundingRect() const
{
        float m = QFontMetricsF(graph->font()).height();
        float f = m/16.0;
    qreal adjust = 2*f;
    /* add in the size of the text */
    qreal realwidth = 40*f;
    
    if (mDeterminedBB)
    {
	realwidth = mBBWidth + adjust;
    }
    if (realwidth < 23*f + adjust)
    {
    	realwidth = 23*f + adjust;
    }

    return QRectF(-10*f - adjust, -10*f - adjust, realwidth, 23*f + adjust);
}

QPainterPath Node::shape() const
{
    QPainterPath path;
    path.addEllipse(-10, -10, 20, 20);
    return path;
}

#if QT_VERSION >= 0x040700
static QColor lightdark(const QColor& col,int l,int d)
{
	int h,s,v ;

   col.getHsv(&h,&s,&v) ;

	v = (int)floor(v * l/(float)d ) ;

	return QColor::fromHsv(h,s,v) ;
}
#endif

void Node::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *)
{
	static QColor type_color[4] = { QColor(255, 230, 0), QColor(0, 255, 130), QColor(0, 220, 255), QColor(180, 180, 180) };

	QColor col0;

	Node* selNode = graph ? graph->selectedNode() : NULL;
	if (selNode == NULL)
		col0 = type_color[_type];
	else if (selNode == this)
		col0 = type_color[0];
	else
	{
		bool found = false;
		foreach (Edge *edge, edgeList)
			if (edge->sourceNode() == selNode || edge->destNode() == selNode)
			{
				col0 = type_color[1];
				found = true;
				break;
			}

		if (!found)
		{
			col0 = type_color[2];
			col0.setAlpha(100); // Semi-transparent for non-related
		}
	}

	painter->setRenderHint(QPainter::Antialiasing);

	// Ombre portée (subtile)
	painter->setPen(Qt::NoPen);
	painter->setBrush(QColor(0, 0, 0, 50));
	int mNodeDrawSize2 = mNodeDrawSize / 2;
	painter->drawEllipse(-mNodeDrawSize2 + 2, -mNodeDrawSize2 + 2, mNodeDrawSize, mNodeDrawSize);

	// Gradient principal pour le nœud
	QRadialGradient gradient(-mNodeDrawSize / 5, -mNodeDrawSize / 5, mNodeDrawSize * 0.7);
	if (option->state & QStyle::State_Sunken)
	{
		gradient.setCenter(mNodeDrawSize / 5, mNodeDrawSize / 5);
		gradient.setFocalPoint(mNodeDrawSize / 5, mNodeDrawSize / 5);
		gradient.setColorAt(0, col0.darker(150));
		gradient.setColorAt(1, col0.darker(100));
	}
	else
	{
		gradient.setColorAt(0, col0.lighter(120));
		gradient.setColorAt(1, col0.darker(110));
	}
	painter->setBrush(gradient);

	QPen pen;
	if (Settings->getSheetName() == ":Standard_Dark") {
		pen = QPen(Qt::white, 1);
	} else {
		pen = QPen(Qt::black, 1);
	}
	if (option->state & QStyle::State_Selected) {
		pen.setWidth(2);
		pen.setColor(Qt::yellow);
	}
	painter->setPen(pen);
	painter->drawEllipse(-mNodeDrawSize2, -mNodeDrawSize2, mNodeDrawSize, mNodeDrawSize);

	// Texte
	QString txt = QString::fromUtf8(_desc_string.c_str());
	QFont font = graph->font();
	font.setBold(true);
	painter->setFont(font);

	float m = QFontMetricsF(font).height();
	float f = m / 16.0;

	// Halo complet (8 directions) pour une lisibilité parfaite
	painter->setPen(Settings->getSheetName() == ":Standard_Dark" ? Qt::black : Qt::white);
	for(int dx=-1; dx<=1; ++dx)
		for(int dy=-1; dy<=1; ++dy)
			if(dx != 0 || dy != 0)
				painter->drawText(-10 + dx, 5 * f + dy, txt);

	painter->setPen(Settings->getSheetName() == ":Standard_Dark" ? Qt::white : Qt::black);
	painter->drawText(-10, 5 * f, txt);

	if (!mDeterminedBB)
	{
		QRect textBox = painter->boundingRect(-10, 5 * f, QFontMetrics_horizontalAdvance(QFontMetricsF(font), txt), 1.5 * m, Qt::AlignVCenter, txt);
		mBBWidth = textBox.width() + 40 * f;
		mDeterminedBB = true;
	}
}

QVariant Node::itemChange(GraphicsItemChange change, const QVariant &value)
{
    switch (change) {
    case ItemPositionHasChanged:
        foreach (Edge *edge, edgeList)
            edge->adjust();
        graph->itemMoved();
        break;
    default:
        break;
    };

    return QGraphicsItem::itemChange(change, value);
}

void Node::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
	if(event->button() == Qt::LeftButton)
	{
		if (graph) graph->setSelectedNode(this);
	}
	update();
	QGraphicsItem::mousePressEvent(event);
}

void Node::peerDetails()
{
#ifdef DEBUG_ELASTIC
	std::cerr << "Calling peer details" << std::endl;
#endif
    PGPKeyDialog::showIt(_gpg_id, PGPKeyDialog::PageDetails);
}
void Node::makeFriend()
{
    PGPKeyDialog::showIt(_gpg_id, PGPKeyDialog::PageDetails);
}
void Node::denyFriend()
{
    PGPKeyDialog::showIt(_gpg_id, PGPKeyDialog::PageDetails);
}

void Node::contextMenuEvent(QGraphicsSceneContextMenuEvent *event) 
{
	QMenu contextMnu ;

	if(_type == GraphWidget::ELASTIC_NODE_TYPE_FRIEND)
        contextMnu.addAction(FilesDefs::getIconFromQtResourcePath(IMAGE_DENIED), QObject::tr( "Deny friend" ), this, SLOT(denyFriend()) );
	else if(_type != GraphWidget::ELASTIC_NODE_TYPE_OWN)
        contextMnu.addAction(FilesDefs::getIconFromQtResourcePath(IMAGE_MAKEFRIEND), QObject::tr( "Make friend" ), this, SLOT(makeFriend()) );

    contextMnu.addAction(FilesDefs::getIconFromQtResourcePath(IMAGE_MAKEFRIEND), QObject::tr( "Peer details" ), this, SLOT(peerDetails()) );

	contextMnu.exec(event->screenPos());
}

void Node::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
	if (graph) {
		graph->setSelectedNode(NULL);
		graph->forceRedraw();
	}

    update();
    QGraphicsItem::mouseReleaseEvent(event);
}

