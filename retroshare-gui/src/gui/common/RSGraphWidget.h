/*******************************************************************************
 * retroshare-gui/src/gui/common/RSGraphWidget.h                               *
 *                                                                             *
 * Copyright (C) 2014 RetroShare Team                                          *
 * Copyright (c) 2006-2007, crypton                                            *
 * Copyright (c) 2006, Matt Edman, Justin Hipple                               *
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

#include <map>
#include <set>
#include <vector>
#include <list>
#include <stdint.h>
#include <cmath>

#include <QApplication>
#include <QFrame>
#include <QTimer>
#include <QObject>
#include <QPainter>

#define GRAPH_BASE     2
#define SCALE_WIDTH   75
#define MINUSER_SCALE 2000
#define SCROLL_STEP   4

#define BACK_COLOR         Qt::white
#define SCALE_COLOR        Qt::black
#define GRID_COLOR         Qt::lightGray
#define BACK_COLOR_DARK    Qt::black
#define SCALE_COLOR_DARK   Qt::green
#define GRID_COLOR_DARK    Qt::darkGreen
#define RSDHT_COLOR        Qt::magenta
#define ALLDHT_COLOR       Qt::yellow

struct ZeroInitFloat
{
    ZeroInitFloat() { v=0; }
    ZeroInitFloat(float f) { v=f; }
    float operator()() const { return v ; }
    float& operator()() { return v ; }
    float v ;
};

class RSGraphSource: public QObject
{
    Q_OBJECT

public:
    RSGraphSource();
    virtual ~RSGraphSource() ;

    void start() ;
    void stop() ;
    void clear() ;
    void reset() ;

    virtual QString unitName() const { return "" ; }
    int n_values() const ;
    virtual QString displayValue(float v) const ;

    virtual void getCurrentValues(std::vector<QPointF>& vals) const ;
    virtual void getCumulatedValues(std::vector<float>& vals) const;

    virtual QString legend(int i, float v, bool show_value=true) const ;
    virtual void getDataPoints(int index, std::vector<QPointF>& pts, float filter_factor=0.0f) const ;
    virtual QString displayName(int index) const ;

    void setCollectionTimeLimit(qint64 msecs) ;
    void setCollectionTimePeriod(qint64 msecs) ;
    void setFiltering(bool b) { _filtering_enabled = b; }
    void setDigits(int d) { _digits = d ;}

protected slots:
    virtual void update() ;
    void updateIfPossible() ;

protected:
    virtual void getValues(std::map<std::string,float>& values) const = 0 ;

    qint64 getTime() const ;

    std::map<std::string, std::list<std::pair<qint64,float> > > _points ;
    std::map<std::string, ZeroInitFloat> _totals ;

    QTimer *_timer ;

    qint64 _time_limit_msecs ;
    qint64 _update_period_msecs ;
    qint64 _time_orig_msecs ;
    int _digits ;
    bool _filtering_enabled ;
};

class RSGraphWidget: public QFrame
{
    Q_OBJECT

public:
    static const uint32_t RSGRAPH_FLAGS_AUTO_SCALE_Y        = 0x0001 ;
    static const uint32_t RSGRAPH_FLAGS_LOG_SCALE_Y         = 0x0002 ;
    static const uint32_t RSGRAPH_FLAGS_ALWAYS_COLLECT      = 0x0004 ;
    static const uint32_t RSGRAPH_FLAGS_PAINT_STYLE_PLAIN   = 0x0008 ;
    static const uint32_t RSGRAPH_FLAGS_SHOW_LEGEND         = 0x0010 ;
    static const uint32_t RSGRAPH_FLAGS_PAINT_STYLE_FLAT    = 0x0020 ;
    static const uint32_t RSGRAPH_FLAGS_LEGEND_CUMULATED    = 0x0040 ;
    static const uint32_t RSGRAPH_FLAGS_PAINT_STYLE_DOTS    = 0x0080 ;
    static const uint32_t RSGRAPH_FLAGS_LEGEND_INTEGER      = 0x0100 ;
    static const uint32_t RSGRAPH_FLAGS_DARK_STYLE          = 0x0200 ;

    enum GraphStyle
    {
        SolidLine = 0,
        AreaGraph = 1
    };

    RSGraphWidget(QWidget *parent = 0);
    ~RSGraphWidget();

    void setTimerPeriod(int miliseconds) ;
    void setSource(RSGraphSource *gs) ;
    void setTimeScale(float pixels_per_second) ;

    void resetGraph();
    void setShowEntry(uint32_t entry, bool show) ;
    void setCurvesOpacity(float f) ;
    void setFiltering(bool b) ;

    uint32_t getFlags() const { return _flags ; }
    void setFlags(uint32_t flag) { _flags |= flag ; }
    void resetFlags(uint32_t flag) { _flags &= ~flag ; }

    void setOverrideMax(qreal max);

    // NOUVEAU : Permet de définir le plancher de l'échelle log (ex: 0.01 pour 10ms)
    void setLogMin(qreal min);

protected:
    void paintEvent(QPaintEvent *event);

protected slots:
    void updateIfPossible() ;
    virtual void wheelEvent(QWheelEvent *e);

private:
    int getNumPoints();

    void paintData();
    void paintTotals();
    void paintLegend();
    void paintScale1();
    void paintScale2();

    QColor getColor(const std::string &name) ;
    QString totalToStr(qreal total);
    void pointsFromData(const std::vector<QPointF>& values, QVector<QPointF> &points ) ;
    void paintLine(const QVector<QPointF>& points, QColor color, Qt::PenStyle lineStyle = Qt::SolidLine);
    void paintDots(const QVector<QPointF>& points, QColor color);
    void paintIntegral(const QVector<QPointF>& points, QColor color, qreal alpha = 1.0);

    QPainter* _painter;
    QRect _rec;
    qreal _maxValue;
    qreal _y_scale ;
    qreal _opacity ;
    qreal _graph_base;

    qreal pixelsToValue(qreal) ;
    qreal valueToPixels(qreal) ;
    int _maxPoints;

    std::set<std::string> _masked_entries ;

    qreal _time_scale ;
    qreal _time_filter ;
    float _linewidthscale ;

    uint32_t _flags ;
    QTimer *_timer ;

    RSGraphSource *_source ;
    qreal _override_max;
    qreal _log_min; // Minimum value for log scale
};
