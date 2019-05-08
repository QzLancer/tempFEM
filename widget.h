#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include "qcustomplot/qcustomplot.h"
//#include "mainwindow.h"

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();
     QCustomPlot *customplot;
};

#endif // WIDGET_H
