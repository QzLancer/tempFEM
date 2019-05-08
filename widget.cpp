#include "widget.h"
#include <QVBoxLayout>
Widget::Widget(QWidget *parent)
    : QWidget(parent),
      customplot(new QCustomPlot(this))
{
    QVBoxLayout *layout = new QVBoxLayout(this);

//    qDebug()<<this->size();

    layout->addWidget(customplot);
    this->resize(800,600);
    customplot->xAxis->setLabel("x");
    customplot->xAxis->setRange(0, 0.09);
    customplot->xAxis->setAutoTickStep(false);
    customplot->xAxis->setTicks(false);
    customplot->yAxis->setLabel("y");
    customplot->yAxis->setRange(-0.09, 0.09);
    customplot->xAxis2->setTicks(false);
//    customplot->yAxis->setScaleRatio(customplot->xAxis, 1.0);

    customplot->yAxis->setAutoTickStep(false);
    customplot->yAxis->setAutoTickLabels(false);
    customplot->yAxis->setTicks(false);
    customplot->yAxis->grid()->setVisible(false);
    customplot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
    customplot->yAxis2->setTicks(false);
    customplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
//    layout->addWidget(customplot);


}

Widget::~Widget()
{

}
