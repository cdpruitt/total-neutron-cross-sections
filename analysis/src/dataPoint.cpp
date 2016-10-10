#include "../include/dataPoint.h"
#include <math.h>

DataPoint::DataPoint()
{
}

DataPoint::DataPoint(double x, double xE,
                     double y, double yE)
{
    xValue = x;
    xError = xE;
    yValue = y;
    yError = yE;
}

double DataPoint::getXValue() const
{
    return xValue;
}

double DataPoint::getXError() const
{
    return xError;
}

double DataPoint::getYValue() const
{
    return yValue;
}

double DataPoint::getYError() const
{
    return yError;
}

void DataPoint::setXValue(double xv)
{
    xValue = xv;
}

void DataPoint::setXError(double xe)
{
    xError = xe;
}

void DataPoint::setYValue(double yv)
{
    yValue = yv;
}

void DataPoint::setYError(double ye)
{
    yError = ye;
}

DataPoint operator-(const DataPoint& minuend, const DataPoint& subtrahend)
{
    DataPoint outputDataPoint(minuend.getXValue()-subtrahend.getXValue(),
                              pow(
                                  pow(minuend.getXError(),2)+
                                  pow(subtrahend.getXError(),2),0.5),
                              minuend.getYValue()-subtrahend.getYValue(),
                              pow(
                                  pow(minuend.getYError(),2)+
                                  pow(subtrahend.getYError(),2),0.5));
    return outputDataPoint;
}

DataPoint DataPoint::mergePoints(DataPoint point2)
{
    DataPoint mergedPoint;
    mergedPoint.setXValue(this->xValue+point2.getXValue());
    mergedPoint.setXError(pow(pow(this->xError,2)+pow(point2.getXError(),2),0.5));
    mergedPoint.setYValue(this->yValue+point2.getYValue());
    mergedPoint.setYError(pow(pow(this->yError,2)+pow(point2.getYError(),2),0.5));
    return mergedPoint;
}
