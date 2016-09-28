#include "../include/dataPoint.h"
#include <math.h>

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
