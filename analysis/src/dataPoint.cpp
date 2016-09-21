#include "../include/dataPoint.h"

DataPoint::DataPoint(double x, double xE,
                     double y, double yE)
{
    xValue = x;
    xError = xE;
    yValue = y;
    yError = yE;
}

double DataPoint::getXValue()
{
    return xValue;
}

double DataPoint::getXError()
{
    return xError;
}

double DataPoint::getYValue()
{
    return yValue;
}

double DataPoint::getYError()
{
    return yError;
}
