#include "../include/dataPoint.h"
#include <math.h>
#include <iostream>

using namespace std;

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

DataPoint operator/(const DataPoint& dividend, const DataPoint& divisor)
{
    if(divisor.getXValue()<=0 || divisor.getYValue()<=0)
    {
        cerr << "Error: cannot divide by <=0 (DataPoint division)" << endl;
        exit(1);
    }

    DataPoint outputDataPoint(dividend.getXValue(),
                              pow(
                                  pow(dividend.getXError()/dividend.getXValue(),2)+
                                  pow(divisor.getXError()/divisor.getXValue(),2),0.5),
                              dividend.getYValue()/divisor.getYValue(),
                              pow(
                                  pow(dividend.getYError()/dividend.getYValue(),2)+
                                  pow(divisor.getYError()/divisor.getYValue(),2),0.5));
    return outputDataPoint;
}

DataPoint DataPoint::mergePoints(DataPoint point2)
{
    DataPoint mergedPoint;
    if(xValue!=point2.getXValue())
    {
        cerr << "Error: tried to merge points with different xValues. Exiting..." << endl;
    }

    mergedPoint.setXValue(xValue);
    mergedPoint.setXError(pow(pow(xError,2)+pow(point2.getXError(),2),0.5));
    mergedPoint.setYValue((yValue+point2.getYValue())/2);
    mergedPoint.setYError(pow(pow(yError,2)+pow(point2.getYError(),2),0.5));

    return mergedPoint;
}
