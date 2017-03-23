#include "../include/dataPoint.h"
#include <cmath>
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

DataPoint::DataPoint(double x, double xE,
                     double y, double yE,
                     long bMC,
                     long tMC,
                     long bDC,
                     long tDC)
{
    xValue = x;
    xError = xE;
    yValue = y;
    yError = yE;

    blankMonitorCounts = bMC;
    targetMonitorCounts = tMC;
    blankDetCounts = bDC;
    targetDetCounts = tDC;
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

long DataPoint::getBlankMonitorCounts() const
{
    return blankMonitorCounts;
}

long DataPoint::getTargetMonitorCounts() const
{
    return targetMonitorCounts;
}

long DataPoint::getBlankDetCounts() const
{
    return blankDetCounts;
}

long DataPoint::getTargetDetCounts() const
{
    return targetDetCounts;
}

void DataPoint::setBlankMonitorCounts(long bMC)
{
    blankMonitorCounts = bMC;
}

void DataPoint::setTargetMonitorCounts(long tMC)
{
    targetMonitorCounts = tMC;
}

void DataPoint::setBlankDetCounts(long bDC)
{
    blankDetCounts = bDC;
}

void DataPoint::setTargetDetCounts(long tDC)
{
    targetDetCounts = tDC;
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

DataPoint operator+(const DataPoint& augend, const DataPoint& addend)
{
    DataPoint outputDataPoint((augend.getXValue()+addend.getXValue())/2,
                              pow(
                                  pow(augend.getXError(),2)+
                                  pow(addend.getXError(),2),0.5),
                              augend.getYValue()+addend.getYValue(),
                              pow(
                                  pow(augend.getYError(),2)+
                                  pow(addend.getYError(),2),0.5));
    return outputDataPoint;
}

DataPoint operator-(const DataPoint& minuend, const DataPoint& subtrahend)
{
    DataPoint outputDataPoint(minuend.getXValue(),
                              minuend.getXError(),
                              minuend.getYValue()-subtrahend.getYValue(),
                              pow(
                                  pow(minuend.getYError(),2)/*+
                                  pow(subtrahend.getYError(),2)*/,0.5));

    outputDataPoint.setBlankMonitorCounts(minuend.getBlankMonitorCounts());
    outputDataPoint.setTargetMonitorCounts(minuend.getTargetMonitorCounts());
    outputDataPoint.setBlankDetCounts(minuend.getBlankDetCounts());
    outputDataPoint.setTargetDetCounts(minuend.getTargetDetCounts());

    return outputDataPoint;
}

DataPoint operator*(const DataPoint& multiplicand, const double multiplier)
{
    return DataPoint(multiplicand.getXValue(), multiplicand.getXError(),
                     multiplicand.getYValue()*multiplier, multiplicand.getYError()*multiplier);
}

DataPoint operator/(const DataPoint& dividend, const double divisor)
{
    return DataPoint(dividend.getXValue(), dividend.getXError(),
                     dividend.getYValue()/divisor, dividend.getYError()/divisor);
}

DataPoint operator/(const DataPoint& dividend, const DataPoint& divisor)
{
    if(divisor.getYValue()==0)
    {
        cerr << "Error: cannot divide a DataPoint y-value by 0 (DataPoint division)" << endl;
        return DataPoint();
    }

    DataPoint outputDataPoint(dividend.getXValue(),
                              dividend.getXError(),
                              dividend.getYValue()/divisor.getYValue(),
                              pow(abs((dividend.getYValue()/divisor.getYValue()))*
                                  pow(dividend.getYError()/dividend.getYValue(),2)+
                                  pow(divisor.getYError()/divisor.getYValue(),2),0.5));
    return outputDataPoint;
}
