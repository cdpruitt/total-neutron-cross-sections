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
    xErrorL = xE;
    xErrorR = xE;
    yValue = y;
    yError = yE;
}

DataPoint::DataPoint(double x, double xEL, double xER,
                     double y, double yE)
{
    xValue = x;
    xErrorL = xEL;
    xErrorR = xER;
    yValue = y;
    yError = yE;
}

DataPoint::DataPoint(double x, double xEL, double xER,
                     double y, double yE, double statE,
                     double sysE)
{
    xValue = x;
    xErrorL = xEL;
    xErrorR = xER;
    yValue = y;
    yError = yE;

    statError = statE;
    sysError = sysE;
}

DataPoint::DataPoint(double x, double xE, double y,
        double yE, double statE, double sysE)
{
    xValue = x;
    xErrorL = xE;
    xErrorR = xE;
    yValue = y;
    yError = yE;

    statError = statE;
    sysError = sysE;
}

DataPoint::DataPoint(double x, double xE,
                     double y, double yE,
                     long bMC,
                     long tMC,
                     long bDC,
                     long tDC)
{
    xValue = x;
    xErrorL = xE;
    xErrorR = xE;
    yValue = y;
    yError = yE;

    blankMonitorCounts = bMC;
    targetMonitorCounts = tMC;
    blankDetCounts = bDC;
    targetDetCounts = tDC;
}

DataPoint::DataPoint(double x, double xE,
                     double y, double yE,
                     double statE, double sysE,
                     long bMC,
                     long tMC,
                     long bDC,
                     long tDC)
{
    xValue = x;
    xErrorL = xE;
    xErrorR = xE;
    yValue = y;
    yError = yE;

    statError = statE;
    sysError = sysE;

    blankMonitorCounts = bMC;
    targetMonitorCounts = tMC;
    blankDetCounts = bDC;
    targetDetCounts = tDC;
}

DataPoint::DataPoint(double x, double xEL, double xER,
                     double y, double yE, double statE,
                     double sysE,
                     long bMC,
                     long tMC,
                     long bDC,
                     long tDC)
{
    xValue = x;
    xErrorL = xEL;
    xErrorR = xER;
    yValue = y;
    yError = yE;
    statError = statE;
    sysError = sysE;

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
    return xErrorL;
}

double DataPoint::getXErrorL() const
{
    return xErrorL;
}

double DataPoint::getXErrorR() const
{
    return xErrorR;
}

double DataPoint::getYValue() const
{
    return yValue;
}

double DataPoint::getYError() const
{
    return yError;
}

double DataPoint::getStatError() const
{
    return statError;
}

double DataPoint::getSysError() const
{
    return sysError;
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
    xErrorL = xe;
    xErrorR = xe;
}

void DataPoint::setXErrorL(double xel)
{
    xErrorL = xel;
}

void DataPoint::setXErrorR(double xer)
{
    xErrorR = xer;
}

void DataPoint::setYValue(double yv)
{
    yValue = yv;
}

void DataPoint::setYError(double ye)
{
    yError = ye;
}

void DataPoint::setStatError(double statE)
{
    statError = statE;
}

void DataPoint::setSysError(double sysE)
{
    sysError = sysE;
}

DataPoint operator+(const DataPoint& augend, const DataPoint& addend)
{
    DataPoint outputDataPoint(
            (augend.getXValue()+addend.getXValue())/2, // x-value
            pow(
                pow(augend.getXErrorL(),2)+
                pow(addend.getXErrorL(),2),0.5), // x-error, left
            pow(
                pow(augend.getXErrorR(),2)+
                pow(addend.getXErrorR(),2),0.5), // x-error, right
            augend.getYValue()+addend.getYValue(), // y-value
            pow(
                pow(augend.getYError(),2)+
                pow(addend.getYError(),2),0.5), // y-error
            pow(
                pow(augend.getStatError(),2)+
                pow(addend.getStatError(),2),0.5), // statistical error
            pow(
                pow(augend.getSysError(),2)+
                pow(addend.getSysError(),2),0.5)); // systematic error

    return outputDataPoint;
}

DataPoint operator-(const DataPoint& minuend, const DataPoint& subtrahend)
{
    DataPoint outputDataPoint(minuend.getXValue(), // x-value
                            pow(
                              pow(minuend.getXErrorL(),2)+ // x-error, left
                              pow(subtrahend.getXErrorL(),2),0.5),
                            pow(
                              pow(minuend.getXErrorR(),2)+ // x-error, right
                              pow(subtrahend.getXErrorR(),2),0.5),
                              minuend.getYValue()-subtrahend.getYValue(),
                              pow(
                                  pow(minuend.getYError(),2)+
                                  pow(subtrahend.getYError(),2),0.5),
                              pow(
                                  pow(minuend.getStatError(),2)+
                                  pow(subtrahend.getStatError(),2),0.5), // statistical error
                              pow(
                                  pow(minuend.getSysError(),2)+
                                  pow(subtrahend.getSysError(),2),0.5)); // systematic error

    outputDataPoint.setBlankMonitorCounts(minuend.getBlankMonitorCounts());
    outputDataPoint.setTargetMonitorCounts(minuend.getTargetMonitorCounts());
    outputDataPoint.setBlankDetCounts(minuend.getBlankDetCounts());
    outputDataPoint.setTargetDetCounts(minuend.getTargetDetCounts());

    return outputDataPoint;
}

DataPoint operator+(const DataPoint& augend, const double addend)
{
    return DataPoint(augend.getXValue(),
                     augend.getXError(),
                     augend.getYValue()+addend,
                     augend.getYError(), augend.getStatError(),
                     augend.getSysError());
}

DataPoint operator*(const DataPoint& multiplicand, const double multiplier)
{
    return DataPoint(multiplicand.getXValue(), multiplicand.getXError(),
                     multiplicand.getYValue()*multiplier, multiplicand.getYError()*multiplier,
                     multiplicand.getStatError()*multiplier, multiplicand.getSysError()*multiplier);
}

DataPoint operator/(const DataPoint& dividend, const double divisor)
{
    return DataPoint(dividend.getXValue(), dividend.getXError(),
                     dividend.getYValue()/divisor, dividend.getYError()/divisor,
                     dividend.getStatError()/divisor, dividend.getSysError()/divisor);
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
                              abs((dividend.getYValue()/divisor.getYValue()))
                              *pow(pow(dividend.getYError()/dividend.getYValue(),2)+
                                   pow(divisor.getYError()/divisor.getYValue(),2),0.5), // y-error
                              abs((dividend.getYValue()/divisor.getYValue()))
                              *pow(pow(dividend.getStatError()/dividend.getYValue(),2)+
                                  pow(divisor.getStatError()/divisor.getYValue(),2),0.5), // statistical error
                              abs((dividend.getYValue()/divisor.getYValue()))
                              *pow(pow(dividend.getSysError()/dividend.getYValue(),2)+
                                  pow(divisor.getSysError()/divisor.getYValue(),2),0.5) // systematic error
                              );
    return outputDataPoint;
}
