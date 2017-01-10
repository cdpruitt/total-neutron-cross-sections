#include <iostream>

#include "../include/target.h"
#include "../include/crossSection.h"
#include "../include/CSPrereqs.h"
#include "../include/dataPoint.h"
#include "../include/dataSet.h"
#include "../include/dataStructures.h"
#include "../include/analysisConstants.h"
#include "../include/physicalConstants.h"
#include "../include/plottingConstants.h"
#include "../include/plots.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <utility>
#include "TH1.h"
#include "TFile.h"
#include "TAxis.h"
#include "TGraphErrors.h"

using namespace std;

CrossSection::CrossSection()
{
}

void CrossSection::addDataSet(DataSet dataSet)
{
    this->dataSet = dataSet;
}

void CrossSection::addDataPoint(DataPoint dataPoint)
{
    dataSet.addPoint(dataPoint);
}

DataSet CrossSection::getDataSet()
{
    return this->dataSet;
}

int CrossSection::getNumberOfPoints() const
{
    return dataSet.getNumberOfPoints();
}

DataPoint CrossSection::getDataPoint(int i) const
{
    if(i>dataSet.getNumberOfPoints())
    {
        cout <<
        "Error: tried to retrieve a non-existent cross section data point" <<
        endl;

        exit(1);
    }

    return dataSet.getPoint(i);
}


vector<double> CrossSection::getEnergyValues() const
{
    vector<double> energyValues;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        energyValues.push_back(dataSet.getPoint(i).getXValue());
    }
    return energyValues;
}

vector<double> CrossSection::getEnergyErrors() const
{
    vector<double> energyErrors;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        energyErrors.push_back(dataSet.getPoint(i).getXError());
    }
    return energyErrors;
}

vector<double> CrossSection::getCrossSectionValues() const
{
    vector<double> crossSectionValues;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        crossSectionValues.push_back(dataSet.getPoint(i).getYValue());
    }
    return crossSectionValues;
}

vector<double> CrossSection::getCrossSectionErrors() const
{
    vector<double> crossSectionErrors;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        crossSectionErrors.push_back(dataSet.getPoint(i).getYError());
    }
    return crossSectionErrors;
}

double CrossSection::getArealDensity() const
{
    return arealDensity;
}

void CrossSection::setArealDensity(double ad)
{
    arealDensity = ad;
}

CrossSection operator+(const CrossSection& augend, const CrossSection& addend)
{
    int n = augend.getNumberOfPoints();
    CrossSection outputCS;

    for(int i=0; i<n; i++)
    {
        outputCS.addDataPoint(augend.getDataPoint(i)+addend.getDataPoint(i));
    }

    return outputCS;
}

CrossSection operator-(const CrossSection& minuend, const CrossSection& subtrahend)
{
    int n = minuend.getNumberOfPoints();
    CrossSection outputCS;

    for(int i=0; i<n; i++)
    {
        outputCS.addDataPoint(minuend.getDataPoint(i)-subtrahend.getDataPoint(i));
    }

    return outputCS;
}

CrossSection operator/(const CrossSection& dividend, const CrossSection& divisor)
{
    int n = dividend.getNumberOfPoints();
    CrossSection outputCS;

    for(int i=0; i<n; i++)
    {
        outputCS.addDataPoint(dividend.getDataPoint(i)/divisor.getDataPoint(i));
    }

    return outputCS;
}

void CrossSection::createCSGraph(string name)
{
    TGraphErrors* t = new TGraphErrors(getNumberOfPoints(),
                                      &getEnergyValues()[0],
                                      &getCrossSectionValues()[0],
                                      &getEnergyErrors()[0],
                                      &getCrossSectionErrors()[0]);
    t->SetNameTitle(name.c_str(),name.c_str());
    t->Write();
}

double CrossSection::calculateRMSError()
{
    double RMSError = 0;

    if(!dataSet.getNumberOfPoints())
    {
        cerr << "Error: no data points found in in data set during RMS error calculation.\
             Returning 0 for RMS Error." << endl;
        return RMSError;
    }

    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        double pointError = dataSet.getPoint(i).getYError();
        RMSError += pow(pointError,2);
    }

    RMSError /= dataSet.getNumberOfPoints();
    RMSError = pow(RMSError,0.5);

    return RMSError;
}

double getPartialError(DataPoint aPoint, DataPoint bPoint, double aArealDensity)
{
    double aYValue = aPoint.getYValue();
    double bYValue = bPoint.getYValue();

    double aBMon = aPoint.getBlankMonitorCounts();
    double aTMon = aPoint.getTargetMonitorCounts();
    double aBDet = aPoint.getBlankDetCounts();
    double aTDet = aPoint.getTargetDetCounts();

    if(aBMon<=0 || aTMon<=0 || aBDet<=0 || aTDet<=0)
    {
        cerr << "Error: could not calculate partial error with non-finite count ratio. " << endl;
        return 0;
    }

    // first calculate ratio of monitor/detector counts for each data set
    double aCountRatio = (aBDet/aTDet)/(aBMon/aTMon);

    // calculate error of aCountRatio
    double errorACountRatio = aCountRatio * pow(1/aBDet+1/aTDet+1/aBMon+1/aTMon,0.5);

    // calculate error of data point from data set a partial derivative
    double leftExpression = errorACountRatio/(aCountRatio*aArealDensity*pow(10,-24)*(aYValue+bYValue));
    double rightExpression = 1-(aYValue/(aYValue+bYValue));

    return leftExpression*rightExpression;
}

CrossSection calculateRelative(CrossSection a, CrossSection b)
{
    CrossSection relative; 

    DataSet aData = a.getDataSet();
    DataSet bData = b.getDataSet();

    if(aData.getNumberOfPoints()!=bData.getNumberOfPoints())
    {
        cerr << "Error: can't calculate relative cross section from "
             << "data sets of different sizes. Returning empty cross section."
             << endl;
        return relative;
    }

    DataSet relativeDataSet;

    // for each point, calculate the relative cross section, including error
    for(int i=0; i<aData.getNumberOfPoints(); i++)
    {
        DataPoint aPoint = aData.getPoint(i);
        DataPoint bPoint = bData.getPoint(i);

        double aXValue = aPoint.getXValue();
        double bXValue = bPoint.getXValue();
        
        if(aXValue != bXValue)
        {
            cerr << "Error: can't calculate relative cross section from "
                 << "data points with different x values. Returning cross section."
                 << endl;
            return relative;
        }

        double aYValue = aPoint.getYValue();
        double bYValue = bPoint.getYValue();

        double yValue = (aYValue-bYValue)/(aYValue+bYValue);

        // calculate cross section error
        double aError = getPartialError(aPoint, bPoint, a.getArealDensity());
        double bError = getPartialError(bPoint, aPoint, b.getArealDensity());
        double totalError = pow(pow(aError,2)+pow(bError,2),0.5);

        relativeDataSet.addPoint(
                DataPoint(aXValue, aPoint.getXError(), yValue, totalError,
                          aPoint.getBlankMonitorCounts()+bPoint.getBlankMonitorCounts(),
                          aPoint.getTargetMonitorCounts()+bPoint.getTargetMonitorCounts(),
                          aPoint.getBlankDetCounts()+bPoint.getBlankDetCounts(),
                          aPoint.getTargetDetCounts()+bPoint.getTargetDetCounts()));
    }

    relative.addDataSet(relativeDataSet);
    return relative;
}

CrossSection CrossSection::subtractCS(string subtrahendFileName, string subtrahendGraphName, double factor)
{
    // get subtrahend graph
    TFile* subtrahendFile = new TFile(subtrahendFileName.c_str(),"READ");
    TGraphErrors* subtrahendGraph = (TGraphErrors*)subtrahendFile->Get(subtrahendGraphName.c_str());
    if(!subtrahendGraph)
    {
        cerr << "Error: failed to find " << subtrahendGraphName << " in " << subtrahendFileName << endl;
        exit(1);
    }

    DataSet subtrahendData = DataSet();
    DataSet rawCSData = this->getDataSet();

    // for each y-value of the raw data set, read the y-value of the subtrahend
    // and the y-error
    for(int i=0; i<rawCSData.getNumberOfPoints(); i++)
    {
        subtrahendData.addPoint(
                DataPoint(rawCSData.getPoint(i).getXValue(),
                    rawCSData.getPoint(i).getXError(),
                    subtrahendGraph->Eval(rawCSData.getPoint(i).getXValue()),
                    subtrahendGraph->GetErrorY(rawCSData.getPoint(i).getXValue()))); 
    }

    // perform the subtraction
    CrossSection differenceCS = CrossSection();
    differenceCS.addDataSet(rawCSData-subtrahendData*factor);

    return differenceCS;
}

CrossSection subtractCS(string rawCSFileName, string rawCSGraphName,
        string subtrahendFileName, string subtrahendGraphName,
        double factor,double divisor)
{
    // get rawCS graph
    TFile* rawCSFile = new TFile(rawCSFileName.c_str(),"UPDATE");
    TGraphErrors* rawCSGraph = (TGraphErrors*)rawCSFile->Get(rawCSGraphName.c_str());
    if(!rawCSGraph)
    {
        cerr << "Error: failed to find " << rawCSGraphName << " in " << rawCSFileName << endl;
        exit(1);
    }

    // get subtrahend graph
    TFile* subtrahendFile = new TFile(subtrahendFileName.c_str(),"READ");
    TGraphErrors* subtrahendGraph = (TGraphErrors*)subtrahendFile->Get(subtrahendGraphName.c_str());
    if(!subtrahendGraph)
    {
        cerr << "Error: failed to find " << subtrahendGraphName << " in " << subtrahendFileName << endl;
        exit(1);
    }

    DataSet subtrahendData = DataSet();
    DataSet rawCSData = DataSet(rawCSGraph, rawCSGraphName);

    // for each y-value of the raw CS graph, read the y-value of the subtrahend
    // and the y-error
    for(int i=0; i<rawCSData.getNumberOfPoints(); i++)
    {
        subtrahendData.addPoint(
                DataPoint(rawCSData.getPoint(i).getXValue(),
                    rawCSData.getPoint(i).getXError(),
                    subtrahendGraph->Eval(rawCSData.getPoint(i).getXValue()),
                    subtrahendGraph->GetErrorY(rawCSData.getPoint(i).getXValue()))); 
    }

    // perform the subtraction
    CrossSection differenceCS = CrossSection();
    differenceCS.addDataSet((rawCSData-subtrahendData*factor)/divisor);

    // create graph of difference
    rawCSFile->cd();
    string differenceName = "(" + rawCSGraphName + " - " +
        to_string(factor) + "*" + subtrahendGraphName + ")/" + to_string(divisor);
    differenceCS.createCSGraph(differenceName);

    return differenceCS;
}


