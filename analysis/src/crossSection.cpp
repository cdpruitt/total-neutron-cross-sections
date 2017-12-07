#include <iostream>

#include "../include/crossSection.h"
#include "../include/CSPrereqs.h"
#include "../include/dataPoint.h"
#include "../include/dataSet.h"
#include "../include/dataStructures.h"
#include "../include/plots.h"
#include "../include/physicalConstants.h"

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

CrossSection::CrossSection() {}

CrossSection::CrossSection(string n) : name(n) {}

void CrossSection::addDataSet(DataSet dS)
{
    dataSet = dS;
}

void CrossSection::addDataPoint(DataPoint dataPoint)
{
    dataSet.addPoint(dataPoint);
}

DataSet CrossSection::getDataSet()
{
    return dataSet;
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

CrossSection operator*(const CrossSection& cs, double factor)
{
    int n = cs.getNumberOfPoints();
    CrossSection outputCS;

    for(int i=0; i<n; i++)
    {
        outputCS.addDataPoint(cs.getDataPoint(i)*factor);
    }

    return outputCS;
}

void CrossSection::createCSGraph(string name, string title)
{
    TGraphErrors* t = new TGraphErrors(getNumberOfPoints(),
                                      &getEnergyValues()[0],
                                      &getCrossSectionValues()[0],
                                      &getEnergyErrors()[0],
                                      &getCrossSectionErrors()[0]);
    t->SetNameTitle(name.c_str(),title.c_str());
    gDirectory->WriteTObject(t);
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
        //cerr << "Error: could not calculate partial error with non-finite count ratio. " << endl;
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

void CrossSection::subtractCS(string subtrahendFileName, string subtrahendGraphName, double factor)
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
    this->addDataSet(rawCSData-subtrahendData*factor);
    subtrahendFile->Close();
}

CrossSection subtractCS(string rawCSFileName, string rawCSGraphName,
        string subtrahendFileName, string subtrahendGraphName,
        double factor, double divisor, string name)
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
    differenceCS.createCSGraph(name, name);

    return differenceCS;
}

CrossSection multiplyCS(string rawCSFileName, string rawCSGraphName,
        double factor, string name)
{
    // get rawCS graph
    TFile* rawCSFile = new TFile(rawCSFileName.c_str(),"UPDATE");
    TGraphErrors* rawCSGraph = (TGraphErrors*)rawCSFile->Get(rawCSGraphName.c_str());
    if(!rawCSGraph)
    {
        cerr << "Error: failed to find " << rawCSGraphName << " in " << rawCSFileName << endl;
        exit(1);
    }

    DataSet rawDS = DataSet(rawCSGraph,rawCSGraphName);

    // perform the multiplication
    DataSet productDS = rawDS*factor;

    CrossSection productCS = CrossSection();
    productCS.addDataSet(productDS);

    // create graph of difference
    rawCSFile->cd();
    productCS.createCSGraph(name, name);

    return productCS;
}

void produceRunningRMS(DataSet firstDS, DataSet secondDS, string name)
{
    DataSet rms;

    for(int i=0; i<firstDS.getNumberOfPoints(); i++)
    {
        double relDiff = 
            (firstDS.getPoint(i).getYValue() -
            secondDS.getPoint(i).getYValue())/
            (firstDS.getPoint(i).getYValue() +
             secondDS.getPoint(i).getYValue());

        if(i==0)
        {
            rms.addPoint(
                    DataPoint(firstDS.getPoint(i).getXValue(),
                        0,
                        pow(pow(relDiff,2),0.5),
                        0)
                    );
        }

        else
        {
            rms.addPoint(DataPoint(firstDS.getPoint(i).getXValue(),
                        0,
                        pow((pow(relDiff,2)+pow(rms.getPoint(i-1).getYValue(),2)*i)/(i+1),0.5),
                        0));
        }
    }

    cout << "Total RMS at " << rms.getPoint(rms.getNumberOfPoints()-1).getXValue()
        << " MeV = " << rms.getPoint(rms.getNumberOfPoints()-1).getYValue() << endl;

    CrossSection rmsPlot = CrossSection();
    rmsPlot.addDataSet(rms);

    string n = name + "rms";
    rmsPlot.createCSGraph(n.c_str(), n.c_str());
}

CrossSection relativeCS(string firstCSFileName, string firstCSGraphName,
        string secondCSFileName, string secondGraphName,
        string name)
{
    // get firstCS graph
    TFile* firstCSFile = new TFile(firstCSFileName.c_str(),"UPDATE");
    TGraphErrors* firstCSGraph = (TGraphErrors*)firstCSFile->Get(firstCSGraphName.c_str());
    if(!firstCSGraph)
    {
        cerr << "Error: failed to find " << firstCSGraphName << " in " << firstCSFileName << endl;
        exit(1);
    }

    // get second graph
    TFile* secondCSFile = new TFile(secondCSFileName.c_str(),"READ");
    TGraphErrors* secondGraph = (TGraphErrors*)secondCSFile->Get(secondGraphName.c_str());
    if(!secondGraph)
    {
        cerr << "Error: failed to find " << secondGraphName << " in " << secondCSFileName << endl;
        exit(1);
    }

    DataSet firstCSData = DataSet(firstCSGraph, firstCSGraphName);
    DataSet secondCSData = DataSet();

    // for each y-value of the first CS graph, read the y-value of the second
    // and the y-error
    for(int i=0; i<firstCSData.getNumberOfPoints(); i++)
    {
        secondCSData.addPoint(
                DataPoint(firstCSData.getPoint(i).getXValue(),
                    firstCSData.getPoint(i).getXError(),
                    secondGraph->Eval(firstCSData.getPoint(i).getXValue()),
                    0
                    /*secondGraph->GetErrorY(firstCSData.getPoint(i).getXValue())*/)); 
    }

    // perform the difference
    CrossSection differenceCS = CrossSection();
    differenceCS.addDataSet(firstCSData-secondCSData);

    // perform the sum
    CrossSection sumCS = CrossSection();
    sumCS.addDataSet(firstCSData+secondCSData);

    // perform the division
    CrossSection relDiffCS = CrossSection();
    relDiffCS.addDataSet(differenceCS.getDataSet()/sumCS.getDataSet());

    // create graph of relative difference
    firstCSFile->cd();
    relDiffCS.createCSGraph(name, name);

    // create running RMS plot
    CrossSection firstCS = CrossSection();
    firstCS.addDataSet(firstCSData);

    CrossSection secondCS = CrossSection();
    secondCS.addDataSet(secondCSData);

    produceRunningRMS(firstCS.getDataSet(), secondCS.getDataSet(), name);

    firstCSFile->Close();

    return relDiffCS;
}

void applyCSCorrectionFactor(string CSCorrectionFileName, string CSCorrectionGraphName, string CSToBeCorrectedFileName, string CSToBeCorrectedGraphName, string outputFileName, string outputGraphName)
{
    cout << "in applyCSCorrectionFactor" << endl;
    // get CS correction graph
    TFile* CSCorrectionFile = new TFile(CSCorrectionFileName.c_str(),"READ");
    TGraphErrors* CSCorrectionGraph = (TGraphErrors*)CSCorrectionFile->Get(CSCorrectionGraphName.c_str());
    if(!CSCorrectionGraph)
    {
        cerr << "Error: failed to find " << CSCorrectionGraphName << " in " << CSCorrectionFileName << endl;
        exit(1);
    }

    // get CSToBeCorrected graph
    TFile* CSToBeCorrectedFile = new TFile(CSToBeCorrectedFileName.c_str(),"READ");
    TGraphErrors* CSToBeCorrectedGraph = (TGraphErrors*)CSToBeCorrectedFile->Get(CSToBeCorrectedGraphName.c_str());
    if(!CSToBeCorrectedGraph)
    {
        cerr << "Error: failed to find " << CSToBeCorrectedGraphName << " in " << CSToBeCorrectedFileName << endl;
        exit(1);
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    DataSet CSToBeCorrectedData = DataSet(CSToBeCorrectedGraph, CSToBeCorrectedGraphName);
    DataSet CSCorrectionData = DataSet();

    // for each y-value of the CSToBeCorrected graph, read the y-value of the CSCorrection graph
    for(int i=0; i<CSToBeCorrectedData.getNumberOfPoints(); i++)
    {
        CSCorrectionData.addPoint(
                DataPoint(CSToBeCorrectedData.getPoint(i).getXValue(),
                    CSToBeCorrectedData.getPoint(i).getXError(),
                    CSCorrectionGraph->Eval(CSToBeCorrectedData.getPoint(i).getXValue()),
                    CSCorrectionGraph->GetErrorY(CSToBeCorrectedData.getPoint(i).getXValue())));
    }

    // perform the correction
    CrossSection correctedCS = CrossSection();
    correctedCS.addDataSet(correctCSUsingControl(CSToBeCorrectedData,CSCorrectionData));

    // create graph of correctedCS
    outputFile->cd();
    correctedCS.createCSGraph(outputGraphName, outputGraphName);

    outputFile->Close();
}

void scaledownCS(string CSToBeCorrectedFileName, string CSToBeCorrectedGraphName, int scaledown, string outputFileName, string outputGraphName)
{
    // get CSToBeCorrected graph
    TFile* CSToBeCorrectedFile = new TFile(CSToBeCorrectedFileName.c_str(),"READ");
    TGraphErrors* CSToBeCorrectedGraph = (TGraphErrors*)CSToBeCorrectedFile->Get(CSToBeCorrectedGraphName.c_str());
    if(!CSToBeCorrectedGraph)
    {
        cerr << "Error: failed to find " << CSToBeCorrectedGraphName << " in " << CSToBeCorrectedFileName << endl;
        exit(1);
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    DataSet CSToBeCorrectedData = DataSet(CSToBeCorrectedGraph, CSToBeCorrectedGraphName);
    DataSet CorrectedCSData = DataSet();

    double rebinnedXValue = 0;
    double rebinnedYValue = 0;
    double rebinnedYError = 0;
    int numberOfPoints = 0;

    // create downscaled CSToBeCorrected graph
    for(int i=0; i<CSToBeCorrectedData.getNumberOfPoints(); i++)
    {
        rebinnedXValue += CSToBeCorrectedData.getPoint(i).getXValue();
        rebinnedYValue += CSToBeCorrectedData.getPoint(i).getYValue();
        rebinnedYError += CSToBeCorrectedData.getPoint(i).getYError();
        numberOfPoints++;

        if(i%scaledown==scaledown-1)
        {
            rebinnedYError /= numberOfPoints;
            CorrectedCSData.addPoint(
                DataPoint(rebinnedXValue/numberOfPoints,
                    CSToBeCorrectedData.getPoint(i).getXError(),
                    rebinnedYValue/numberOfPoints,
                    rebinnedYError));

            rebinnedXValue = 0;
            rebinnedYValue = 0;
            numberOfPoints = 0;
        }
    }

    // perform the correction
    CrossSection correctedCS = CrossSection();
    correctedCS.addDataSet(CorrectedCSData);

    // create graph of correctedCS
    outputFile->cd();
    correctedCS.createCSGraph(outputGraphName, outputGraphName);

    outputFile->Close();
}
