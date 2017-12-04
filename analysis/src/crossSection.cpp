#include <iostream>

#include "../include/target.h"
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
                DataPoint(aXValue, aPoint.getXError(), 100*yValue, 100*totalError, /* convert to % */
                          aPoint.getBlankMonitorCounts()+bPoint.getBlankMonitorCounts(),
                          aPoint.getTargetMonitorCounts()+bPoint.getTargetMonitorCounts(),
                          aPoint.getBlankDetCounts()+bPoint.getBlankDetCounts(),
                          aPoint.getTargetDetCounts()+bPoint.getTargetDetCounts()));
    }

    relative.addDataSet(relativeDataSet);
    return relative;
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

CrossSection calculateCS(const CSPrereqs& targetData, const CSPrereqs& blankData)
{
    // make sure the monitor recorded counts for both the blank and the target
    // of interest so we can normalize the flux between them
    if(targetData.monitorCounts == 0 || blankData.monitorCounts == 0)
    {
        cerr << "Error - didn't find any monitor counts for target while trying to calculate cross sections. Returning empty cross section..." << endl;
        return CrossSection();
    }

    // define variables to hold cross section information
    CrossSection crossSection;
    double energyValue;
    double energyError;
    double crossSectionValue;
    double crossSectionError;

    int numberOfBins = targetData.energyHisto->GetNbinsX();

    // calculate the ratio of target/blank monitor counts (normalize
    // flux/macropulse)
    double tMon = targetData.monitorCounts;
    double bMon = blankData.monitorCounts;
    double monitorRatio = tMon/bMon;

    // calculate the ratio of target/blank good macropulse ratio (normalize
    // macropulse number)
    double tGMN = targetData.goodMacroNumber;
    double bTMN = blankData.totalMacroNumber;
    double macroNumberRatio = (targetData.goodMacroNumber/targetData.totalMacroNumber)
        /(blankData.goodMacroNumber/blankData.totalMacroNumber);
    //macroNumberRatio = 1;

    // calculate number of atoms in this target
    long double numberOfAtoms =
        (targetData.target.getMass()/targetData.target.getMolarMass())*AVOGADROS_NUMBER;

    // calculate areal density (atoms/cm^2) in target
    double arealDensity =
        numberOfAtoms/(pow(targetData.target.getDiameter()/2,2)*M_PI); // area of cylinder end

    // calculate volume density (atoms/cm^3) in target
    double volumeDensity =
        arealDensity/targetData.target.getLength();

    // save this areal density for later error propagation
    crossSection.setArealDensity(arealDensity);

    // loop through each bin in the energy histo, calculating a cross section
    // for each bin
    for(int i=1; i<=numberOfBins; i++) // skip the overflow and underflow bins
    {
        // read data from detector histograms for target and blank
        TH1D* bCounts = blankData.energyHisto;
        TH1D* tCounts = targetData.energyHisto;

        energyValue = tCounts->GetBinCenter(i);
        energyError = tCounts->GetBinWidth(i)/2;

        long bDet = bCounts->GetBinContent(i);
        long tDet = tCounts->GetBinContent(i);

        // if any essential values are 0, return an empty DataPoint
        if(bMon <=0 || tMon <=0 || bDet <=0 || tDet <=0)
        {
            crossSectionValue = 0;
            crossSectionError = 0;
            crossSection.addDataPoint(
                DataPoint(energyValue, energyError, crossSectionValue, crossSectionError,
                          bMon, tMon, bDet, tDet));
            continue;
        }

        // calculate the ratio of target/blank counts in the detector
        double detectorRatio = (double)tDet/bDet;

        crossSectionValue =
            -log(detectorRatio/(macroNumberRatio*monitorRatio))/arealDensity; // in cm^2

        crossSectionValue *= pow(10,24); // in barns 
            
        // calculate the statistical error
        crossSectionError =
            pow(1/tDet+1/bDet+1/bMon+1/tMon,0.5)/arealDensity; // in cm^2

        crossSectionError *= pow(10,24); // in barns

        crossSection.addDataPoint(
                DataPoint(energyValue, energyError, crossSectionValue, crossSectionError,
                    bMon, tMon, bDet, tDet));
    }

    string name = targetData.target.getName();
    crossSection.createCSGraph(name, name);

    //CrossSection corrected = correctForBlank(crossSection, volumeDensity, expName, "literatureData.root");
    //name += "blankCorrected";

    //corrected.createCSGraph(name, name);

    //return corrected;
    return crossSection;
}

CrossSection correctForBlank(CrossSection rawCS, double targetNumberDensity, string expName, string graphFileName)
{
    string blankDataLocation = "../" + expName + "/targetData/" + "blank.txt";
    ifstream blankDataFile(blankDataLocation.c_str());
    if(!blankDataFile.is_open())
    {
        std::cout << "Attempted to read blank data, but failed to find data in " << blankDataLocation << std::endl;
        exit(1);
    }

    vector<Target> blankComposition; // for holding multiple elements that comprise the blank

    string str;

    while(getline(blankDataFile,str))
    {
        // parse into tokens

        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        if(tokens[0]=="Name:")
        {
            blankComposition.push_back(Target());
            blankComposition.back().setName(tokens[1]);
        }

        else if(tokens[0]=="Length:")
        {
            blankComposition.back().setLength(atof(tokens[1].c_str()));
        }

        else if(tokens[0]=="Diameter:")
        {
            blankComposition.back().setDiameter(atof(tokens[1].c_str()));
        }

        else if(tokens[0]=="Mass:")
        {
            blankComposition.back().setMass(atof(tokens[1].c_str()));
        }

        else if(tokens[0]=="Molar")
        {
            blankComposition.back().setMolarMass(atof(tokens[2].c_str()));
        }

        else
        {
            cerr << "Error - couldn't parse a line in a targetData text file" << endl;
            exit(1);
        }
    }

    CrossSection correctedCS = rawCS;

    for(Target t : blankComposition)
    {
        double factor = (t.getMass()/t.getMolarMass())*AVOGADROS_NUMBER/(t.getLength()*M_PI*pow((t.getDiameter()/2),2));

        if(targetNumberDensity==0)
        {
            continue;
        }

        factor /= targetNumberDensity; // ratio of atoms of this element in blank, compared to target
        factor *= -1; // the correction should be additive, not subtractive

        string graphFileLocation = "../" + expName + "/literatureData/" + graphFileName;
        string graphFileName = t.getName() + "(n,tot)";
        correctedCS.subtractCS(graphFileLocation, graphFileName, factor);
    }

    return correctedCS;
}

int producePlots(string dataLocation, const vector<CSPrereqs>& allCSPrereqs)
{
    string outFileName = dataLocation + "/total.root";
    TFile* outFile = new TFile(outFileName.c_str(), "UPDATE");

    vector<CrossSection> crossSections;

    cout << endl << "Total statistics over all runs: " << endl << endl;

    CSPrereqs blank;
    for(auto& p : allCSPrereqs)
    {
        if(p.target.getName()=="blank")
        {
            blank = p;
            break;
        }
    }

    for(auto& p : allCSPrereqs)
    {
        cout << "all CS Prereqs name = " << p.target.getName() << endl;

        long totalCounts = 0;
        for(int i=0; i<p.energyHisto->GetNbinsX(); i++)
        {
            long tempCounts = p.energyHisto->GetBinContent(i);
            if(tempCounts < 0)
            {
                continue;
            }

            totalCounts += tempCounts;
        }

        cout << p.target.getName() << ": total events in energy histo = "
            << totalCounts << ", total monitor events = "
            << p.monitorCounts << ", good macro number = "
            << p.goodMacroNumber << ", total macro number = "
            << p.totalMacroNumber << endl;
        crossSections.push_back(calculateCS(p, blank));
        cout << "Target " << crossSections.back().getDataSet().getReference() <<
            " RMS error: " << crossSections.back().calculateRMSError() << endl << endl;

        p.energyHisto->SetDirectory(outFile);
        p.energyHisto->Write();

        string name = p.target.getName() + "TOF";
        p.TOFHisto->SetNameTitle(name.c_str(),name.c_str());
        p.TOFHisto->SetDirectory(outFile);

        p.TOFHisto->Write();
    }

    outFile->Close();

    /**************************************************************************
      Create relative cross section plots
     **************************************************************************/

    /*string relativeFileName = dataLocation + "/relative.root";
    TFile* relativeFile = new TFile(relativeFileName.c_str(), "UPDATE");

    // read which relative cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> relativeTargets = getRelativePlotNames(expName,"relativePlots.txt");

    for(pair<string,string> p : relativeTargets)
    {
        int largerTarget = -1;
        int smallerTarget = -1;
        for(int i=0; (size_t)i<config.target.TARGET_ORDER.size(); i++)
        {
            if(config.target.TARGET_ORDER[i]==p.first)
            {
                largerTarget = i;
            }
            else if(config.target.TARGET_ORDER[i]==p.second)
            {
                smallerTarget = i;
            }
        }

        if(largerTarget>=0 && smallerTarget>=0)
        {
            // found cross section plots for both individual targets,
            // so use them to make a relative cross section plot
            cout << "Producing relative cross section plot of " << config.target.TARGET_ORDER[largerTarget] << " and " << config.target.TARGET_ORDER[smallerTarget] << endl;

            //CrossSection sum = crossSections[largerTarget]+crossSections[smallerTarget];
            //cout << "sum plot " << sum.getDataSet().getReference() <<
            //        " RMS error: " << sum.calculateRMSError() << endl << endl;

            //CrossSection difference = crossSections[largerTarget]-crossSections[smallerTarget];
            //cout << "difference plot " << difference.getDataSet().getReference() <<
            //        " RMS error: " << difference.calculateRMSError() << endl << endl;

            CrossSection relative = calculateRelative(crossSections[largerTarget],crossSections[smallerTarget]);
            string relativeTitle = "#frac{#sigma_{" + allCSPrereqs[largerTarget].target.getName() +
                "}-#sigma_{" + allCSPrereqs[smallerTarget].target.getName() + 
                "}}{#sigma_{" + allCSPrereqs[largerTarget].target.getName() +
                "}+#sigma_{" + allCSPrereqs[smallerTarget].target.getName() + "}}";
            string relativeName = "relative" + allCSPrereqs[largerTarget].target.getName() + allCSPrereqs[smallerTarget].target.getName();
            relative.createCSGraph(relativeName.c_str(), relativeTitle.c_str());
            cout << "Relative plot " << relative.getDataSet().getReference() <<
                " RMS error: " << relative.calculateRMSError() << endl;
        }

        else
        {
            cout << "Failed to find cross section plot for either " << config.target.TARGET_ORDER[largerTarget] << " or " << config.target.TARGET_ORDER[smallerTarget] << endl;
        }
    }

    relativeFile->Close();
    */

    /**************************************************************************
      Create subtracted cross section plots
     **************************************************************************/
    /*
       string subtractedFileName = dataLocation + "/subtracted.root";
       TFile* subtractedFile = new TFile(subtractedFileName.c_str(), "RECREATE");

    // read which subtracted cross section plots to make from the experimental
    // directory
    vector<pair<string,string>> subtractedTargets = getSubtractedPlotNames(expName,"subtractedPlots.txt");

    for(pair<string,string> p : subtractedTargets)
    {
    int largerTarget = -1;
    int smallerTarget = -1;
    for(int i=0; (size_t)i<config.target.TARGET_ORDER.size(); i++)
    {
    if(config.target.TARGET_ORDER[i]==p.first)
    {
    largerTarget = i;
    }
    else if(config.target.TARGET_ORDER[i]==p.second)
    {
    smallerTarget = i;
    }
    }

    if(largerTarget>=0)
    {
    // found experimental cross section plot for this target
    cout << "Producing subtracted cross section plot of " << config.target.TARGET_ORDER[largerTarget] << " and " << config.target.TARGET_ORDER[smallerTarget] << endl;

    CrossSection litData = 

    CrossSection subtracted = calculateSubtracted(crossSections[largerTarget],litData);
    string subtractedName = "#sigma_{" + allCSPrereqs[largerTarget].target.getName() + "}-" +
    "#sigma_{" + allCSPrereqs[smallerTarget].target.getName() + "}";

    subtracted.createCSGraph(relativeName.c_str());
    cout << "subtracted plot " << relative.getDataSet().getReference() <<
    " RMS error: " << subtracted.calculateRMSError() << endl;
    }

    else
    {
    cout << "Failed to find cross section plot for either " << config.target.TARGET_ORDER[largerTarget] << " or " << config.target.TARGET_ORDER[smallerTarget] << endl;
    }
    }

    subtractedFile->Close();
    */

    return 0;
}
