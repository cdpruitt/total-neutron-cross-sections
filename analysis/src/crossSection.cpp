#include <iostream>

#include "../include/config.h"
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
#include "TF1.h"
#include "TFile.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;

extern Config config;

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

vector<double> CrossSection::getEnergyErrorsL() const
{
    vector<double> energyErrorsL;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        energyErrorsL.push_back(dataSet.getPoint(i).getXErrorL());
    }

    return energyErrorsL;
}

vector<double> CrossSection::getEnergyErrorsR() const
{
    vector<double> energyErrorsR;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        energyErrorsR.push_back(dataSet.getPoint(i).getXErrorR());
    }
    return energyErrorsR;
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

void CrossSection::createGraph(string name, string title)
{
    TGraphAsymmErrors* t = new TGraphAsymmErrors(getNumberOfPoints(),
                                      &getEnergyValues()[0],
                                      &getCrossSectionValues()[0],
                                      &getEnergyErrorsL()[0],
                                      &getEnergyErrorsR()[0],
                                      &getCrossSectionErrors()[0],
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

double calculateTOFSigma(TH1D* TOFHisto)
{
    // define gamma times
    const double GAMMA_TIME = pow(10,7)*config.facility.FLIGHT_DISTANCE/C;
    const double GAMMA_WINDOW_WIDTH = config.time.GAMMA_WINDOW_SIZE/2;

    // calculate width of gamma peak
    TF1* gammaPeakFit = new TF1("gammaPeakFit","gaus",
            GAMMA_TIME-GAMMA_WINDOW_WIDTH, GAMMA_TIME+GAMMA_WINDOW_WIDTH);
    TOFHisto->Fit("gammaPeakFit","BQ0", "",
            GAMMA_TIME-GAMMA_WINDOW_WIDTH, GAMMA_TIME+GAMMA_WINDOW_WIDTH);
    gammaPeakFit = TOFHisto->GetFunction("gammaPeakFit");
    double tofSigma = gammaPeakFit->GetParameter(2);

    /*cout << "Target gamma peak FWHM = "
        << 2.355*tofSigma
        << " ns (fit range: " << GAMMA_TIME-GAMMA_WINDOW_WIDTH
        << " - " << GAMMA_TIME+GAMMA_WINDOW_WIDTH << " ns)." << endl;
        */

    return tofSigma;
}

void CrossSection::calculateCS(const CSPrereqs& targetData, const CSPrereqs& blankData)
{
    // define variables to hold cross section information
    int numberOfBins = targetData.energyHisto->GetNbinsX();

    // calculate the ratio of target/blank monitor counts (normalize
    // flux/macropulse)
    double tMon = targetData.monitorCounts;
    double bMon = blankData.monitorCounts;

    double monitorRatio = tMon/bMon;

    // read data from detector histograms for target and blank
    TH1D* bEnergyHisto = blankData.energyHisto;
    TH1D* tEnergyHisto = targetData.energyHisto;

    // calculate the ratio of target/blank good macropulse ratio (normalize
    // macropulse number)
    double goodMacroRatio = targetData.goodMacroNumber/blankData.goodMacroNumber;
    double totalMacroRatio = targetData.totalMacroNumber/blankData.totalMacroNumber;

    double avgFluxRatio = monitorRatio*goodMacroRatio;

    // calculate number of atoms in this target
    double numberOfAtoms =
        (targetData.target.getMass()/targetData.target.getMolarMass())*AVOGADROS_NUMBER;

    // calculate areal density (atoms/cm^2) in target
    double arealDensity =
        numberOfAtoms/(pow(targetData.target.getDiameter()/2,2)*M_PI); // area of cylinder end

    cout << "arealDensity for " << targetData.target.getName() << " = " << arealDensity << endl;

    double tofSigma = calculateTOFSigma(targetData.TOFHisto);

    // loop through each bin in the energy histo, calculating a cross section
    // for each bin
    for(int i=1; i<=numberOfBins; i++) // skip the overflow and underflow bins
    {
        double energyValue = tEnergyHisto->GetBinCenter(i);
        double energyErrorL = calculateEnergyErrorL(tEnergyHisto->GetBinCenter(i), tofSigma);
        double energyErrorR = calculateEnergyErrorR(tEnergyHisto->GetBinCenter(i), tofSigma);

        double tCounts = tEnergyHisto->GetBinContent(i);
        double bCounts = bEnergyHisto->GetBinContent(i);

        // calculate the ratio of target/blank counts in the detector
        double detectorRatio = tCounts/bCounts;

        // if any essential values are 0, return an empty DataPoint
        if(detectorRatio <=0 || monitorRatio <=0
                || arealDensity <=0)
        {
            //addDataPoint(
            //    DataPoint(energyValue, energyErrorL, energyErrorR, 0,
            //              bMon, tMon, bCounts, tCounts));
            cerr << "Error: tried to produce a cross section point at bin " << i << ", but an cross section prerequisite was 0." << endl;
            continue;
        }

        double crossSectionValue =
            -log(detectorRatio/(monitorRatio))/arealDensity; // in cm^2

        crossSectionValue *= pow(10,24); // in barns 
            
        // calculate the statistical error
        double crossSectionError =
            pow((1/tCounts+1/bCounts+1/bMon+1/tMon),0.5)/arealDensity; // in cm^2

        crossSectionError *= pow(10,24); // in barns

        addDataPoint(
                DataPoint(energyValue, energyErrorL, energyErrorR, crossSectionValue, crossSectionError,
                    bMon, tMon, bCounts, tCounts));
    }

    name = targetData.target.getName();
}
