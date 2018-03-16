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
#include "TGraphAsymmErrors.h"
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

vector<double> CrossSection::getStatErrors() const
{
    vector<double> statErrors;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        statErrors.push_back(dataSet.getPoint(i).getStatError());
    }
    return statErrors;
}

vector<double> CrossSection::getSysErrors() const
{
    vector<double> sysErrors;
    for(int i=0; i<dataSet.getNumberOfPoints(); i++)
    {
        sysErrors.push_back(dataSet.getPoint(i).getSysError());
    }
    return sysErrors;
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

void CrossSection::createStatErrorsGraph(string name, string title)
{
    TGraphAsymmErrors* t = new TGraphAsymmErrors(getNumberOfPoints(),
                                      &getEnergyValues()[0],
                                      &getCrossSectionValues()[0],
                                      &getEnergyErrorsL()[0],
                                      &getEnergyErrorsR()[0],
                                      &getStatErrors()[0],
                                      &getStatErrors()[0]);
    t->SetNameTitle(name.c_str(),title.c_str());
    gDirectory->WriteTObject(t);
}
void CrossSection::createSysErrorsGraph(string name, string title)
{
    TGraphAsymmErrors* t = new TGraphAsymmErrors(getNumberOfPoints(),
                                      &getEnergyValues()[0],
                                      &getCrossSectionValues()[0],
                                      &getEnergyErrorsL()[0],
                                      &getEnergyErrorsR()[0],
                                      &getSysErrors()[0],
                                      &getSysErrors()[0]);
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

    double massError =
        targetData.target.getMassUncertainty()/targetData.target.getMass();

    double molarMassError =
        targetData.target.getMolarMassUncertainty()/targetData.target.getMolarMass();

    double diameterError =
        targetData.target.getDiameterUncertainty()/targetData.target.getDiameter();

    double arealDensityError =
        pow(pow(massError,2) + pow(molarMassError,2) + pow(diameterError,2),0.5); // as percent

    //cout << "arealDensity for " << targetData.target.getName() << " = " << arealDensity << endl;
    //cout << "arealDensityError for " << targetData.target.getName() << " = " << 100*arealDensityError << "%" << endl;

    double tofSigma = 1; //calculateTOFSigma(targetData.TOFHisto);

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
        double statisticalError = 
            pow((1/tCounts+1/bCounts+1/bMon+1/tMon),0.5)/(arealDensity/pow(10,24)); // as percent

        double crossSectionError =
            pow(pow(statisticalError,2)+pow(arealDensityError,2),0.5); // as percent

        crossSectionError *= crossSectionValue; // in barns

        addDataPoint(
                DataPoint(energyValue, energyErrorL, energyErrorR, crossSectionValue, crossSectionError,
                    crossSectionValue*statisticalError, crossSectionValue*arealDensityError, bMon, tMon, bCounts, tCounts));
    }

    name = targetData.target.getName();
}

/*void CrossSection::calculateRelDiffCS(const CSPrereqs& target1Data, const CSPrereqs& target2Data, const CSPrereqs& blankData)
{
    CrossSection target1CS = calculateCS(target1Data, blankData);
    CrossSection target2CS = calculateCS(target1Data, blankData);

    double bMon = blankData.monitorCounts;
    TH1D* bEnergyHisto = blankData.energyHisto;

    // calculate number of atoms in target1
    double set1NumberOfAtoms =
        (target1Data.target.getMass()/target1Data.target.getMolarMass())*AVOGADROS_NUMBER;

    // calculate areal density (atoms/barn) in target 1
    double set1ArealDensity =
        set1NumberOfAtoms*pow(10,-24)/(pow(target1Data.target.getDiameter()/2,2)*M_PI); // area of cylinder end

    // calculate number of atoms in target2
    double set2NumberOfAtoms =
        (target1Data.target.getMass()/target1Data.target.getMolarMass())*AVOGADROS_NUMBER;

    // calculate areal density (atoms/barn) in target 2
    double set2ArealDensity =
        set2NumberOfAtoms*pow(10,-24)/(pow(target1Data.target.getDiameter()/2,2)*M_PI); // area of cylinder end

    double arealDensityCorrectionFactor =
        1/(set2ArealDensity*set1ArealDensity)
      - 1/(set1ArealDensity*set1ArealDensity)
      - 1/(set2ArealDensity*set2ArealDensity);

    // create summed cross section (blank counts and blank monitor counts are
    // the same for both targets' cross sections, so don't double-count in error
    // propagation)
    DataSet summed = target1CS.getDataSet()+target2CS.getDataSet();

    for(int i=0; i<summed.getNumberOfPoints(); i++)
    {
        double blankCounts = bEnergyHisto->GetBinContent(i+1);

        summed.getPoint(i).setYError(
                summed.getPoint(i).getYError(i)
              + arealDensityCorrectionFactor*(1/blankCounts + 1/blankMonitorCounts));

        cout << "uncorrected Y error for point " << i << " = "
            << summed.getPoint(i).getYError(i)
            << " Y error correction is "
            << arealDensityCorrectionFactor*(1/blankCounts + 1/blankMonitorCounts)
            << endl;
    }

    // calculate first term of difference
    DataSet firstTermOfD = target1CS.getDataSet()/summed;

    for(int i=0; i<firstTermOfD.getNumberOfPoints(); i++)
    {
        double blankCounts = bEnergyHisto->GetBinContent(i+1);

        firstTermOfD.getPoint(i).setYError(
                firstTermOfD.getPoint(i).getYError(i)
              + arealDensityCorrectionFactor*(1/blankCounts + 1/blankMonitorCounts));

        cout << "uncorrected Y error for point " << i << " = "
            << firstTermOfD.getPoint(i).getYError(i)
            << " Y error correction is "
            << arealDensityCorrectionFactor*(1/blankCounts + 1/blankMonitorCounts)
            << endl;
    }

    addDataSet(difference/summed);

    name = "RelDiff" + target1Data.target.getName() + target2Data.target.getName();
}*/
