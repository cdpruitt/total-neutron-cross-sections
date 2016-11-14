#include <iostream>

#include "../include/target.h"
#include "../include/targetConstants.h"
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
#include "TRandom3.h"
#include "TGraphErrors.h"

using namespace std;

const double BLANK_MON_SCALING = 1;

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

void correctForDeadtime(string histoFileName, string deadtimeFileName, string directory)
{
    TFile* deadtimeFile = new TFile(deadtimeFileName.c_str(),"READ");
    TFile* histoFile = new TFile(histoFileName.c_str(),"UPDATE");
    gDirectory->cd("/");
    gDirectory->cd(directory.c_str());

    vector<Plots*> uncorrectedPlots;
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        string name = positionNames[i];
        uncorrectedPlots.push_back(new Plots(name,histoFile,directory));
    }

    vector<Plots*> correctedPlots;
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        string name = positionNames[i] + "Corrected";
        correctedPlots.push_back(new Plots(name));
    }

    vector<Plots*> deadtimePlots;
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        string name = positionNames[i];
        deadtimePlots.push_back(new Plots(name, deadtimeFile, directory));
    }

    // extract deadtime from waveform-mode fit

    TRandom3 *randomizeBin = new TRandom3();

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        // "deadtimeFraction" records the fraction of time that the detector is dead, for
        // neutrons of a certain energy.

        vector<double> deadtimeFraction;

        //string temp;
        //temp = "deadtime" + t.getName() + "Waveform";
        //plots.waveformDeadtimes.push_back((TH1I*)deadtimeFile->Get(temp.c_str()));

        /*if(!t.getDeadtime.back())
        {
            cerr << "Error: couldn't find waveform deadtime histograms." << endl;
            exit(1);
        }*/

        TH1I* deadtimeHisto = deadtimePlots[i]->getDeadtimeHisto();
        if(!deadtimeHisto)
        {
            cout << "Couldn't find deadtimeHisto for target " << i << endl;
            break;
        }

        int deadtimeBins = deadtimeHisto->GetNbinsX();

        for(int j=0; j<deadtimeBins; j++)
        {
            deadtimeFraction.push_back(deadtimeHisto->GetBinContent(j)/(double)pow(10,3));
        }

        // create deadtime-corrected histograms

        deadtimeHisto->Write();

        //vector<vector<double>> eventsPerBinPerMicro(6,vector<double>(0));

        //const double FULL_DEADTIME = 183; // total amount of time after firing when
        // detector is at least partially dead to
        // incoming pulses (in ns)
        //const double PARTIAL_DEADTIME = 9; // amount of time after the end of
        // FULL_DEADTIME when detector is
        // becoming live again, depending on
        // amplitude (in ns)

        /*************************************************************************/
        // Perform deadtime correction
        /*************************************************************************/

        // loop through all TOF histos

        TH1I* tof = uncorrectedPlots[i]->getTOFHisto();
        //TH1I* en = uncorrectedPlots[i]->getEnergyHisto();

        TH1I* tofC = correctedPlots[i]->getTOFHisto();
        TH1I* enC = correctedPlots[i]->getEnergyHisto();

        int tofBins = tofC->GetNbinsX();

        // apply deadtime correction to TOF histos
        for(int j=0; j<tofBins; j++)
        {
            if(deadtimeFraction[j] > 0)
            {
                tofC->SetBinContent(j,(tof->GetBinContent(j)/(1-deadtimeFraction[j])));
            }

            // convert microTime into neutron velocity based on flight path distance
            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/(tofC->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)(2*TOF_BINS),TOF_RANGE/(double)(2*TOF_BINS))); // in meters/sec 

            // convert velocity to relativistic kinetic energy
            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            enC->Fill(rKE,tofC->GetBinContent(j));
            tofC->SetBinError(j,pow(tofC->GetBinContent(j),0.5));
            enC->SetBinError(j,pow(enC->GetBinContent(j),0.5));
        }

    }
    histoFile->Write();
    histoFile->Close();
}

int calculateCS(string CSFileName, CSPrereqs& targetData, CSPrereqs& blankData)
{
    // make sure the monitor recorded counts for both the blank and the target
    // of interest so we can normalize the flux between them
    if(targetData.monitorCounts == 0 || blankData.monitorCounts == 0)
    {
        cerr << "Error - didn't find any monitor counts for target while trying to calculate cross sections. Exiting..." << endl;
        return 1;
    }

    // define variables to hold cross section information
    CrossSection crossSection;
    double crossSectionValue;
    double crossSectionError;
    double energyValue;
    double energyError;

    // loop through each bin in the energy histo, calculating a cross section
    // for each bin
    int numberOfBins = targetData.energyHisto->GetNbinsX();

    for(int j=1; j<=numberOfBins; j++) // start j at 1 to skip the underflow bin
    {
        energyValue = targetData.energyHisto->GetBinCenter(j);
        energyError = targetData.energyHisto->GetBinError(j);

        // avoid "divide by 0" and "log of 0" errors
        if(blankData.energyHisto->GetBinContent(j) <= 0 || targetData.energyHisto->GetBinContent(j) <= 0)
        {
            crossSectionValue = 0;
            crossSectionError = 0;
            crossSection.addDataPoint(
                DataPoint(energyValue, energyError, crossSectionValue, crossSectionError));
            continue;
        }

        // calculate the cross section
        crossSectionValue =
            -log(
                    ((double)targetData.energyHisto->GetBinContent(j) // counts in target
                     /blankData.energyHisto->GetBinContent(j))// counts in blank
                    *((blankData.monitorCounts*BLANK_MON_SCALING)/(double)targetData.monitorCounts) // scale by monitor counts
                    )
            /
            (
             targetData.target.getMass()
             *AVOGADROS_NUMBER
             *pow(10.,-24) // convert cm^2 to barns 
             /
             (pow(targetData.target.getDiameter()/2,2)*M_PI // area of cylinder end
              *targetData.target.getMolMass())
            );

        // calculate the statistical error
        crossSectionError =
            pow((1/(double)targetData.energyHisto->GetBinContent(j) 
                        +1/(double)blankData.energyHisto->GetBinContent(j)
                        +1/(double)blankData.monitorCounts
                        +1/(double)targetData.monitorCounts
                ),0.5)
            /(targetData.target.getMass()
                    *AVOGADROS_NUMBER
                    *pow(10.,-24) // convert cm^2 to barns 
                    *crossSectionValue // error of log(x) ~ (errorOfX)/x
                    /
                    ((pow(targetData.target.getDiameter()/2,2)*M_PI // area of cylinder end
                      *targetData.target.getMolMass())));

        crossSection.addDataPoint(
                DataPoint(energyValue, energyError, crossSectionValue, crossSectionError));
    }

    // get hydrogen data
    TFile* hydrogenFile = new TFile("../tin2/literatureData/literatureData.root","READ");
    TGraphErrors* hydrogenGraph = (TGraphErrors*)hydrogenFile->Get("Hydrogen (n,tot)\r");
    if(!hydrogenGraph)
    {
        cerr << "Error: failed to find hydrogen graph in literature data." << endl;
        exit(1);
    }

    DataSet hydrogenData = DataSet();

    DataSet rawData = crossSection.getDataSet();

    for(int j=0; j<rawData.getNumberOfPoints(); j++)
    {

        hydrogenData.addPoint(
                DataPoint(rawData.getPoint(j).getXValue(),
                    rawData.getPoint(j).getXError(),
                    hydrogenGraph->Eval(rawData.getPoint(j).getXValue()),
                    0)
                ); 
    }

    crossSection.addDataSet(crossSection.getDataSet()-hydrogenData*2);

    TFile* CSFile = new TFile(CSFileName.c_str(),"UPDATE");
    crossSection.createCSGraph(targetData.target.getName());

    CSFile->Write();
    CSFile->Close();

    return 0;
}
