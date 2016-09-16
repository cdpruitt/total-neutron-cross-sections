#include <iostream>

#include "../include/target.h"
#include "../include/targetConstants.h"
#include "../include/crossSection.h"
#include "../include/DataPoint.h"
#include "../include/analysisConstants.h"
#include "../include/physicalConstants.h"
#include "../include/plottingConstants.h"
#include "../include/plots.h"

#include <iostream>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "TAxis.h"

using namespace std;

void getMonitorCounts(vector<long>& monitorCounts, TFile*& histoFile)
{
    histoFile->cd(dirs[1].c_str());

    for(int i=0; i<6; i++)
    {
        monitorCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
        cout << "target position " << i << " monitor counts = " << monitorCounts.back() << endl;
    }
}

void calculateCS(const vector<string>& targetOrder, string histoFileName, string CSFileName)
{
    // Find number of events in the monitor for each target to use in scaling
    // cross-sections
    vector<long> monitorCounts;

    TFile* histoFile = new TFile(histoFileName.c_str(),"READ");

    // to normalize flux between all channels, we'll need to keep track of the
    // number of counts that recorded by the monitor paddle for each target
    getMonitorCounts(monitorCounts, histoFile);

    histoFile->cd("/");

    vector<TH1I*> energyHistos;
    for(int i = 0; i<NUMBER_OF_TARGETS; i++)
    {
        string name = positionNames[i] + "CorrectedEnergy";
        energyHistos.push_back((TH1I*)histoFile->Get(name.c_str()));
    }

    // Loop through the relativistic kinetic energy histograms and use them
    // to populate cross-section for each target

    // calculate cross sections for each bin of each target
    double crossSectionValue;
    double crossSectionError;
    double energyValue;
    double energyError;

    vector<Target*> targets;
    for(string s : targetOrder)
    {
        targets.push_back(
                new Target(TARGET_DATA_FILE_PATH + s + TARGET_DATA_FILE_EXTENSION));
    }

    Target* blank = targets[0];
    TH1I* blankEnergy = energyHistos[0];

    TFile* CSFile = new TFile(CSFileName.c_str(),"RECREATE");

    for(int i=0; (size_t)i<targets.size(); i++)
    {
        Target* t = targets[i];
        TH1I* energy = energyHistos[i];

        CrossSection crossSection;

        int numberOfBins = blankEnergy->GetNbinsX();

        long targetMonCounts = monitorCounts[i];
        long blankMonCounts = monitorCounts[0];
        if(targetMonCounts == 0 || blankMonCounts == 0)
        {
            cout << "Error - didn't find any monitor counts for target while trying to calculate cross sections. Exiting..." << endl;
            exit(1);
        }

        for(int j=1; j<numberOfBins; j++)
        {
            energyValue = energy->GetBinCenter(j);
            energyError = 0;

            // avoid "divide by 0" and "log of 0" errors
            if(blankEnergy->GetBinContent(j) <= 0 || energy->GetBinContent(j) <= 0)
            {
                crossSectionValue = 0;
                crossSectionError = 0;
            }

            else
            {
                // calculate the cross section
                crossSectionValue =
                    -log(
                            ((double)energy->GetBinContent(j) // counts in target
                             /blankEnergy->GetBinContent(j))// counts in blank
                            *(blankMonCounts/(double)targetMonCounts) // scale by monitor counts
                        )
                    /
                        (
                            t->getMass()
                            *AVOGADROS_NUMBER
                            *pow(10.,-24) // convert cm^2 to barns 
                            /
                                (pow(t->getDiameter()/2,2)*M_PI // area of cylinder end
                                *t->getMolMass())
                        );

                // calculate the statistical error
                crossSectionError =
                    pow((1/(double)energy->GetBinContent(j) 
                        +1/(double)blankEnergy->GetBinContent(j)
                        +1/(double)blankMonCounts
                        +1/(double)targetMonCounts
                        ),0.5)
                        /(t->getMass()
                         *AVOGADROS_NUMBER
                         *pow(10.,-24) // convert cm^2 to barns 
                         *crossSectionValue // error of log(x) ~ (errorOfX)/x
                         /
                         ((pow(t->getDiameter()/2,2)*M_PI // area of cylinder end
                           *t->getMolMass())));
            }

            crossSection.addDataPoint(
                    DataPoint(energyValue, energyError, crossSectionValue, crossSectionError));

        }

        crossSection.createCSGraph();
    }

    histoFile->Close();
    CSFile->Close();
}

CrossSection::CrossSection()
{
}

void CrossSection::addDataPoint(DataPoint dataPoint)
{
    data.push_back(dataPoint);
}

DataPoint CrossSection::getDataPoint(int i)
{
    if((size_t)i>data.size())
    {
        cout <<
        "Error: tried to retrieve a non-existent cross section data point" <<
        endl;

        exit(1);
    }

    return data[i];
}

int CrossSection::getNumberOfPoints()
{
    return data.size();
}

vector<double> CrossSection::getEnergyValues()
{
    vector<double> energyValues;
    for(DataPoint d : data)
    {
        energyValues.push_back(d.getXValue());
    }
    return energyValues;
}

vector<double> CrossSection::getEnergyErrors()
{
    vector<double> energyErrors;
    for(DataPoint d : data)
    {
        energyErrors.push_back(d.getXError());
    }
    return energyErrors;
}

vector<double> CrossSection::getCrossSectionValues()
{
    vector<double> crossSectionValues;
    for(DataPoint d : data)
    {
        crossSectionValues.push_back(d.getYValue());
    }
    return crossSectionValues;
}

vector<double> CrossSection::getCrossSectionErrors()
{
    vector<double> crossSectionErrors;
    for(DataPoint d : data)
    {
        crossSectionErrors.push_back(d.getYError());
    }
    return crossSectionErrors;
}

void CrossSection::createCSGraph()
{
    TGraphErrors* t = new TGraphErrors(getNumberOfPoints(),
                                      &getEnergyValues()[0],
                                      &getCrossSectionValues()[0],
                                      &getEnergyErrors()[0],
                                      &getCrossSectionErrors()[0]);
    t->Write();
}
