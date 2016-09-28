#include <iostream>

#include "../include/target.h"
#include "../include/targetConstants.h"
#include "../include/crossSection.h"
#include "../include/dataPoint.h"
#include "../include/analysisConstants.h"
#include "../include/physicalConstants.h"
#include "../include/plottingConstants.h"
#include "../include/plots.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "TH1.h"
#include "TFile.h"
#include "TAxis.h"

using namespace std;

vector<Target*> getTargetOrder(string expName, int runNumber)
{
    expName = "../" + expName + "/targetOrder.txt";
    ifstream dataFile(expName.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target order data in " << expName << std::endl;
        exit(1);
    }

    string str;
    vector<string> targetOrder;

    while(getline(dataFile,str))
    {
        // ignore comments in data file
        string delimiter = "-";
        string token = str.substr(0,str.find(delimiter));
        if(!atoi(token.c_str()))
        {
            // This line starts with a non-integer and is thus a comment; ignore
            continue;
        }

        // parse data lines into space-delimited tokens
        vector<string> tokens;
        istringstream iss(str);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        // extract run numbers from first token
        string lowRun = tokens[0].substr(0,tokens[0].find(delimiter));
        tokens[0] = tokens[0].erase(0,tokens[0].find(delimiter) + delimiter.length());

        delimiter = "\n";
        string highRun = tokens[0].substr(0,tokens[0].find(delimiter));
        
        if(atoi(lowRun.c_str()) <= runNumber && runNumber <= atoi(highRun.c_str()))
        {
            for(int i=1; (size_t)i<tokens.size(); i++)
            {
                targetOrder.push_back(tokens[i]);
            }
            break;
        }
    }

    vector<Target*> targets;

    for(string s : targetOrder)
    {
        targets.push_back(
                new Target(TARGET_DATA_FILE_PATH + s + TARGET_DATA_FILE_EXTENSION));
    }
    return targets;
}

void getMonitorCounts(vector<long>& monitorCounts, TFile*& histoFile)
{
    histoFile->cd(dirs[1].c_str());

    for(int i=0; i<6; i++)
    {
        monitorCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
        cout << "target position " << i << " monitor counts = " << monitorCounts.back() << endl;
    }
}

void calculateCS(string histoFileName, string CSFileName, int runNumber)
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

    vector<Target*> targets = getTargetOrder("tin",runNumber);

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

        crossSection.createCSGraph(t->getName());
    }

    histoFile->Close();
    CSFile->Close();
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

CrossSection::CrossSection()
{
}

void CrossSection::addDataPoint(DataPoint dataPoint)
{
    data.push_back(dataPoint);
}

DataPoint CrossSection::getDataPoint(int i) const
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

int CrossSection::getNumberOfPoints() const
{
    return data.size();
}

vector<double> CrossSection::getEnergyValues() const
{
    vector<double> energyValues;
    for(DataPoint d : data)
    {
        energyValues.push_back(d.getXValue());
    }
    return energyValues;
}

vector<double> CrossSection::getEnergyErrors() const
{
    vector<double> energyErrors;
    for(DataPoint d : data)
    {
        energyErrors.push_back(d.getXError());
    }
    return energyErrors;
}

vector<double> CrossSection::getCrossSectionValues() const
{
    vector<double> crossSectionValues;
    for(DataPoint d : data)
    {
        crossSectionValues.push_back(d.getYValue());
    }
    return crossSectionValues;
}

vector<double> CrossSection::getCrossSectionErrors() const
{
    vector<double> crossSectionErrors;
    for(DataPoint d : data)
    {
        crossSectionErrors.push_back(d.getYError());
    }
    return crossSectionErrors;
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
