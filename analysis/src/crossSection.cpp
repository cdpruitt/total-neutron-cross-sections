#include <iostream>

#include "../include/target.h"
#include "../include/targetConstants.h"
#include "../include/crossSection.h"
#include "../include/dataPoint.h"
#include "../include/dataSet.h"
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

vector<Target*> getTargetOrder(string expName, int runNumber)
{
    string targetOrderLocation = "../" + expName + "/targetOrder.txt";
    ifstream dataFile(targetOrderLocation.c_str());
    if(!dataFile.is_open())
    {
        std::cout << "Failed to find target order data in " << targetOrderLocation << std::endl;
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
        string targetDataLocation = "../" + expName + "/targetData/" + s + ".txt";
        targets.push_back(
                new Target(targetDataLocation));
    }
    return targets;
}

void getMonitorCounts(vector<long>& monitorCounts, TFile*& histoFile)
{
    histoFile->cd(dirs[1].c_str());

    for(int i=0; i<6; i++)
    {
        monitorCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
    }
}

void calculateCS(string histoFileName, string directory, string CSFileName, string expName, int runNumber)
{
    // Find number of events in the monitor for each target to use in scaling
    // cross-sections
    vector<long> monitorCounts;

    TFile* histoFile = new TFile(histoFileName.c_str(),"READ");

    // to normalize flux between all channels, we'll need to keep track of the
    // number of counts that recorded by the monitor paddle for each target
    getMonitorCounts(monitorCounts, histoFile);

    histoFile->cd("/");
    gDirectory->cd(directory.c_str());

    vector<TH1I*> energyHistos;
    for(int i = 0; i<NUMBER_OF_TARGETS; i++)
    {
        string name = positionNames[i] + "CorrectedEnergy";
        energyHistos.push_back((TH1I*)gDirectory->Get(name.c_str()));
    }

    // Loop through the relativistic kinetic energy histograms and use them
    // to populate cross-section for each target

    // calculate cross sections for each bin of each target
    double crossSectionValue;
    double crossSectionError;
    double energyValue;
    double energyError;

    vector<Target*> targets = getTargetOrder(expName,runNumber);

    TH1I* blankEnergy = energyHistos[0];

    TFile* CSFile = new TFile(CSFileName.c_str(),"RECREATE");

    for(int i=0; (size_t)i<targets.size(); i++)
    {
        Target* t = targets[i];
        TH1I* energy = energyHistos[i];

        CrossSection crossSection;

        int numberOfBins = blankEnergy->GetNbinsX();

        long targetMonCounts = monitorCounts[i];
        long blankMonCounts = BLANK_MON_SCALING*monitorCounts[0];
        if(targetMonCounts == 0 || blankMonCounts == 0)
        {
            cerr << "Error - didn't find any monitor counts for target while trying to calculate cross sections. Exiting..." << endl;
            exit(1);
        }

        for(int j=1; j<=numberOfBins; j++) // start j at 1 to skip the underflow bin
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
    CSFile->Write();
    CSFile->Close();
}
