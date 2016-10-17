#include "../include/plottingConstants.h"
#include "../include/plots.h"
#include "../include/physicalConstants.h"
#include <iostream>


using namespace std;

void scaleBins(vector<double> inputBins, vector<double>& outputBins, int scaledown)
{
    outputBins.resize((int)floor(inputBins.size()/scaledown));
    for(int i=0; (size_t)i<outputBins.size(); i++)
    {
        outputBins[i] = inputBins[i*scaledown+scaledown/2];
    }

    /*for(int i=0; i<nInputBins/scaledown; i++)
    {
        outputBins[i] /= scaledown;
    }*/
}

double tofToRKE(double TOF)
{
    double velocity = pow(10.,7.)*FLIGHT_DISTANCE/TOF; // in meters/sec 

    if (velocity>C)
    {
        return -1;
    }

    // convert velocity to relativistic kinetic energy
    double RKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV
    if(RKE<0)
    {
        return -1;
    }
    return RKE;
}

TH1I* timeBinsToRKEBins(TH1I *inputHisto, string name)
{
    // extract the total number of bins in the input Histo (minus the
    // overflow and underflow bins)
    int nOldBins = inputHisto->GetSize()-2;

    double minimumTime = (((TAxis*)inputHisto->GetXaxis())->GetXmin());
    int minimumBin = 0;

    for(int i=0; i<nOldBins; i++)
    {
        if(tofToRKE(minimumTime)>0 && tofToRKE(minimumTime)<ENERGY_UPPER_BOUND)
        {
            break;
        }

        minimumTime = inputHisto->GetBinLowEdge(i);
        minimumBin = i;
    }

    double tentativeEnergy = tofToRKE(minimumTime);
    if(tentativeEnergy==-1)
    {
        cerr << "Error: energy of old min time " << minimumTime << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    //double newXMax = tentativeEnergy;

    double maximumTime = (((TAxis*)inputHisto->GetXaxis())->GetXmax());
    int maximumBin = nOldBins;

    for(int i=nOldBins; i>0; i--)
    {
        if(tofToRKE(maximumTime)>ENERGY_LOWER_BOUND)
        {
            break;
        }
        maximumTime = inputHisto->GetBinLowEdge(i);
        maximumBin = i;
    }

    tentativeEnergy = tofToRKE(maximumTime);
    if(tentativeEnergy==-1)
    {
        cerr << "Error: energy of old maximum time " << maximumTime << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    //double newXMin = tentativeEnergy;

    TAxis* oldAxis = inputHisto->GetXaxis();

    // Remap bins from old histo to new histo
    int nUnscaledEnergyBins = maximumBin-minimumBin;
    vector<double> unscaledEnergyBins;

    // Reorder bins to go from lowest energy (shortest time) to highest energy (longest time)
    // n bins are defined n+1 points (like fence sections and fence posts)
    for(int i=0; i<nUnscaledEnergyBins+1; i++)
    {
        double newBin = tofToRKE(oldAxis->GetBinLowEdge(maximumBin-i));
        if(newBin<=0)
        {
            continue;
        }
        unscaledEnergyBins.push_back(newBin);
    }

    // Downscale bins to desired granularity
    vector<double> scaledEnergyBins;
    scaleBins(unscaledEnergyBins, scaledEnergyBins, nUnscaledEnergyBins/NUMBER_ENERGY_BINS);

    TH1I* outputHisto = new TH1I(name.c_str(),
            name.c_str(),
            scaledEnergyBins.size()-1,
            &scaledEnergyBins[0]);
            //newXMin,
            //newXMax);

    // Assign the remapped bins to the new histo
    //TH1* outputHistoNonZero = outputHisto->Rebin(scaledEnergyBins.size()-2,"outputHistoNonZero",&scaledEnergyBins[0]);

    //double test = outputHistoNonZero->GetXaxis()->GetBinLowEdge(scaledEnergyBins.size()-2);
    //double test2 = outputHistoNonZero->GetXaxis()->GetBinLowEdge(0);

    //return outputHistoNonZero;
    return outputHisto;
}

// create new plots
Plots::Plots(string name)
{
    string tofName = name + "TOF";
    string energyName = name + "Energy";
    string deadtimeName = name + "Deadtime";

    TOFHisto = new TH1I(tofName.c_str(),tofName.c_str(),TOF_BINS,0,TOF_RANGE);
    energyHisto = timeBinsToRKEBins(TOFHisto,energyName);
    deadtimeHisto = new TH1I(deadtimeName.c_str(),deadtimeName.c_str(),TOF_BINS,0,TOF_RANGE);
}

// reconnect to old plots
Plots::Plots(string name, TFile*& inputFile)
{
    string tofName = name + "TOF";
    string energyName = name + "Energy";
    string deadtimeName = name + "Deadtime";

    if(!inputFile->IsOpen())
    {
        // can't find histograms - exit with error
        cout << "error - couldn't open histograms" << endl;
    }

    else
    {
        TOFHisto = (TH1I*)inputFile->Get(tofName.c_str());
        energyHisto = (TH1I*)inputFile->Get(energyName.c_str());
        deadtimeHisto = (TH1I*)inputFile->Get(deadtimeName.c_str());
    }
}

TH1I* Plots::getTOFHisto()
{
    return TOFHisto;
}

TH1I* Plots::getEnergyHisto()
{
    return energyHisto;
}

TH1I* Plots::getDeadtimeHisto()
{
    return deadtimeHisto;
}


