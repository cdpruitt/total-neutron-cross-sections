#include "../include/plottingConstants.h"
#include "../include/plots.h"
#include "../include/physicalConstants.h"
#include <iostream>


using namespace std;

vector<double> scaleBins(vector<double> inputBins, int scaledown)
{
    vector<double> outputBins;

    /*if(((int)inputBins.size())%scaledown!=0)
    {
        cerr << "Error: cannot scale down bins to non-integral bin sizes." << endl;
        exit(1);
    }*/

    for(size_t i=0; i<inputBins.size()/scaledown; i++)
    {
        outputBins.push_back(inputBins[scaledown*i]);
    }

    return outputBins;
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

double RKEToTOF(double RKE)
{
    // convert relativistic kinetic energy to velocity
    double velocity = pow(1-pow((1/((RKE/NEUTRON_MASS)+1)),2),0.5)*C;

    if(velocity<0 || velocity>C)
    {
        return -1;
    }

    double TOF = pow(10.,7.)*FLIGHT_DISTANCE/velocity; // in meters/sec 

    return TOF; // in ns
}

TH1D* timeBinsToRKEBins(TH1D* inputHisto, string name)
{
    // extract the total number of bins in the input Histo (minus the
    // overflow and underflow bins)
    int nOldBins = inputHisto->GetSize()-2;

    double minimumTime = (((TAxis*)inputHisto->GetXaxis())->GetXmin());
    int minimumBin = 0;

    for(int i=1; i<nOldBins+1; i++)
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

        maximumTime = inputHisto->GetBinLowEdge(i) + inputHisto->GetBinWidth(i);
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
    for(int i=0; i<nUnscaledEnergyBins; i++)
    {
        double newBin = tofToRKE(oldAxis->GetBinLowEdge(maximumBin-i)+oldAxis->GetBinWidth(maximumBin-i));
        if(newBin<=0)
        {
            cerr << "Error: tried to make negative energy bin." << endl;
            exit(1);
        }

        unscaledEnergyBins.push_back(newBin);
    }

    unscaledEnergyBins.push_back(tofToRKE(oldAxis->GetBinLowEdge(minimumBin)));
    
    // Downscale bins to desired granularity
    vector<double> scaledEnergyBins = scaleBins(unscaledEnergyBins, nUnscaledEnergyBins/NUMBER_ENERGY_BINS);

    TH1D* outputHisto = new TH1D(name.c_str(),
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

TH1D* RKEBinsToTimeBins(TH1D *inputHisto, string name)
{
    // extract the total number of bins in the input Histo (minus the
    // overflow and underflow bins)
    int nOldBins = inputHisto->GetSize()-2;

    double minimumEnergy = (((TAxis*)inputHisto->GetXaxis())->GetXmin());
    int minimumBin = 0;

    for(int i=0; i<nOldBins; i++)
    {
        if(RKEToTOF(minimumEnergy)>0 && RKEToTOF(minimumEnergy)<TOF_UPPER_BOUND)
        {
            break;
        }

        minimumEnergy = inputHisto->GetBinLowEdge(i);
        minimumBin = i;
    }

    double tentativeTime = RKEToTOF(minimumEnergy);
    if(tentativeTime==-1)
    {
        cerr << "Error: time of old min energy " << minimumEnergy << " was not finite: " << tentativeTime << " (ns)" << endl;
        exit(1);
    }

    //double newXMax = tentativeEnergy;

    double maximumEnergy = (((TAxis*)inputHisto->GetXaxis())->GetXmax());
    int maximumBin = nOldBins;

    for(int i=nOldBins; i>0; i--)
    {
        if(RKEToTOF(maximumEnergy)>TOF_LOWER_BOUND)
        {
            break;
        }
        maximumEnergy = inputHisto->GetBinLowEdge(i);
        maximumBin = i;
    }

    tentativeTime = RKEToTOF(maximumEnergy);
    if(tentativeTime==-1)
    {
        cerr << "Error: time of old maximum energy " << maximumEnergy << " was not finite: " << tentativeTime << " (ns)" << endl;
        exit(1);
    }

    //double newXMin = tentativeEnergy;

    TAxis* oldAxis = inputHisto->GetXaxis();

    // Remap bins from old histo to new histo
    int nUnscaledTimeBins = maximumBin-minimumBin;
    vector<double> unscaledTimeBins;

    // Reorder bins to go from lowest time (highest energy) to highest time (lowest energy)
    // n bins are defined n+1 points (like fence sections and fence posts)
    for(int i=0; i<nUnscaledTimeBins+1; i++)
    {
        double newBin = RKEToTOF(oldAxis->GetBinLowEdge(maximumBin-i));
        if(newBin<=0)
        {
            continue;
        }
        unscaledTimeBins.push_back(newBin);
    }

    // Downscale bins to desired granularity
    vector<double> scaledTimeBins = unscaledTimeBins;
    //scaleBins(unscaledTimeBins, scaledTimeBins, nUnscaledTimeBins/NUMBER_TOF_BINS);

    TH1D* outputHisto = new TH1D(name.c_str(),
            name.c_str(),
            scaledTimeBins.size()-1,
            &scaledTimeBins[0]);
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
    string rawTOFName = name + "rawTOF";
    string energyName = name + "Energy";
    string deadtimeName = name + "Deadtime";

    TOFHisto = new TH1D(tofName.c_str(),tofName.c_str(),TOF_BINS,TOF_LOWER_BOUND,TOF_UPPER_BOUND);
    rawTOFHisto = new TH1D(rawTOFName.c_str(),rawTOFName.c_str(),TOF_BINS,TOF_LOWER_BOUND,TOF_UPPER_BOUND);
    deadtimeHisto = new TH1D(deadtimeName.c_str(),deadtimeName.c_str(),TOF_BINS,TOF_LOWER_BOUND,TOF_UPPER_BOUND);

    energyHisto = timeBinsToRKEBins(TOFHisto,energyName);

    //energyHisto = new TH1D(energyName.c_str(),energyName.c_str(),ENERGY_BINS,0,ENERGY_RANGE);
    //TOFHisto = RKEBinsToTimeBins(energyHisto,tofName);
    //deadtimeHisto = RKEBinsToTimeBins(energyHisto,deadtimeName);
}

// reconnect to old plots
Plots::Plots(string name, TFile*& inputFile, string directory)
{
    string tofName = name + "TOF";
    string rawTOFName = name + "rawTOF";
    string energyName = name + "Energy";
    string deadtimeName = name + "Deadtime";

    if(!inputFile->IsOpen())
    {
        // can't find histograms - exit with error
        cout << "error - couldn't open histograms" << endl;
    }

    else
    {
        gDirectory->cd("/");
        gDirectory->GetDirectory(directory.c_str())->cd();

        TOFHisto = (TH1D*)gDirectory->Get(tofName.c_str());
        rawTOFHisto = (TH1D*)gDirectory->Get(rawTOFName.c_str());
        energyHisto = (TH1D*)gDirectory->Get(energyName.c_str());
        deadtimeHisto = (TH1D*)gDirectory->Get(deadtimeName.c_str());
    }
}

TH1D* Plots::getTOFHisto()
{
    return TOFHisto;
}

TH1D* Plots::getRawTOFHisto()
{
    return rawTOFHisto;
}

TH1D* Plots::getEnergyHisto()
{
    return energyHisto;
}

TH1D* Plots::getDeadtimeHisto()
{
    return deadtimeHisto;
}
