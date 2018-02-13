#include "TRandom3.h"

#include "../include/target.h"
#include "../include/plots.h"
#include "../include/physicalConstants.h"
#include "../include/config.h"
#include "../include/CSPrereqs.h"
#include "../include/CSUtilities.h"
#include "../include/crossSection.h"
#include "../include/experiment.h"

#include <iostream>

using namespace std;

TH1D* convertTOFtoEnergy(TH1D* tof, string name)
{
    TH1D* energy = timeBinsToRKEBins(tof, name); 

    if(!tof)
    {
        cerr << "Error: cannot convert empty TOF histogram to energy units in convertTOFtoEnergy()" << endl;
        return energy;
    }

    int tofBins = tof->GetNbinsX();

    for(int j=1; j<=tofBins; j++)
    {
        // convert time into neutron velocity based on flight path distance
        double velocity = pow(10.,7.)*(config.facility.FLIGHT_DISTANCE)
            /tof->GetBinCenter(j); // in meters/sec 

        // convert velocity to relativistic kinetic energy
        double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

        energy->Fill(rKE,tof->GetBinContent(j));
    }

    // calculate energy error
    int energyBins = energy->GetNbinsX();
    for(int j=1; j<=energyBins; j++)
    {
        energy->SetBinError(j,pow(energy->GetBinContent(j),0.5));
    }

    return energy;
}

vector<double> scaleBins(vector<double> inputBins, double scaledown)
{
    vector<double> outputBins;

    for(double i=0; i<inputBins.size(); i+=scaledown)
    {
        outputBins.push_back(inputBins[((int)i)]);
    }

    // add the max edge of the input bins
    if(outputBins.back() != inputBins.back())
    {
        outputBins.push_back(inputBins.back());
    }

    return outputBins;
}

double tofToRKE(double TOF)
{
    double velocity = pow(10.,7.)*config.facility.FLIGHT_DISTANCE/TOF; // in meters/sec 

    if(velocity>C)
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

    double TOF = pow(10.,7.)*config.facility.FLIGHT_DISTANCE/velocity; // in meters/sec 

    return TOF; // in ns
}

TH1D* timeBinsToRKEBins(TH1D* inputHisto, string name)
{
    if(!inputHisto)
    {
        cerr << "Error: tried to convert time bins to energy bins for histo " << name << ", but histo pointer was null." << endl;
    }

    int nOldBins = inputHisto->GetNbinsX();
    TAxis* oldAxis = inputHisto->GetXaxis();

    double minimumTime;
    int minimumTimeBinEdgeNumber;

    for(int i=1; i<=nOldBins; i++)
    {
        minimumTime = inputHisto->GetBinLowEdge(i);
        minimumTimeBinEdgeNumber = i;

        double rke =  tofToRKE(minimumTime);
        if(rke>0 && rke<config.plot.ENERGY_UPPER_BOUND)
        {
            break;
        }
    }

    if(tofToRKE(minimumTime)==-1)
    {
        cerr << "Error: energy of old min time " << minimumTime << " was not finite." << endl;
        exit(1);
    }

    double maximumTime;
    int maximumTimeBinEdgeNumber;

    for(int i=nOldBins; i>=1; i--)
    {
        maximumTime = inputHisto->GetBinLowEdge(i) + inputHisto->GetBinWidth(i);
        maximumTimeBinEdgeNumber = i;

        if(tofToRKE(maximumTime)>config.plot.ENERGY_LOWER_BOUND)
        {
            break;
        }
    }

    if(tofToRKE(maximumTime)==-1)
    {
        cerr << "Error: energy of old maximum time " << maximumTime << " was not finite." << endl;
        exit(1);
    }

    // Remap bins from old histo to new histo
    int numberEnergyBins = maximumTimeBinEdgeNumber-minimumTimeBinEdgeNumber+1;
    vector<double> unscaledEnergyBinEdges;

    // Reorder bins to go from lowest energy (shortest time) to highest energy (longest time)
    // n points are defined for n+1 bin edges (like fence sections and fence posts)
    for(int i=0; i<numberEnergyBins-1; i++)
    {
        double newBinEdge = tofToRKE(oldAxis->GetBinLowEdge(maximumTimeBinEdgeNumber-i)+oldAxis->GetBinWidth(maximumTimeBinEdgeNumber-i));

        if(newBinEdge<=0)
        {
            cerr << "Error: tried to make negative energy bin." << endl;
            exit(1);
        }

        unscaledEnergyBinEdges.push_back(newBinEdge);
    }

    unscaledEnergyBinEdges.push_back(tofToRKE(oldAxis->GetBinLowEdge(minimumTimeBinEdgeNumber)));

    // Downscale bins to desired granularity
    double scaledown = ((double)unscaledEnergyBinEdges.size())/config.plot.NUMBER_ENERGY_BINS;

    vector<double> scaledEnergyBinEdges = scaleBins(unscaledEnergyBinEdges, scaledown);

    TH1D* outputHisto = new TH1D(name.c_str(),
            name.c_str(),
            scaledEnergyBinEdges.size()-1,
            &scaledEnergyBinEdges[0]);

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
        if(RKEToTOF(minimumEnergy)>0 && RKEToTOF(minimumEnergy)<config.plot.TOF_UPPER_BOUND)
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
        if(RKEToTOF(maximumEnergy)>config.plot.TOF_LOWER_BOUND)
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

    TOFHisto = new TH1D(tofName.c_str(),tofName.c_str(),config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);
    rawTOFHisto = new TH1D(rawTOFName.c_str(),rawTOFName.c_str(),config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);
    deadtimeHisto = new TH1D(deadtimeName.c_str(),deadtimeName.c_str(),config.plot.TOF_BINS,config.plot.TOF_LOWER_BOUND,config.plot.TOF_UPPER_BOUND);

    energyHisto = timeBinsToRKEBins(TOFHisto,energyName);
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
