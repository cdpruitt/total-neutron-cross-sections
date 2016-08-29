#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TMath.h"
#include "TRandom3.h"

using namespace std;

ofstream monitorCounts;

/* Experimental constants */

const double FLIGHT_DISTANCE = 2672; // detector distance from source, in cm

// Period of micropulses
const double MICRO_PERIOD = 1788.814; // in ns

// Physical constants
const double C = 299792458; // speed of light in m/s

const double NEUTRON_MASS = 939.56536; // in MeV/c^2
const double AVOGADROS_NUMBER = 6.02214*pow(10.,23.); // in atoms/mol
const double PI = 3.14159265;

// experimentally-determined  digitizer deadtime
const int DEADTIME_PERIOD = 183;

/* Target data */

const int NUMBER_OF_TARGETS = 6;

const vector<string> positionNames = {"blank", "target1", "target2", "target3", "target4", "target5"}; 

const vector<string> targetNames = {"blank", "shortCarbon", "longCarbon", "Sn112", "NatSn", "Sn124"}; 

// physical target data, listed in order of targetNames:
const double targetLength[NUMBER_OF_TARGETS] =    {0,     1.366,  2.737,  1.369,  1.373,  1.371}; // in cm
const double targetDiameter[NUMBER_OF_TARGETS] =  {0.827, 0.826,  0.827,  0.825,  0.827,  0.827}; // in cm
const double targetMass[NUMBER_OF_TARGETS] =      {0,     1.2363, 2.4680, 4.9749, 5.3286, 5.5505}; // in grams
const double targetMolMass[NUMBER_OF_TARGETS] =   {0,     12.01,  12.01,  112,    118.7,  124}; // in g/mol

/* Plotting data*/

// number of bins in the raw energy histograms and in the cross-section
// histograms
const int NUMBER_ENERGY_BINS = 200;

const int TOF_RANGE = 1800; // in ns
const int TOF_BINS = 18000;

const double ENERGY_LOWER_BOUND = 4; // cross-section plots' lower bound, in MeV
const double ENERGY_UPPER_BOUND = 600; // cross-section plots' upper bound, in MeV

/* event variables */

// variables for holding tree event data to be histogrammed
unsigned int evtNo, macroNo, targetPos;
unsigned int sgQ, lgQ;
double completeTime, macroTime;

vector<int> *waveform; // for holding one event's waveform data


/* ROOT and organizational variables */

// ROOT file directory structure 
const string dirs[4] = {"targetChanger","monitor","detS","scavenger"};

vector<TTree*> orchard; // holds DPP channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

vector<TTree*> orchardW; // holds waveform-mode channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

TH1I* waveformH; // holds the current histogram being filled with waveform data

// a small subset of events occur in both channels 4 and 6 in the microTime
// window of ~275 ns to ~325 ns after a micropulse start. We'd like to plot
// these on the same axis to make sure they look the same.
TH1I* waveformCh4;
TH1I* waveformCh6;

TDirectory *waveformsDir;

struct Plots
{
    vector<TH1I*> TOFHistos;
    vector<TH1I*> TOFHistosCorrected;
    vector<TH1I*> TOFHistosFirstInMicro;

    vector<TH1I*> energyHistos;
    vector<TH1I*> energyHistosCorrected;
    vector<TH1I*> energyHistosUngated;
    vector<TH1I*> energyHistosNoGamma;

    vector<TH1I*> gateRatios;
    vector<TH1I*> waveformDeadtimes;
    vector<TH1I*> deadtimeHistos;
    vector<TH1I*> orderInMicro;

    vector<TGraph*> CSGraphs;
    vector<TGraph*> CSGraphsScaledToLit;

} plots;

struct Output
{
    // declare vectors to hold the cross sections of each target
    // i = target position (0 to 5), j = bin (0 to numberOfBins)
    vector<vector<double>*> crossSections;
    vector<vector<double>*> crossSectionsError;
    vector<vector<double>*> crossSectionsScaledToLit;
    vector<vector<double>*> crossSectionsScaledToLitError;

    // declare vector to hold the energy bins for cross sections
    vector<double> energyBins;
} output;


// Re-link to an already-existing tree's data so we can read the tree
void setBranches(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("macroTime",&macroTime);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("targetPos",&targetPos);
    tree->SetBranchAddress("sgQ",&sgQ);
    tree->SetBranchAddress("lgQ",&lgQ);
    tree->SetBranchAddress("waveform",&waveform);
}

// Re-link to an already-existing tree's data so we can read the tree
void setBranchesW(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("targetPos",&targetPos);
    tree->SetBranchAddress("waveform",&waveform);
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

void scaleBins(vector<double> inputBins, int nInputBins, int scaledown, vector<double>& outputBins)
{
    outputBins.resize((int)floor(nInputBins/scaledown));
    for(int i=0; i<outputBins.size(); i++)
    {
        outputBins[i] = inputBins[i*scaledown];
    }

    /*for(int i=0; i<nInputBins/scaledown; i++)
    {
        outputBins[i] /= scaledown;
    }*/
}

// Map an input histogram with bins in the time domain to equivalent bins in the relativistic
// kinetic energy domain
TH1* timeBinsToRKEBins(TH1 *inputHisto, string name)
{
    TAxis* oldAxis = inputHisto->GetXaxis();
    int nOldBins = oldAxis->GetNbins();

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
        cout << "Error: energy of old min time " << minimumTime << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    double newXMax = tentativeEnergy;

    double maximumTime = (((TAxis*)inputHisto->GetXaxis())->GetXmax());
    int maximumBin = nOldBins;

    for(int i=nOldBins; i>0; i--)
    {
        if(tofToRKE(maximumTime)>ENERGY_LOWER_BOUND)
        {
            break;
        }
        maximumTime = inputHisto->GetBinLowEdge(i)+inputHisto->GetBinWidth(i);
        maximumBin = i;
    }

    tentativeEnergy = tofToRKE(maximumTime);
    if(tentativeEnergy==-1)
    {
        cout << "Error: energy of old maximum time " << maximumTime << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    double newXMin = tentativeEnergy;

    // Remap bins from old histo to new histo
    int nUnscaledEnergyBins = maximumBin-minimumBin;
    vector<double> unscaledEnergyBins(nUnscaledEnergyBins);

    // Reorder bins to go from lowest energy (shortest time) to highest energy (longest time)
    for(int i=0; i<nUnscaledEnergyBins; i++)
    {
        unscaledEnergyBins[i] = tofToRKE(oldAxis->GetBinLowEdge(maximumBin-(i+1))+oldAxis->GetBinWidth(maximumBin-(i+1)));
    }

    // Downscale bins to desired granularity
    vector<double> scaledEnergyBins;
    scaleBins(unscaledEnergyBins, nUnscaledEnergyBins, nUnscaledEnergyBins/NUMBER_ENERGY_BINS, scaledEnergyBins);
    // n bins are defined n+1 points (like fence sections and fence posts)
    scaledEnergyBins.push_back(newXMax);

    TH1* outputHisto = new TH1D(name.c_str(),
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

// Populate advanced histograms (TOFs, cross-section, etc) calculated using
// data from ch4 and ch6 trees
void fillAdvancedHistos(int detIndex, string waveformFileName, TFile *histoFile)
{
    /*************************************************************************/
    // Prepare histograms

    // navigate to the correct directory for channel 2*i 
    gDirectory->cd("/");
    gDirectory->GetDirectory(dirs[detIndex].c_str())->cd();

    // Initialize histograms for this channel (to be filled by events after they
    // pass through various energy, time, and charge filters below

    // diagnostic histograms
    TH1I *TOF = new TH1I("TOF","Summed-detector time of flight",TOF_BINS,0,TOF_RANGE);
    TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",TOF_RANGE,0,MICRO_PERIOD+1,2048,0,65536);
    TH2I *triangleRKE = new TH2I("triangleRKE","Pulse integral vs. relativistic KE",NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,2048,0,65536);
    TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
    TH1I *QRatio = new TH1I("QRatio","short gate Q/long gate Q",1000,0,1);
    TH2I *rKElgQ = new TH2I("lgQrKE","relativistic KE vs. long gate Q",500,ENERGY_LOWER_BOUND,100,2048,0,65536);
    TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
    microNoH->GetXaxis()->SetTitle("micropulse number of each event");

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        string temp = positionNames[i] + "TOF";
        plots.TOFHistos.push_back(new TH1I(temp.c_str(),temp.c_str(),TOF_BINS,0,TOF_RANGE));

        temp = positionNames[i] + "Energy";
        plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(plots.TOFHistos[i],temp.c_str())); 

        temp = positionNames[i] + "EnergyUngated";
        plots.energyHistosUngated.push_back((TH1I*)plots.energyHistos[i]->Clone(temp.c_str()));

        temp = positionNames[i] + "EnergyNoGamma";
        plots.energyHistosNoGamma.push_back((TH1I*)plots.energyHistos[i]->Clone(temp.c_str()));

        temp = positionNames[i] + "FirstInMicro";
        plots.TOFHistosFirstInMicro.push_back((TH1I*)plots.TOFHistos[i]->Clone(temp.c_str()));
    }

    //double test = plots.energyHistos[0]->GetXaxis()->GetBinLowEdge(plots.energyHistos[0]->GetXaxis()->GetNbins());
    //double test2 = plots.energyHistosCorrected[0]->GetXaxis()->GetBinLowEdge(plots.energyHistosCorrected[0]->GetXaxis()->GetNbins());

    // create neutron energy plots using only micropulses with no gammas
    
    // create TOF plots for events that come first, second, or third in their
    // micro
    for(int i=0; i<3; i++)
    {
        stringstream temp;
        temp << "order" << i << "InMicro";
        plots.orderInMicro.push_back(new TH1I(temp.str().c_str(),temp.str().c_str(),TOF_BINS,0,TOF_RANGE));
    }

        /*************************************************************************/
    // Prepare variables used to fill histograms (TOF, order of event in micro,
    // etc.)

    // create TIME VARIABLES used for filling histograms
    double microTime;
    int microNo, prevMicroNo;
    long totalMicros = 0;
    vector<long> microsPerTarget(6,0);

    // create variables for DESCRIBING MICROPULSE TYPE (i.e., whether a
    // micropulse has a gamma at the start and keep track of the influence of
    // early micropulse events on later ones)
    int orderInMicro = 0;
    bool isGamma = false;
    bool gammaInMicro = false;
    /*************************************************************************/

    /*************************************************************************/
    // Prepare GAMMA GATE for filtering events depending on gamma presence

    // channel-dependent time offset relative to the target changer's macropulse
    // start time
    int gammaGate[2];

    // adjust time parameters based on channel identity
    switch(detIndex)
    {
        case 1:
            // monitor
            gammaGate[0] = 25;
            gammaGate[1] = 40;
            break;
        case 2:
            // summed detector
            gammaGate[0] = 85;
            gammaGate[1] = 95;
            break;
        case 3:
            // scavenger
            gammaGate[0] = 85;
            gammaGate[1] = 95;
            break;
    }
    /*************************************************************************/

    /*************************************************************************/
    // Loop through sorted trees to calculate advanced histogram variables

    // point at correct tree in preparation for reading data
    setBranches(orchard[detIndex]);

    int totalEntries = orchard[detIndex]->GetEntries();

    // reduce # of total entries for testing purposes
    //totalEntries /= 2;

    cout << "Populating advanced histograms for channel " << 2*detIndex << endl;

    // MAIN LOOP for sorting through channel-specific events
    for(int j=0; j<totalEntries; j++)
    {
        orchard[detIndex]->GetEntry(j);

        // calculate time since start of macro (includes time offsets)
        double timeDiff = completeTime-macroTime;

        // Apply gates:
        if (timeDiff < 650000 && timeDiff > 0           // require events to come during the macropulse's beam-on period
                && targetPos != 0                              // discard events during target-changer movement
                && lgQ<65500 && sgQ<32750                      // discard events with unphysically-large integrated charges
                && lgQ>sgQ                                     // discard events with short gate charge is larger than long gate charge
                //&& (sgQ/(double)lgQ<0.25 || sgQ/(double)lgQ>0.35)  // discard events outside accepted range of sgQ/lgQ
                //&& lgQ>100                                     // discard gammas at lowest range of energy
           )
        {
            /*****************************************************************/
            // Calculate event properties

            // find which micropulse the event is in and the time since
            // the start of the micropulse

            // first, save previous event's micropulse number (we'll need this
            // to calculate event ordering in each micropulse)
            prevMicroNo = microNo;

            microNo = floor(timeDiff/MICRO_PERIOD);
            microTime = fmod(timeDiff,MICRO_PERIOD);

            // convert microTime into neutron velocity based on flight path distance
            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/microTime; // in meters/sec 

            // convert velocity to relativistic kinetic energy
            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            // check to see whether this event is a gamma
            if (microTime>gammaGate[0] && microTime<gammaGate[1])
            {
                // this event IS a gamma => indicate as such
                isGamma = true;
                gammaInMicro = true;
            }

            else
            {
                isGamma = false;
            }
            /*****************************************************************/

            // GATE: discard events with too low of an integrated charge for their energy
            //if (lgQ>50)
            //if (lgQ>500*exp(-(microTime-100)/87))
            //if (lgQ<30*rKE)
            {
                // tag this event by its order in the micropulse
                if (microNo==prevMicroNo)
                {
                    // still in same micropulse => increment order counter
                    orderInMicro++;
                }

                else
                {
                    // new micropulse => return order counter to 1
                    orderInMicro = 1;

                    // new micropulse => return gamma indicator to false
                    gammaInMicro = false;
                    totalMicros++;
                }

                /*****************************************************************/
                // Fill troubleshooting plots with event variables (rKE, microtime, etc.)
                TOF->Fill(microTime);
                triangle->Fill(microTime,lgQ);
                sgQlgQ->Fill(sgQ,lgQ);
                QRatio->Fill(sgQ/(double)lgQ);
                rKElgQ->Fill(rKE,lgQ);
                triangleRKE->Fill(microTime,rKE);
                microNoH->Fill(microNo);

                // populate only first three orderInMicro plots
                if(orderInMicro<4)
                {
                    plots.orderInMicro[orderInMicro-1]->Fill(microTime);
                }
                
                /*****************************************************************/

                /*****************************************************************/
                // Fill target-specific plots

                if (targetPos>0 && targetPos<=NUMBER_OF_TARGETS)
                {
                    plots.TOFHistos[targetPos-1]->Fill(microTime);
                    plots.energyHistos[targetPos-1]->Fill(rKE);

                    if(!gammaInMicro)
                    {
                        plots.energyHistosNoGamma[targetPos-1]->Fill(rKE);
                    }

                    if(orderInMicro==1)
                    {
                        plots.TOFHistosFirstInMicro[targetPos-1]->Fill(microTime);
                        microsPerTarget[targetPos-1]++;
                    }
                }

                // end of energy, time, cross-section gates on events
                /*****************************************************************/
            }

            // fill ungated plots
            if (targetPos>0 && targetPos<=NUMBER_OF_TARGETS)
            {
                plots.energyHistosUngated[targetPos-1]->Fill(rKE);
            }

            // end of main event loop
            /*****************************************************************/

            /*if(macroNo>0)
              {
              break;
              }*/

            if(j%1000==0)
            {
                cout << "Processed " << j << " events...\r";
                fflush(stdout);
            }
        }
    }

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        string temp = positionNames[i] + "GateRatio";
        plots.gateRatios.push_back((TH1I*)plots.energyHistos[i]->Clone(temp.c_str()));
        plots.gateRatios[i]->Divide(plots.energyHistosUngated[i]);
    }

    // Find number of macropulses for each target to use in error calculation
    gDirectory->cd("/");
    gDirectory->GetDirectory("targetChanger")->cd();

    vector<long> tarCounts;

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
    }
    
    /*
    // Scale perMacroHistos histogram by number of macropulses in that target
    for(int i=0; i<6; i++)
    {
    if(tarCounts[i]==0)
    {
    continue;
    }
    perMacroHistos[i]->Scale(pow(10,3)/tarCounts[i]);
    }
    */

       /*************************************************************************/
}

void calculateDeadtime(string waveformFileName, bool skippedHistoFilling)
{
    // save reference to the histo file for later use
    TFile* histoFile = gDirectory->GetFile();

    // extract deadtime from waveform-mode fit

    // "deadtimeFraction" records the fraction of time that the detector is dead, for
    // neutrons of a certain energy. It is target-dependent.
    vector<vector<double>> deadtimeFraction(6, vector<double>(0));

    TFile *waveformFile = new TFile(waveformFileName.c_str(),"READ");

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        string temp;
        temp = "deadtime" + targetNames[i] + "Waveform";
        plots.waveformDeadtimes.push_back((TH1I*)waveformFile->Get(temp.c_str()));
        if(!plots.waveformDeadtimes.back())
        {
            cout << "Error: couldn't find waveform deadtime histograms." << endl;
            exit(1);
        }

        for(int j=0; j<plots.waveformDeadtimes.back()->GetNbinsX(); j++)
        {
            deadtimeFraction[i].push_back(plots.waveformDeadtimes.back()->GetBinContent(j)/(double)pow(10,6));
        }
    }

    waveformFile->Close();

    // create deadtime-corrected histograms
    histoFile->cd();
    gDirectory->cd(dirs[2].c_str());

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        // test to see if plot vectors have already been linked to existing plots (i.e.,
        // fillHistos and fillAdvancedHistos were run).
        if(skippedHistoFilling)
        {
            string temp = positionNames[i] + "TOF";
            plots.TOFHistos.push_back((TH1I*)gDirectory->Get(temp.c_str()));

            temp = positionNames[i] + "Energy";
            plots.energyHistos.push_back((TH1I*)gDirectory->Get(temp.c_str()));
        }

        // prepare deadtime-corrected plots
        string temp = positionNames[i] + "TOFCorrected";
        plots.TOFHistosCorrected.push_back((TH1I*)plots.TOFHistos[i]->Clone(temp.c_str()));

        temp = positionNames[i] + "EnergyCorrected";
        plots.energyHistosCorrected.push_back((TH1I*)plots.energyHistos[i]->Clone(temp.c_str()));
        plots.energyHistosCorrected.back()->Reset();
    }

    //vector<vector<double>> eventsPerBinPerMicro(6,vector<double>(0));

    //const double FULL_DEADTIME = 183; // total amount of time after firing when
    // detector is at least partially dead to
    // incoming pulses (in ns)
    //const double PARTIAL_DEADTIME = 9; // amount of time after the end of
    // FULL_DEADTIME when detector is
    // becoming live again, depending on
    // amplitude (in ns)

    TRandom3 *randomizeBin = new TRandom3();

    /*************************************************************************/
    // Perform deadtime correction
    /*************************************************************************/
    
    // loop through all TOF histos
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
/*        cout << "microsPerTarget[i] = " << microsPerTarget[i] << endl;

        //plots.energyHistosCorrected[i]->Sumw2();

        // loop through all bins
        for(int j=0; j<TOFCorrectedHistos[i]->GetNbinsX(); j++)
        {
            if(microsPerTarget[i] > 0)
            {
                eventsPerBinPerMicro[i].push_back(TOFCorrectedHistos[i]->GetBinContent(j)/(double)microsPerTarget[i]);
            }

            else
            {
                eventsPerBinPerMicro[i].push_back(0);
            }

            deadtimeFraction[i].push_back(0);
        }

        // find the fraction of the time that the detector is dead for each bin in the micropulse
        for(int j=0; (size_t)j<eventsPerBinPerMicro[i].size(); j++)
        {
            *//*int k = j-(FULL_DEADTIME+PARTIAL_DEADTIME)*(TOF_BINS/TOF_RANGE);
            while(k<j)
            {
                if(k>=0)
                {
                    // partially-dead region
                    if((j-k)>=(FULL_DEADTIME)*(TOF_BINS/TOF_RANGE))
                    {
                        deadtimeFraction[i][j] +=
                            eventsPerBinPerMicro[i][k]*
                            (k+(FULL_DEADTIME+PARTIAL_DEADTIME)*(TOF_BINS/TOF_RANGE)-j)/(PARTIAL_DEADTIME*(TOF_BINS/TOF_RANGE));
                    }

                    else
                    {
                        deadtimeFraction[i][j] += eventsPerBinPerMicro[i][k];
                    } 
                }
                k++;
            }*/
            
            /*deadtimeFraction[i][j] = deadtimeFraction[i][j-1]+eventsPerBinPerMicro[i][j];
            if(j>(TOF_BINS/TOF_RANGE)*DEADTIME_PERIOD)
            {
                deadtimeFraction[i][j] -= eventsPerBinPerMicro[i][j-(TOF_BINS/TOF_RANGE)*DEADTIME_PERIOD];
            }
        }
        */

        // apply deadtime correction to TOF histos
        for(int j=0; j<plots.TOFHistosCorrected[i]->GetNbinsX(); j++)
        {
            if(deadtimeFraction[i][j] > 0)
            {
                plots.TOFHistosCorrected[i]->SetBinContent(j,(plots.TOFHistosCorrected[i]->GetBinContent(j)/(1-deadtimeFraction[i][j])));
            }

            // convert microTime into neutron velocity based on flight path distance
            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/(plots.TOFHistosCorrected[i]->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)(2*TOF_BINS),TOF_RANGE/(double)(2*TOF_BINS))); // in meters/sec 

            // convert velocity to relativistic kinetic energy
            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            plots.energyHistosCorrected[i]->Fill(rKE,plots.TOFHistosCorrected[i]->GetBinContent(j));
            plots.TOFHistosCorrected[i]->SetBinError(j,pow(plots.TOFHistosCorrected[i]->GetBinContent(j),0.5));
            plots.energyHistosCorrected[i]->SetBinError(j,pow(plots.energyHistosCorrected[i]->GetBinContent(j),0.5));

            /*if(j%100==0 && i==3)
              {
              cout << "bin spacing randomized = " << (plots.TOFHistosCorrected[i]->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)TOF_BINS,TOF_RANGE/(double)TOF_BINS))
              -(plots.TOFHistosCorrected[i]->GetBinCenter(j-1)+randomizeBin->Uniform(-TOF_RANGE/(double)TOF_BINS,TOF_RANGE/(double)TOF_BINS))
              << endl;
              cout << "bin spacing normal = " << (plots.TOFHistosCorrected[i]->GetBinCenter(j))
              -(plots.TOFHistosCorrected[i]->GetBinCenter(j-1))
              << endl;
              }*/
        }

        plots.TOFHistosCorrected[i]->Write();
        plots.energyHistosCorrected[i]->Write();

        /*for(int j=0; j<plots.energyHistosCorrected[i]->GetNbinsX(); j++)
          {
          csCorrection[i].push_back(0);
          csCorrection[i][j]

          }*/

        //cout << deadtimeCorrection[i][150] << endl;

        // Plot the calculated dead time fractions (for debugging)
        string temp;
        temp = "deadtime" + targetNames[i];
        plots.deadtimeHistos.push_back(new TH1I(temp.c_str(),temp.c_str(),TOF_BINS,0,TOF_RANGE));
        for(int j=0; (size_t)j<deadtimeFraction[i].size(); j++)
        {
            plots.deadtimeHistos.back()->SetBinContent(j,1000000*deadtimeFraction[i][j]);
        }
        plots.deadtimeHistos.back()->Write();
    }
}

void calculateCS(vector<int> targetOrder, TFile* histoFile)
{
    // Find number of events in the monitor for each target to use in scaling
    // cross-sections

    // switch to the monitor directory
    //histoFile->cd();
    histoFile->cd("/");
    histoFile->cd(dirs[1].c_str());

    // to normalize flux between all channels, we'll need to keep track of the
    // number of counts that recorded by the monitor paddle for each target
    vector<long> monCounts;

    for(int i=2; i<8; i++)
    {
        monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i));
        cout << "target position " << i << " monitor counts = " << monCounts[i-2] << endl;
        monitorCounts << monCounts[i-2] << endl;
    }

    // switch to the detector directory
    histoFile->cd("/");
    histoFile->cd(dirs[2].c_str());

    // Loop through the relativistic kinetic energy histograms and use them
    // to populate cross-section for each target

    // Link up to the energy histograms if they haven't already been accessed
    if(!plots.energyHistos[0] || !plots.energyHistosCorrected[0]
    || !monCounts[0])
    {
        cout << "Error: failed to find necessary histos for cross section"
            << "calculation in calculateCS(). Exiting..." << endl;
        histoFile->Close();
        exit(1);
    }

    int numberOfBins = ((TAxis*)plots.energyHistos[0]->GetXaxis())->GetNbins();

    // initialize vectors for storing cross section data
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        output.crossSections.push_back(new vector<double>);
        output.crossSectionsError.push_back(new vector<double>);
        output.crossSectionsScaledToLit.push_back(new vector<double>);
        output.crossSectionsScaledToLitError.push_back(new vector<double>);
    }

    // calculate cross sections for each target (i) and each bin (j)
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        for(int j=1; j<numberOfBins; j++)
        {
            // avoid "divide by 0" and "log of 0" errors
            if(plots.energyHistosCorrected[0]->GetBinContent(j) <= 0 || plots.energyHistosCorrected[i]->GetBinContent(j) <= 0)
            {
                output.crossSections[targetOrder[i]]->push_back(0);
                output.crossSectionsError[targetOrder[i]]->push_back(0);
            }

            else
            {
                // calculate the cross section
                output.crossSections[targetOrder[i]]->push_back(
                        -log(
                            ((double)plots.energyHistosCorrected[i]->GetBinContent(j) // counts in target
                                    /plots.energyHistosCorrected[0]->GetBinContent(j))// counts in blank
                            *(monCounts[0]/(double)monCounts[i]) // scale by monitor counts
                            )
                        /(targetMass[targetOrder[i]]
                            *AVOGADROS_NUMBER
                            *pow(10.,-24) // convert cm^2 to barns 
                            /
                            ((pow(targetDiameter[targetOrder[i]]/2,2)*PI // area of cylinder end
                              *(double)targetMolMass[targetOrder[i]])))
                        );

                // calculate the statistical error
                output.crossSectionsError[targetOrder[i]]->push_back(
                    pow((1/(double)plots.energyHistosCorrected[i]->GetBinContent(j) 
                        +1/(double)plots.energyHistosCorrected[0]->GetBinContent(j)
                        +1/(double)monCounts[0]
                        +1/(double)monCounts[i]
                        ),0.5)
                        /(targetMass[targetOrder[i]]
                         *AVOGADROS_NUMBER
                         *pow(10.,-24) // convert cm^2 to barns 
                         *output.crossSections[targetOrder[i]]->back() // error of log(x) ~ (errorOfX)/x
                         /
                         ((pow(targetDiameter[targetOrder[i]]/2,2)*PI // area of cylinder end
                           *(double)targetMolMass[targetOrder[i]])))
                        );
            }

            // record the energy bin for each cross section
            if(i==0)
            {
                output.energyBins.push_back(plots.energyHistosCorrected[0]->GetBinCenter(j));
            }
        }
    }
}

void fillCSGraphs(string CSFileName, vector<int> targetOrder)
{
    // create a file to hold calculated cross sections
    TFile *CSfile = new TFile(CSFileName.c_str(),"RECREATE");

    // remake the cross-section histogram file
    vector<vector<double>*> xError;
    vector<vector<double>*> targetCounts;

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        xError.push_back(new vector<double>);
        targetCounts.push_back(new vector<double>);

        for(int j=0; j<output.energyBins.size(); j++)
        {
            xError[i]->push_back(0);
            targetCounts[i]->push_back(plots.energyHistos[i]->GetBinContent(j));
        }

        plots.CSGraphs.push_back(new TGraphErrors(output.energyBins.size(),&output.energyBins[0],&(output.crossSections[targetOrder[i]]->at(0)),&xError[targetOrder[i]]->at(0),&output.crossSectionsError[targetOrder[i]]->at(0)));
        plots.CSGraphs[i]->SetNameTitle(targetNames[targetOrder[i]].c_str(),targetNames[targetOrder[i]].c_str());
        plots.CSGraphs[i]->Write();
    }

    /*plots.CSGraphsScaledToLit.push_back(new TGraphErrors(output.energyBins.size(),&output.energyBins[0],&(output.crossSectionsScaledToLit[3]->at(0)),&xError[3]->at(0),&output.crossSectionsError[3]->at(0)));
    string temp = targetNames[3] + "Scaled";
    plots.CSGraphsScaledToLit.back()->SetNameTitle(temp.c_str(),temp.c_str());
    plots.CSGraphsScaledToLit.back()->Write();

    temp = "";
    temp = targetNames[5] + "Scaled";
    plots.CSGraphsScaledToLit.push_back(new TGraphErrors(output.energyBins.size(),&output.energyBins[0],&(output.crossSectionsScaledToLit[5]->at(0)),&xError[5]->at(0),&output.crossSectionsError[5]->at(0)));
    plots.CSGraphsScaledToLit.back()->SetNameTitle(temp.c_str(),temp.c_str());
    plots.CSGraphsScaledToLit.back()->Write();
    */

    CSfile->Write();
    CSfile->Close();

    // read literature data for natural Sn
    TFile *litData = new TFile("/data2/analysis/literatureData.root","READ");
    TGraphErrors *SnNatLitData = (TGraphErrors*)litData->Get("Natural Sn (n,tot)");

    // scale Sn112 and Sn124 cross sections using the literature value for natural Sn
    for(int i=0; i<output.energyBins.size(); i++)
    {
        // avoid "divide by 0"
        if(output.crossSections[4]->at(i) != 0)
        {
            double litValue = SnNatLitData->Eval(output.energyBins[i]);
            output.crossSectionsScaledToLit[3]->push_back(output.crossSections[3]->at(i)*(litValue/output.crossSections[4]->at(i)));
            output.crossSectionsScaledToLit[5]->push_back(output.crossSections[5]->at(i)*(litValue/output.crossSections[4]->at(i)));
        }
    }

    litData->Close();
}

// Loop through all trees (one per channel) and populate their data into basic
// histograms. Then, calculate the TOF, cross-section, etc using the channel 4
// data and produce histograms of these calculated quantities.
void fillHistos()
{
    // fill basic histograms for DPP mode in each channel

    // first loop through all channel-specific DPP-mode trees
    for(int i=0; (size_t)i<orchard.size(); i++)
    {
        // create a channel-specific directory for putting histograms inside
        gDirectory->cd("/");
        gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());
        gDirectory->GetDirectory(dirs[i].c_str())->cd();

        // instantiate DPP-mode histograms
        TH1I* macroNoH = new TH1I("macroNoH","macroNo",150000,0,150000);
        macroNoH->GetXaxis()->SetTitle("macropulse number of each event");

        TH1I* evtNoH = new TH1I("evtNoH","evtNo",2500,0,2500);
        evtNoH->GetXaxis()->SetTitle("event number of each event");

        TH1I* macroTimeH = new TH1I("macroTimeH","macroTime",6000,0,6000000000);
        macroTimeH->GetXaxis()->SetTitle("macropulse time zero for each event");

        //TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",6000,0,6000000000);
        //completeTimeH->GetXaxis()->SetTitle("complete time for each event");

        //TH1I* microTimeH = new TH1I("microTimeH","microTime",2000,0,2000);
        //microTimeH->GetXaxis()->SetTitle("time since start of micro for each event");

        TH1I* targetPosH = new TH1I("targetPosH","targetPos",7,0,7);
        targetPosH->GetXaxis()->SetTitle("target position of each event");

        TH1I* sgQH = new TH1I("sgQH","sgQ",3500,0,35000);
        sgQH->GetXaxis()->SetTitle("short gate integrated charge for each event");

        TH1I* lgQH = new TH1I("lgQH","lgQ",7000,0,70000);
        lgQH->GetXaxis()->SetTitle("long gate integrated charge for each event");

        // create a subdirectory for holding DPP-mode waveform data
        gDirectory->mkdir("waveformsDir","raw DPP waveforms");

        // reattach to the channel-specific tree for reading out data
        setBranches(orchard[i]);
        int totalEntries = orchard[i]->GetEntries();

        // reduce entries to sort for diagnostic purposes
        //totalEntries /= 2;
        cout << "Populating " << dirs[i] << " histograms..." << endl;

        waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
        waveformsDir->cd();

        // loop through the channel-specific tree and populate histos
        for(int j=0; j<totalEntries; j++)
        {
            orchard[i]->GetEntry(j);

            if(lgQ!=65535 && sgQ!=32767)
            {
                macroNoH->Fill(macroNo);
                evtNoH->Fill(evtNo);
                macroTimeH->Fill(macroTime);
                //completeTimeH->Fill(completeTime);
                targetPosH->Fill(targetPos);
                sgQH->Fill(sgQ);
                lgQH->Fill(lgQ);
            }

            // if waveform data for this event exist, we want to populate
            // a histogram to display it

            // only plot 1 out of 10000 waveforms to save space and processing
            // time

            if(waveform->size() > 0 && j%10000 == 0)
            {
                stringstream temp;
                temp << "macroNo " << macroNo << ", evtNo " << evtNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

                // loop through waveform data and fill histo
                for(int k=0; (size_t)k<waveform->size(); k++)
                {
                    waveformH->SetBinContent(k,waveform->at(k));
                }
            }

            /*cout << "completeTime = " << completeTime << endl;

              if(macroNo>0)
              {
              break;
              }*/
            if(j%1000==0)
            {
                cout << "Processed " << j << " events...\r";
                fflush(stdout);
            }
        }
        cout << "Processed " << totalEntries << " in " << dirs[i] << " histograms." << endl;
    }

    // fill basic histograms for waveform mode in each channel

    // loop through all channel-specific waveform-mode trees
    for(int i=0; (size_t)i<orchardW.size(); i++)
    {
        // create a channel-specific directory for each tree
        gDirectory->cd("/");
        string tempName = dirs[i] + "WaveformMode";
        gDirectory->mkdir(tempName.c_str(),tempName.c_str());

        gDirectory->GetDirectory(tempName.c_str())->cd();

        // instantiate histograms inside the channel-specific directory
        TH1I* macroNoH = new TH1I("macroNoH","macroNo",100000,0,100000);
        macroNoH->GetXaxis()->SetTitle("macropulse number of each event");

        TH1I* evtNoH = new TH1I("evtNoH","evtNo",2500,0,2500);
        evtNoH->GetXaxis()->SetTitle("event number of each event");

        TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",6000,0,60000000);
        completeTimeH->GetXaxis()->SetTitle("complete time for each event");

        // create subdirectory for holding waveform-mode waveform data
        gDirectory->mkdir("waveformsDir","raw waveform-mode waveforms");
        waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
        waveformsDir->cd();

        // fill waveform mode histograms
        setBranchesW(orchardW[i]);
        int totalEntries = orchardW[i]->GetEntries();
        cout << totalEntries << " total waveforms detected in channel " << i*2 << endl;

        // we need label the number of waveform-mode macropulses to make
        // uniquely named histograms
        int waveformNo = 0;

        for(int j=0; j<totalEntries; j++)
        {
            //cout << "looping thru waveform mode events, #" << j << endl;
            //fflush(stdout);

            orchardW[i]->GetEntry(j);

            macroNoH->Fill(macroNo);
            evtNoH->Fill(evtNo);
            completeTimeH->Fill(completeTime);

            if(evtNo==0)
            {
                // new macropulse in waveform mode
                // create a new plot to hold the new waveforms of this macropulse
                stringstream temp;
                temp << "full waveform " << waveformNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),360000,0,720000);

                // set the start of the macropulse to the first event timer
                waveformNo++;
            }

            for(int k=0; (size_t)k<waveform->size(); k++)
            {
                waveformH->SetBinContent(k+floor((fmod(waveform->size()*2,MICRO_PERIOD)+1)*evtNo*MICRO_PERIOD/(double)2),waveform->at(k));
            }
        }
    }
}

/*void matchWaveforms()
  {
// we want to pull events from the overlap window where both ch4 and ch6
// record the same events (the signal should be identical in both channels)
int ch6Entries = orchard[3]->GetEntries();
int ch4Entries = orchard[2]->GetEntries();

TH1I* ch4ch6timeDiff = new TH1I("ch4ch6timeDiff","ch6-ch4 times for same event",300,-1,1);
TH2I *timeDiffQ = new TH2I("timeDiffQ","Delta-t between ch6-ch4 vs. lgQ",60,-1.2,1.2,256,0,65536);

cout << "Searching through orchard[3] for events in the overlap region" << endl;

gDirectory->cd("/");

// keep track of how many histograms have been drawn; stop looping after a
// set number have been filled (i.e. first 50)
int plots = 0;

setBranches(orchard[3]); // pointed at the ch6Tree now

// keep track of how far we've looped through the ch4 tree so we don't have
// to loop from the beginning each time we're looking for a parallel event
int currentIndex = 0;

// voltage offset (in ADC units) between channel 4 and channel 6
// calculated by averaging all waveform values for each matched
// event and taking the difference between ch4 & ch6
//int offset = 0;

// loop through the channel-specific tree and populate histos
for(int j=0; j<ch6Entries && plots<=10000 plot only first 50; j++)
{
orchard[3]->GetEntry(j);

double ch6timeDiff = completeTime-macroTime;
double ch6microTime = fmod(ch6timeDiff,MICRO_PERIOD);
int ch6microNo = floor(ch6timeDiff/MICRO_PERIOD);

// keep track of this event's macroNo to compare it with events in ch4;
// this will help us identify the parallel event in ch4 with the same
// timestamp
unsigned int ch6MacroNo = macroNo;

if(ch6microTime > 275 && ch6microTime < 325)
{
// used to calculate the channel 6 average value of the waveform and
// find the baseline offset between ch6 and ch4
//int ch6Average = 0; 
stringstream temp;
temp << "macroNo " << macroNo << ", scavenger event number " << evtNo;
waveformCh6= new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

// loop through waveform data and fill histo
for(int k=0; (size_t)k<waveform->size(); k++)
{
waveformCh6->SetBinContent(k,waveform->at(k)+84);
//ch6Average += waveform->at(k);
}

// loop through ch4 tree to find the same event (same timestamp and
// macropulse as the ch6 event)

//cout << "Found a ch6 event in the window. Moving to ch4 events..." << endl;

setBranches(orchard[2]); // pointed at the ch6Tree now
int microNo = 0;

for(int i=currentIndex; i<ch4Entries; i++)
{
orchard[2]->GetEntry(i);

if(macroNo == ch6MacroNo)
{
    // in the correct macropulse - loop through ch4 to find the
    // the event whose timestamp matches the one in channel 6
    double ch4timeDiff = completeTime-macroTime;
    double ch4microTime = fmod(ch4timeDiff,MICRO_PERIOD);
    microNo = floor(ch4timeDiff/MICRO_PERIOD);

    //cout << "Found same macro; ch4timeDiff = " << ch4timeDiff << endl;

    if(microNo == ch6microNo)
    {
        //cout << "Found same micro; ch4microTime = " << ch4microTime << endl;
        if(ch4microTime >= ch6microTime-5 && ch4microTime <= ch6microTime+5)
        {
            // found the parallel event in ch4 - same time, same
            // micropulse
            //
            // plot both the ch6 and ch4 events so we can plot the
            // histograms together later

            ch4ch6timeDiff->Fill(ch6microTime-ch4microTime);
            timeDiffQ->Fill(ch6microTime-ch4microTime, lgQ);

            //cout << "ch4microTime = " << ch4microTime << ", ch4microNo = " << microNo << ", ch6microTime = " << ch6microTime << ", ch6microNo = " << ch6microNo << endl;

            // used to calculate the channel 4 average value of the waveform and
            // find the baseline offset between ch6 and ch4
            //int ch4Average = 0;

            // clear the histo title stringstream
            temp.str("");
            temp << "macroNo " << macroNo << ", ch4 event number " << evtNo << " w/ lgQ=65535";
            waveformCh4 = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

            for(int k=0; (size_t)k<waveform->size(); k++)
            {
                waveformCh4->SetBinContent(k,waveform->at(k));
                //   ch4Average += waveform->at(k);
            }

            //offset += (ch4Average-ch6Average)/(double)waveform->size();

            // Draw both histos on same axes
            //waveformCh6->Draw();
            //waveformCh4->Draw("same");

            plots++;

            break;
        }
    }

    else if (microNo > ch6microNo)
    {
        // went past ch4 micropulse without finding the event
        // parallel to ch6
        //cout << "Failed to find parallel event. ch6microTime = " << ch6microTime << ", ch6microNo = " << ch6microNo << endl;
        //cout << "currentIndex = " << i << endl;
        break;
    }
}

else if (macroNo > ch6MacroNo)
{
    break;
}
currentIndex = i;
}
// point back to ch6 events in preparation for next loop through
// ch6Tree
setBranches(orchard[3]);
}
}
//offset /= (double)plots;
//cout << offset << endl;
}
*/

int main(int argc, char *argv[])
{
    cout << endl << "Entering ./histos..." << endl;

    if(TOF_BINS%NUMBER_ENERGY_BINS!=0)
    {
        cout << "Error: number of TOF bins must be a multiple of the number of energy bins." << endl;
        exit(1);
    }

    string sortedDataFileName = argv[1];
    sortedDataFileName += "resort.root";

    string histoFileName = argv[1];
    histoFileName += "histos.root";

    string CSFileName = argv[1];
    CSFileName += "cross-sections.root";

    string waveformFileName = argv[1];
    waveformFileName += "waveform.root";

    string runNumber = argv[2];

    // report the number of counts in the monitor paddle for each target (diagnostic)
    /*stringstream monitorCountsName;
      monitorCountsName << outpath <<  "/analysis/run" << runNumber << "/" << treeName.str() << "_monitorCounts.log";
      monitorCounts.open(monitorCountsName.str());
      monitorCounts.precision(13);
      */

    TFile* sortedDataFile = new TFile(sortedDataFileName.c_str(),"READ");
    if(!sortedDataFile->IsOpen())
    {
        cout << "Error: failed to open resort.root" << endl;
        exit(1);
    }

    TTree* ch0Tree = (TTree*)sortedDataFile->Get("ch0ProcessedTree");
    TTree* ch2Tree = (TTree*)sortedDataFile->Get("ch2ProcessedTree");
    TTree* ch4Tree = (TTree*)sortedDataFile->Get("ch4ProcessedTree");
    //TTree* ch6Tree = (TTree*)sortedDataFile->Get("ch6ProcessedTree");
    TTree* ch0TreeW = (TTree*)sortedDataFile->Get("ch0ProcessedTreeW");
    TTree* ch2TreeW = (TTree*)sortedDataFile->Get("ch2ProcessedTreeW");
    TTree* ch4TreeW = (TTree*)sortedDataFile->Get("ch4ProcessedTreeW");

    orchard.push_back(ch0Tree);
    orchard.push_back(ch2Tree);
    orchard.push_back(ch4Tree);
    // uncomment to include scavenger data
    //orchard.push_back(ch6Tree);

    orchardW.push_back(ch0TreeW);
    orchardW.push_back(ch2TreeW);
    orchardW.push_back(ch4TreeW);

    // Target order changes between runs. So use the run number to map the
    // correct target to where it was during that run in the target changer

    // Targets are labeled by number as follows:
    // blank = 0, sc = 1, lc = 2, Sn112 = 3, NatSn = 4, Sn124 = 5

    vector<int> targetOrder;

    if(stoi(runNumber)<=5)
    {
        cout << "Neon run - stopping sort." << endl;
        exit(0);
    }

    if(stoi(runNumber)<=151)
    {
        // blank, short carbon, long carbon, Sn112, NatSn, Sn124
        targetOrder.push_back(0);
        targetOrder.push_back(1);
        targetOrder.push_back(2);
        targetOrder.push_back(3);
        targetOrder.push_back(4);
        targetOrder.push_back(5);
    }

    else if(stoi(runNumber)==152)
    {
        // blank, Sn112, NatSn, Sn124, short carbon, long carbon
        targetOrder.push_back(0);
        targetOrder.push_back(3);
        targetOrder.push_back(4);
        targetOrder.push_back(5);
        targetOrder.push_back(1);
        targetOrder.push_back(2);
    }

    else if(stoi(runNumber)>=153 && stoi(runNumber)<=168)
    {
        // blank, Sn112, NatSn, Sn124
        targetOrder.push_back(0);
        targetOrder.push_back(3);
        targetOrder.push_back(4);
        targetOrder.push_back(5);
    }

    else if(stoi(runNumber)>=169 && stoi(runNumber)<=180)
    {
        // blank, Sn112, NatSn, Sn124, short carbon
        targetOrder.push_back(0);
        targetOrder.push_back(3);
        targetOrder.push_back(4);
        targetOrder.push_back(5);
        targetOrder.push_back(1);
    }

    // increase precision to handle outputted times (for troubleshooting)
    cout.precision(13);

    bool skippedHistoFilling = true;

    // open output file to contain histos
    TFile* histoFile;
    histoFile = new TFile(histoFileName.c_str());
    if(!histoFile->IsOpen())
    {
        // No histogram file - need to create it and fill it before moving on to
        // cross-section histos
        histoFile = new TFile(histoFileName.c_str(),"RECREATE");

        // prepare the root file with 4 directories, one for each channel
        // these directories will hold basic variable histograms showing the
        // raw data in each tree, plus TOF, x-sections, etc histograms
        fillHistos();
        // fill TOF, cross-section, etc. histos for channels 4, 6
        fillAdvancedHistos(2, waveformFileName, histoFile);
        //fillAdvancedHistos(3);

        // indicate that plot vectors have already been linked to histograms,
        // so there's no need to relink later
        skippedHistoFilling = false;

        sortedDataFile->Close();

        histoFile->Write();
    }
    else
    {
        histoFile->Close();
        histoFile = new TFile(histoFileName.c_str(),"UPDATE");
        histoFile->cd(dirs[2].c_str());
        gDirectory->Delete("deadtime*;*");
        gDirectory->Delete("*Corrected;*");
    }

    // Calculate deadtime using waveform-mode data, and apply correction to
    // DPP-mode data
    calculateDeadtime(waveformFileName, skippedHistoFilling);

    // Calculate cross-sections using target data and corrected energy histograms
    calculateCS(targetOrder, histoFile);

    // Fill cross-section histograms using calculated cross-section data
    fillCSGraphs(CSFileName, targetOrder);

    // Modify plots
    /*for(int i=0; i<NUMBER_OF_TARGETS; i++)
      {
      plots.energyHistos[i]->GetXaxis()->SetRangeUser(0,700);
      plots.energyHistosCorrected[i]->GetXaxis()->SetRangeUser(0,700);
      }*/

    // print out waveforms for first 50 events that occur in both ch4 and ch6
    //matchWaveforms();

    sortedDataFile->Close();
    histoFile->Close();
}
