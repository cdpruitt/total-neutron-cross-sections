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

// scavenger event stream
ofstream scavengerEvents;

// summedDet event stream
ofstream summedDetEvents;

// fileIn/fileOut names to be accessed to open files
stringstream fileInName, fileOutName, fileCSName;



/* Experimental constants */

const double FLIGHT_DISTANCE = 2672; // detector distance from source, in cm
// 2080 for Neon  

// Period of micropulses
const double MICRO_PERIOD = 1788.814; // in ns

// Physical constants
const double C = 299792458; // speed of light in m/s

const double NEUTRON_MASS = 939.56536; // in MeV/c^2

double AVOGADROS_NUMBER = 6.022*pow(10.,23.); // in atoms/mol



/* Target data */

const int NUM_TARGETS = 6;

const string targetNames[NUM_TARGETS] = {"blank", "shortCarbon", "longCarbon", "Sn112", "NatSn", "Sn124"}; 

// physical target data, listed in order of target names above:
const double targetlength[NUM_TARGETS] =  {0, 1.37,  2.74,  1.365, 1.370, 1.370}; // cm
const double targetMolMass[NUM_TARGETS] = {0, 12.01, 12.01, 112,   118.7, 124}; // g/mol
const double targetdensity[NUM_TARGETS] = {0, 2.2,   2.2,   6.89,  7.31,  7.63}; // g/cm^3


/* Plotting data*/

// number of bins in the raw energy histograms and in the cross-section
// histograms
const int NUM_ENERGY_BINS = 360;

// declare arrays to hold the scaled cross-sections of each target; i = target #, j = bin
vector<vector<double>*> sigma;
vector<vector<double>*> sigmaLog;

int energyCSBinRatio = NUM_ENERGY_BINS/NUM_ENERGY_BINS;

vector<double> sigmaXAxis;
vector<double> sigmaXAxisLog;

const int TOF_RANGE = 1800; // in ns
const int TOF_BINS = 18000;

const double ENERGY_LOWER_BOUND = 1; // cross-section plots' lower bound, in MeV
const double ENERGY_UPPER_BOUND = 500; // cross-section plots' upper bound, in MeV

/* event variables */

// variables for holding tree event data to be histogrammed
unsigned int evtNo, macroNo, targetPos;
unsigned int sgQ, lgQ;
double completeTime, macroTime;

vector<int> *waveform; // for holding one event's waveform data


/* ROOT and organizational variables */

// ROOT file directory structure 
string dirs[4] = {"targetChanger","monitor","detS","scavenger"};

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

// keep track of which order the targets are in, based on which run number we're
// sorting
vector<int> order;

struct Plots
{
    vector<TH1I*> energyHistos;
    vector<TH1I*> energyHistosLog;
    vector<TH1I*> correctedEnergyHistos;
    vector<TH1I*> correctedEnergyHistosLog;

    vector<TGraph*> CSGraphs;
    vector<TGraph*> CSGraphsLog;

} plots;

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

// Create a copy of an input histogram and re-bin it to the log scale
TH1* logBins(TH1 *inputHisto)
{
    string newName;
    newName = inputHisto->GetName();
    newName += "Log";

    double newXMin = (((TAxis*)inputHisto->GetXaxis())->GetXmin());
    if (newXMin <= 0)
    {
        cout << "Error: can't take log of negative energy on cross-section plot" << endl;
        exit(1);
    }

    newXMin = TMath::Log10(newXMin);

    // Pull bin data from input histo, and map to the log scale:
    TAxis* axis = inputHisto->GetXaxis();
    int nBins = axis->GetNbins();

    TH1* outputHisto = new TH1D(newName.c_str(),newName.c_str(),nBins, newXMin,
            TMath::Log10(((TAxis*)inputHisto->GetXaxis())->GetXmax()));

    TAxis* newAxis = outputHisto->GetXaxis();

    double xMin = newAxis->GetXmin();
    double xMax = newAxis->GetXmax();

    double binWidth = (xMax-xMin)/nBins;
    double *newBins = new double[nBins+1];

    for(int i=0; i<=nBins; i++)
    {
        newBins[i] = TMath::Power(10, xMin+i*binWidth);
    }

    // Assign the log-scale bins to the new histo
    ((TAxis*)outputHisto->GetXaxis())->Set(nBins,newBins);
    delete newBins;

    return outputHisto;
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

vector<double> scaleBins(vector<double> inputBins, int nInputBins, int scaledown)
{
    vector<double> outputBins((int)floor(nInputBins/scaledown));
    for(int i=0; i<nInputBins; i++)
    {
        if(i>floor(nInputBins/scaledown)*scaledown)
        {
            break;
        }
        outputBins[(int)floor(i/scaledown)] += inputBins[i];
    }

    for(int i=0; i<nInputBins/scaledown; i++)
    {
        outputBins[i] /= scaledown;
    }

    return outputBins;
}

// Map an input histogram with bins in the time domain to equivalent bins in the relativistic
// kinetic energy domain
TH1* timeBinsToRKEBins(TH1 *inputHisto)
{
    string newName;
    newName = inputHisto->GetName();
    newName += "RKE";

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

        minimumTime = inputHisto->GetBinCenter(i);
        minimumBin = i;
    }

    int nUnscaledEnergyBins = nOldBins-minimumBin;
    double tentativeEnergy = tofToRKE(minimumTime);
    if(tentativeEnergy==-1)
    {
        cout << "Error: energy of old min time " << minimumTime << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    double newXMax = tentativeEnergy;

    double oldXMax = (((TAxis*)inputHisto->GetXaxis())->GetXmax());

    tentativeEnergy = tofToRKE(oldXMax);
    if(tentativeEnergy==-1)
    {
        cout << "Error: energy of old max time " << oldXMax << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    double newXMin = tentativeEnergy;

    // Remap bins from old histo to new histo
    vector<double> unscaledEnergyBins(nUnscaledEnergyBins);

    // Reorder bins to go from lowest energy (shortest time) to highest energy (longest time)
    for(int i=0; i<unscaledEnergyBins.size(); i++)
    {
        unscaledEnergyBins[i] = tofToRKE(oldAxis->GetBinCenter(nOldBins-i-1));
    }

    // Downscale bins to desired granularity
    vector<double> scaledEnergyBins = scaleBins(unscaledEnergyBins, nUnscaledEnergyBins, TOF_BINS/NUM_ENERGY_BINS);

    TH1* outputHisto = new TH1D(newName.c_str(),
                                newName.c_str(),
                                floor(nUnscaledEnergyBins/(TOF_BINS/NUM_ENERGY_BINS)),
                                newXMin,
                                newXMax);

    // Assign the remapped bins to the new histo
    ((TAxis*)outputHisto->GetXaxis())->Set(floor(nUnscaledEnergyBins/(TOF_BINS/NUM_ENERGY_BINS)),
&scaledEnergyBins[0]);

    return outputHisto;
}

// Populate advanced histograms (TOFs, cross-section, etc) calculated using
// data from ch4 and ch6 trees
void fillAdvancedHistos(int detIndex)
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
    TH2I *triangleRKE = new TH2I("triangleRKE","Pulse integral vs. relativistic KE",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,2048,0,65536);
    TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
    TH2I *rKElgQ = new TH2I("lgQrKE","relativistic KE vs. long gate Q",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,2048,0,65536);
    TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
    microNoH->GetXaxis()->SetTitle("micropulse number of each event");

    TH1I *blankTOF = new TH1I("blankTOF","blank TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target1TOF = new TH1I("target1TOF","target 1 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target2TOF = new TH1I("target2TOF","target 2 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target3TOF = new TH1I("target3TOF","target 3 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target4TOF = new TH1I("target4TOF","target 4 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target5TOF = new TH1I("target5TOF","target 5 TOF",TOF_BINS,0,TOF_RANGE);

    // create raw (unnormalized) neutron energy plots
    plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(blankTOF)); 
    plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(target1TOF)); 
    plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(target2TOF));
    plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(target3TOF));
    plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(target4TOF));
    plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(target5TOF));

    // create raw log-scaled neutron energy plots
    plots.energyHistosLog.push_back((TH1I*)logBins(plots.energyHistos[0]));
    plots.energyHistosLog.push_back((TH1I*)logBins(plots.energyHistos[1]));  
    plots.energyHistosLog.push_back((TH1I*)logBins(plots.energyHistos[2]));
    plots.energyHistosLog.push_back((TH1I*)logBins(plots.energyHistos[3]));
    plots.energyHistosLog.push_back((TH1I*)logBins(plots.energyHistos[4]));
    plots.energyHistosLog.push_back((TH1I*)logBins(plots.energyHistos[5]));

    // create corrected (but still flux unnormalized) neutron energy plots
    plots.correctedEnergyHistos.push_back((TH1I*)plots.energyHistos[0]->Clone("blankCorrected"));
    plots.correctedEnergyHistos.push_back((TH1I*)plots.energyHistos[1]->Clone("target1Corrected")); 
    plots.correctedEnergyHistos.push_back((TH1I*)plots.energyHistos[2]->Clone("target2Corrected"));
    plots.correctedEnergyHistos.push_back((TH1I*)plots.energyHistos[3]->Clone("target3Corrected"));
    plots.correctedEnergyHistos.push_back((TH1I*)plots.energyHistos[4]->Clone("target4Corrected"));
    plots.correctedEnergyHistos.push_back((TH1I*)plots.energyHistos[5]->Clone("target5Corrected"));

    // create corrected (but still flux unnormalized) log-scaled neutron energy plots
    plots.correctedEnergyHistosLog.push_back((TH1I*)plots.energyHistosLog[0]->Clone("blankCorrectedLog"));
    plots.correctedEnergyHistosLog.push_back((TH1I*)plots.energyHistosLog[1]->Clone("target1CorrectedLog"));
    plots.correctedEnergyHistosLog.push_back((TH1I*)plots.energyHistosLog[2]->Clone("target2CorrectedLog"));
    plots.correctedEnergyHistosLog.push_back((TH1I*)plots.energyHistosLog[3]->Clone("target3CorrectedLog"));
    plots.correctedEnergyHistosLog.push_back((TH1I*)plots.energyHistosLog[4]->Clone("target4CorrectedLog"));
    plots.correctedEnergyHistosLog.push_back((TH1I*)plots.energyHistosLog[5]->Clone("target5CorrectedLog"));

    // create neutron energy plots using only micropulses with no gammas
    TH1I *noGBlank = new TH1I("noGBlank","no gamma in micro, blank",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);
    TH1I *noGTarget1 = new TH1I("noGTarget1","no gamma in micro, target 1",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);
    TH1I *noGTarget2 = new TH1I("noGTarget2","no gamma in micro, target 2",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);
    TH1I *noGTarget3 = new TH1I("noGTarget3","no gamma in micro, target 3",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);
    TH1I *noGTarget4 = new TH1I("noGTarget4","no gamma in micro, target 4",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);
    TH1I *noGTarget5 = new TH1I("noGTarget5","no gamma in micro, target 5",NUM_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);

    // create log-scaled neutron energy plots using only micropulses with no gammas
    TH1I *noGBlankLog = new TH1I("noGBlankLog","no gamma in micro, log E, blank",NUM_ENERGY_BINS,TMath::Log10(ENERGY_LOWER_BOUND),TMath::Log10(ENERGY_UPPER_BOUND));
    TH1I *noGTarget1Log = new TH1I("noGTarget1Log","no gamma in micro, log E, target 1",NUM_ENERGY_BINS,TMath::Log10(ENERGY_LOWER_BOUND),TMath::Log10(ENERGY_UPPER_BOUND));
    TH1I *noGTarget2Log = new TH1I("noGTarget2Log","no gamma in micro, log E, target 2",NUM_ENERGY_BINS,TMath::Log10(ENERGY_LOWER_BOUND),TMath::Log10(ENERGY_UPPER_BOUND));
    TH1I *noGTarget3Log = new TH1I("noGTarget3Log","no gamma in micro, log E, target 3",NUM_ENERGY_BINS,TMath::Log10(ENERGY_LOWER_BOUND),TMath::Log10(ENERGY_UPPER_BOUND));
    TH1I *noGTarget4Log = new TH1I("noGTarget4Log","no gamma in micro, log E, target 4",NUM_ENERGY_BINS,TMath::Log10(ENERGY_LOWER_BOUND),TMath::Log10(ENERGY_UPPER_BOUND));
    TH1I *noGTarget5Log = new TH1I("noGTarget5Log","no gamma in micro, log E, target 5",NUM_ENERGY_BINS,TMath::Log10(ENERGY_LOWER_BOUND),TMath::Log10(ENERGY_UPPER_BOUND));

    // create TOF plots for events that come first, second, or third in their
    // micro
    TH1I *firstInMicro = new TH1I("firstInMicro","first in micro time of flight",TOF_BINS,0,TOF_RANGE);
    TH1I *secondInMicro = new TH1I("secondInMicro","second in micro time of flight",TOF_BINS,0,TOF_RANGE);
    TH1I *thirdInMicro = new TH1I("thirdInMicro","third in micro time of flight",TOF_BINS,0,TOF_RANGE);

    // create TOF plots for events that come first in their micro, split by
    // target
    TH1I *fimBlank = new TH1I("fimBlank","first in micro, blank",TOF_BINS,0,TOF_RANGE);
    TH1I *fimTarget1 = new TH1I("fimTarget1","first in micro, target 1",TOF_BINS,0,TOF_RANGE);
    TH1I *fimTarget2 = new TH1I("fimTarget2","first in micro, target 2",TOF_BINS,0,TOF_RANGE);
    TH1I *fimTarget3 = new TH1I("fimTarget3","first in micro, target 3",TOF_BINS,0,TOF_RANGE);
    TH1I *fimTarget4 = new TH1I("fimTarget4","first in micro, target 4",TOF_BINS,0,TOF_RANGE);
    TH1I *fimTarget5 = new TH1I("fimTarget5","first in micro, target 5",TOF_BINS,0,TOF_RANGE);

    TH1I *fimBlankLog = new TH1I("fimBlankLog","first in micro, blank",TOF_BINS,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget1Log = new TH1I("fimTarget1Log","first in micro, target 1",TOF_BINS,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget2Log = new TH1I("fimTarget2Log","first in micro, target 2",TOF_BINS,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget3Log = new TH1I("fimTarget3Log","first in micro, target 3",TOF_BINS,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget4Log = new TH1I("fimTarget4Log","first in micro, target 4",TOF_BINS,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget5Log = new TH1I("fimTarget5Log","first in micro, target 5",TOF_BINS,0,TMath::Log10(TOF_RANGE));
    /*************************************************************************/

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

        // GATE: require events to come during the macropulse's beam-on period
        // GATE: discard events during target-changer movement
        // GATE: discard events with unphysical integrated charges
        // GATE: discard events with lgQ>=65500 or sgQ>=32750
        //      (i.e., short gate charge is larger than long gate charge)
        if (timeDiff < 650000 && timeDiff > 0 && targetPos != 0 && lgQ<65500 && sgQ<32750 /*&& lgQ/(double)sgQ<2.1 && lgQ/(double)sgQ>1.5*/)
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

            // GATE: discard events with too low of an integrated charge for their energy
            //if (lgQ>500*exp(-(microTime-100)/87))
            //if (lgQ>10*rKE)
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

                /*****************************************************************/
                // Fill troubleshooting plots with event variables (rKE, microtime, etc.)
                TOF->Fill(microTime);
                triangle->Fill(microTime,lgQ);
                sgQlgQ->Fill(sgQ,lgQ);
                rKElgQ->Fill(rKE,lgQ);
                triangleRKE->Fill(microTime,rKE);
                microNoH->Fill(microNo);

                switch(orderInMicro)
                {
                    case 1:
                        firstInMicro->Fill(microTime);
                        break;
                    case 2:
                        secondInMicro->Fill(microTime);
                        break;
                    case 3:
                        thirdInMicro->Fill(microTime);
                        break;
                }
                /*****************************************************************/

                /*****************************************************************/
                // Fill target-specific plots

                if (targetPos>0 && targetPos<NUM_TARGETS)
                {
                    plots.energyHistos[targetPos-1]->Fill(rKE);
                    plots.energyHistosLog[targetPos-1]->Fill(rKE);
                }

                switch (targetPos)
                {
                    case 1:
                        // BLANK

                        if(orderInMicro==1)
                        {
                            fimBlank->Fill(microTime);
                            fimBlankLog->Fill(TMath::Log10(microTime));

                            microsPerTarget[targetPos-1]++;
                        }

                        if(!gammaInMicro)
                        {
                            noGBlank->Fill(rKE);
                            noGBlankLog->Fill(TMath::Log10(rKE));
                        }

                        blankTOF->Fill(microTime);
                        break;

                    case 2:
                        // TARGET 1

                        if(orderInMicro==1)
                        {
                            fimTarget1->Fill(microTime);
                            fimTarget1Log->Fill(TMath::Log10(microTime));

                            microsPerTarget[targetPos-1]++;
                        }

                        if(!gammaInMicro)
                        {
                            noGTarget1->Fill(rKE);
                            noGTarget1Log->Fill(TMath::Log10(rKE));
                        }

                        target1TOF->Fill(microTime);
                        break;

                    case 3:
                        // TARGET 2

                        if(orderInMicro==1)
                        {
                            fimTarget2->Fill(microTime);
                            fimTarget2Log->Fill(TMath::Log10(microTime));

                            microsPerTarget[targetPos-1]++;
                        }

                        if(!gammaInMicro)
                        {
                            noGTarget2->Fill(rKE);
                            noGTarget2Log->Fill(TMath::Log10(rKE));
                        }

                        target2TOF->Fill(microTime);
                        break;

                    case 4:
                        // TARGET 3

                        if(orderInMicro==1)
                        {
                            fimTarget3->Fill(microTime);
                            fimTarget3Log->Fill(TMath::Log10(microTime));

                            microsPerTarget[targetPos-1]++;
                        }

                        if(!gammaInMicro)
                        {
                            noGTarget3->Fill(rKE);
                            noGTarget3Log->Fill(TMath::Log10(rKE));
                        }

                        target3TOF->Fill(microTime);
                        break;

                    case 5:
                        // TARGET 4

                        if(orderInMicro==1)
                        {
                            fimTarget4->Fill(microTime);
                            fimTarget4Log->Fill(TMath::Log10(microTime));

                            microsPerTarget[targetPos-1]++;
                        }

                        if(!gammaInMicro)
                        {
                            noGTarget4->Fill(rKE);
                            noGTarget4Log->Fill(TMath::Log10(rKE));
                        }

                        target4TOF->Fill(microTime);
                        break;

                    case 6:
                        // TARGET 5

                        if(orderInMicro==1)
                        {
                            fimTarget5->Fill(microTime);
                            fimTarget5Log->Fill(TMath::Log10(microTime));

                            microsPerTarget[targetPos-1]++;
                        }

                        if(!gammaInMicro)
                        {
                            noGTarget5->Fill(rKE);
                            noGTarget5Log->Fill(TMath::Log10(rKE));
                        }

                        target5TOF->Fill(microTime);
                        break;

                    default:
                        break;
                }

                // end of energy, time, cross-section gates on events
                /*****************************************************************/
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

    // create raw log-scaled neutron energy plots
    TH1I *blankRawLogPerMacro = (TH1I*)plots.energyHistosLog[0]->Clone("blankLogPerMacro");
    TH1I *target1RawLogPerMacro = (TH1I*)plots.energyHistosLog[1]->Clone("target1LogPerMacro");
    TH1I *target2RawLogPerMacro = (TH1I*)plots.energyHistosLog[2]->Clone("target2LogPerMacro");
    TH1I *target3RawLogPerMacro = (TH1I*)plots.energyHistosLog[3]->Clone("target3LogPerMacro");
    TH1I *target4RawLogPerMacro = (TH1I*)plots.energyHistosLog[4]->Clone("target4LogPerMacro");
    TH1I *target5RawLogPerMacro = (TH1I*)plots.energyHistosLog[5]->Clone("target5LogPerMacro");

    TH1I* perMacroHistos[6] = {blankRawLogPerMacro, target1RawLogPerMacro,
                               target2RawLogPerMacro, target3RawLogPerMacro,
                               target4RawLogPerMacro, target5RawLogPerMacro};

    // Find number of macropulses for each target to use in error calculation
    gDirectory->cd("/");
    gDirectory->GetDirectory("targetChanger")->cd();

    vector<long> tarCounts;

    tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2));
    tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(3));
    tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(4));
    tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(5));
    tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(6));
    tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(7));

    // Scale perMacroHistos histogram by number of macropulses in that target
    for(int i=0; i<6; i++)
    {
        if(tarCounts[i]==0)
        {
            continue;
        }
        perMacroHistos[i]->Scale(pow(10,3)/tarCounts[i]);
    }

    gDirectory->cd("/");
    gDirectory->GetDirectory("detS")->cd();

    // create raw log-scaled neutron energy plots
    //TH1I *blankRawLogPerMicro = (TH1I*)blankRawLog->Clone("blankLogPerMicro");
    //TH1I *target1RawLogPerMicro = (TH1I*)target1RawLog->Clone("target1LogPerMicro");
    //TH1I *target2RawLogPerMicro = (TH1I*)target2RawLog->Clone("target2LogPerMicro");
    //TH1I *target3RawLogPerMicro = (TH1I*)target3RawLog->Clone("target3LogPerMicro");
    //TH1I *target4RawLogPerMicro = (TH1I*)target4RawLog->Clone("target4LogPerMicro");
    //TH1I *target5RawLogPerMicro = (TH1I*)target5RawLog->Clone("target5LogPerMicro");

    /*TH1I* perMicroHistos[6] = {blankRawLogPerMicro, target1RawLogPerMicro,
                               target2RawLogPerMicro, target3RawLogPerMicro,
                               target4RawLogPerMicro, target5RawLogPerMicro};
    */

    TH1I *blankTOFCorrected = (TH1I*)blankTOF->Clone("blankTOFCorrected");
    TH1I *target1TOFCorrected = (TH1I*)target1TOF->Clone("target1TOFCorrected");
    TH1I *target2TOFCorrected = (TH1I*)target2TOF->Clone("target2TOFCorrected");
    TH1I *target3TOFCorrected = (TH1I*)target3TOF->Clone("target3TOFCorrected");
    TH1I *target4TOFCorrected = (TH1I*)target4TOF->Clone("target4TOFCorrected");
    TH1I *target5TOFCorrected = (TH1I*)target5TOF->Clone("target5TOFCorrected");

    TH1I* TOFCorrectedHistos[6] = {blankTOFCorrected, target1TOFCorrected,
                                  target2TOFCorrected, target3TOFCorrected,
                                  target4TOFCorrected, target5TOFCorrected};

    vector<vector<double>> eventsPerBinPerMicro(6,vector<double>(0));

    // "deadtimeFraction" records the fraction of time that the detector is dead, for
    // neutrons of a certain energy. It is target-dependent.
    vector<vector<double>> deadtimeFraction(6, vector<double>(0));

    vector<vector<double>> csCorrection(6, vector<double>(0));

    const double FULL_DEADTIME = 183; // total amount of time after firing when
                                       // detector is at least partially dead to
                                       // incoming pulses (in ns)
    const double PARTIAL_DEADTIME = 9; // amount of time after the end of
                                        // FULL_DEADTIME when detector is
                                        // becoming live again, depending on
                                        // amplitude (in ns)

    TRandom3 *randomizeBin = new TRandom3();

    /*************************************************************************/
    // Perform 'manual' deadtime correction
    /*************************************************************************/

    // loop through all TOF histos
    for(int i=0; i<6; i++)
    {
        cout << "microsPerTarget[i] = " << microsPerTarget[i] << endl;

        //plots.correctedEnergyHistos[i]->Sumw2();
        //plots.correctedEnergyHistosLog[i]->Sumw2();

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
            int k = j-(FULL_DEADTIME+PARTIAL_DEADTIME)*(TOF_BINS/TOF_RANGE);
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
            }
        }

        for(int j=0; j<TOFCorrectedHistos[i]->GetNbinsX(); j++)
        {
            if(deadtimeFraction[i][j] > 0)
            {
                TOFCorrectedHistos[i]->SetBinContent(j,(TOFCorrectedHistos[i]->GetBinContent(j)/(1-deadtimeFraction[i][j])));
            }

            // convert microTime into neutron velocity based on flight path distance
            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/(TOFCorrectedHistos[i]->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)(2*TOF_BINS),TOF_RANGE/(double)(2*TOF_BINS))); // in meters/sec 

            // convert velocity to relativistic kinetic energy
            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            plots.correctedEnergyHistos[i]->Fill(rKE,TOFCorrectedHistos[i]->GetBinContent(j));
            plots.correctedEnergyHistosLog[i]->Fill(rKE,TOFCorrectedHistos[i]->GetBinContent(j));

            TOFCorrectedHistos[i]->SetBinError(j,pow(TOFCorrectedHistos[i]->GetBinContent(j),0.5));
            plots.correctedEnergyHistos[i]->SetBinError(j,pow(plots.correctedEnergyHistos[i]->GetBinContent(j),0.5));
            plots.correctedEnergyHistosLog[i]->SetBinError(j,pow(plots.correctedEnergyHistosLog[i]->GetBinContent(j),0.5));

            /*if(j%100==0 && i==3)
            {
                cout << "bin spacing randomized = " << (TOFCorrectedHistos[i]->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)TOF_BINS,TOF_RANGE/(double)TOF_BINS))
                                                      -(TOFCorrectedHistos[i]->GetBinCenter(j-1)+randomizeBin->Uniform(-TOF_RANGE/(double)TOF_BINS,TOF_RANGE/(double)TOF_BINS))
                                                    << endl;
                cout << "bin spacing normal = " << (TOFCorrectedHistos[i]->GetBinCenter(j))
                                                  -(TOFCorrectedHistos[i]->GetBinCenter(j-1))
                                                    << endl;
            }*/
        }

        /*for(int j=0; j<plots.correctedEnergyHistos[i]->GetNbinsX(); j++)
        {
            csCorrection[i].push_back(0);
            csCorrection[i][j]

        }*/

        //cout << deadtimeCorrection[i][150] << endl;
    }

    // Plot the calculated dead time fractions (for debugging)
    TH1I* deadtimeFractionH = new TH1I("deadtimeFractionH","deadtimeFractionH",TOF_BINS,0,TOF_RANGE);

    for(int j=0; (size_t)j<deadtimeFraction[0].size(); j++)
    {
        deadtimeFractionH->SetBinContent(j,1000000*deadtimeFraction[0][j]);
    }

    /*************************************************************************/
}

void calculateCS()
{
    // Find number of events in the monitor for each target to use in scaling
    // cross-sections

    // switch to the monitor directory
    gDirectory->cd("/");
    gDirectory->GetDirectory("monitor")->cd();

    // to normalize flux between all channels, we'll need to keep track of the
    // number of counts that came in through the monitor during that target's
    // beam time
    vector<long> monCounts;

    for(int i=2; i<8; i++)
    {
        monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i));
        cout << "target position " << i << " monitor counts = " << monCounts[i-2] << endl;
    }

    // switch to the detector directory
    gDirectory->cd("/");
    gDirectory->GetDirectory("detS")->cd();

    // holds the raw target-specific energy histograms in preparation for populating
    // cross-sections
    
    // switch to the scavenger directory
    /*gDirectory->cd("/");
    gDirectory->GetDirectory("scavenger")->cd();

    // create vectors to hold the scavenger data from the dead-time region
    vector<TH1I*> scavHistos;
    vector<TH1I*> scavLogHistos;

    scavLogHistos.push_back((TH1I*)gDirectory->Get("blankLog"));
    scavLogHistos.push_back((TH1I*)gDirectory->Get("target1Log"));
    scavLogHistos.push_back((TH1I*)gDirectory->Get("target2Log"));
    scavLogHistos.push_back((TH1I*)gDirectory->Get("target3Log"));
    scavLogHistos.push_back((TH1I*)gDirectory->Get("target4Log"));
    //scavLogHistos.push_back((TH1I*)gDirectory->Get("target5Log"));
*/
    // Loop through the relativistic kinetic energy histograms and use them
    // to populate cross-section for each target

    int numberOfBins = ((TAxis*)plots.energyHistos[0]->GetXaxis())->GetNbins();

    for(int i=0; (size_t)i<order.size(); i++)
    {
        sigma.push_back(new vector<double>);
        sigmaLog.push_back(new vector<double>);
    }

    for(int i=0; (size_t)i<order.size(); i++)
    {
        for(int j=0; j<numberOfBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(plots.energyHistos[0]->GetBinContent(j) <= 0 || plots.energyHistos[i]->GetBinContent(j) <= 0)
            {
                sigma[order[i]]->push_back(0);
            }

            else
            {
                // calculate the cross-section
                sigma[order[i]]->push_back(-log((plots.energyHistos[i]->GetBinContent(j)/(double)plots.energyHistos[0]->GetBinContent(j))*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)AVOGADROS_NUMBER*pow(10.,-24)/(double)targetMolMass[order[i]])); // in barns
            }

            if(plots.energyHistosLog[0]->GetBinContent(j) <= 0 || plots.energyHistosLog[i]->GetBinContent(j) <= 0)
            {
                sigmaLog[i]->push_back(0);
            }

            else
            {
                // calculate the cross-section
                sigmaLog[order[i]]->push_back(-log((plots.energyHistosLog[i]->GetBinContent(j))/((double)plots.energyHistosLog[0]->GetBinContent(j))*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)AVOGADROS_NUMBER*pow(10.,-24)/(double)targetMolMass[order[i]])); // in barns
            }

            if(i==0)
            {
                sigmaXAxis.push_back(plots.energyHistos[0]->GetBinCenter(j));
                sigmaXAxisLog.push_back(plots.energyHistosLog[0]->GetBinCenter(j));
            }
        }

        /*
        // downscale the sigma array into the correct number of cross-section bins
        double binWidth = (ENERGY_UPPER_BOUND-ENERGY_LOWER_BOUND)/NUM_ENERGY_BINS;
        double binWidthLog = (TMath::Log10(ENERGY_UPPER_BOUND)-TMath::Log10(ENERGY_LOWER_BOUND))/NUM_ENERGY_BINS;

        for(int k=0; k<NUM_ENERGY_BINS; k++)
        {
            sigmaXAxis[k] = ENERGY_LOWER_BOUND+binWidth*(k);
            sigmaLogXAxis[k] = TMath::Power(10,TMath::Log10(ENERGY_LOWER_BOUND)+binWidthLog*k);

            for(int l=0; l<energyCSBinRatio; l++)
            {
                sigma[order[i]][k] += sigma[order[i]][k*energyCSBinRatio+l];
                sigmaLog[order[i]][k] += sigmaLog[order[i]][k*energyCSBinRatio+l];
            }
            sigma[order[i]][k] /= energyCSBinRatio;
            sigmaLog[order[i]][k] /= energyCSBinRatio;
        }
        */
    }
}

void fillCSGraphs()
{
    // remake the cross-section histogram file
    TFile *CSfile = new TFile(fileCSName.str().c_str(),"RECREATE");

    for(int i=0; i<order.size(); i++)
    {
        plots.CSGraphs.push_back(new TGraphErrors(sigmaXAxis.size(),&sigmaXAxis[0],&(sigma[order[i]]->at(0))));
        plots.CSGraphs[i]->SetNameTitle(targetNames[order[i]].c_str(),targetNames[order[i]].c_str());
        plots.CSGraphs[i]->Write();

        plots.CSGraphsLog.push_back(new TGraphErrors(sigmaXAxisLog.size(),&sigmaXAxisLog[0],&(sigmaLog[order[i]]->at(0))));
        plots.CSGraphsLog[i]->SetNameTitle((targetNames[order[i]]+"Log").c_str(),(targetNames[order[i]]+"Log").c_str());
        plots.CSGraphsLog[i]->Write();

        //plots.CSGraphsLog.push_back(logBins(plots.CSGraphs[i]new TGraphErrors(NUM_ENERGY_BINS,&sigma[i],&sigmaXAxis));

        /*for(int j=0; j<NUM_ENERGY_BINS; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0

            if(sigma[i][j] == 0)
            {
                continue;
            }

            else
            {
                sigma[i][j] /= energyCSBinRatio;
                plots.CSGraphs[i]->SetBinContent(j,sigma[i][j]);
            }
        }

        for(int j=0; j<NUM_ENERGY_BINS; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0

            if(sigmaLog[i][j] == 0)
            {
                continue;
            }

            else
            {
                sigmaLog[i][j] /= energyCSBinRatio;
                plots.CSGraphsLog[i]->SetBinContent(j,sigmaLog[i][j]);
            }
        }*/
    }

    // create relative cross-section plot for 112Sn/124Sn
    /*TH1D *Sn124_plus_Sn112CSLog = (TH1D*)plots.CSGraphsLog[5]->Clone("Sn124_plus_Sn112CSLog");
    Sn124_plus_Sn112CSLog->Add(plots.CSGraphsLog[3]);

    TH1D *Sn124_minus_Sn112CSLog = (TH1D*)plots.CSGraphsLog[5]->Clone("Sn124_minus_Sn112CSLog");
    Sn124_minus_Sn112CSLog->Add(plots.CSGraphsLog[3],-1);

    // rebin and scale these dummy cross-section plots (124+112 and 124-112)
    Sn124_plus_Sn112CSLog->Rebin(NUM_ENERGY_BINS/(double)20);
    Sn124_minus_Sn112CSLog->Rebin(NUM_ENERGY_BINS/(double)20);

    Sn124_plus_Sn112CSLog->Scale(1/(NUM_ENERGY_BINS/(double)20));
    Sn124_minus_Sn112CSLog->Scale(1/(NUM_ENERGY_BINS/(double)20));

    // Divide 124-112 by 124+112 for relative cross-section
    TH1D *relativeSnCSLog = (TH1D*)Sn124_minus_Sn112CSLog->Clone("relativeSnCSLog");
    relativeSnCSLog->SetTitle("#frac{#sigma_{^{124}Sn}-#sigma_{^{112}Sn}}{#sigma_{^{124}Sn}+#sigma_{^{112}Sn}}");
    relativeSnCSLog->Divide(Sn124_plus_Sn112CSLog);
*/
    ifstream SnData("SnNatData.dat");
    if(!SnData.is_open())
    {
        cout << "No Previous Data..." << endl;
        return;
    }

    char dummy[200];

    SnData.getline(dummy,200);
    SnData.getline(dummy,200);

    vector<float> energy;
    vector<float> xsection;
    vector<float> error;

    float dum,dum2,dum3;

    while(!SnData.eof())
    {
        SnData >> dum >> dum2 >> dum3;

        energy.push_back(dum);
        xsection.push_back(dum2);
        error.push_back(dum3);
    }

    TGraphErrors *SnLitLog = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
    //SnLitLog->Draw("AP");

    SnLitLog->GetXaxis()->SetTitle("Energy [MeV]");
    SnLitLog->GetXaxis()->CenterTitle();
    SnLitLog->GetXaxis()->SetRangeUser(ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);

    SnLitLog->GetYaxis()->SetTitle("sigma [b]");
    SnLitLog->GetYaxis()->CenterTitle();
    //SnLitLog->GetYaxis()->SetTitleOffSet(1.5);
    SnLitLog->Write();

    // carbon literature data
    ifstream carbonData("CarbonData.dat");
    if(!carbonData.is_open())
    {
        cout << "No Previous Data..." << endl;
        return;
    }

    carbonData.getline(dummy,200);
    carbonData.getline(dummy,200);

    energy.clear();
    xsection.clear();
    error.clear();

    while(!carbonData.eof())
    {
        carbonData >> dum >> dum2 >> dum3;

        energy.push_back(dum);
        xsection.push_back(dum2);
        error.push_back(dum3);
    }

    TGraphErrors *carbonLitLog = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
    //carbonLitLog->Draw("AP");

    carbonLitLog->GetXaxis()->SetTitle("Energy [MeV]");
    carbonLitLog->GetXaxis()->CenterTitle();
    carbonLitLog->GetXaxis()->SetRangeUser(ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);

    carbonLitLog->GetYaxis()->SetTitle("sigma [b]");
    carbonLitLog->GetYaxis()->CenterTitle();
    //carbonLitLog->GetYaxis()->SetTitleOffSet(1.5);
    carbonLitLog->Write();

    CSfile->Write();
    CSfile->Close();
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

            /*// print out first 5000 scavenger events inside the correct time
            // window (just after the gamma dead time)
            if(i==2 && j <10000)
            {
                summedDetEvents << "summed det event " << j << " completeTime = " << completeTime << endl;
            }

            // print out first 100 scavenger events inside the correct time
            // window (gamma dead time)
            else if(i==3 && microTime<350 && microTime>200)
            {
                scavengerEvents << "scavenger event " << j << " completeTime = " << completeTime << endl;
            }

            else if (j>10000)
            {
                break;
            }*/

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

    // fill TOF, cross-section, etc. histos for channels 4, 6
    fillAdvancedHistos(2);
    //fillAdvancedHistos(3);

    // fill basic histograms for waveform mode in each channel

    // first loop through all channel-specific waveform-mode trees
    for(int i=0; (size_t)i<orchardW.size(); i++)
    {
        // create a channel-specific directory for each tree
        gDirectory->cd("/");
        dirs[i] = dirs[i] + "WaveformMode";
        gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());

        gDirectory->GetDirectory(dirs[i].c_str())->cd();

        // instantiate histograms inside the channel-specific directory
        TH1I* macroNoH = new TH1I("macroNoH","macroNo",100000,0,100000);
        macroNoH->GetXaxis()->SetTitle("macropulse number of each event");

        TH1I* evtNoH = new TH1I("evtNoH","evtNo",2500,0,2500);
        evtNoH->GetXaxis()->SetTitle("event number of each event");

        TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",6000,0,6000000000);
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

void matchWaveforms()
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
    for(int j=0; j<ch6Entries && plots<=10000 /*plot only first 50 */; j++)
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

int main(int argc, char *argv[])
{
    cout << endl << "Entering ./histos..." << endl;
    // needed to avoid ROOT error for header files being incorrectly brought in
    // in both resort.cpp and histos.cpp
    // Look online for more info (I'm not really sure why it's necessary)
    //TApplication app("app",&argc,argv);

    if(TOF_BINS%NUM_ENERGY_BINS!=0)
    {
        cout << "Error: number of TOF bins must be a multiple of the number of energy bins." << endl;
        exit(1);
    }

    // read in the raw file name
    string runDir = argv[1];
    string runNo = argv[2];
    string outpath = argv[3];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    stringstream treeName;
    //stringstream scavengerEventsName;
    //stringstream summedDetEventsName;

    treeName << "run" << runDir << "-" << runNo; 

    //scavengerEventsName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_scavenger.csv";
    //summedDetEventsName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_summedDet.csv";

    //scavengerEvents.open(scavengerEventsName.str());
    //summedDetEvents.open(summedDetEventsName.str());

    //scavengerEvents.precision(10);
    //summedDetEvents.precision(10);

    fileInName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_sorted.root";
    fileOutName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_histos.root";
    fileCSName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_cross-sections.root";

    TFile* file = new TFile(fileInName.str().c_str(),"READ");

    if(!file->IsOpen())
    {
        cout << "Error: failed to open resort.root" << endl;
        exit(1);
    }

    TTree* ch0Tree = (TTree*)file->Get("ch0ProcessedTree");
    TTree* ch2Tree = (TTree*)file->Get("ch2ProcessedTree");
    TTree* ch4Tree = (TTree*)file->Get("ch4ProcessedTree");
    //TTree* ch6Tree = (TTree*)file->Get("ch6ProcessedTree");
    TTree* ch0TreeW = (TTree*)file->Get("ch0ProcessedTreeW");
    TTree* ch2TreeW = (TTree*)file->Get("ch2ProcessedTreeW");
    TTree* ch4TreeW = (TTree*)file->Get("ch4ProcessedTreeW");

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

    if(stoi(runDir)<=5)
    {
        cout << "Neon run - stopping sort." << endl;
        exit(0);
    }

    if(stoi(runDir)<=151)
    {
        // blank, short carbon, long carbon, Sn112, NatSn, Sn124
        order.push_back(0);
        order.push_back(1);
        order.push_back(2);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
    }

    else if(stoi(runDir)==152)
    {
        // blank, Sn112, NatSn, Sn124, short carbon, long carbon
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
        order.push_back(1);
        order.push_back(2);
    }

    else if(stoi(runDir)>=153 && stoi(runDir)<=168)
    {
        // blank, Sn112, NatSn, Sn124
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
    }

    else if(stoi(runDir)>=169 && stoi(runDir)<=180)
    {
        // blank, Sn112, NatSn, Sn124, short carbon
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
        order.push_back(1);
    }

    // increase precision to handle outputted times (for troubleshooting)
    cout.precision(13);
    
    // open output file to contain histos
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"READ");
    if(!fileOut->IsOpen())
    {
        // No histogram file - need to create it and fill it before moving on to
        // cross-section histos
        TFile* fileOut = new TFile(fileOutName.str().c_str(),"CREATE");

        // prepare the root file with 4 directories, one for each channel
        // these directories will hold basic variable histograms showing the
        // raw data in each tree, plus TOF, x-sections, etc histograms
        fillHistos();
        fileOut->Write();
        file->Close();
    }

    // Calculate cross-sections using channels 2, 4, and 6
    calculateCS();

    // Fill cross-section histograms using calculated cross-section data
    fillCSGraphs();

    // Modify plots
    for(int i=0; i<NUM_TARGETS; i++)
    {
        plots.energyHistos[i]->GetXaxis()->SetRangeUser(0,700);
        plots.energyHistosLog[i]->GetXaxis()->SetRangeUser(0,700);
        plots.correctedEnergyHistos[i]->GetXaxis()->SetRangeUser(0,700);
        plots.correctedEnergyHistosLog[i]->GetXaxis()->SetRangeUser(0,700);
    }

    // print out waveforms for first 50 events that occur in both ch4 and ch6
    //matchWaveforms();
        
    fileOut->Close();
}
