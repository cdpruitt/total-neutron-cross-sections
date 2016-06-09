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

double avo = 6.022*pow(10.,23.); // Avogadro's number, in atoms/mol



/* Target data */

// physical target data, listed in order:
// {blank, Sn112, Natural Sn, Sn124, short carbon, long carbon} 

// lengths of each target:
double targetlength[6] = {0,1.37,2.74,1.365,1.370,1.370}; //cm

// molar mass of each target:
double targetMolMass[6] = {0,12.01,12.01,112,118.7,124}; //g/mol

// density of each target:
double targetdensity[6] = {0,2.2,2.2,6.89,7.31,7.63}; //g/cm^3



/* Plotting data*/

// number of bins in the raw energy histograms and in the cross-section
// histograms
const int noBins = 200;

// declare arrays to hold the scaled cross-sections of each target; i = target #, j = bin
double sigma[6][noBins] = {0};
double sigmaLog[6][noBins] = {0};

const int TOF_RANGE = 1800;
const int TOF_BINS = 18000;

const double CS_LOWER_BOUND = 1; // cross-section plots' lower bound, in MeV
const double CS_UPPER_BOUND = 700; // cross-section plots' upper bound, in MeV

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

TH1* logBins(TH1 *inputHisto)
{
    string newName;
    newName = inputHisto->GetName();
    newName += "Log";

    double newXMin = TMath::Log10(((TAxis*)inputHisto->GetXaxis())->GetXmin());
    if (newXMin <= 0)
    {
        newXMin = 1;
    }
    TH1* outputHisto = new TH1D(newName.c_str(),newName.c_str(),noBins, newXMin,
            TMath::Log10(((TAxis*)inputHisto->GetXaxis())->GetXmax()));

    TAxis* axis = outputHisto->GetXaxis();
    int nBins = axis->GetNbins();

    double xMin = axis->GetXmin();
    double xMax = axis->GetXmax();

    double binWidth = (xMax-xMin)/nBins;
    double *newBins = new double[nBins+1];

    for(int i=0; i<=nBins; i++)
    {
        newBins[i] = TMath::Power(10, xMin+i*binWidth);
    }

    ((TAxis*)outputHisto->GetXaxis())->Set(nBins,newBins);
    delete newBins;

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
    TH1I *TOF = new TH1I("TOF","Summed-detector time of flight",1800,0,TOF_RANGE);
    TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",TOF_RANGE,0,MICRO_PERIOD+1,2048,0,65536);
    TH2I *triangleRKE = new TH2I("triangleRKE","Pulse integral vs. relativistic KE",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND,2048,0,65536);
    TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
    TH2I *rKElgQ = new TH2I("lgQrKE","relativistic KE vs. long gate Q",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND,2048,0,65536);
    TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
    microNoH->GetXaxis()->SetTitle("micropulse number of each event");

    TH1I *blankTOF = new TH1I("blankTOF","blank TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target1TOF = new TH1I("target1TOF","target 1 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target2TOF = new TH1I("target2TOF","target 2 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target3TOF = new TH1I("target3TOF","target 3 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target4TOF = new TH1I("target4TOF","target 4 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target5TOF = new TH1I("target5TOF","target 5 TOF",TOF_BINS,0,TOF_RANGE);

    // create raw (unnormalized) neutron energy plots
    TH1I *blankRaw = new TH1I("blank","blank",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target1Raw = new TH1I("target1","target1",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target2Raw = new TH1I("target2","target2",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target3Raw = new TH1I("target3","target3",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target4Raw = new TH1I("target4","target4",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target5Raw = new TH1I("target5","target5",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    // create raw log-scaled neutron energy plots
    TH1I *blankRawLog = new TH1I("blankLog","blank",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target1RawLog = new TH1I("target1Log","target1",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target2RawLog = new TH1I("target2Log","target2",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target3RawLog = new TH1I("target3Log","target3",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target4RawLog = new TH1I("target4Log","target4",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target5RawLog = new TH1I("target5Log","target5",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));

    // create corrected (but still flux unnormalized) neutron energy plots
    TH1I *blankCorrected = new TH1I("blankCorrected","blankCorrected",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target1Corrected = new TH1I("target1Corrected","target1Corrected",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target2Corrected = new TH1I("target2Corrected","target2Corrected",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target3Corrected = new TH1I("target3Corrected","target3Corrected",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target4Corrected = new TH1I("target4Corrected","target4Corrected",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *target5Corrected = new TH1I("target5Corrected","target5Corrected",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    // create corrected (but still flux unnormalized) neutron energy plots
    TH1I *blankCorrectedLog = new TH1I("blankCorrectedLog","blankCorrectedLog",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target1CorrectedLog = new TH1I("target1CorrectedLog","target1CorrectedLog",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target2CorrectedLog = new TH1I("target2CorrectedLog","target2CorrectedLog",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target3CorrectedLog = new TH1I("target3CorrectedLog","target3CorrectedLog",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target4CorrectedLog = new TH1I("target4CorrectedLog","target4CorrectedLog",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *target5CorrectedLog = new TH1I("target5CorrectedLog","target5CorrectedLog",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));

    // create neutron energy plots using only micropulses with no gammas
    TH1I *noGBlank = new TH1I("noGBlank","no gamma in micro, blank",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *noGTarget1 = new TH1I("noGTarget1","no gamma in micro, target 1",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *noGTarget2 = new TH1I("noGTarget2","no gamma in micro, target 2",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *noGTarget3 = new TH1I("noGTarget3","no gamma in micro, target 3",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *noGTarget4 = new TH1I("noGTarget4","no gamma in micro, target 4",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1I *noGTarget5 = new TH1I("noGTarget5","no gamma in micro, target 5",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    // create log-scaled neutron energy plots using only micropulses with no gammas
    TH1I *noGBlankLog = new TH1I("noGBlankLog","no gamma in micro, log E, blank",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *noGTarget1Log = new TH1I("noGTarget1Log","no gamma in micro, log E, target 1",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *noGTarget2Log = new TH1I("noGTarget2Log","no gamma in micro, log E, target 2",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *noGTarget3Log = new TH1I("noGTarget3Log","no gamma in micro, log E, target 3",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *noGTarget4Log = new TH1I("noGTarget4Log","no gamma in micro, log E, target 4",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1I *noGTarget5Log = new TH1I("noGTarget5Log","no gamma in micro, log E, target 5",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));

    // create TOF plots for events that come first, second, or third in their
    // micro
    TH1I *firstInMicro = new TH1I("firstInMicro","first in micro time of flight",1800,0,TOF_RANGE);
    TH1I *secondInMicro = new TH1I("secondInMicro","second in micro time of flight",1800,0,TOF_RANGE);
    TH1I *thirdInMicro = new TH1I("thirdInMicro","third in micro time of flight",1800,0,TOF_RANGE);

    // create TOF plots for events that come first in their micro, split by
    // target
    TH1I *fimBlank = new TH1I("fimBlank","first in micro, blank",1800,0,TOF_RANGE);
    TH1I *fimTarget1 = new TH1I("fimTarget1","first in micro, target 1",1800,0,TOF_RANGE);
    TH1I *fimTarget2 = new TH1I("fimTarget2","first in micro, target 2",1800,0,TOF_RANGE);
    TH1I *fimTarget3 = new TH1I("fimTarget3","first in micro, target 3",1800,0,TOF_RANGE);
    TH1I *fimTarget4 = new TH1I("fimTarget4","first in micro, target 4",1800,0,TOF_RANGE);
    TH1I *fimTarget5 = new TH1I("fimTarget5","first in micro, target 5",1800,0,TOF_RANGE);

    TH1I *fimBlankLog = new TH1I("fimBlankLog","first in micro, blank",1800,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget1Log = new TH1I("fimTarget1Log","first in micro, target 1",1800,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget2Log = new TH1I("fimTarget2Log","first in micro, target 2",1800,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget3Log = new TH1I("fimTarget3Log","first in micro, target 3",1800,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget4Log = new TH1I("fimTarget4Log","first in micro, target 4",1800,0,TMath::Log10(TOF_RANGE));
    TH1I *fimTarget5Log = new TH1I("fimTarget5Log","first in micro, target 5",1800,0,TMath::Log10(TOF_RANGE));

    /*TH1I *fimBlankLog = (TH1I*)logBins(fimBlank);
    TH1I *fimTarget1Log = (TH1I*)logBins(fimTarget1);
    TH1I *fimTarget2Log = (TH1I*)logBins(fimTarget2);
    TH1I *fimTarget3Log = (TH1I*)logBins(fimTarget3);
    TH1I *fimTarget4Log = (TH1I*)logBins(fimTarget4);
    TH1I *fimTarget5Log = (TH1I*)logBins(fimTarget5);
    */

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
        //      (i.e., short gate charge is larger than long gate charge)
        if (timeDiff < 650000 && timeDiff > 0 && targetPos != 0 /*&& lgQ/(double)sgQ<2.1 && lgQ/(double)sgQ>1.5*/)
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

            switch (targetPos)
            {
                case 1:
                    // BLANK
                    blankRaw->Fill(rKE);
                    blankRawLog->Fill(TMath::Log10(rKE));

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

                    target1Raw->Fill(rKE);
                    target1RawLog->Fill(TMath::Log10(rKE));

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

                    target2Raw->Fill(rKE);
                    target2RawLog->Fill(TMath::Log10(rKE));

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

                    target3Raw->Fill(rKE);
                    target3RawLog->Fill(TMath::Log10(rKE));

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

                    target4Raw->Fill(rKE);
                    target4RawLog->Fill(TMath::Log10(rKE));

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

                    target5Raw->Fill(rKE);
                    target5RawLog->Fill(TMath::Log10(rKE));

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

    // create raw log-scaled neutron energy plots
    TH1I *blankRawLogPerMacro = (TH1I*)blankRawLog->Clone("blankLogPerMacro");
    TH1I *target1RawLogPerMacro = (TH1I*)target1RawLog->Clone("target1LogPerMacro");
    TH1I *target2RawLogPerMacro = (TH1I*)target2RawLog->Clone("target2LogPerMacro");
    TH1I *target3RawLogPerMacro = (TH1I*)target3RawLog->Clone("target3LogPerMacro");
    TH1I *target4RawLogPerMacro = (TH1I*)target4RawLog->Clone("target4LogPerMacro");
    TH1I *target5RawLogPerMacro = (TH1I*)target5RawLog->Clone("target5LogPerMacro");

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
    TH1I *blankRawLogPerMicro = (TH1I*)blankRawLog->Clone("blankLogPerMicro");
    TH1I *target1RawLogPerMicro = (TH1I*)target1RawLog->Clone("target1LogPerMicro");
    TH1I *target2RawLogPerMicro = (TH1I*)target2RawLog->Clone("target2LogPerMicro");
    TH1I *target3RawLogPerMicro = (TH1I*)target3RawLog->Clone("target3LogPerMicro");
    TH1I *target4RawLogPerMicro = (TH1I*)target4RawLog->Clone("target4LogPerMicro");
    TH1I *target5RawLogPerMicro = (TH1I*)target5RawLog->Clone("target5LogPerMicro");

    TH1I* perMicroHistos[6] = {blankRawLogPerMicro, target1RawLogPerMicro,
                               target2RawLogPerMicro, target3RawLogPerMicro,
                               target4RawLogPerMicro, target5RawLogPerMicro};

    TH1I *blankTOFCorrected = (TH1I*)blankTOF->Clone("blankTOFCorrected");
    TH1I *target1TOFCorrected = (TH1I*)target1TOF->Clone("target1TOFCorrected");
    TH1I *target2TOFCorrected = (TH1I*)target2TOF->Clone("target2TOFCorrected");
    TH1I *target3TOFCorrected = (TH1I*)target3TOF->Clone("target3TOFCorrected");
    TH1I *target4TOFCorrected = (TH1I*)target4TOF->Clone("target4TOFCorrected");
    TH1I *target5TOFCorrected = (TH1I*)target5TOF->Clone("target5TOFCorrected");

    TH1I* TOFCorrectedHistos[6] = {blankTOFCorrected, target1TOFCorrected,
                                  target2TOFCorrected, target3TOFCorrected,
                                  target4TOFCorrected, target5TOFCorrected};

    TH1I* energyCorrectedHistos[6] = {blankCorrected, target1Corrected,
                                  target2Corrected, target3Corrected,
                                  target4Corrected, target5Corrected};

    TH1I* energyCorrectedLogHistos[6] = {blankCorrectedLog, target1CorrectedLog,
                                  target2CorrectedLog, target3CorrectedLog,
                                  target4CorrectedLog, target5CorrectedLog};

    vector<vector<double>> eventsPerBinPerMicro(6,vector<double>(0));
    vector<vector<double>> deadtimeCorrection(6, vector<double>(0));
    vector<vector<double>> csCorrection(6, vector<double>(0));

    double deadtime = 186; // in ns
    TRandom3 *randomizeBin = new TRandom3();

    /*************************************************************************/
    // Perform 'manual' deadtime correction
    /*************************************************************************/

    // loop through all TOF histos
    for(int i=0; i<6; i++)
    {
        cout << "microsPerTarget[i] = " << microsPerTarget[i] << endl;

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

            deadtimeCorrection[i].push_back(0);
        }

        // find the fraction of the time that the detector is dead for each bin in the micropulse
        for(int j=0; j<eventsPerBinPerMicro[i].size(); j++)
        {
            int k = j-deadtime*(TOF_BINS/TOF_RANGE);
            while(k<j)
            {
                if(k>=0)
                {
                    deadtimeCorrection[i][j] += eventsPerBinPerMicro[i][k];
                }
                k++;
            }
        }

        for(int j=0; j<TOFCorrectedHistos[i]->GetNbinsX(); j++)
        {
            if(deadtimeCorrection[i][j] > 0)
            {
                TOFCorrectedHistos[i]->SetBinContent(j,(TOFCorrectedHistos[i]->GetBinContent(j)/(1-deadtimeCorrection[i][j])));
            }

            // convert microTime into neutron velocity based on flight path distance
            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/(TOFCorrectedHistos[i]->GetBinCenter(j)+randomizeBin->Uniform(-0.5,0.5)); // in meters/sec 

            // convert velocity to relativistic kinetic energy
            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            energyCorrectedHistos[i]->Fill(rKE,TOFCorrectedHistos[i]->GetBinContent(j));
            energyCorrectedLogHistos[i]->Fill(TMath::Log10(rKE),TOFCorrectedHistos[i]->GetBinContent(j));

            //cout << "rKE = " << rKE << ", weight = " << TOFCorrectedHistos[i]->GetBinContent(j) << endl;
        }

        /*for(int j=0; j<energyCorrectedHistos[i]->GetNbinsX(); j++)
        {
            csCorrection[i].push_back(0);
            csCorrection[i][j]

        }*/

        //cout << deadtimeCorrection[i][150] << endl;
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

    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(3));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(4));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(5));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(6));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(7));

    for(int k = 0; k<monCounts.size(); k++)
    {
        cout << "target position " << k << " monitor counts = " << monCounts[k] << endl;
    }

    // switch to the detector directory
    gDirectory->cd("/");
    gDirectory->GetDirectory("detS")->cd();

    // holds the raw target-specific energy histograms in preparation for populating
    // cross-sections
    vector<TH1I*> rawHistos;
    vector<TH1I*> rawLogHistos;

    if (gDirectory->Get("blankCorrected"))
    {
        rawHistos.push_back((TH1I*)gDirectory->Get("blankCorrected"));
    }

    if (gDirectory->Get("target1Corrected"))
    {
        rawHistos.push_back((TH1I*)gDirectory->Get("target1Corrected"));
    }

    if (gDirectory->Get("target2Corrected"))
    {
        rawHistos.push_back((TH1I*)gDirectory->Get("target2Corrected"));
    }

    if (gDirectory->Get("target3Corrected"))
    {
        rawHistos.push_back((TH1I*)gDirectory->Get("target3Corrected"));
    }

    if (gDirectory->Get("target4Corrected"))
    {
        rawHistos.push_back((TH1I*)gDirectory->Get("target4Corrected"));
    }

    if (gDirectory->Get("target5Corrected"))
    {
        rawHistos.push_back((TH1I*)gDirectory->Get("target5Corrected"));
    }

    
    /*for(int k = 0; k<rawHistos.size(); k++)
    {
        cout << "noG counts " << k << " = " << rawHistos[k]->GetEntries() << endl;
    }*/

    if (gDirectory->Get("blankCorrectedLog"))
    {
        rawLogHistos.push_back((TH1I*)gDirectory->Get("blankCorrectedLog"));
    }

    if (gDirectory->Get("target1CorrectedLog"))
    {
        rawLogHistos.push_back((TH1I*)gDirectory->Get("target1CorrectedLog"));
    }

    if (gDirectory->Get("target2CorrectedLog"))
    {
        rawLogHistos.push_back((TH1I*)gDirectory->Get("target2CorrectedLog"));
    }

    if (gDirectory->Get("target3CorrectedLog"))
    {
        rawLogHistos.push_back((TH1I*)gDirectory->Get("target3CorrectedLog"));
    }

    if (gDirectory->Get("target4CorrectedLog"))
    {
        rawLogHistos.push_back((TH1I*)gDirectory->Get("target4CorrectedLog"));
    }

    if (gDirectory->Get("target5CorrectedLog"))
    {
        rawLogHistos.push_back((TH1I*)gDirectory->Get("target5CorrectedLog"));
    }

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
    for(int i=0; i<order.size(); i++)
    {
        for(int j=0; j<noBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(rawHistos[0]->GetBinContent(j) <= 0 || rawHistos[i]->GetBinContent(j) <= 0)
            {
                sigma[i][j] = 0;
            }

            else
            {
                // calculate the cross-section
                sigma[order[i]][j] = -log((rawHistos[i]->GetBinContent(j)/(double)rawHistos[0]->GetBinContent(j))*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)avo*pow(10.,-24)/(double)targetMolMass[order[i]]); // in barns

                //sigma[order[i]][j] = -log((rawHistos[i]->GetBinContent(j)/(double)rawHistos[0]->GetBinContent(j))*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)avo*pow(10.,-24)/(double)targetMolMass[order[i]]); // in barns
                //cout << "sigma = " << i << ", bin content at " << j << " = " << sigma[i][j] << endl;
            }

            if(rawLogHistos[0]->GetBinContent(j) <= 0 || rawLogHistos[i]->GetBinContent(j) <= 0)
            {
                sigmaLog[i][j] = 0;
            }

            else
            {
                // calculate the cross-section
                sigmaLog[order[i]][j] = -log((rawLogHistos[i]->GetBinContent(j)/*+scavLogHistos[i]->GetBinContent(j)*/)/((double)rawLogHistos[0]->GetBinContent(j)/*+scavLogHistos[0]->GetBinContent(j)*/)*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)avo*pow(10.,-24)/(double)targetMolMass[order[i]]); // in barns
            }
        }
    }
}

void fillCShistos()
{
    // remake the cross-section histogram file
    TFile *CSfile = new TFile(fileCSName.str().c_str(),"RECREATE");

    // declare the cross-section histograms to be filled
    TH1D *blankCS = new TH1D("blankCS","blank cross-section",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonSCS = new TH1D("carbonSCS","short carbon cross-section",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonLCS = new TH1D("carbonLCS","long carbon cross-section",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn112CS = new TH1D("Sn112CS","Sn112 cross-section",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *SnNatCS = new TH1D("SnNatCS","SnNat cross-section",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn124CS = new TH1D("Sn124CS","Sn124 cross-section",noBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    // use holder for cross-section histograms to make looping through
    // histograms easier when we calculate cross-sections below
    vector<TH1D*> csHistos;

    csHistos.push_back(blankCS);
    csHistos.push_back(carbonSCS);
    csHistos.push_back(carbonLCS);
    csHistos.push_back(Sn112CS);
    csHistos.push_back(SnNatCS);
    csHistos.push_back(Sn124CS);

    for(int i=0; i<csHistos.size(); i++)
    {
        for(int j=1; j<noBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(sigma[i][j] == 0)
            {
                continue;
            }

            else
            {
                //cout << "i = " << i << ", j = " << j << ", sigma[i][j] = " << sigma[i][j] << endl;
                csHistos[i]->SetBinContent(j,sigma[i][j]);
            }
        }
    }

    // declare the cross-section histograms to be filled
    TH1D *blankCSLog = (TH1D*)logBins(blankCS);
    TH1D *carbonSCSLog = (TH1D*)logBins(carbonSCS);
    TH1D *carbonLCSLog = (TH1D*)logBins(carbonLCS);
    TH1D *Sn112CSLog = (TH1D*)logBins(Sn112CS);
    TH1D *SnNatCSLog = (TH1D*)logBins(SnNatCS);
    TH1D *Sn124CSLog = (TH1D*)logBins(Sn124CS);

    /*TH1D *blankCSLog = new TH1D("blankCSLog","blank cross-section",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *carbonSCSLog = new TH1D("carbonSCSLog","short carbon cross-section",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *carbonLCSLog = new TH1D("carbonLCSLog","long carbon cross-section",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *Sn112CSLog = new TH1D("Sn112CSLog","Sn112 cross-section",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *SnNatCSLog = new TH1D("SnNatCSLog","SnNat cross-section",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *Sn124CSLog = new TH1D("Sn124CSLog","Sn124 cross-section",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    */
    
    // create log-scaled cross-section plots with scavenger added back in
    /*
    TH1D *blankScavLog = new TH1D("blankScavLog","blank",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *target1ScavLog = new TH1D("target1ScavLog","target1",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *target2ScavLog = new TH1D("target2ScavLog","target2",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *target3ScavLog = new TH1D("target3ScavLog","target3",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *target4ScavLog = new TH1D("target4ScavLog","target4",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    TH1D *target5ScavLog = new TH1D("target5ScavLog","target5",noBins,CS_LOWER_BOUND,TMath::Log10(CS_UPPER_BOUND));
    */
    // use holder for cross-section histograms to make looping through
    // histograms easier when we calculate cross-sections below
    vector<TH1D*> csLogHistos;

    csLogHistos.push_back(blankCSLog);
    csLogHistos.push_back(carbonSCSLog);
    csLogHistos.push_back(carbonLCSLog);
    csLogHistos.push_back(Sn112CSLog);
    csLogHistos.push_back(SnNatCSLog);
    csLogHistos.push_back(Sn124CSLog);

    for(int i=0; i<csLogHistos.size(); i++)
    {
        for(int j=1; j<noBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(sigmaLog[i][j] == 0)
            {
                continue;
            }

            else
            {
                csLogHistos[i]->SetBinContent(j,sigmaLog[i][j]);
            }
        }
    }

    // create relative cross-section plot for 112Sn/124Sn
    /*TH1D *Sn124_plus_Sn112CSLog = (TH1D*)csLogHistos[5]->Clone("Sn124_plus_Sn112CSLog");
    Sn124_plus_Sn112CSLog->Add(csLogHistos[3]);

    TH1D *Sn124_minus_Sn112CSLog = (TH1D*)csLogHistos[5]->Clone("Sn124_minus_Sn112CSLog");
    Sn124_minus_Sn112CSLog->Add(csLogHistos[3],-1);

    // rebin and scale these dummy cross-section plots (124+112 and 124-112)
    Sn124_plus_Sn112CSLog->Rebin(noBins/(double)20);
    Sn124_minus_Sn112CSLog->Rebin(noBins/(double)20);

    Sn124_plus_Sn112CSLog->Scale(1/(noBins/(double)20));
    Sn124_minus_Sn112CSLog->Scale(1/(noBins/(double)20));

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
    SnLitLog->GetXaxis()->SetRangeUser(CS_LOWER_BOUND,CS_UPPER_BOUND);

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
    carbonLitLog->GetXaxis()->SetRangeUser(CS_LOWER_BOUND,CS_UPPER_BOUND);

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
    for(int i=0; i<orchard.size(); i++)
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
                for(int k=0; k<waveform->size(); k++)
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
    for(int i=0; i<orchardW.size(); i++)
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

            for(int k=0; k<waveform->size(); k++)
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
    int offset = 0;

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
        int ch6MacroNo = macroNo;

        if(ch6microTime > 275 && ch6microTime < 325)
        {
            // used to calculate the channel 6 average value of the waveform and
            // find the baseline offset between ch6 and ch4
            int ch6Average = 0; 
            stringstream temp;
            temp << "macroNo " << macroNo << ", scavenger event number " << evtNo;
            waveformCh6= new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

            // loop through waveform data and fill histo
            for(int k=0; k<waveform->size(); k++)
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

                            for(int k=0; k<waveform->size(); k++)
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
    TTree* ch6Tree = (TTree*)file->Get("ch6ProcessedTree");
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
    fillCShistos();

    // print out waveforms for first 50 events that occur in both ch4 and ch6
    //matchWaveforms();
        
    fileOut->Close();
}
