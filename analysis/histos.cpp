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

using namespace std;

// scavenger event stream
ofstream scavengerEvents;

// summedDet event stream
ofstream summedDetEvents;

// fileIn/fileOut names to be accessed to open files
stringstream fileInName, fileOutName, fileCSName;

// path to sorted tree data
const string analysispath =  "/media/Drive3/";



/* Experimental constants */

const double FLIGHT_DISTANCE = 2672; // detector distance from source, in cm
// 2080 for Neon  

// Period of micropulses
const double MICRO_PERIOD = 1788.820; // in ns

// Physical constants
const double C = 299792458; // speed of light in m/s

const double NEUTRON_MASS = 939.56536; // in MeV/c^2

// target ordering:
// {blank, Sn112, Nat. Sn, Sn124, short carbon, long carbon} 

const int noTargets = 4;

// lengths of each target, respectively:
double targetlength[6] = {0,1.365,1.370,1.370,1.37,2.74}; //cm

// molar mass of each target, respectively
double targetMolMass[6] = {0,112,118.7,124,12.01,12.01}; //g/mol

// density of each target, respectively
double targetdensity[6] = {0,6.89,7.31,7.63,2.3,2.3}; //g/cm^3

// Avogadro's number
double avo = 6.022*pow(10.,23.); //atoms/mol

// number of bins in the raw energy histograms and in the cross-section
// histograms
const int noBins = 1000;

// holder for the scaled cross-section of each target; i = target #, j = bin
double sigma[6][noBins] = {0};
double sigmaLog[6][noBins] = {0};



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



// Re-link to an already-existing tree's data so we can read the tree
void setBranches(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    //tree->SetBranchAddress("microNo",&microNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("macroTime",&macroTime);
    tree->SetBranchAddress("completeTime",&completeTime);
    //tree->SetBranchAddress("microTime",&microTime);
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
    tree->SetBranchAddress("waveform",&waveform);
}

// Populate advanced histograms (TOFs, cross-section, etc) calculated using
// data from ch4 and ch6 trees
void fillAdvancedHistos(int i)
{
    gDirectory->cd("/");
    gDirectory->GetDirectory(dirs[i].c_str())->cd();

    /*    TH1I *firstSTOFblank = new TH1I("firstSTOFblank","first in micro time of flight, blank",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFcs = new TH1I("firstSTOFcs","first in micro time of flight, carbon s",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFcl = new TH1I("firstSTOFcl","first in micro time of flight, carbon l",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFsn112 = new TH1I("firstSTOFsn112","first in micro time of flight, sn112",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFsnnat = new TH1I("firstSTOFsnnat","first in micro time of flight, snnat",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *firstSTOFsn124 = new TH1I("firstSTOFsn124","first in micro time of flight, sn124",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    */

    TH1I *TOF = new TH1I("TOF","Summed-detector time of flight",1800,0,MICRO_PERIOD*1.05);
    TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",1800,0,MICRO_PERIOD+1,2048,0,65536);

    TH1I* microNoH = new TH1I("microNoH","microNo",360,0,360);
    microNoH->GetXaxis()->SetTitle("micropulse number of each event");

    TH1I *firstInMicro = new TH1I("firstInMicro","first in micro time of flight",1800,0,MICRO_PERIOD*1.05);
    TH1I *secondInMicro = new TH1I("secondInMicro","second in micro time of flight",1800,0,MICRO_PERIOD*1.05);
    TH1I *thirdInMicro = new TH1I("thirdInMicro","third in micro time of flight",1800,0,MICRO_PERIOD*1.05);

    TH1I *fimBlank = new TH1I("fimBlank","first in micro, blank",1800,0,MICRO_PERIOD*1.05);
    TH1I *fimTarget1 = new TH1I("fimTarget1","first in micro, target 1",1800,0,MICRO_PERIOD*1.05);
    TH1I *fimTarget2 = new TH1I("fimTarget2","first in micro, target 2",1800,0,MICRO_PERIOD*1.05);
    TH1I *fimTarget3 = new TH1I("fimTarget3","first in micro, target 3",1800,0,MICRO_PERIOD*1.05);
    TH1I *fimTarget4 = new TH1I("fimTarget4","first in micro, target 4",1800,0,MICRO_PERIOD*1.05);
    TH1I *fimTarget5 = new TH1I("fimTarget5","first in micro, target 5",1800,0,MICRO_PERIOD*1.05);

    TH1I *noGBlank = new TH1I("noGBlank","no gamma in micro, blank",noBins,0,700);
    TH1I *noGTarget1 = new TH1I("noGTarget1","no gamma in micro, target 1",noBins,0,700);
    TH1I *noGTarget2 = new TH1I("noGTarget2","no gamma in micro, target 2",noBins,0,700);
    TH1I *noGTarget3 = new TH1I("noGTarget3","no gamma in micro, target 3",noBins,0,700);
    TH1I *noGTarget4 = new TH1I("noGTarget4","no gamma in micro, target 4",noBins,0,700);
    TH1I *noGTarget5 = new TH1I("noGTarget5","no gamma in micro, target 5",noBins,0,700);

    TH1I *noGBlankLog = new TH1I("noGBlankLog","no gamma in micro, log E, blank",noBins,0,TMath::Log10(700));
    TH1I *noGTarget1Log = new TH1I("noGTarget1Log","no gamma in micro, log E, target 1",noBins,0,TMath::Log10(700));
    TH1I *noGTarget2Log = new TH1I("noGTarget2Log","no gamma in micro, log E, target 2",noBins,0,TMath::Log10(700));
    TH1I *noGTarget3Log = new TH1I("noGTarget3Log","no gamma in micro, log E, target 3",noBins,0,TMath::Log10(700));
    TH1I *noGTarget4Log = new TH1I("noGTarget4Log","no gamma in micro, log E, target 4",noBins,0,TMath::Log10(700));
    TH1I *noGTarget5Log = new TH1I("noGTarget5Log","no gamma in micro, log E, target 5",noBins,0,TMath::Log10(700));

    // create cross-section plots
    TH1I *blankRaw = new TH1I("blank","blank",noBins,0,700);
    TH1I *target1Raw = new TH1I("target1","target1",noBins,0,700);
    TH1I *target2Raw = new TH1I("target2","target2",noBins,0,700);
    TH1I *target3Raw = new TH1I("target3","target3",noBins,0,700);
    TH1I *target4Raw = new TH1I("target4","target4",noBins,0,700);
    TH1I *target5Raw = new TH1I("target5","target5",noBins,0,700);
    TH1I *totalRaw = new TH1I("total","total",noBins,0,700);

    // create log-scaled cross-section plots
    TH1I *blankRawLog = new TH1I("blankLog","blank",noBins,0,TMath::Log10(700));
    TH1I *target1RawLog = new TH1I("target1Log","target1",noBins,0,TMath::Log10(700));
    TH1I *target2RawLog = new TH1I("target2Log","target2",noBins,0,TMath::Log10(700));
    TH1I *target3RawLog = new TH1I("target3Log","target3",noBins,0,TMath::Log10(700));
    TH1I *target4RawLog = new TH1I("target4Log","target4",noBins,0,TMath::Log10(700));
    TH1I *target5RawLog = new TH1I("target5Log","target5",noBins,0,TMath::Log10(700));
    TH1I *totalRawLog = new TH1I("totalLog","total",noBins,0,TMath::Log10(700));

    double microTime;
    int microNo, prevMicroNo;

    // create variables to gate on whether a micropulse has a gamma at the start
    // and keep track of the influence of early micropulse events on later ones
    int orderInMicro = 0;
    bool isGamma = false;
    bool gammaInMicro = false;

    // point at correct tree in preparation for reading data
    setBranches(orchard[i]);

    // channel-dependent time offset relative to the target changer's macropulse
    // start time
    int gammaGate[2];

    // adjust time parameters based on channel identity
    switch(i)
    {
        case 1:
            // monitor
            gammaGate[0] = 25;
            gammaGate[1] = 40;
            break;
        case 2:
            // summed detector
            gammaGate[0] = 84;
            gammaGate[1] = 95;
            break;
        case 3:
            // scavenger
            gammaGate[0] = 84;
            gammaGate[1] = 95;
            break;
    }

    int totalEntries = orchard[i]->GetEntries();
    cout << "Populating advanced histograms for channel " << 2*i << endl;

    for(int j=0; j<totalEntries; j++)
    {
        orchard[i]->GetEntry(j);

        double timeDiff = completeTime-macroTime;

        // only plot events within the beam-on period for each macropulse, and
        // throw away events when the macropulse is between target changer
        // positions
        if (timeDiff < 650000 && timeDiff > 0 && targetPos != 0)
        {
            // valid event within the beam-on period of each macropulse
            // calculate the correct micropulse-referenced time and micropulse
            // number

            // if a monitor event, fill monitor counter to allow for relative
            // scaling of cross-sections
            /*if(i==1)
              {
            // in the monitor channel; fill a counter for monitor events
            // so we can scale the cross-sections relative to each other
            //monCounts[0] is blank, monCounts[1] is target 1, etc.
            monCounts[targetPos-1]++;
            }*/

            microTime = fmod(timeDiff,MICRO_PERIOD);

            TOF->Fill(microTime);
            triangle->Fill(microTime,lgQ);

            // convert microTime into neutron velocity based on flight path distance
            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/microTime; // in meters/sec 

            // convert velocity to relativistic kinetic energy
            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            // assign a micropulse number to the event
            prevMicroNo = microNo;
            microNo = floor(timeDiff/MICRO_PERIOD);

            // assign an ordering of events in each micropulse
            if (microNo==prevMicroNo)
            {
                // still in same micropulse
                orderInMicro++;
            }

            else
            {
                // new micropulse
                orderInMicro = 1;

                // reset gamma indicator
                gammaInMicro = false;
            }

            // now check to see whether event is a gamma
            if (microTime>gammaGate[0] && microTime<gammaGate[1])
            {
                isGamma = true;
                // indicate that this micro started with a gamma
                gammaInMicro = true;
            }

            else
            {
                isGamma = false;
            }

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

            //if (!gammaInMicro) // gate disallowing gammas in detector channels
            //{
            // to scale the target's cross-sections properly, we need to
            // know how many events we're throwing out by using the
            // gammaInMicro gate.
            // So we should keep track of how many gammaInMicro events there
            // are on a target-by-target basis, and divide each
            // cross-section bin with them (just as we do with monitor
            // counts), to scale bins properly

            switch (targetPos)
            {
                case 1:
                    // BLANK
                    blankRaw->Fill(rKE);
                    blankRawLog->Fill(TMath::Log10(rKE));

                    if(orderInMicro==1)
                    {
                        fimBlank->Fill(microTime);
                    }

                    if(!gammaInMicro)
                    {
                        noGBlank->Fill(rKE);
                        noGBlankLog->Fill(TMath::Log10(rKE));
                    }
                    break;

                case 2:
                    // TARGET 1
                    target1Raw->Fill(rKE);
                    target1RawLog->Fill(TMath::Log10(rKE));

                    if(orderInMicro==1)
                    {
                        fimTarget1->Fill(microTime);
                    }

                    if(!gammaInMicro)
                    {
                        noGTarget1->Fill(rKE);
                        noGTarget1Log->Fill(TMath::Log10(rKE));
                    }
                    break;

                case 3:
                    // TARGET 2
                    target2Raw->Fill(rKE);
                    target2RawLog->Fill(TMath::Log10(rKE));

                    if(orderInMicro==1)
                    {
                        fimTarget2->Fill(microTime);
                    }

                    if(!gammaInMicro)
                    {
                        noGTarget2->Fill(rKE);
                        noGTarget2Log->Fill(TMath::Log10(rKE));
                    }
                    break;

                case 4:
                    // TARGET 3
                    target3Raw->Fill(rKE);
                    target3RawLog->Fill(TMath::Log10(rKE));

                    if(orderInMicro==1)
                    {
                        fimTarget3->Fill(microTime);
                    }

                    if(!gammaInMicro)
                    {
                        noGTarget3->Fill(rKE);
                        noGTarget3Log->Fill(TMath::Log10(rKE));
                    }
                    break;

                case 5:
                    // TARGET 4
                    target4Raw->Fill(rKE);
                    target4RawLog->Fill(TMath::Log10(rKE));

                    if(orderInMicro==1)
                    {
                        fimTarget4->Fill(microTime);
                    }

                    if(!gammaInMicro)
                    {
                        noGTarget4->Fill(rKE);
                        noGTarget4Log->Fill(TMath::Log10(rKE));
                    }
                    break;

                case 6:
                    // TARGET 5
                    target5Raw->Fill(rKE);
                    target5RawLog->Fill(TMath::Log10(rKE));

                    if(orderInMicro==1)
                    {
                        fimTarget5->Fill(microTime);
                    }

                    if(!gammaInMicro)
                    {
                        noGTarget5->Fill(rKE);
                        noGTarget5Log->Fill(TMath::Log10(rKE));
                    }
                    break;

                default:
                    break;
            }

            totalRaw->Fill(rKE);
            totalRawLog->Fill(TMath::Log10(rKE));
            //}
        }
    }
}

void calculateCS()
{
    // Find number of events in the monitor for each target to use in scaling
    // cross-sections
    gDirectory->cd("/");

    // switch to the monitor directory
    gDirectory->GetDirectory("monitor")->cd();

    // to normalize flux between all channels, we'll need to keep track of the
    // number of counts that came in through the monitor during that target's
    // beam time
    vector<int> monCounts;

    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(3));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(4));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(5));
    monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(6));

    for(int k = 0; k<monCounts.size(); k++)
    {
        cout << "mon counts " << k << " = " << monCounts[k] << endl;
    }

    gDirectory->cd("/");

    // switch to the detector directory
    gDirectory->GetDirectory("detS")->cd();

    // holds the raw target-specific energy histograms in preparation for populating
    // cross-sections
    vector<TH1I*> rawHistos;
    vector<TH1I*> rawLogHistos;

    rawHistos.push_back((TH1I*)gDirectory->Get("noGBlank"));
    rawHistos.push_back((TH1I*)gDirectory->Get("noGTarget1"));
    rawHistos.push_back((TH1I*)gDirectory->Get("noGTarget2"));
    rawHistos.push_back((TH1I*)gDirectory->Get("noGTarget3"));
    rawHistos.push_back((TH1I*)gDirectory->Get("noGTarget4"));
    //rawHistos.push_back(TH1I*)gDirectory->Get("target5"));

    for(int k = 0; k<rawHistos.size(); k++)
    {
        cout << "noG counts " << k << " = " << rawHistos[k]->GetEntries() << endl;
    }

    rawLogHistos.push_back((TH1I*)gDirectory->Get("noGBlankLog"));
    rawLogHistos.push_back((TH1I*)gDirectory->Get("noGTarget1Log"));
    rawLogHistos.push_back((TH1I*)gDirectory->Get("noGTarget2Log"));
    rawLogHistos.push_back((TH1I*)gDirectory->Get("noGTarget3Log"));
    rawLogHistos.push_back((TH1I*)gDirectory->Get("noGTarget4Log"));
    //rawLogHistos.push_back((TH1I*)gDirectory->Get("target5Log"));

    for(int i=0; i<=noTargets; i++)
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
                // we must have found positive-definite values found for the raw
                // histogram bins in questions; calculate the cross-section and
                // fill the relevant csHisto

                sigma[i][j] = -log((rawHistos[i]->GetBinContent(j)/(double)rawHistos[0]->GetBinContent(j))*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[i]*(double)targetdensity[i]*(double)avo*pow(10.,-24)/(double)targetMolMass[i]); // in barns
                //cout << "sigma = " << i << ", bin content at " << j << " = " << sigma[i][j] << endl;
            }
        }
    }

    for(int i=0; i<=noTargets; i++)
    {
        for(int j=0; j<noBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(rawLogHistos[0]->GetBinContent(j) <= 0 || rawLogHistos[i]->GetBinContent(j) <= 0)
            {
                sigmaLog[i][j] = 0;
            }

            else
            {
                // we must have found positive-definite values found for the rawLog
                // histogram bins in questions; calculate the cross-section and
                // fill the relevant csHisto

                sigmaLog[i][j] = -log((rawLogHistos[i]->GetBinContent(j)/(double)rawLogHistos[0]->GetBinContent(j))*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[i]*(double)targetdensity[i]*(double)avo*pow(10.,-24)/(double)targetMolMass[i]); // in barns
            }
            //cout << "sigmaLog = " << i << ", bin content at " << j << " = " << sigmaLog[i][j] << endl;
        }
    }
}

void fillCShistos()
{
    // remake the cross-section histogram file
    TFile *CSfile = new TFile(fileCSName.str().c_str(),"RECREATE");

    // declare the cross-section histograms to be filled
    TH1D *blankcs = new TH1D("blankcs","blank cross-section",noBins,0,700);
    TH1D *target1cs = new TH1D("target1cs","target 1 cross-section",noBins,0,700);
    TH1D *target2cs = new TH1D("target2cs","target 2 cross-section",noBins,0,700);
    TH1D *target3cs = new TH1D("target3cs","target 3 cross-section",noBins,0,700);
    TH1D *target4cs = new TH1D("target4cs","target 4 cross-section",noBins,0,700);
    //TH1D *target5cs = new TH1D("target5cs","target 5 cross-section",noBins,0,TMath::Log(10700);

    // use holder for cross-section histograms to make looping through
    // histograms easier when we calculate cross-sections below
    vector<TH1D*> csHistos;

    csHistos.push_back(blankcs);
    csHistos.push_back(target1cs);
    csHistos.push_back(target2cs);
    csHistos.push_back(target3cs);
    csHistos.push_back(target4cs);
    //csHistos.push_back(target5cs);

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
    TH1D *blankcsLog = new TH1D("blankcsLog","blank cross-section",noBins,0,TMath::Log10(700));
    TH1D *target1csLog = new TH1D("target1csLog","target 1 cross-section",noBins,0,TMath::Log10(700));
    TH1D *target2csLog = new TH1D("target2csLog","target 2 cross-section",noBins,0,TMath::Log10(700));
    TH1D *target3csLog = new TH1D("target3csLog","target 3 cross-section",noBins,0,TMath::Log10(700));
    TH1D *target4csLog = new TH1D("target4csLog","target 4 cross-section",noBins,0,TMath::Log10(700));
    //TH1D *target5csLog = new TH1D("target5csLog","target 5 cross-section",noBins,0,TMath::Log10(700));

    // use holder for cross-section histograms to make looping through
    // histograms easier when we calculate cross-sections below
    vector<TH1D*> csLogHistos;

    csLogHistos.push_back(blankcsLog);
    csLogHistos.push_back(target1csLog);
    csLogHistos.push_back(target2csLog);
    csLogHistos.push_back(target3csLog);
    csLogHistos.push_back(target4csLog);
    //csLogHistos.push_back(target5csLog);

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

    ifstream SnData("/home/wudaq/WashUDAQ/analysis/SnNatData.dat");
    if(!SnData.is_open())
    {
        cout << "No Previous Data..." << endl;
        return;
    }

    char dummy[200];
    SnData.getline(dummy,200);

    vector<float> energy;
    vector<float> xsection;
    vector<float> error;

    float dum,dum2,dum3;

    while(!SnData.eof())
    {
        SnData >> dum >> dum2 >> dum3;

        energy.push_back(TMath::Log10(dum));
        xsection.push_back(dum2);
        error.push_back(dum3);
    }

    TGraphErrors *SnLitLog = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
    //SnLitLog->Draw("AP");

    SnLitLog->GetXaxis()->SetTitle("Energy [MeV] (log10)");
    SnLitLog->GetXaxis()->CenterTitle();
    SnLitLog->GetXaxis()->SetRangeUser(0,TMath::Log10(700.));

    SnLitLog->GetYaxis()->SetTitle("sigma [b]");
    SnLitLog->GetYaxis()->CenterTitle();
    //SnLitLog->GetYaxis()->SetTitleOffSet(1.5);
    SnLitLog->Write();

    // carbon literature data
    ifstream carbonData("/home/wudaq/WashUDAQ/analysis/CarbonData.dat");
    if(!carbonData.is_open())
    {
        cout << "No Previous Data..." << endl;
        return;
    }

    carbonData.getline(dummy,200);

    energy.clear();
    xsection.clear();
    error.clear();

    while(!carbonData.eof())
    {
        carbonData >> dum >> dum2 >> dum3;

        energy.push_back(TMath::Log10(dum));
        xsection.push_back(dum2);
        error.push_back(dum3);
    }

    TGraphErrors *carbonLitLog = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
    carbonLitLog->Draw("AP");

    carbonLitLog->GetXaxis()->SetTitle("Energy [MeV] (log10)");
    carbonLitLog->GetXaxis()->CenterTitle();
    carbonLitLog->GetXaxis()->SetRangeUser(0,TMath::Log10(700.));

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

        TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",6000,0,6000000000);
        completeTimeH->GetXaxis()->SetTitle("complete time for each event");

        //TH1I* microTimeH = new TH1I("microTimeH","microTime",2000,0,2000);
        //microTimeH->GetXaxis()->SetTitle("time since start of micro for each event");

        TH1I* targetPosH = new TH1I("targetPosH","targetPos",6,0,6);
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
        cout << "Populating " << dirs[i] << " histograms..." << endl;

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

            macroNoH->Fill(macroNo);
            evtNoH->Fill(evtNo);
            macroTimeH->Fill(macroTime);
            completeTimeH->Fill(completeTime);
            targetPosH->Fill(targetPos);
            sgQH->Fill(sgQ);
            lgQH->Fill(lgQ);

            // if waveform data for this event exist, we want to populate
            // a histogram to display it

            // only plot 1 out of 10000 waveforms to save space and processing
            // time
            if(waveform->size() > 0 && j%10000 == 0)
            {
                waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
                waveformsDir->cd();

                stringstream temp;
                temp << "macroNo " << macroNo << ", evtNo " << evtNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

                // loop through waveform data and fill histo
                for(int k=0; k<waveform->size(); k++)
                {
                    waveformH->SetBinContent(k,waveform->at(k));
                }
                gDirectory->cd("..");
            }
        }
    }

    // fill TOF, cross-section, etc. histos for channels 4, 6
    fillAdvancedHistos(2);
    fillAdvancedHistos(3);

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
        gDirectory->mkdir("waveformsDir","raw DPP waveforms");
        waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
        waveformsDir->cd();

        // fill waveform mode histograms
        setBranchesW(orchardW[i]);
        int totalEntries = orchardW[i]->GetEntries();
        cout << "Populating " << dirs[i] << " histograms..." << endl;

        // create a holder for the time-zero of each macropulse waveform
        // so we can plot all 11 waveform chunks on the same histogram
        // initialize to -650000 to ensure that the first event in the waveform
        // trees is always considered the start of a macropulse waveform
        double waveformStart = -650000;

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

            if(completeTime >= waveformStart+650000)
            {
                // new macropulse in waveform mode because no other waveform
                // data came before it in the last macropulse period
                // create a new plot to stitch 11 waveforms together
                stringstream temp;
                temp << "full waveform " << waveformNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),350000,0,700000);

                // set the start of the macropulse to the first event timer
                waveformStart = completeTime;
                waveformNo++;
            }

            for(int k=0; k<waveform->size(); k++)
            {
                waveformH->SetBinContent(k+(completeTime-waveformStart)/2,waveform->at(k));
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

int main(int argc, char* argv[])
{
    cout.precision(13);
    // needed to avoid ROOT error for header files being incorrectly brought in
    // in both resort.cpp and histos.cpp
    // Look online for more info (I'm not really sure why it's necessary)
    TApplication app("app",&argc,argv);

    // read in the raw file name
    string runDir = argv[1];
    string runNo = argv[2];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    stringstream treeName;
    stringstream scavengerEventsName;
    stringstream summedDetEventsName;

    treeName << runDir << "-" << runNo; 

    scavengerEventsName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_scavenger.csv";
    summedDetEventsName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_summedDet.csv";

    scavengerEvents.open(scavengerEventsName.str());
    summedDetEvents.open(summedDetEventsName.str());

    scavengerEvents.precision(10);
    summedDetEvents.precision(10);

    fileInName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_sorted.root";
    fileOutName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_histos.root";
    fileCSName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << "_cross-sections.root";

    TFile* file = new TFile(fileInName.str().c_str(),"READ");

    if(file->Get("ch0Tree"))
    {
        cout << "Located channel-specific trees in " << fileInName << "." << endl;
    }

    TTree* ch0Tree = (TTree*)file->Get("ch0Tree");
    TTree* ch2Tree = (TTree*)file->Get("ch2Tree");
    TTree* ch4Tree = (TTree*)file->Get("ch4Tree");
    TTree* ch6Tree = (TTree*)file->Get("ch6Tree");
    TTree* ch0TreeW = (TTree*)file->Get("ch0TreeW");
    TTree* ch2TreeW = (TTree*)file->Get("ch2TreeW");
    TTree* ch4TreeW = (TTree*)file->Get("ch4TreeW");

    orchard.push_back(ch0Tree);
    orchard.push_back(ch2Tree);
    orchard.push_back(ch4Tree);
    orchard.push_back(ch6Tree);

    orchardW.push_back(ch0TreeW);
    orchardW.push_back(ch2TreeW);
    orchardW.push_back(ch4TreeW);
    
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
