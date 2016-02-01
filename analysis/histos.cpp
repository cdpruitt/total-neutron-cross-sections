#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TMath.h"

using namespace std;

// scavenger event stream
ofstream scavengerEvents;

// summedDet event stream
ofstream summedDetEvents;


/* Experimental constants */

const double FLIGHT_DISTANCE = 2580; // detector distance from source, in cm
// 2080 for Neon  

// Target changer time delay after the macropulse start time as calculated using
// the summed detector signals' gamma peaks
const double MACROPULSE_OFFSET = 840; // in ns
// 814 for Neon

// Time correction to MACROPULSE_OFFSET for the monitor; takes into account the
// additional ~10 meters between the summed detector and the extra cable length
const double MONITOR_OFFSET = 0; // in ns

// Time correction to scavenger times, relative to summed det
const double SCAVENGER_OFFSET = -12.5; // in ns

// Period of micropulses
const double MICRO_PERIOD = 1788.820; // in ns

// Physical constants
const double C = 299792458; // speed of light in m/s

const double NEUTRON_MASS = 939.56536; // in MeV/c^2

const string analysispath =  "/media/Drive3/";


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

vector<vector<TH1I*>> histos; // holds all histograms for the run
// split into sub-vectors on a per-channel basis

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
// data from ch4Tree (summed detector signal)
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

    TH1I *secondSTOFblank = new TH1I("secondSTOFblank","second in micro time of flight, blank",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFcs = new TH1I("secondSTOFcs","second in micro time of flight, carbon s",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFcl = new TH1I("secondSTOFcl","second in micro time of flight, carbon l",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFsn112 = new TH1I("secondSTOFsn112","second in micro time of flight, sn112",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFsnnat = new TH1I("secondSTOFsnnat","second in micro time of flight, snnat",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOFsn124 = new TH1I("secondSTOFsn124","second in micro time of flight, sn124",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);

    TH1I *thirdSTOFblank = new TH1I("thirdSTOFblank","third in micro time of flight, blank",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFcs = new TH1I("thirdSTOFcs","third in micro time of flight, carbon s",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFcl = new TH1I("thirdSTOFcl","third in micro time of flight, carbon l",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFsn112 = new TH1I("thirdSTOFsn112","third in micro time of flight, sn112",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFsnnat = new TH1I("thirdSTOFsnnat","third in micro time of flight, snnat",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOFsn124 = new TH1I("thirdSTOFsn124","third in micro time of flight, sn124",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    */

    TH1I *TOF = new TH1I("TOF","Summed-detector time of flight",1800,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",1800,0,MICRO_PERIOD+1,2048,0,65536);

    TH1I* microNoH = new TH1I("microNoH","microNo",360,0,360);
    microNoH->GetXaxis()->SetTitle("micropulse number of each event");

    TH1I *firstInMicro = new TH1I("firstInMicro","first in micro time of flight",1800,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondInMicro = new TH1I("secondInMicro","second in micro time of flight",1800,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdInMicro = new TH1I("thirdInMicro","third in micro time of flight",1800,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);

    // create cross-section plots
    TH1I *blankRaw = new TH1I("blank","blank",500,0,700);
    TH1I *carbonSRaw = new TH1I("carbonS","carbonS",500,0,700);
    TH1I *carbonLRaw = new TH1I("carbonL","carbonL",500,0,700);
    TH1I *Sn112Raw = new TH1I("Sn112","Sn112",500,0,700);
    TH1I *NatSnRaw = new TH1I("NatSn","NatSn",500,0,700);
    TH1I *Sn124Raw = new TH1I("Sn124","Sn124",500,0,700);
    TH1I *totalRaw = new TH1I("total","total",500,0,700);

    // create log-scaled cross-section plots
    TH1I *blankRawLog = new TH1I("blankLog","blank",500,0,TMath::Log10(700));
    TH1I *carbonSRawLog = new TH1I("carbonSLog","carbonS",500,0,TMath::Log10(700));
    TH1I *carbonLRawLog = new TH1I("carbonLLog","carbonL",500,0,TMath::Log10(700));
    TH1I *Sn112RawLog = new TH1I("Sn112Log","Sn112",500,0,TMath::Log10(700));
    TH1I *NatSnRawLog = new TH1I("NatSnLog","NatSn",500,0,TMath::Log10(700));
    TH1I *Sn124RawLog = new TH1I("Sn124Log","Sn124",500,0,TMath::Log10(700));
    TH1I *totalRawLog = new TH1I("totalLog","total",500,0,TMath::Log10(700));

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
    double TIME_OFFSET;
    int gammaGate[2];

    // adjust time parameters based on channel identity
    switch(i)
    {
        case 1:
            // monitor
            TIME_OFFSET = MACROPULSE_OFFSET+MONITOR_OFFSET;
            gammaGate[1] = 40;
            gammaGate[2] = 50;
            break;
        case 2:
            // summed detector
            TIME_OFFSET = MACROPULSE_OFFSET;
            gammaGate[1] = 80;
            gammaGate[2] = 92;
            break;
        case 3:
            // scavenger
            TIME_OFFSET = MACROPULSE_OFFSET+SCAVENGER_OFFSET;
            gammaGate[1] = 80;
            gammaGate[2] = 92;
            break;
    }

    int totalEntries = orchard[i]->GetEntries();
    cout << "Populating advanced histograms for channel " << 2*i << endl;

    for(int j=0; j<totalEntries; j++)
    {
        orchard[i]->GetEntry(j);

        double timeDiff = completeTime-macroTime+TIME_OFFSET;

        // only plot events within the beam-on period for each macropulse, and
        // throw away events when the macropulse is between target changer
        // positions
        if (timeDiff < 650000 && timeDiff > 0 && targetPos != 0)
        {
            // valid event within the beam-on period of each macropulse
            // calculate the correct micropulse-referenced time and micropulse
            // number

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
            if (microTime < gammaGate[1] && microTime > gammaGate[0])
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

            if (!gammaInMicro) // gate disallowing gammas
            {
                switch (targetPos)
                {
                    case 1:
                        // BLANK
                        blankRaw->Fill(rKE);
                        //blankRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 2:
                        // SHORT CARBON
                        carbonSRaw->Fill(rKE);
                        //carbonSRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 3:
                        // LONG CARBON
                        carbonLRaw->Fill(rKE);
                        //carbonLRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 4:
                        // Sn112
                        Sn112Raw->Fill(rKE);
                        //Sn112RawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 5:
                        // Natural Sn
                        NatSnRaw->Fill(rKE);
                        //NatSnRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 6:
                        // Sn124
                        Sn124Raw->Fill(rKE);
                        //Sn124RawLog->Fill(TMath::Log10(rKE));
                        break;

                    default:
                        break;
                }

                totalRaw->Fill(rKE);
                //totalRawLog->Fill(TMath::Log10(rKE));
            }
        }
    }
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
            //microNoH->Fill(microNo);
            evtNoH->Fill(evtNo);
            macroTimeH->Fill(macroTime);
            completeTimeH->Fill(completeTime);
            //microTimeH->Fill(microTime);
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

    // fill TOF, cross-section, etc. histos for channels 2, 4, 6
    for(int i=1; i<4; i++)
    {
        fillAdvancedHistos(i);
    }

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

    cout << "Searching through orchard[3] for events in the overlap region" << endl;

    gDirectory->cd("/");

    // keep track of how many histograms have been drawn; stop looping after a
    // set number have been filled (i.e. first 50)
    int plots = 0;

    setBranches(orchard[3]); // pointed at the ch6Tree now

    // keep track of how far we've looped through the ch4 tree so we don't have
    // to loop from the beginning each time we're looking for a parallel event
    int currentIndex = 0;

    // loop through the channel-specific tree and populate histos
    for(int j=0; j<ch6Entries && plots<=1000 /*plot only first 50 */; j++)
    {
        orchard[3]->GetEntry(j);

        double ch6timeDiff = completeTime-macroTime+MACROPULSE_OFFSET+SCAVENGER_OFFSET;
        double ch6microTime = fmod(ch6timeDiff,MICRO_PERIOD);
        int ch6microNo = floor(ch6timeDiff/MICRO_PERIOD);

        // keep track of this event's macroNo to compare it with events in ch4;
        // this will help us identify the parallel event in ch4 with the same
        // timestamp
        int ch6MacroNo = macroNo;

        if(ch6microTime > 275 && ch6microTime < 325)
        {
            stringstream temp;
            temp << "macroNo " << macroNo << ", scavenger event number " << evtNo;
            waveformCh6= new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

            // loop through waveform data and fill histo
            for(int k=0; k<waveform->size(); k++)
            {
                waveformCh6->SetBinContent(k,waveform->at(k));
            }

            // loop through ch4 tree to find the same event (same timestamp and
            // macropulse as the ch6 event)

            setBranches(orchard[2]); // pointed at the ch6Tree now

            for(int i=currentIndex; i<ch4Entries; i++)
            {
                orchard[2]->GetEntry(i);

                if(macroNo == ch6MacroNo)
                {
                    // in the correct macropulse - loop through ch4 to find the
                    // the event whose timestamp matches the one in channel 6
                    double ch4timeDiff = completeTime-macroTime+MACROPULSE_OFFSET;
                    double ch4microTime = fmod(ch4timeDiff,MICRO_PERIOD);
                    int microNo = floor(ch4timeDiff/MICRO_PERIOD);

                    if(microNo == ch6microNo && (ch4microTime >= ch6microTime-5 && ch4microTime <= ch6microTime+5))
                    {
                        // found the parallel event in ch4 - same time, same
                        // micropulse
                        //
                        // plot both the ch6 and ch4 events so we can plot the
                        // histograms together later

                        ch4ch6timeDiff->Fill(ch6microTime-ch4microTime);
                        cout << "ch4microTime = " << ch4microTime << ", ch4microNo = " << microNo << ", ch6microTime = " << ch6microTime << ", ch6microNo = " << ch6microNo << endl;

                        // clear the histo title stringstream
                        temp.str("");
                        temp << "macroNo " << macroNo << ", ch4 event number " << evtNo;
                        waveformCh4 = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

                        for(int k=0; k<waveform->size(); k++)
                        {
                            waveformCh4->SetBinContent(k,waveform->at(k));
                        }

                        // Draw both histos on same axes
                        //waveformCh6->Draw();
                        //waveformCh4->Draw("same");

                        plots++;

                        break;
                    }
                }

                else if (macroNo > ch6MacroNo)
                {
                    // went past ch4 macropulse without finding the event
                    // parallel to ch6
                    plots++;
                    break;
                }

                currentIndex = i;
            }
            // point back to ch6 events in preparation for next loop through
            // ch6Tree
            setBranches(orchard[3]);
        }
    }
    cout << endl;
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
    stringstream fileInName;
    stringstream fileOutName;
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
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"RECREATE");

    // print out waveforms for first 50 events that occur in both ch4 and ch6
    matchWaveforms();

    // prepare the root file with 4 directories, one for each channel
    // these directories will hold basic variable histograms showing the
    // raw data in each tree, plus TOF, x-sections, etc histograms
//    fillHistos();

    fileOut->Write();

    fileOut->Close();
}
