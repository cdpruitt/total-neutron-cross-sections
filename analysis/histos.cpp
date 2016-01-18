#include <iostream>
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


/* Experimental constants */

const double FLIGHT_DISTANCE = 2580; // detector distance from source, in cm
// 2080 for Neon  

// Time delay of target changer coarse time after real macropulse start time
const double TIME_OFFSET = 836; // in ns
// 814 for Neon

// Period of micropulses
const double MICRO_PERIOD = 1788.820; // in ns

// Physical constants
const double C = 299792458; // speed of light in m/s

const double NEUTRON_MASS = 939.56536; // in MeV/c^2

const string analysispath =  "/media/Drive3/";


/* event variables */

// variables for holding tree event data to be histogrammed
unsigned int evtNo, macroNo, microNo, targetPos;
unsigned int sgQ, lgQ;
double completeTime, macroTime, trueTime;

vector<int> *waveform; // for holding one event's waveform data


/* ROOT and organizational variables */

// ROOT file directory structure 
string dirs[4] = {"targetChanger","monitor","detS","scavenger"};

vector<TTree*> orchard; // holds all channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

vector<vector<TH1I*>> histos; // holds all histograms for the run
// split into sub-vectors on a per-channel basis

TH1I* DPPWaveform;

TDirectory *DPPWaveformsDir;
TDirectory *WaveWaveformsDir;


// Re-link to an already-existing tree's data so we can read the tree
void setBranches(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("microNo",&microNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("macroTime",&macroTime);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("trueTime",&trueTime);
    tree->SetBranchAddress("targetPos",&targetPos);
    tree->SetBranchAddress("sgQ",&sgQ);
    tree->SetBranchAddress("lgQ",&lgQ);
    tree->SetBranchAddress("waveform",&waveform);
}

// Populate advanced histograms (TOFs, cross-section, etc) calculated using
// data from ch4Tree (summed detector signal)
void fillAdvanced(int i)
{
    gDirectory->cd("/");
    gDirectory->GetDirectory(dirs[i].c_str())->cd();

    TH1I *blankRaw = new TH1I("blank","blank",20000,0,700);
    TH1I *carbonSRaw = new TH1I("carbonS","carbonS",20000,0,700);
    TH1I *carbonLRaw = new TH1I("carbonL","carbonL",20000,0,700);
    TH1I *Sn112Raw = new TH1I("Sn112","Sn112",20000,0,700);
    TH1I *NatSnRaw = new TH1I("NatSn","NatSn",20000,0,700);
    TH1I *Sn124Raw = new TH1I("Sn124","Sn124",20000,0,700);
    TH1I *totalRaw = new TH1I("total","total",20000,0,700);

    TH1I *blankRawLog = new TH1I("blankLog","blank",20000,0,TMath::Log10(700));
    TH1I *carbonSRawLog = new TH1I("carbonSLog","carbonS",20000,0,TMath::Log10(700));
    TH1I *carbonLRawLog = new TH1I("carbonLLog","carbonL",20000,0,TMath::Log10(700));
    TH1I *Sn112RawLog = new TH1I("Sn112Log","Sn112",20000,0,TMath::Log10(700));
    TH1I *NatSnRawLog = new TH1I("NatSnLog","NatSn",20000,0,TMath::Log10(700));
    TH1I *Sn124RawLog = new TH1I("Sn124Log","Sn124",20000,0,TMath::Log10(700));
    TH1I *totalRawLog = new TH1I("totalLog","total",20000,0,TMath::Log10(700));

    TH1I *TOF = new TH1I("TOF","Summed-detector time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",500,0,MICRO_PERIOD+1,2000,0,65536);

    // point at the ch4Tree in preparation for reading data
    setBranches(orchard[i]);

    int totalEntries = orchard[i]->GetEntries();
    cout << "Populating TOF, cross-section, etc. histograms..." << endl;

    for(int j=0; j<totalEntries; j++)
    {
        orchard[i]->GetEntry(j);

        double timeDiff = completeTime-macroTime+TIME_OFFSET;

        if (timeDiff < 650000 && timeDiff > 0)
        {
            TOF->Fill(trueTime);
            triangle->Fill(trueTime,lgQ);

            // convert trueTime into neutron velocity based on flight path distance
            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/trueTime; // in meters/sec 

            // convert velocity to relativistic kinetic energy
            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            if (trueTime > 90 /* && !gammaInMicro*/) // gate disallowing gammas
            {
                switch (targetPos)
                {
                    case 1:
                        // BLANK
                        blankRaw->Fill(rKE);
                        blankRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 2:
                        // SHORT CARBON
                        carbonSRaw->Fill(rKE);
                        carbonSRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 3:
                        // LONG CARBON
                        carbonLRaw->Fill(rKE);
                        carbonLRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 4:
                        // Sn112
                        Sn112Raw->Fill(rKE);
                        Sn112RawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 5:
                        // Natural Sn
                        NatSnRaw->Fill(rKE);
                        NatSnRawLog->Fill(TMath::Log10(rKE));
                        break;

                    case 6:
                        // Sn124
                        Sn124Raw->Fill(rKE);
                        Sn124RawLog->Fill(TMath::Log10(rKE));
                        break;

                    default:
                        break;
                }

                totalRaw->Fill(rKE);
                totalRawLog->Fill(TMath::Log10(rKE));
            }
        }
    }

    //TH1I *TOF = new TH1I("TTOF","All events time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    
    /*TH1I *firstSTOF = new TH1I("firstSTOF","first in micro time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *secondSTOF = new TH1I("secondSTOF","second in micro time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
    TH1I *thirdSTOF = new TH1I("thirdSTOF","third in micro time of flight",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);

    TH1I *firstSTOFblank = new TH1I("firstSTOFblank","first in micro time of flight, blank",10000,-MICRO_PERIOD*1.1,MICRO_PERIOD*1.1);
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
    

    //TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",500,0,MICRO_PERIOD+1,2000,0,65536);

    bool firstInMicro = false;
    bool secondInMicro = false;
    bool thirdInMicro = false;
    bool gammaInMicro = false;

    double trueTime = -1;
    int microNo = -1;
    int microNoPrev = -1;
    int microNo2Prev = -1;
    int microNo3Prev = -1;
    */

    /*for(int j = 0; j<channelList.size(); j++)
    {
        // loop through tree once per channel number 
        cout << "Populating " << dirs[j] << " histograms..." << endl;
        nE = 0;

        gDirectory->cd("/");
        gDirectory->GetDirectory(dirs[j].c_str())->cd();

        DPPWaveformsDir = (TDirectory*)gDirectory->Get("DPPWaveformsDir");
        WaveWaveformsDir = (TDirectory*)gDirectory->Get("WaveWaveformsDir");

        int targetCounter = 0;
        targetTime = get<1>(targetTimeList[targetCounter]);

        tree->SetEntryList(channelList[j]);
        totalEntries = tree->GetEntries();

        rng = new TRandom3();

        double fullTimeP = 0;

        for (int i=0; i<totalEntries; i++)

        {
            tree->GetEntry(i);

            if (chNo == chNoList[j])
            {
                // prepare for filling basic histos
                fullTime = timetag+pow(2,32)*extTime;

                if (chNo == 4 || chNo == 6 || chNo == 7)
                {
                    fullTime += fineTime*(2./1024.);
                }

                if (evtType == 1)
                {
*/
                    /*if (fullTime < targetTime-1000000)
                      {
                      cout << fullTime << " " << targetTime << endl;
                      }*/
/*
                    while (fullTime-get<1>(targetTimeList[targetCounter+1])+TIME_OFFSET > 0)
                    {
                        // if it's been too long since the last target changer event,
                        // step to the next target changer event - provided
                        // we haven't reset the time because of a recent switch
                        // to waveform mode

                        if ((get<1>(targetTimeList[targetCounter]) < get<1>(targetTimeList[targetCounter+1])) || fullTimeP > fullTime)
                        {
                            targetCounter++;
                            targetTime = get<1>(targetTimeList[targetCounter]);
                            targetPos = get<2>(targetTimeList[targetCounter]);
                            targetType = get<3>(targetTimeList[targetCounter]);
                            fill(evtNo.begin(),evtNo.end(),0);

                            fill(waveformStart.begin(),waveformStart.end(),0); // prepare for next waveform mode

                            fullTimeP = fullTime; // update the time of the last event

                        }

                        else
                        {
                            break;
                        }
                    }
                    */

                    /*if (chNo==4 && targetCounter > 23388)
                      {
                      cout << "fullTime " << fullTime << " targetTime " << targetTime << " fineTime " << fineTime << "target counter " << targetCounter << endl;
                      abort();
                      }*/

                    // if event has associated target changer event, fill DPP histo
        /*
                    if (fullTime-targetTime+TIME_OFFSET < 650000 && fullTime-targetTime+TIME_OFFSET > 0) 
                    {
                        // within macropulse window; fill histos
                        outMacro = (TH1I*)(gDirectory->Get("outMacro"));
                        outMacro->Fill(get<0>(targetTimeList[targetCounter]));

                        outEvt = (TH1I*)(gDirectory->Get("outEvt"));
                        outEvt->Fill(evtNo[chNo]);

                        outExtTime = (TH1I*)(gDirectory->Get("outExtTime"));
                        outExtTime->Fill(extTime);

                        outTime = (TH1I*)(gDirectory->Get("outTime"));
                        outTime->Fill(timetag);

                        outSGQ = (TH1I*)(gDirectory->Get("outSGQ"));
                        outSGQ->Fill(sgQ);

                        outLGQ = (TH1I*)(gDirectory->Get("outLGQ"));
                        outLGQ->Fill(lgQ);

                        outFT = (TH1I*)(gDirectory->Get("outFT"));
                        outFT->Fill(fineTime + 16*rng->Rndm());
                        

                        if (dummyWaveform->size() > 0)
                        {
                            DPPWaveformsDir->cd();

                            stringstream temp;
                            temp << "macroNo " << get<0>(targetTimeList[targetCounter]) << "evtNo " << evtNo[chNo];
                            DPPWaveform = new TH1I(temp.str().c_str(),temp.str().c_str(),dummyWaveform->size(),0,dummyWaveform->size()*2);

                            for(int i=0;i<dummyWaveform->size();i++)
                            {
                                DPPWaveform->SetBinContent(i,dummyWaveform->at(i));
                            }

                            gDirectory->cd("/");
                            gDirectory->GetDirectory(dirs[j].c_str())->cd();

                        }

                        if (chNo==2 || chNo==4 || chNo==6 || chNo==7)
                        {
                            microNo3Prev = microNo2Prev;
                            microNo2Prev = microNoPrev;
                            microNoPrev = microNo;
                            trueTime = fmod(((double)extTime*pow(2,32)+(double)timetag+((double)fineTime*2./1024.)-targetTime+TIME_OFFSET),MICRO_PERIOD);
                            microNo = floor(((double)extTime*pow(2,32)+(double)timetag+((double)fineTime*2./1024.)-targetTime+TIME_OFFSET)/MICRO_PERIOD);

                            firstInMicro = false;
                            secondInMicro = false;
                            thirdInMicro = false;
                            gammaInMicro = false;

                            if (microNo != microNoPrev)
                            {
                                firstInMicro = true;
                            }

                            if (microNo != microNo2Prev && microNo==microNoPrev)
                            {
                                secondInMicro = true;
                            }

                            if (microNo !=microNo3Prev && microNo==microNoPrev && microNoPrev==microNo2Prev)
                            {
                                thirdInMicro = true;
                            }

                            if (trueTime < 100)
                            {
                                gammaInMicro = true;
                            }

                            timeDiff << trueTime << " " << timetag << " " << targetTime << endl;

                            // convert trueTime into neutron velocity based on flight path distance
                            double velocity = pow(10.,7.)*FLIGHT_DISTANCE/trueTime; // in meters/sec 

                            // convert velocity to relativistic kinetic energy
                            double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

                            if (trueTime > 100 && chNo != 2 *//* && !gammaInMicro*//*) // gate disallowing gammas and monitors
                            {
                                switch (targetPos)
                                {
                                    case 1:
                                        // BLANK
                                        blankRaw->Fill(rKE);
                                        blankRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 2:
                                        // SHORT CARBON
                                        carbonSRaw->Fill(rKE);
                                        carbonSRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 3:
                                        // LONG CARBON
                                        carbonLRaw->Fill(rKE);
                                        carbonLRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 4:
                                        // Sn112
                                        Sn112Raw->Fill(rKE);
                                        Sn112RawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 5:
                                        // Natural Sn
                                        NatSnRaw->Fill(rKE);
                                        NatSnRawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    case 6:
                                        // Sn124
                                        Sn124Raw->Fill(rKE);
                                        Sn124RawLog->Fill(TMath::Log10(rKE));
                                        break;

                                    default:
                                        break;
                                }

                                totalRaw->Fill(rKE);
                                totalRawLog->Fill(TMath::Log10(rKE));
                            }

                            if (trueTime > 50 && chNo == 2) // gate disallowing gammas and non-monitors
                            {
                                switch (targetPos)
                                {
                                    case 1:
                                        // BLANK
                                        monBlankRaw->Fill(rKE);
                                        break;
                                    case 2:
                                        // SHORT CARBON
                                        monCarbonSRaw->Fill(rKE);
                                        break;
                                    case 3:
                                        // LONG CARBON
                                        monCarbonLRaw->Fill(rKE);
                                        break;
                                    case 4:
                                        // Sn112
                                        monSn112Raw->Fill(rKE);
                                        break;
                                    case 5:
                                        // Natural Sn
                                        monNatSnRaw->Fill(rKE);
                                        break;
                                    case 6:
                                        // Sn124
                                        monSn124Raw->Fill(rKE);
                                        break;
                                    default:
                                        break;
                                }
                                monTotalRaw->Fill(rKE);
                            }

                            if (targetPos > 0 *//* && gammaInMicro*//*)
                            {
                                switch (chNo)
                                {
                                    case 2:
                                        MTOF->Fill(trueTime);
                                        break;
                                    case 4:
                                        STOF->Fill(trueTime);

                                        if (firstInMicro)
                                        {

                                            firstSTOF->Fill(trueTime);
                                            switch (targetPos)
                                            {
                                                case 1:
                                                    firstSTOFblank->Fill(trueTime);
                                                    break;
                                                case 2:
                                                    firstSTOFcs->Fill(trueTime);
                                                    break;
                                                case 3:
                                                    firstSTOFcl->Fill(trueTime);
                                                    break;
                                                case 4:
                                                    firstSTOFsn112->Fill(trueTime);
                                                    break;
                                                case 5:
                                                    firstSTOFsnnat->Fill(trueTime);
                                                    break;
                                                case 6:
                                                    firstSTOFsn124->Fill(trueTime);
                                                    break;
                                            }
                                        }

                                        else if (secondInMicro)
                                        {
                                            secondSTOF->Fill(trueTime);
                                            switch (targetPos)
                                            {
                                                case 1:
                                                    secondSTOFblank->Fill(trueTime);
                                                    break;
                                                case 2:
                                                    secondSTOFcs->Fill(trueTime);
                                                    break;
                                                case 3:
                                                    secondSTOFcl->Fill(trueTime);
                                                    break;
                                                case 4:
                                                    secondSTOFsn112->Fill(trueTime);
                                                    break;
                                                case 5:
                                                    secondSTOFsnnat->Fill(trueTime);
                                                    break;
                                                case 6:
                                                    secondSTOFsn124->Fill(trueTime);
                                                    break;
                                            }
                                        }

                                        else if (thirdInMicro)
                                        {
                                            thirdSTOF->Fill(trueTime);
                                            switch (targetPos)
                                            {
                                                case 1:
                                                    thirdSTOFblank->Fill(trueTime);
                                                    break;
                                                case 2:
                                                    thirdSTOFcs->Fill(trueTime);
                                                    break;
                                                case 3:
                                                    thirdSTOFcl->Fill(trueTime);
                                                    break;
                                                case 4:
                                                    thirdSTOFsn112->Fill(trueTime);
                                                    break;
                                                case 5:
                                                    thirdSTOFsnnat->Fill(trueTime);
                                                    break;
                                                case 6:
                                                    thirdSTOFsn124->Fill(trueTime);
                                                    break;
                                            }
                                        }

                                        triangle->Fill(trueTime,lgQ);
                                        break;

                                    case 6:
                                        LTOF->Fill(trueTime);
                                        break;
                                    case 7:
                                        RTOF->Fill(trueTime);
                                }
                            } 
                        }
                    }

                    prevTarget=1;
                }

                else if (evtType == 2)
                {
                    WaveWaveformsDir->cd();

                    TH1I* waveformHolder;

                    if (fullTime >= waveformStart[chNo]+650000 || prevTarget==1)
                    {
                        // new macropulse in waveform mode - create new plot

                        stringstream temp;
                        //cout << "waveform mode targetCounter" << targetCounter << " and supposed macropulse number " << get<0>(targetTimeList[targetCounter]) << endl;
                        temp << "macropulse " << get<0>(targetTimeList[targetCounter]) << " event no " << evtNo[chNo];
                        waveformStart[chNo] = fullTime;

                        switch (chNo)
                        {
                            case 0:
                                WaveWaveform0 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform0;
                                break;

                            case 2:
                                WaveWaveform2 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform2;
                                break;

                            case 4:
                                WaveWaveform4 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform4;
                                break;

                            case 6:
                                WaveWaveform6 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform6;
                                break;

                            case 7:
                                WaveWaveform7 = new TH1I(temp.str().c_str(),temp.str().c_str(),400000,0,800000);
                                waveformHolder = WaveWaveform7;
                                break;
                        }
                    }

                    for(int i=0;i<dummyWaveform->size();i++)
                    {
                        waveformHolder->SetBinContent(i+(fullTime-waveformStart[chNo])/2,dummyWaveform->at(i));
                    }

                    gDirectory->cd("/");
                    gDirectory->GetDirectory(dirs[j].c_str())->cd();

                    prevTarget=2;
                }

                nE++;
                evtNo[chNo]++;

                if(nE%100==0)
                {
                    cout << nE << " events\r";
                    fflush(stdout);
                }
            }
        }
        cout << endl;
    }
    */
}


// Loop through all trees (one per channel) and populate their data into basic
// histograms. Then, calculate the TOF, cross-section, etc using the channel 4
// data and produce histograms of these calculated quantities.
void fillHistos()
{
    for(int i=0; i<orchard.size(); i++)
    {
        gDirectory->cd("/");
        gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());

        gDirectory->GetDirectory(dirs[i].c_str())->cd();

        // instantiate histograms

        TH1I* macroNoH = new TH1I("macroNoH","macroNo",100000,0,100000);
        macroNoH->GetXaxis()->SetTitle("macropulse number of each event");

        TH1I* microNoH = new TH1I("microNoH","microNo",360,0,360);
        microNoH->GetXaxis()->SetTitle("micropulse number of each event");

        TH1I* evtNoH = new TH1I("evtNoH","evtNo",2500,0,2500);
        evtNoH->GetXaxis()->SetTitle("event number of each event");

        TH1I* macroTimeH = new TH1I("macroTimeH","macroTime",6000,0,6000000000);
        macroTimeH->GetXaxis()->SetTitle("macropulse time zero for each event");

        TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",6000,0,6000000000);
        completeTimeH->GetXaxis()->SetTitle("complete time for each event");

        TH1I* trueTimeH = new TH1I("trueTimeH","trueTime",2000,0,2000);
        trueTimeH->GetXaxis()->SetTitle("time since start of micro for each event");

        TH1I* targetPosH = new TH1I("targetPosH","targetPos",6,0,6);
        targetPosH->GetXaxis()->SetTitle("target position of each event");

        TH1I* sgQH = new TH1I("sgQH","sgQ",3500,0,35000);
        sgQH->GetXaxis()->SetTitle("short gate integrated charge for each event");

        TH1I* lgQH = new TH1I("lgQH","lgQ",7000,0,70000);
        lgQH->GetXaxis()->SetTitle("long gate integrated charge for each event");

        gDirectory->mkdir("DPPWaveformsDir","raw DPP waveforms");
        gDirectory->mkdir("WaveWaveformsDir","concatenated waveform waveforms");

        // set up tree access
        setBranches(orchard[i]);

        int totalEntries = orchard[i]->GetEntries();
        cout << "Populating " << dirs[i] << " histograms..." << endl;

        for(int j=0; j<totalEntries; j++)
        {
            orchard[i]->GetEntry(j);
            if(completeTime-macroTime+TIME_OFFSET < 650000 && completeTime-macroTime+TIME_OFFSET > 0)
            {
                macroNoH->Fill(macroNo);
                microNoH->Fill(microNo);
                evtNoH->Fill(evtNo);
                macroTimeH->Fill(macroTime);
                completeTimeH->Fill(completeTime);
                trueTimeH->Fill(trueTime);
                targetPosH->Fill(targetPos);
                sgQH->Fill(sgQ);
                lgQH->Fill(lgQ);

                //cout << "waveform first element is " << waveform->at(0) << endl;

                if(waveform->size() > 0)
                {
                    DPPWaveformsDir = (TDirectory*)gDirectory->Get("DPPWaveformsDir");
                    DPPWaveformsDir->cd();

                    stringstream temp;
                    temp << "macroNo " << macroNo << ", evtNo " << evtNo;
                    DPPWaveform = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),0,waveform->size()*2);

                    for(int k=0; k<waveform->size(); k++)
                    {
                        DPPWaveform->SetBinContent(k,waveform->at(k));
                    }
                    gDirectory->cd("..");
                }
            }
        }
    }

    // fill ch4 advanced histograms
    fillAdvanced(2);

    // fill ch6 advanced histograms
    fillAdvanced(3);
}

int main(int argc, char* argv[])
{
    TApplication app("app",&argc,argv);

    // read in the raw file name
    string runDir = argv[1];
    string runNo = argv[2];

    // Open the raw tree from the initial sort. If it doesn't exist, exit.
 
    stringstream treeName;
    stringstream fileName;

    treeName << runDir << "-" << runNo; 
    fileName << analysispath <<"analysis/" << runDir << "/" << treeName.str() << ".root";

    TFile* file = new TFile(fileName.str().c_str(),"UPDATE");

    TTree* ch0Tree = (TTree*)file->Get("ch0Tree");
    TTree* ch2Tree = (TTree*)file->Get("ch2Tree");
    TTree* ch4Tree = (TTree*)file->Get("ch4Tree");
    TTree* ch6Tree = (TTree*)file->Get("ch6Tree");
    TTree* ch0TreeW = (TTree*)file->Get("ch0TreeW");
    TTree* ch2TreeW = (TTree*)file->Get("ch2TreeW");
    TTree* ch4TreeW = (TTree*)file->Get("ch4TreeW");
    TTree* ch6TreeW = (TTree*)file->Get("ch6TreeW");

    orchard.push_back(ch0Tree);
    orchard.push_back(ch2Tree);
    orchard.push_back(ch4Tree);
    orchard.push_back(ch6Tree);
    
    // prepare the root file with 4 directories, one for each channel
    // these directories will hold basic variable histograms showing the
    // raw data in each tree, plus TOF, x-sections, etc histograms

    fillHistos();

    file->Write();

    file->Close();
}
