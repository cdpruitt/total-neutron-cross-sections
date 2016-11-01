#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TMath.h"
#include "../include/physicalConstants.h"
#include "../include/analysisConstants.h"
#include "../include/helperFunctions.h"
#include "../include/dataStructures.h"
#include "../include/plottingConstants.h"
#include "../include/targetConstants.h"
#include "../include/branches.h"
#include "../include/plots.h"
#include "../include/histos.h"
#include "../include/waveform.h"

using namespace std;

long totalMicros = 0;

extern ProcessedEvent procEvent;

vector<TTree*> orchard; // holds DPP channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

vector<TTree*> orchardW; // holds waveform-mode channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

// a small subset of events occur in both channels 4 and 6 in the microTime
// window of ~275 ns to ~325 ns after a micropulse start. We'd like to plot
// these on the same axis to make sure they look the same.
TH1I* waveformCh4;
TH1I* waveformCh6;

TDirectory *waveformsDir;

// Loop through all trees (one per channel) and populate their data into basic
// histograms. Then, calculate the TOF, cross-section, etc using the channel 4
// data and produce histograms of these calculated quantities.
void fillHistos(vector<Plots*>& plots)
{
    cout << endl << "Entering ./histos..." << endl;
    TH1* waveformH;
    // fill basic histograms for DPP mode in each channel

    // first loop through all channel-specific DPP-mode trees
    for(int i=0; (size_t)i<orchard.size(); i++)
    {
        // create a channel-specific directory for putting histograms inside
        gDirectory->cd("/");
        gDirectory->mkdir(dirs[i].c_str(),dirs[i].c_str());
        gDirectory->GetDirectory(dirs[i].c_str())->cd();

        // instantiate DPP-mode histograms
        TH1I* macroNoH = new TH1I("macroNoH","macroNo",500000,0,500000);
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

        // diagnostic histograms
        TH1I *TOF = new TH1I("TOF","Summed-detector time of flight",TOF_BINS,0,TOF_RANGE);
        TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",TOF_RANGE,0,MICRO_LENGTH+1,2048,0,65536);
        TH2I *triangleRKE = new TH2I("triangleRKE","Pulse integral vs. relativistic KE",NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,2048,0,65536);
        TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
        TH1I *QRatio = new TH1I("QRatio","short gate Q/long gate Q",1000,0,1);
        TH2I *rKElgQ = new TH2I("lgQrKE","relativistic KE vs. long gate Q",500,ENERGY_LOWER_BOUND,100,2048,0,65536);
        TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
        microNoH->GetXaxis()->SetTitle("micropulse number of each event");

        // create a subdirectory for holding DPP-mode waveform data
        gDirectory->mkdir("waveformsDir","raw DPP waveforms");

        // reattach to the channel-specific tree for reading out data
        if(i==0)
        {
            setBranchesTC(orchard[i]);
        }
        else
        {
            setBranchesHistos(orchard[i]);
        }

        // create TIME VARIABLES used for filling histograms
        double microTime;
        int microNo, prevMicroNo;

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
        switch(i)
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

        cout << "Populating advanced histograms for channel " << 2*i << endl;

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

            if(i==0)
            {
                macroNoH->Fill(procEvent.macroNo);
                macroTimeH->Fill(procEvent.macroTime);
                targetPosH->Fill(procEvent.targetPos);
                continue;
            }
            if(procEvent.lgQ!=65535 && procEvent.sgQ!=32767)
            {
                macroNoH->Fill(procEvent.macroNo);
                evtNoH->Fill(procEvent.evtNo);
                macroTimeH->Fill(procEvent.macroTime);
                //completeTimeH->Fill(procEvent.completeTime);
                targetPosH->Fill(procEvent.targetPos);
                sgQH->Fill(procEvent.sgQ);
                lgQH->Fill(procEvent.lgQ);
            }

            // if waveform data for this event exist, we want to populate
            // a histogram to display it

            // only plot 1 out of 10000 waveforms to save space and processing
            // time

            if(procEvent.waveform->size() > 0 && j%100000 == 0)
            {
                stringstream temp;
                temp << "macroNo " << procEvent.macroNo << ", evtNo " << procEvent.evtNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,procEvent.waveform->size()*2);

                // loop through waveform data and fill histo
                for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
                {
                    waveformH->SetBinContent(k,procEvent.waveform->at(k));
                }
            }

            /*cout << "completeTime = " << procEvent.completeTime << endl;

              if(procEvent.macroNo>0)
              {
              break;
              }*/
        
            if(i==2)
            {
                // calculate time since start of macro (includes time offsets)
                double timeDiff = procEvent.completeTime-procEvent.macroTime;

                // Apply gates:
                if (timeDiff < 650000 && timeDiff > 0           // require events to come during the macropulse's beam-on period
                        && procEvent.targetPos != 0                              // discard events during target-changer movement
                        && procEvent.lgQ<65500 && procEvent.sgQ<32750                      // discard events with unphysically-large integrated charges
                        && procEvent.lgQ>procEvent.sgQ                                     // discard events with short gate charge is larger than long gate charge
                        //&& (procEvent.sgQ/(double)procEvent.lgQ<0.25 || procEvent.sgQ/(double)procEvent.lgQ>0.35)  // discard events outside accepted range of sgQ/lgQ
                        //&& procEvent.lgQ>100                                     // discard gammas at lowest range of energy
                   )
                {
                    /*****************************************************************/
                    // Calculate event properties

                    // find which micropulse the event is in and the time since
                    // the start of the micropulse

                    // first, save previous event's micropulse number (we'll need this
                    // to calculate event ordering in each micropulse)
                    prevMicroNo = microNo;

                    microNo = floor(timeDiff/MICRO_LENGTH);
                    microTime = fmod(timeDiff,MICRO_LENGTH);

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
                    //if (procEvent.lgQ>50)
                    //if (procEvent.lgQ>500*exp(-(microTime-100)/87))
                    //if (procEvent.lgQ<30*rKE)

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
                    }

                    /*****************************************************************/
                    // Fill troubleshooting plots with event variables (rKE, microtime, etc.)
                    TOF->Fill(microTime);
                    triangle->Fill(microTime,procEvent.lgQ);
                    sgQlgQ->Fill(procEvent.sgQ,procEvent.lgQ);
                    QRatio->Fill(procEvent.sgQ/(double)procEvent.lgQ);
                    rKElgQ->Fill(rKE,procEvent.lgQ);
                    triangleRKE->Fill(microTime,rKE);
                    microNoH->Fill(microNo);

                    // populate only first three orderInMicro plots
                    /*if(orderInMicro<4)
                      {
                      plots.orderInMicro[orderInMicro-1]->Fill(microTime);
                      }*/

                    /*****************************************************************/

                    /*****************************************************************/
                    // Fill target-specific plots

                    if (procEvent.targetPos>0 && procEvent.targetPos<=NUMBER_OF_TARGETS)
                    {
                        plots[procEvent.targetPos-1]->getTOFHisto()->Fill(microTime);
                        plots[procEvent.targetPos-1]->getEnergyHisto()->Fill(rKE);

                        if(!gammaInMicro)
                        {
                            //plots.energyHistosNoGamma[procEvent.targetPos-1]->Fill(rKE);
                        }

                        if(orderInMicro==1)
                        {
                            //plots.TOFHistosFirstInMicro[procEvent.targetPos-1]->Fill(microTime);
                        }

                        // end of energy, time, cross-section gates on events
                        /*****************************************************************/
                    }

                    // fill ungated plots
                    /*if (procEvent.targetPos>0 && procEvent.targetPos<=NUMBER_OF_TARGETS)
                      {
                      plots.energyHistosUngated[procEvent.targetPos-1]->Fill(rKE);
                      }*/

                    // end of main event loop
                    /*****************************************************************/

                    /*if(procEvent.macroNo>0)
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

        TH1I* targetPosH = new TH1I("targetPosH","targetPos",7,0,7);
        targetPosH->GetXaxis()->SetTitle("target position of each event");

        // create subdirectory for holding waveform-mode waveform data
        gDirectory->mkdir("waveformsDir","raw waveform-mode waveforms");
        waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
        waveformsDir->cd();

        // fill waveform mode histograms
        setBranchesHistosW(orchardW[i]);

        int totalEntries = orchardW[i]->GetEntries();

        // we need label the number of waveform-mode macropulses to make
        // uniquely named histograms
        int waveformNo = 0;

        for(int j=0; j<totalEntries; j++)
        {
            //cout << "looping thru waveform mode events, #" << j << endl;
            //fflush(stdout);

            orchardW[i]->GetEntry(j);

            if(i==0)
            {
                macroNoH->Fill(procEvent.macroNo);
                targetPosH->Fill(procEvent.targetPos);
                continue;
            }

            macroNoH->Fill(procEvent.macroNo);
            evtNoH->Fill(procEvent.evtNo);
            completeTimeH->Fill(procEvent.completeTime);

            if(procEvent.evtNo==0)
            {
                // new macropulse in waveform mode
                // create a new plot to hold the new waveforms of this macropulse
                stringstream temp;
                temp << "full waveform " << waveformNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),360000,0,720000);

                // set the start of the macropulse to the first event timer
                waveformNo++;
            }

            for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
            {
                waveformH->SetBinContent(k+floor((fmod(procEvent.waveform->size()*2,MICRO_LENGTH)+1)*procEvent.evtNo*MICRO_LENGTH/(double)2),procEvent.waveform->at(k));
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

double ch6timeDiff = procEvent.completeTime-procEvent.macroTime;
double ch6microTime = fmod(ch6timeDiff,MICRO_LENGTH);
int ch6microNo = floor(ch6timeDiff/MICRO_LENGTH);

// keep track of this event's macroNo to compare it with events in ch4;
// this will help us identify the parallel event in ch4 with the same
// timestamp
unsigned int ch6MacroNo = procEvent.macroNo;

if(ch6microTime > 275 && ch6microTime < 325)
{
// used to calculate the channel 6 average value of the waveform and
// find the baseline offset between ch6 and ch4
//int ch6Average = 0; 
stringstream temp;
temp << "macroNo " << procEvent.macroNo << ", scavenger event number " << procEvent.evtNo;
waveformCh6= new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,procEvent.waveform->size()*2);

// loop through waveform data and fill histo
for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
{
waveformCh6->SetBinContent(k,procEvent.waveform->at(k)+84);
//ch6Average += procEvent.waveform->at(k);
}

// loop through ch4 tree to find the same event (same timestamp and
// macropulse as the ch6 event)

//cout << "Found a ch6 event in the window. Moving to ch4 events..." << endl;

setBranches(orchard[2]); // pointed at the ch6Tree now
int microNo = 0;

for(int i=currentIndex; i<ch4Entries; i++)
{
orchard[2]->GetEntry(i);

if(procEvent.macroNo == ch6MacroNo)
{
    // in the correct macropulse - loop through ch4 to find the
    // the event whose timestamp matches the one in channel 6
    double ch4timeDiff = procEvent.completeTime-procEvent.macroTime;
    double ch4microTime = fmod(ch4timeDiff,MICRO_LENGTH);
    microNo = floor(ch4timeDiff/MICRO_LENGTH);

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
            timeDiffQ->Fill(ch6microTime-ch4microTime, procEvent.lgQ);

            //cout << "ch4microTime = " << ch4microTime << ", ch4microNo = " << microNo << ", ch6microTime = " << ch6microTime << ", ch6microNo = " << ch6microNo << endl;

            // used to calculate the channel 4 average value of the waveform and
            // find the baseline offset between ch6 and ch4
            //int ch4Average = 0;

            // clear the histo title stringstream
            temp.str("");
            temp << "macroNo " << procEvent.macroNo << ", ch4 event number " << procEvent.evtNo << " w/ lgQ=65535";
            waveformCh4 = new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,procEvent.waveform->size()*2);

            for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
            {
                waveformCh4->SetBinContent(k,procEvent.waveform->at(k));
                //   ch4Average += procEvent.waveform->at(k);
            }

            //offset += (ch4Average-ch6Average)/(double)procEvent.waveform->size();

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

else if (procEvent.macroNo > ch6MacroNo)
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

int histos(string sortedFileName, string histoFileName)
{
    TFile* sortedFile = new TFile(sortedFileName.c_str(),"READ");
    if(!sortedFile->IsOpen())
    {
        cerr << "Error: failed to open resort.root" << endl;
        exit(1);
    }

    TTree* ch0Tree = (TTree*)sortedFile->Get("targetChangerTree");
    TTree* ch2Tree = (TTree*)sortedFile->Get("ch2ProcessedTree");
    TTree* ch4Tree = (TTree*)sortedFile->Get("ch4ProcessedTree");
    //TTree* ch6Tree = (TTree*)sortedFile->Get("ch6ProcessedTree");
    TTree* ch0TreeW = (TTree*)sortedFile->Get("ch0ProcessedTreeW");
    TTree* ch2TreeW = (TTree*)sortedFile->Get("ch2ProcessedTreeW");
    TTree* ch4TreeW = (TTree*)sortedFile->Get("ch4ProcessedTreeW");

    orchard.push_back(ch0Tree);
    orchard.push_back(ch2Tree);
    orchard.push_back(ch4Tree);

    orchardW.push_back(ch0TreeW);
    orchardW.push_back(ch2TreeW);
    orchardW.push_back(ch4TreeW);

    // increase precision to handle outputted times (for troubleshooting)
    cout.precision(13);

    vector<Plots*> plots;

    // open output file to contain histos
    TFile* histoFile = new TFile(histoFileName.c_str(),"READ");
    if(!histoFile->IsOpen())
    {
        // No histogram file - create and fill
        histoFile = new TFile(histoFileName.c_str(),"CREATE");

        for(int i=0; i<NUMBER_OF_TARGETS; i++)
        {
            plots.push_back(new Plots(positionNames[i]));
        }

        // prepare the root file with 4 directories, one for each channel
        // these directories will hold basic variable histograms showing the
        // raw data in each tree, plus TOF, x-sections, etc histograms
        fillHistos(plots);
        // fill TOF, cross-section, etc. histos for channels 4, 6

        sortedFile->Close();
        histoFile->Write();
    }

    else
    {
        for(int i=0; i<NUMBER_OF_TARGETS; i++)
        {
            plots.push_back(new Plots(positionNames[i], histoFile));
        }
    }


    std::vector<long> macrosPerTarget;
    std::vector<long> microsPerTarget;

    long totalMacros = 0;

    gDirectory->cd("/");
    gDirectory->cd(dirs[0].c_str());

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        macrosPerTarget.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
        totalMacros+=macrosPerTarget.back();

        microsPerTarget.push_back(macrosPerTarget.back()*(MACRO_LENGTH/MICRO_LENGTH));
        totalMicros+=microsPerTarget.back();

        cout << "Micropulses on target " << i+1 << ": " << microsPerTarget.back() << endl;
    }

    int numberTotalTriggers = 0;
    for(Plots* p : plots)
    {
        numberTotalTriggers += p->getTOFHisto()->GetEntries();
    }

    cout << "Total number of triggers: " << numberTotalTriggers << endl;
    cout << "Triggers/micropulse: " << numberTotalTriggers/(double)totalMicros << endl;

    calculateDeadtime(microsPerTarget, plots);
    histoFile->Write();

    // Modify plots
    /*for(int i=0; i<NUMBER_OF_TARGETS; i++)
      {
      plots.energyHistos[i]->GetXaxis()->SetRangeUser(0,700);
      plots.energyHistosCorrected[i]->GetXaxis()->SetRangeUser(0,700);
      }*/

    // print out waveforms for first 50 events that occur in both ch4 and ch6
    //matchWaveforms();

    //sortedFile->Close();
    //histoFile->Close();

    return 0;
}
