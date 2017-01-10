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
#include "TRandom3.h"

#include "../include/physicalConstants.h"
#include "../include/analysisConstants.h"
#include "../include/dataStructures.h"
#include "../include/plottingConstants.h"
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

TDirectory *waveformsDir;

// Loop through all trees (one per channel) and populate their data into basic
// histograms. Then, calculate the TOF, cross-section, etc using the channel 4
// data and produce histograms of these calculated quantities.

void fillVetoedHistos(TFile* vetoFile, TFile* histoFile)
{
    for(string s : detectorNames)
    {
        cout << "Filling vetoed histos for " << s << endl;

        vetoFile->cd();
        string treeName = s + "Clean";
        TTree* tree = (TTree*)vetoFile->Get(treeName.c_str());

        // move to detector directory in histo file
        histoFile->cd("/");
        histoFile->cd(s.c_str());

        // create diagnostic histograms
        TH2I *triangle = new TH2I("triangle","Pulse integral vs. TOF",TOF_RANGE,0,MICRO_LENGTH+1,2048,0,65536);
        TH2I *triangleRKE = new TH2I("triangleRKE","Pulse integral vs. relativistic KE",NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,2048,0,65536);
        TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
        TH1I *QRatio = new TH1I("QRatio","short gate Q/long gate Q",1000,0,1);
        TH2I *rKElgQ = new TH2I("lgQrKE","relativistic KE vs. long gate Q",500,ENERGY_LOWER_BOUND,100,2048,0,65536);
        TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
        microNoH->GetXaxis()->SetTitle("micropulse number of each event");

        // make target-specific plots
        vector<Plots*> plots;
        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            plots.push_back(new Plots(positionNames[i]));
        }

        // reattach to the channel-specific tree for reading out data
        setBranchesHistos(tree);

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
        // Prepare GAMMA GATE for filtering events depending on gamma presence

        // channel-dependent time offset relative to the target changer's macropulse
        // start time
        int gammaGate[2] = {80,90};

        long totalEntries = tree->GetEntries();
        for(long i=0; i<totalEntries; i++)
        {
            tree->GetEntry(i);

            // calculate time since start of macro (includes time offsets)
            double timeDiff = procEvent.completeTime-procEvent.macroTime;

            // Apply gates:
            if (timeDiff < MACRO_LENGTH && timeDiff > 0 // omit events outside beam-on period
                    && procEvent.targetPos != 0                // discard events during target-changer movement
                    && procEvent.lgQ<65500 && procEvent.sgQ<32750  // discard events with unphysical integrated charges
                    //&& (procEvent.sgQ/(double)procEvent.lgQ<0.25 || procEvent.sgQ/(double)procEvent.lgQ>0.35)  // discard events outside accepted range of sgQ/lgQ
                    //&& procEvent.lgQ>100    // discard gammas at lowest range of energy
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
                triangle->Fill(microTime,procEvent.lgQ);
                sgQlgQ->Fill(procEvent.sgQ,procEvent.lgQ);
                QRatio->Fill(procEvent.sgQ/(double)procEvent.lgQ);
                rKElgQ->Fill(rKE,procEvent.lgQ);
                triangleRKE->Fill(microTime,rKE);
                microNoH->Fill(microNo);

                /*****************************************************************/
                // Fill target-specific plots

                if (procEvent.targetPos>0 && procEvent.targetPos<=tarGates.size())
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
                }
            }

            if(i%10000==0)
            {
                cout << "Processed " << i << " events...\r";
            }
        }

        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            plots[i]->getTOFHisto()->Write();
            plots[i]->getEnergyHisto()->Write();
        }

        std::vector<long> macrosPerTarget;
        std::vector<long> microsPerTarget;

        long totalMacros = 0;

        histoFile->cd("/targetChanger");

        for(unsigned int i=0; i<tarGates.size(); i++)
        {
            macrosPerTarget.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
            totalMacros+=macrosPerTarget.back();

            microsPerTarget.push_back(macrosPerTarget.back()*(MACRO_LENGTH/MICRO_LENGTH));
            totalMicros+=microsPerTarget.back();
        }

        histoFile->cd("/");
        histoFile->cd(s.c_str());

        calculateDeadtime(microsPerTarget, plots);
    }
}

void fillBasicHistos(TFile* histoFile)
{
    // create diagnostic graphs for each channel's DPP mode data
    for(TTree* t : orchard)
    {
        // create a channel-specific directory for putting histograms inside
        histoFile->cd("/");
        histoFile->mkdir(t->GetName(),t->GetName());
        histoFile->GetDirectory(t->GetName())->cd();

        // instantiate DPP-mode histograms
        TH1I* macroNoH = new TH1I("macroNoH","macroNo",200000,0,200000);
        macroNoH->GetXaxis()->SetTitle("macropulse number of each event");

        TH1I* evtNoH = new TH1I("evtNoH","evtNo",150,0,150);
        evtNoH->GetXaxis()->SetTitle("event number of each event");

        TH1I* macroTimeH = new TH1I("macroTimeH","macroTime",6000,0,6000000000);
        macroTimeH->GetXaxis()->SetTitle("macropulse time zero for each event");

        TH1I* targetPosH = new TH1I("targetPosH","targetPos",7,0,7);
        targetPosH->GetXaxis()->SetTitle("target position of each event");

        TH1I* sgQH = new TH1I("sgQH","sgQ",3500,0,35000);
        sgQH->GetXaxis()->SetTitle("short gate integrated charge for each event");

        TH1I* lgQH = new TH1I("lgQH","lgQ",7000,0,70000);
        lgQH->GetXaxis()->SetTitle("long gate integrated charge for each event");

        // create a subdirectory for holding DPP-mode waveform data
        gDirectory->mkdir("waveformsDir","raw DPP waveforms");

        /*************************************************************************/
        // Loop through sorted trees to calculate advanced histogram variables

        string name = t->GetName();
        if(name=="targetChanger")
        {
            setBranchesTC(t);
        }

        else
        {
            setBranchesHistos(t);
        }

        long totalEntries = t->GetEntries();

        cout << "Populating " << t->GetName() << " histograms..." << endl;

        waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
        waveformsDir->cd();

        // loop through the channel-specific tree and populate histos
        for(long j=0; j<totalEntries; j++)
        {
            t->GetEntry(j);

            macroNoH->Fill(procEvent.macroNo);
            targetPosH->Fill(procEvent.targetPos);
            macroTimeH->Fill(procEvent.macroTime);

            if(name=="targetChanger")
            {
                continue;
            }

            evtNoH->Fill(procEvent.evtNo);
            sgQH->Fill(procEvent.sgQ);
            lgQH->Fill(procEvent.lgQ);

            // if waveform data for this event exist, we want to populate
            // a histogram to display it

            // only plot 1 out of 10000 waveforms to save space and processing
            // time

            if(procEvent.waveform->size() > 0 && j%100000 == 0)
            {
                stringstream temp;
                temp << "macroNo " << procEvent.macroNo << ", evtNo " << procEvent.evtNo;
                TH1I* waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,procEvent.waveform->size()*SAMPLE_PERIOD);

                // loop through waveform data and fill histo
                for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
                {
                    waveformH->SetBinContent(k,procEvent.waveform->at(k));
                }

                waveformH->Write();
            }

            if(j%1000==0)
            {
                cout << "Processed " << j << " events...\r";
                fflush(stdout);
            }
        }

    }

    // fill basic histograms for waveform mode in each channel

    // loop through all channel-specific waveform-mode trees
    /*for(TTree* t : orchardW)
    {
        // create a channel-specific directory for each tree
        gDirectory->cd("/");
        string tempName = t->GetName();
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

        setBranchesHistosW(t);

        long totalEntries = t->GetEntries();

        // we need label the number of waveform-mode macropulses to make
        // uniquely named histograms
        int waveformNo = 0;

        TH1I* waveformH;

        for(long j=0; j<totalEntries; j++)
        {
            t->GetEntry(j);

            macroNoH->Fill(procEvent.macroNo);
            targetPosH->Fill(procEvent.targetPos);
            evtNoH->Fill(procEvent.evtNo);
            completeTimeH->Fill(procEvent.completeTime);

            if(procEvent.evtNo==1)
            {
                // new macropulse in waveform mode
                // create a new plot to hold the new waveforms of this macropulse
                stringstream temp;
                temp << "full waveform " << waveformNo;
                waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),360000,0,720000);
                waveformH->Write();
                // set the start of the macropulse to the first event timer
                waveformNo++;
            }

            for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
            {
                waveformH->SetBinContent(k+floor((fmod(procEvent.waveform->size()*SAMPLE_PERIOD,MICRO_LENGTH)+1)*procEvent.evtNo*MICRO_LENGTH/(double)2),procEvent.waveform->at(k));
                waveformH->Write();
            }
        }
    }*/

    histoFile->Write();
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

void correctForDeadtime(string histoFileName, string deadtimeFileName, vector<string> detectorChannels)
{
    TFile* deadtimeFile = new TFile(deadtimeFileName.c_str(),"READ");
    TFile* histoFile = new TFile(histoFileName.c_str(),"UPDATE");

    for(string directory : detectorChannels)
    {
        gDirectory->cd("/");
        gDirectory->cd(directory.c_str());

        vector<Plots*> uncorrectedPlots;
        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            string name = positionNames[i];
            uncorrectedPlots.push_back(new Plots(name,histoFile,directory));
        }

        vector<Plots*> correctedPlots;
        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            string name = positionNames[i] + "Corrected";
            correctedPlots.push_back(new Plots(name));
        }

        vector<Plots*> deadtimePlots;
        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            string name = positionNames[i];
            deadtimePlots.push_back(new Plots(name, deadtimeFile, directory));
        }

        // extract deadtime from waveform-mode fit

        TRandom3 *randomizeBin = new TRandom3();

        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            // "deadtimeFraction" records the fraction of time that the detector is dead, for
            // neutrons of a certain energy.

            vector<double> deadtimeFraction;

            //string temp;
            //temp = "deadtime" + t.getName() + "Waveform";
            //plots.waveformDeadtimes.push_back((TH1I*)deadtimeFile->Get(temp.c_str()));

            /*if(!t.getDeadtime.back())
              {
              cerr << "Error: couldn't find waveform deadtime histograms." << endl;
              exit(1);
              }*/

            TH1I* deadtimeHisto = deadtimePlots[i]->getDeadtimeHisto();
            if(!deadtimeHisto)
            {
                cout << "Couldn't find deadtimeHisto for target " << i << endl;
                continue;
            }

            int deadtimeBins = deadtimeHisto->GetNbinsX();

            for(int j=0; j<deadtimeBins; j++)
            {
                deadtimeFraction.push_back(deadtimeHisto->GetBinContent(j)/(double)pow(10,3));
            }

            // create deadtime-corrected histograms

            deadtimeHisto->Write();

            //vector<vector<double>> eventsPerBinPerMicro(6,vector<double>(0));

            //const double FULL_DEADTIME = 183; // total amount of time after firing when
            // detector is at least partially dead to
            // incoming pulses (in ns)
            //const double PARTIAL_DEADTIME = 9; // amount of time after the end of
            // FULL_DEADTIME when detector is
            // becoming live again, depending on
            // amplitude (in ns)

            /*************************************************************************/
            // Perform deadtime correction
            /*************************************************************************/

            // loop through all TOF histos

            TH1I* tof = uncorrectedPlots[i]->getTOFHisto();
            //TH1I* en = uncorrectedPlots[i]->getEnergyHisto();

            TH1I* tofC = correctedPlots[i]->getTOFHisto();
            TH1I* enC = correctedPlots[i]->getEnergyHisto();

            int tofBins = tofC->GetNbinsX();

            // apply deadtime correction to TOF histos
            for(int j=0; j<tofBins; j++)
            {
                if(deadtimeFraction[j] > 0)
                {
                    tofC->SetBinContent(j,(tof->GetBinContent(j)/(1-deadtimeFraction[j])));
                }

                // convert microTime into neutron velocity based on flight path distance
                double velocity = pow(10.,7.)*FLIGHT_DISTANCE/(tofC->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)(2*TOF_BINS),TOF_RANGE/(double)(2*TOF_BINS))); // in meters/sec 

                // convert velocity to relativistic kinetic energy
                double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

                enC->Fill(rKE,tofC->GetBinContent(j));
                tofC->SetBinError(j,pow(tofC->GetBinContent(j),0.5));
                enC->SetBinError(j,pow(enC->GetBinContent(j),0.5));
            }
        }
    }

    histoFile->Write();
    histoFile->Close();
}

int histos(string sortedFileName, string vetoedFileName, string histoFileName, vector<string> channelMap)
{
    cout << "Entering ./histos..." << endl;

    TFile* sortedFile = new TFile(sortedFileName.c_str(),"READ");
    TFile* vetoedFile = new TFile(vetoedFileName.c_str(),"READ");
    if(!sortedFile->IsOpen() || !vetoedFile->IsOpen())
    {
        cerr << "Error: failed to open " << sortedFileName << " or " << vetoedFileName << endl;
        exit(1);
    }

    for(string s : channelMap)
    {
        if(s=="-")
        {
            continue;
        }
        string treeName = s;
        orchard.push_back((TTree*)sortedFile->Get(treeName.c_str()));
        treeName += "W";
        orchardW.push_back((TTree*)sortedFile->Get(treeName.c_str()));
    }

    // increase precision to handle outputted times (for troubleshooting)
    cout.precision(13);

    TFile* histoFile = new TFile(histoFileName.c_str(),"RECREATE");

    // prepare the root file with 4 directories, one for each channel
    // these directories will hold basic variable histograms showing the
    // raw data in each tree, plus TOF, x-sections, etc histograms
    fillBasicHistos(histoFile);
    fillVetoedHistos(vetoedFile, histoFile);
    // fill TOF, cross-section, etc. histos for channels 4, 6

    vetoedFile->Close();
    sortedFile->Close();

    histoFile->Close();

    // print out waveforms for first 50 events that occur in both ch4 and ch6
    //matchWaveforms();

    return 0;
}
