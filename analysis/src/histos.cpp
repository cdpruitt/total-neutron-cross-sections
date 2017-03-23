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

TH1I* convertTOFtoEn(TH1I* tof, string name)
{
    if(!tof)
    {
        cerr << "Error: cannot convert empty TOF histogram to energy units in convertTOFtoEn()" << endl;
        exit(1);
    }

    TRandom3 *randomizeBin = new TRandom3();

    TH1I* en = timeBinsToRKEBins(tof, name); 

    int tofBins = tof->GetNbinsX();
    for(int j=0; j<tofBins-1; j++)
    {
        // convert time into neutron velocity based on flight path distance
        double velocity = pow(10.,7.)*FLIGHT_DISTANCE/(tof->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)(2*TOF_BINS),TOF_RANGE/(double)(2*TOF_BINS))); // in meters/sec 

        // convert velocity to relativistic kinetic energy
        double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

        en->Fill(rKE,tof->GetBinContent(j));
        en->SetBinError(j,pow(en->GetBinContent(j),0.5));
    }

    return en;
}

// recursive procedure for deadtime-correcting TOF plots
double correctForDeadtime(TH1I* tof, TH1I* oldTOF, TH1I*& newTOF, long totalNumberOfMicros)
{
    vector<double> eventsPerMicroPerBin;
    vector<double> deadtimesPerBin;

    if(totalNumberOfMicros==0)
    {
        cerr << "Error: tried to create deadtimes, but totalNumberOfMicros = 0" << endl;
        return 0;
    }

    const int deadtimeBins = (TOF_BINS/TOF_RANGE)*DEADTIME_PERIOD;

    TH1I* diffTOF = (TH1I*)tof->Clone();
    diffTOF->Add(oldTOF, -1);

    for(int j=0; j<=TOF_BINS; j++)
    {
        eventsPerMicroPerBin.push_back(diffTOF->GetBinContent(j)/((double)totalNumberOfMicros));
        deadtimesPerBin.push_back(0);
    }

    // find the fraction of the time that the detector is dead for each bin in the micropulse
    // set deadtime fraction base case

    // use deadtime base case to calculate deadtime for remaining bins
    for(int j=0; (size_t)j<TOF_BINS; j++)
    {
        for(int k=j-deadtimeBins; k<j; k++)
        {
            if(k<0)
            {
                deadtimesPerBin[j] += eventsPerMicroPerBin[k+TOF_BINS]*(1-deadtimesPerBin[j]);
                continue;
            }

            deadtimesPerBin[j] += eventsPerMicroPerBin[k]*(1-deadtimesPerBin[j]);
        }

        deadtimesPerBin[j] += eventsPerMicroPerBin[j]/2;
    }

    double averageDeadtime = 0;

    string name = tof->GetName();
    name = name + "deadtimeH";
    TH1I* deadtimeH = (TH1I*)tof->Clone(name.c_str());
    newTOF = (TH1I*)tof->Clone();

    for(int j=0; (size_t)j<deadtimesPerBin.size(); j++)
    {
        newTOF->SetBinContent(j,tof->GetBinContent(j)/(1-deadtimesPerBin[j]));

        deadtimeH->SetBinContent(j,pow(10,6)*deadtimesPerBin[j]);
        averageDeadtime += deadtimesPerBin[j];
    }

    deadtimeH->Write();
    newTOF->Write();
    
    return averageDeadtime/deadtimesPerBin.size();
}

void fillVetoedHistos(TFile* vetoFile, TFile* histoFile)
{
    for(string s : detectorNames)
    {
        vetoFile->cd();
        //string treeName = s;

        //Uncomment to use vetoed tree
        cout << "Filling vetoed histos for " << s << endl;
        string treeName = s + "Clean";

        TTree* tree = (TTree*)vetoFile->Get(treeName.c_str());

        // move to detector directory in histo file
        histoFile->cd("/");
        histoFile->cd(s.c_str());

        // create diagnostic histograms
        TH2I *triangle = new TH2I("triangle","TOF vs. lgQ",TOF_RANGE,TOF_LOWER_BOUND,TOF_UPPER_BOUND,2048,0,65536);
        TH2I *triangleRKE = new TH2I("triangleRKE","relativistic KE vs. lgQ",NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,2048,0,65536);
        TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
        TH1I *QRatio = new TH1I("QRatio","short gate Q/long gate Q",1000,0,1);

        TH1I* timeDiffHisto = new TH1I("time since last event","time since last event",TOF_RANGE,0,TOF_RANGE);

        TH2I *timeDiffVEnergy1 = new TH2I("time difference vs. energy of first","time difference vs. energy of first",TOF_RANGE,0,TOF_RANGE,500,2,700);
        TH2I *timeDiffVEnergy2 = new TH2I("time difference vs. energy of second","time difference vs. energy of second",TOF_RANGE,0,TOF_RANGE,NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);
        TH2I *energy1VEnergy2 = new TH2I("energy of first vs. energy of second","energy of first vs. energy of second",10*NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,10*NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND);
        TH2I *time1Vtime2 = new TH2I("time of first vs. time of second","time of first vs. time of second",TOF_RANGE,0,TOF_RANGE,TOF_RANGE,0,TOF_RANGE);
        TH2I *lgQ1VlgQ2 = new TH2I("lgQ of first vs. lgQ of second","lgQ of first vs. lgQ of second",1024,0,65536,1024,0,65536);
        TH2I *lgQ1VlgQ2_background = new TH2I("lgQ of first vs. lgQ of second_background","lgQ of first vs. lgQ of second_background",1024,0,65536,1024,0,65536);

        TH1I *doublePeakH = new TH1I("likelihood of double peak in wavelet","likelihood of double peak in wavelet",TOF_BINS,0,TOF_RANGE);

        TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
        microNoH->GetXaxis()->SetTitle("micropulse number of each event");

        TH1I *orderInMicroH = new TH1I("orderInMicroH","orderInMicro",8,0,8);
        orderInMicroH->GetXaxis()->SetTitle("order in micro");

        TH1I *gammaOffsetHBlank = new TH1I("gammaOffsetHBlank","gammaOffsetHBlank",1000,-5,5);
        TH1I *gammaOffsetHTarget1 = new TH1I("gammaOffsetHTarget1","gammaOffsetHTarget1",1000,-5,5);
        TH1I *gammaOffsetHTarget2 = new TH1I("gammaOffsetHTarget2","gammaOffsetHTarget2",1000,-5,5);
        TH1I *gammaOffsetHTarget3 = new TH1I("gammaOffsetHTarget3","gammaOffsetHTarget3",1000,-5,5);
        TH1I *gammaOffsetHTarget4 = new TH1I("gammaOffsetHTarget4","gammaOffsetHTarget4",1000,-5,5);
        TH1I *gammaOffsetHTarget5 = new TH1I("gammaOffsetHTarget5","gammaOffsetHTarget5",1000,-5,5);

        TH1I *gammaOffsetDiffH = new TH1I("gammaOffsetDiffH","gammaOffsetDiffH",10000,-5,5);

        TH2I* gammaOffsetByGammaNumber = new TH2I("gammaOffsetByGammaNumber","gammaOffsetByGammaNumber",12,0,12,1000,-5,5);

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
        long runningMicroNo = 0;

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
        int gammaGate[2] = {82,89};

        double eventTimeDiff = 0;
        double prevCompleteTime = 0;
        double prevMicroTime = 0;
        double prevRKE = 0;
        double prevlgQ = 0;
        double prevsgQ = 0;
        vector<int> prevWaveform;

        const int WAVEFORM_THRESHOLD = 15500;

        double rKE = 0;
        bool brokeThreshold = false;
        int subThreshold = 0;

        unsigned int numberOfGammas = 0;
        double gammaTimes = 0;
        vector<double> gammaOffsets;
        vector<pair<double,int>> gammaOffsetsByTarget;

        vector<int> macrosByTarget;

        long totalEntries = tree->GetEntries();

        long prevMacroNo = 0;
        int prevTargetPos = 0;

        TH1I* histoToFill;

        // calculate time offset using gammas
        for(long i=0; i<totalEntries; i++)
        {
            tree->GetEntry(i);

            /*if(procEvent.macroNo!=prevMacroNo)
            {
                macrosByTarget.push_back(prevTargetPos);

                // calculate average gamma offset for previous macro
                if(numberOfGammas>=1)
                {
                    gammaOffsets.push_back(gammaTimes/numberOfGammas-(FLIGHT_DISTANCE/C)*pow(10,7));
                }

                else
                {
                    //gammaOffsets.push_back(0);
                }

                gammaOffsetByGammaNumber->Fill(numberOfGammas,gammaOffsets.back());

                numberOfGammas=0;
                gammaTimes=0;
            }*/

            double timeDiff = procEvent.completeTime-procEvent.macroTime;
            microTime = fmod(timeDiff,MICRO_LENGTH);

            // test if gamma
            if(microTime<gammaGate[1]&&microTime>gammaGate[0] && procEvent.targetPos!=0)
            {
                /*numberOfGammas++;
                gammaTimes += microTime;
                */
                //gammasPerTarget[procEvent.targetPos-1]++;
                gammaOffsetsByTarget.push_back(make_pair(microTime-(FLIGHT_DISTANCE/C)*pow(10,7),procEvent.targetPos));

                switch(procEvent.targetPos)
                {
                    case 1:
                        histoToFill = gammaOffsetHBlank;
                        break;
                    case 2:
                        histoToFill = gammaOffsetHTarget1;
                        break;
                    case 3:
                        histoToFill = gammaOffsetHTarget2;
                        break;
                    case 4:
                        histoToFill = gammaOffsetHTarget3;
                        break;
                    case 5:
                        histoToFill = gammaOffsetHTarget4;
                        break;
                    case 6:
                        histoToFill = gammaOffsetHTarget5;
                        break;
                    default:
                        break;
                }

                histoToFill->Fill(gammaOffsetsByTarget.back().first);

                gammaOffsetDiffH->Fill(gammaOffsetsByTarget.back().first);

                //prevMacroNo = procEvent.macroNo;
                //prevTargetPos = procEvent.targetPos;
            }
        }

        gammaOffsetHBlank->Write();
        gammaOffsetHTarget1->Write();
        gammaOffsetHTarget2->Write();
        gammaOffsetHTarget3->Write();
        gammaOffsetHTarget4->Write();
        gammaOffsetHTarget5->Write();

        gammaOffsetDiffH->Write();

        gammaOffsetByGammaNumber->Write();

        // add the very last macro's gammas into the gamma offsets
        //gammaOffsets.push_back(gammaTimes/numberOfGammas-(FLIGHT_DISTANCE/C)*pow(10,7));

        // add the very last macro's target position into counter
        //macrosByTarget.push_back(procEvent.targetPos);

        //vector<int> totalMacrosInTarget(positionNames.size(),0);

        /*for(int i=0; i<macrosByTarget.size(); i++)
        {
            if(macrosByTarget[i]==0)
            {
                continue;
            }

            gammaOffsetsByTarget[macrosByTarget[i]-1] += gammaOffsets[i];
            totalMacrosInTarget[macrosByTarget[i]-1]++;
        }*/

        vector<double> gammaOffsetAverageByTarget(positionNames.size(),0);
        vector<long> gammasByTarget(positionNames.size(),0);

        for(int i=0; i<gammaOffsetsByTarget.size(); i++)
        {
            int targetPos = gammaOffsetsByTarget[i].second;
            gammaOffsetAverageByTarget[targetPos-1] += gammaOffsetsByTarget[i].first;
            gammasByTarget[targetPos-1]++;
        }

        for(int i=0; i<gammaOffsetAverageByTarget.size(); i++)
        {
            gammaOffsetAverageByTarget[i] = gammaOffsetAverageByTarget[i]/(double)gammasByTarget[i];
        }

        for(long i=0; i<totalEntries; i++)
        {
            tree->GetEntry(i);

            /*if(i>10000)
            {
                break;
            }*/

            // calculate time since start of macro (includes time offsets)
            /*if(procEvent.targetPos==1)
            {
                procEvent.completeTime = procEvent.completeTime + 0.5;
            }*/

            //procEvent.completeTime -= gammaOffsets[procEvent.macroNo];

            if(procEvent.targetPos!=0)
            {
                //procEvent.completeTime -= gammaOffsetAverageByTarget[procEvent.targetPos-1];
            }

            double timeDiff = procEvent.completeTime-procEvent.macroTime;
            eventTimeDiff = procEvent.completeTime-prevCompleteTime;

            brokeThreshold = false;

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
            rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

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

                runningMicroNo++;
            }

            // GATE: discard events with too low of an integrated charge for their energy
            //if (procEvent.lgQ>50)
            //if (procEvent.lgQ>500*exp(-(microTime-100)/87))
            /*if (procEvent.lgQ>100*rKE || procEvent.lgQ<10*rKE)
              {
              continue;
              }*/

            /*for(int j=0; j<procEvent.waveform->size(); j++)
            {
                if(procEvent.waveform->at(j) < WAVEFORM_THRESHOLD)
                {
                    brokeThreshold = true;
                    break;
                }
            }

            if(!brokeThreshold)
            {
                subThreshold++;
            }*/

            // Apply gates:
            if (timeDiff < MACRO_LENGTH && timeDiff > 0 // omit events outside beam-on period
                    && procEvent.targetPos != 0         // discard events during target-changer movement
                    /*&& procEvent.sgQ/(double)procEvent.lgQ > 0.9
                    && procEvent.sgQ/(double)procEvent.lgQ < 0.99
                    */
                    && (procEvent.sgQ/(double)procEvent.lgQ < 0.498
                    || procEvent.sgQ/(double)procEvent.lgQ > 0.502)
                    /*&& (procEvent.sgQ/(double)procEvent.lgQ > 0.498
                    && procEvent.sgQ/(double)procEvent.lgQ < 0.502)
                    */
                    && procEvent.lgQ < 65500
                    && procEvent.sgQ < 32750
                    && procEvent.lgQ > 50
                    && procEvent.sgQ > 50
                    //&& brokeThreshold

                    //&& procEvent.lgQ>200
                    //&& prevlgQ>200
                    //&& procEvent.lgQ<65500 && procEvent.sgQ<32750  // discard events with unphysical integrated charges
                    //&& (procEvent.sgQ/(double)procEvent.lgQ<0.25 || procEvent.sgQ/(double)procEvent.lgQ>0.35)  // discard events outside accepted range of sgQ/lgQ
                    //&& procEvent.lgQ>100    // discard gammas at lowest range of energy
               )
            {

                /*****************************************************************/
                // Fill troubleshooting plots with event variables (rKE, microtime, etc.)

                triangle->Fill(microTime,procEvent.lgQ);
                sgQlgQ->Fill(procEvent.sgQ,procEvent.lgQ);
                QRatio->Fill(procEvent.sgQ/(double)procEvent.lgQ);
                triangleRKE->Fill(rKE,procEvent.lgQ);
                orderInMicroH->Fill(orderInMicro);

                /*if(abs(eventTimeDiff)>300 && abs(eventTimeDiff)<400)
                {
                    lgQ1VlgQ2_background->Fill(prevlgQ,procEvent.lgQ);
                }*/

                //if(abs(eventTimeDiff)>300 && abs(eventTimeDiff)<675)
                {
                    timeDiffHisto->Fill(eventTimeDiff);
                    timeDiffVEnergy1->Fill(eventTimeDiff,prevRKE);
                    timeDiffVEnergy2->Fill(eventTimeDiff,rKE);
                    energy1VEnergy2->Fill(prevRKE,rKE);
                    time1Vtime2->Fill(prevMicroTime,microTime);
                    lgQ1VlgQ2->Fill(prevlgQ,procEvent.lgQ);
                }

                microNoH->Fill(microNo);

                /*****************************************************************/
                // Fill target-specific plots

                ((TH1I*)plots[procEvent.targetPos-1]->getTOFHisto())->Fill(microTime);
                ((TH1I*)plots[procEvent.targetPos-1]->getEnergyHisto())->Fill(rKE);

                /*if(!gammaInMicro)
                  {
                  plots.energyHistosNoGamma[procEvent.targetPos-1]->Fill(rKE);
                  }

                  if(orderInMicro==1)
                  {
                  plots.TOFHistosFirstInMicro[procEvent.targetPos-1]->Fill(microTime);
                  }
                  */

                /*if(microTime>100 && microTime<110 && procEvent.waveform->size() > 0 && i%1000==0)
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
                }*/
            }

            //prevWaveform = *procEvent.waveform;
            prevlgQ = procEvent.lgQ;
            prevsgQ = procEvent.sgQ;
            prevMicroTime = microTime;
            prevCompleteTime = procEvent.completeTime;
            prevRKE = rKE;

            if(i%10000==0)
            {
                cout << "Processed " << i << " events...\r";
            }
        }

        triangle->Write();
        sgQlgQ->Write();
        QRatio->Write();
        triangleRKE->Write();
        microNoH->Write();
        orderInMicroH->Write();

        timeDiffHisto->Write();
        timeDiffVEnergy1->Write();
        timeDiffVEnergy2->Write();
        energy1VEnergy2->Write();
        time1Vtime2->Write();
        lgQ1VlgQ2->Write();
        lgQ1VlgQ2_background->Write();

        //cout << endl << "Subthreshold = " << subThreshold << endl;

        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            ((TH1I*)plots[i]->getTOFHisto())->Write();
            ((TH1I*)plots[i]->getEnergyHisto())->Write();
        }

        std::vector<long> macrosPerTarget;
        std::vector<long> microsPerTarget;

        long totalMacros = 0;

        histoFile->cd("/targetChanger");

        for(unsigned int i=0; i<tarGates.size()-1; i++)
        {
            macrosPerTarget.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
            totalMacros+=macrosPerTarget.back();

            microsPerTarget.push_back(macrosPerTarget.back()*(MACRO_LENGTH/MICRO_LENGTH));
            totalMicros+=microsPerTarget.back();
        }

        histoFile->cd("/");
        histoFile->cd(s.c_str());

        // perform iterative deadtime correction, until average deadtime changes
        // <0.1%
        for(unsigned int i=0; (unsigned int)i<microsPerTarget.size(); i++)
        {
            double averageDeadtimeDiff = 0;
            TH1I* tof = ((TH1I*)plots[i]->getTOFHisto());
            TH1I* emptyTOF = (TH1I*)tof->Clone();
            emptyTOF->Reset();

            TH1I* newTOF;

            averageDeadtimeDiff = correctForDeadtime(tof, emptyTOF, newTOF, microsPerTarget[i]) - averageDeadtimeDiff;

            while(averageDeadtimeDiff>0.001)
            {
                averageDeadtimeDiff = correctForDeadtime(newTOF, tof, newTOF, microsPerTarget[i]) - averageDeadtimeDiff;
            }
        }

        // calculate likelihood of double peak in wavelet, for each TOF bin
        if(microsPerTarget[0]==0)
        {
            cerr << "Error: can't calculate double peak likelihood with 0 total micros" << endl;
            continue;
        }

        const int waveletBins = (TOF_BINS/TOF_RANGE)*WAVELET_PERIOD;

        vector<double> eventsPerMicroPerBin(TOF_BINS);
        vector<double> doublePeakOddsPerBin(TOF_BINS);

        TH1I* tof = (TH1I*)plots[0]->getTOFHisto();

        for(int i=0; i<eventsPerMicroPerBin.size(); i++)
        {
            eventsPerMicroPerBin[i] = tof->GetBinContent(i)/(double)microsPerTarget[0];
        }

        for(int i=0; i<TOF_BINS; i++)
        {
            for(int j=i; j<i+waveletBins/2; j++)
            {
                if(j>TOF_BINS)
                {
                    doublePeakOddsPerBin[i] += eventsPerMicroPerBin[j-TOF_BINS];
                    continue;
                }

                doublePeakOddsPerBin[i] += eventsPerMicroPerBin[j];
            }
        }

        for(int i=0; i<doublePeakOddsPerBin.size(); i++)
        {
            doublePeakH->SetBinContent(i,pow(10,3)*doublePeakOddsPerBin[i]);
        }

        doublePeakH->Write();
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

        TH1I* macroNoHBlank = new TH1I("macroNoHBlank","macroNoBlank",200000,0,200000);
        macroNoHBlank->GetXaxis()->SetTitle("macropulse number of each event, blank target");

        TH1I* macroNoHTarget1 = new TH1I("macroNoHTarget1","macroNoTarget1",200000,0,200000);
        macroNoHTarget1->GetXaxis()->SetTitle("macropulse number of each event, target 1");

        TH1I* macroNoHTarget2 = new TH1I("macroNoHTarget2","macroNoTarget2",200000,0,200000);
        macroNoHTarget2->GetXaxis()->SetTitle("macropulse number of each event, target 2");

        TH1I* macroNoHTarget3 = new TH1I("macroNoHTarget3","macroNoTarget3",200000,0,200000);
        macroNoHTarget3->GetXaxis()->SetTitle("macropulse number of each event, target 3");

        TH1I* macroNoHTarget4 = new TH1I("macroNoHTarget4","macroNoTarget4",200000,0,200000);
        macroNoHTarget4->GetXaxis()->SetTitle("macropulse number of each event, target 4");

        TH1I* macroNoHTarget5 = new TH1I("macroNoHTarget5","macroNoTarget5",200000,0,200000);
        macroNoHTarget5->GetXaxis()->SetTitle("macropulse number of each event, target 5");

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

        TH1I* macroNoDiffH = new TH1I("macroNoDiffH","macroNoDiff",100,0,100);
        macroNoDiffH->GetXaxis()->SetTitle("difference between macropulse numbers of consecutive events");

        TH1I* completeTimeH = new TH1I("completeTimeH","completeTime",pow(2,20),0,pow(2,32));
        completeTimeH->GetXaxis()->SetTitle("complete time of event");

        TH1I* diffCompleteTimeH = new TH1I("diffCompleteTimeH","diffCompleteTime",pow(2,20),0,pow(2,26));
        diffCompleteTimeH->GetXaxis()->SetTitle("difference between complete time of consecutive events");

        TH1I* diffMacroCompleteTimesH = new TH1I("diffMacroCompleteTimesH","diffMacroCompleteTime",pow(2,20),0,pow(2,20));
        diffMacroCompleteTimesH->GetXaxis()->SetTitle("difference between complete time of event and its macrotime");

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

        double timeDiff = 0;
        int prevMacroNo = 0;
        double prevCompleteTime = 0;

        // loop through the channel-specific tree and populate histos
        for(long j=0; j<totalEntries; j++)
        {
            t->GetEntry(j);

            if (procEvent.targetPos == 0)         // discard events during target-changer movement
            {
                continue;
            }

            if (name!="targetChanger" &&
               ((procEvent.sgQ/(double)procEvent.lgQ > 0.498
                        && procEvent.sgQ/(double)procEvent.lgQ < 0.502)
               || procEvent.completeTime > procEvent.macroTime + MACRO_LENGTH
               || procEvent.completeTime < procEvent.macroTime))
            {
                    continue;
            }

            /*&& (prevsgQ/(double)prevlgQ > 0.498
              && prevsgQ/(double)prevlgQ < 0.502)
              && procEvent.lgQ < 65500
              && procEvent.sgQ < 32750
              && procEvent.lgQ > 100
              && procEvent.sgQ > 100
              && brokeThreshold
              */

            //&& eventTimeDiff > 425
            //&& eventTimeDiff < 500
            //&& procEvent.lgQ>200
            //&& prevlgQ>200
            //&& procEvent.lgQ<65500 && procEvent.sgQ<32750  // discard events with unphysical integrated charges
            //&& (procEvent.sgQ/(double)procEvent.lgQ<0.25 || procEvent.sgQ/(double)procEvent.lgQ>0.35)  // discard events outside accepted range of sgQ/lgQ
            //&& procEvent.lgQ>100    // discard gammas at lowest range of energy

            switch(procEvent.targetPos)
            {
                case 1:
                    macroNoHBlank->Fill(procEvent.macroNo);
                    break;
                case 2:
                    macroNoHTarget1->Fill(procEvent.macroNo);
                    break;
                case 3:
                    macroNoHTarget2->Fill(procEvent.macroNo);
                    break;
                case 4:
                    macroNoHTarget3->Fill(procEvent.macroNo);
                    break;
                case 5:
                    macroNoHTarget4->Fill(procEvent.macroNo);
                    break;
                case 6:
                    macroNoHTarget5->Fill(procEvent.macroNo);
                    break;
                default:
                    break;
            }

            macroNoH->Fill(procEvent.macroNo);
            targetPosH->Fill(procEvent.targetPos);
            macroTimeH->Fill(procEvent.macroTime);

            macroNoDiffH->Fill(procEvent.macroNo-prevMacroNo);

            completeTimeH->Fill(procEvent.completeTime);
            diffCompleteTimeH->Fill(procEvent.completeTime-prevCompleteTime);

            diffMacroCompleteTimesH->Fill(procEvent.completeTime-procEvent.macroTime);

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

            /*if(j%1000==0)
              {
              cout << "Processed " << j << " events...\r";
              fflush(stdout);
              }*/

            prevMacroNo = procEvent.macroNo;
            prevCompleteTime = procEvent.completeTime;

            //cout << procEvent.completeTime <<  " " << procEvent.macroNo << endl;
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
