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

extern ProcessedEvent procEvent;

vector<TTree*> orchard; // holds DPP channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

vector<TTree*> orchardW; // holds waveform-mode channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree

TH1D* convertTOFtoEnergy(TH1D* tof, string name)
{
    if(!tof)
    {
        cerr << "Error: cannot convert empty TOF histogram to energy units in convertTOFtoEnergy()" << endl;
        exit(1);
    }

    TRandom3 *randomizeBin = new TRandom3();

    TH1D* energy = timeBinsToRKEBins(tof, name); 

    unsigned int tofBins = tof->GetNbinsX();

    for(unsigned int j=0; j<tofBins-1; j++)
    {
        // convert time into neutron velocity based on flight path distance
        double velocity = pow(10.,7.)*FLIGHT_DISTANCE/(tof->GetBinCenter(j)+randomizeBin->Uniform(-TOF_RANGE/(double)(2*TOF_BINS),TOF_RANGE/(double)(2*TOF_BINS))); // in meters/sec 

        // convert velocity to relativistic kinetic energy
        double rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

        energy->Fill(rKE,tof->GetBinContent(j));
        energy->SetBinError(j,pow(energy->GetBinContent(j),0.5));
    }

    return energy;
}

double applyDeadtimeCorrection(TH1D*& rawTOF, vector<double> deadtimesPerBin)
{
    // produce histo showing deadtime for each bin
    //string name = correctedTOF->GetName();
    //name = name + "deadtimeH";
    //TH1D* deadtimeH = (TH1D*)correctedTOF->Clone(name.c_str());

    double sumOfDeadtimes = 0; // for computing average deadtime per bin

    string name = rawTOF->GetName();
    name += "Corrected";
    TH1D* correctedTOFHisto = (TH1D*)rawTOF->Clone(name.c_str());

    for(int i=0; (size_t)i<deadtimesPerBin.size(); i++)
    {
        if(deadtimesPerBin[i]>=1)
        {
            cerr << "Error: attempted to correct for deadtime, but encountered deadtime >100% (deadtime was "
                << deadtimesPerBin[i] << "). Exiting." << endl;
            exit(1);
        }

        correctedTOFHisto->SetBinContent(i+1,rawTOF->GetBinContent(i+1)/(1-deadtimesPerBin[i]));
        sumOfDeadtimes += deadtimesPerBin[i];
    }

    correctedTOFHisto->Write();
    
    return sumOfDeadtimes/deadtimesPerBin.size();
}

vector<double> generateDeadtimeCorrection(TH1D* tof, long totalNumberOfMicros)
{
    vector<double> eventsPerMicroPerBin(TOF_BINS, 0);
    vector<double> deadtimesPerBin(TOF_BINS,0);

    string name = tof->GetName();
    string eventsPerMicroName = name + "EventsPerMicro";
    TH1D* eventsPerMicroPerBinH = new TH1D(eventsPerMicroName.c_str(), eventsPerMicroName.c_str(),TOF_BINS,TOF_LOWER_BOUND,TOF_UPPER_BOUND);

    string deadtimeName = name + "Deadtime";
    TH1D* deadtimeHisto = new TH1D(deadtimeName.c_str(), deadtimeName.c_str(),TOF_BINS,TOF_LOWER_BOUND,TOF_UPPER_BOUND);

    if(totalNumberOfMicros==0)
    {
        cerr << "Error: tried to create deadtimes, but totalNumberOfMicros = 0" << endl;
        return deadtimesPerBin;
    }

    const int deadtimeBins = ((double)TOF_BINS/TOF_RANGE)*DEADTIME_PERIOD;
    const int deadtimeTransitionBins = ((double)TOF_BINS/TOF_RANGE)*DEADTIME_TRANSITION_PERIOD;

    for(int j=0; j<eventsPerMicroPerBin.size(); j++)
    {
        eventsPerMicroPerBin[j] = tof->GetBinContent(j+1)/((double)totalNumberOfMicros);
        eventsPerMicroPerBinH->SetBinContent(j+1,eventsPerMicroPerBin[j]);
    }

    eventsPerMicroPerBinH->Write();

    // find the fraction of the time that the detector is dead for each bin in the micropulse
    // set deadtime fraction base case

    // use deadtime base case to calculate deadtime for remaining bins

    for(int j=0; j<TOF_BINS; j++)
    {
        for(int k=j-(deadtimeBins+deadtimeTransitionBins); k<j; k++)
        {
            if(k<(j-deadtimeBins))
            {
                if(k<0)
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k+TOF_BINS]/*(double)(1-deadtimesPerBin[j])*/*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }

                else
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k]/*(double)(1-deadtimesPerBin[j])*/*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if(k<0)
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k+TOF_BINS]/*(double)(1-deadtimesPerBin[j])*/;
                }

                else
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k]/*(double)(1-deadtimesPerBin[j])*/;
                }
            }
        }

        deadtimesPerBin[j] += eventsPerMicroPerBin[j]/2/*(1-deadtimesPerBin[j])*/; // last bin contributes 1/2 its value

        deadtimeHisto->SetBinContent(j+1,deadtimesPerBin[j]);

    }

    deadtimeHisto->Write();

    return deadtimesPerBin;
}

void fillVetoedHistos(TFile* vetoFile, TFile* histoFile)
{
    for(string s : detectorNames)
    {
        vetoFile->cd();
        //string treeName = s;

        //Uncomment to use vetoed tree
        cout << "Filling vetoed histos for " << s << endl;
        string treeName = s;// + "Clean";

        TTree* tree = (TTree*)vetoFile->Get(treeName.c_str());

        // move to detector directory in histo file
        histoFile->cd("/");
        histoFile->cd(s.c_str());

        // create diagnostic histograms
        TH2I *triangle = new TH2I("triangle","TOF vs. lgQ",TOF_RANGE,TOF_LOWER_BOUND,TOF_UPPER_BOUND/10,8192,0,65536);
        TH2I *triangleRKE = new TH2I("triangleRKE","relativistic KE vs. lgQ",NUMBER_ENERGY_BINS,ENERGY_LOWER_BOUND,ENERGY_UPPER_BOUND,2048,0,65536);
        TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
        TH1I *QRatio = new TH1I("QRatio","short gate Q/long gate Q",1000,0,1);

        TH1I* timeDiffHistoBlank = new TH1I("time since last event, blank","time since last event",TOF_RANGE,0,TOF_RANGE);
        TH1I* timeDiffHistoTarget1 = new TH1I("time since last event, target 1","time since last event",TOF_RANGE,0,TOF_RANGE);
        TH1I* timeDiffHistoTarget2 = new TH1I("time since last event, target 2","time since last event",TOF_RANGE,0,TOF_RANGE);
        TH1I* timeDiffHistoTarget3 = new TH1I("time since last event, target 3","time since last event",TOF_RANGE,0,TOF_RANGE);
        TH1I* timeDiffHistoTarget4 = new TH1I("time since last event, target 4","time since last event",TOF_RANGE,0,TOF_RANGE);
        TH1I* timeDiffHistoTarget5 = new TH1I("time since last event, target 5","time since last event",TOF_RANGE,0,TOF_RANGE);

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

        TH1D *gammaToGammaTimeH = new TH1D("gammaToGammaTimeH","time between consecutive gammas", 1000, -5 , 5);

        TH1I *gammaOffsetHBlank = new TH1I("gammaOffsetHBlank","gammaOffsetHBlank",1000,-5,5);
        TH1I *gammaOffsetHTarget1 = new TH1I("gammaOffsetHTarget1","gammaOffsetHTarget1",1000,-5,5);
        TH1I *gammaOffsetHTarget2 = new TH1I("gammaOffsetHTarget2","gammaOffsetHTarget2",1000,-5,5);
        TH1I *gammaOffsetHTarget3 = new TH1I("gammaOffsetHTarget3","gammaOffsetHTarget3",1000,-5,5);
        TH1I *gammaOffsetHTarget4 = new TH1I("gammaOffsetHTarget4","gammaOffsetHTarget4",1000,-5,5);
        TH1I *gammaOffsetHTarget5 = new TH1I("gammaOffsetHTarget5","gammaOffsetHTarget5",1000,-5,5);

        TH1I *gammaOffsetDiffH = new TH1I("gammaOffsetDiffH","gammaOffsetDiffH",10000,-5,5);

        TH2I* gammaOffsetByGammaNumber = new TH2I("gammaOffsetByGammaNumber","gammaOffsetByGammaNumber",12,0,12,1000,-5,5);

        TH1I* numberOfGammasH = new TH1I("numberOfGammasH", "number of gammas in each macropulse", 20, 0, 20);

        // make target-specific plots
        //vector<Plots*> plots;

        vector<TH1D*> TOFHistos;
        vector<TH1D*> energyHistos;

        for(unsigned int i=0; i<positionNames.size(); i++)
        {
            string TOFName = positionNames[i] + "TOF";
            TOFHistos.push_back(new TH1D(TOFName.c_str(),TOFName.c_str(),TOF_BINS,TOF_LOWER_BOUND,TOF_UPPER_BOUND));

            string energyName =  positionNames[i] + "Energy";
            energyHistos.push_back(timeBinsToRKEBins(TOFHistos.back(),energyName));
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
        const unsigned int gammaGate[2] = {83,88};

        double eventTimeDiff = 0;
        double prevCompleteTime = 0;
        double prevMicroTime = 0;
        double prevRKE = 0;
        double prevlgQ = 0;
        double prevsgQ = 0;
        vector<int> prevWaveform;

        double rKE = 0;
        int subThreshold = 0;

        unsigned int numberOfGammas = 0;
        double gammaTimesSum = 0;

        vector<pair<double,long>> gammaOffsets;
        //vector<pair<double,int>> gammaOffsetsByTarget;

        vector<int> macrosByTarget;

        long totalEntries = tree->GetEntries();

        long prevMacroNo = 0;
        int prevTargetPos = 0;
        double prevGammaTime = 0;

        //TH1I* histoToFill;

        double timeOffset = 0;

        // calculate time offset using gammas
        /*for(long i=0; i<totalEntries; i++)
        {
            tree->GetEntry(i);

            if(procEvent.macroNo!=prevMacroNo)
            {
                //macrosByTarget.push_back(prevTargetPos);

                numberOfGammasH->Fill(numberOfGammas);

                if(numberOfGammas>=1)
                {
                    // calculate average gamma offset for previous macro
                    gammaOffsets.push_back(make_pair(gammaTimesSum/(double)numberOfGammas-(FLIGHT_DISTANCE/C)*pow(10,7),prevMacroNo));
                }

                else
                {
                    if(gammaOffsets.size())
                    {
                        gammaOffsets.push_back(make_pair(gammaOffsets.back().first, prevMacroNo));
                    }

                    else
                    {
                        gammaOffsets.push_back(make_pair(0, prevMacroNo));
                    }
                }

                gammaOffsetByGammaNumber->Fill(numberOfGammas,gammaOffsets.back().first);

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

                histoToFill->Fill(gammaOffsets.back().first);

                numberOfGammas=0;
                gammaTimesSum=0;
            }

            double timeDiff = procEvent.completeTime-procEvent.macroTime;
            microTime = fmod(timeDiff,MICRO_LENGTH);

            // test if gamma
            if(microTime<gammaGate[1] && microTime>gammaGate[0] && procEvent.lgQ > 300)
            {
                numberOfGammas++;
                gammaTimesSum += microTime;

                gammaToGammaTimeH->Fill(microTime-prevGammaTime);

                prevGammaTime = microTime;
            }

            prevMacroNo = procEvent.macroNo;

            if(i%10000==0)
            {
                cout << "Processed " << i << " events through gamma offset correction...\r";
            }
        }*/

        //gammaOffsetHBlank->Write();
        //gammaOffsetHTarget1->Write();
        //gammaOffsetHTarget2->Write();
        //gammaOffsetHTarget3->Write();
        //gammaOffsetHTarget4->Write();
        //gammaOffsetHTarget5->Write();

        //gammaOffsetByGammaNumber->Write();

        //numberOfGammasH->Write();

        //gammaToGammaTimeH->Write();

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

        for(long i=0; i<totalEntries; i++)
        {
            tree->GetEntry(i);

            /*if(i>10000)
            {
                break;
            }*/

            /*if(procEvent.targetPos!=0)
            {
                procEvent.completeTime -= gammaOffsets[procEvent.macroNo].first;
            }*/

            double timeDiff = procEvent.completeTime-procEvent.macroTime;

            eventTimeDiff = procEvent.completeTime-prevCompleteTime;

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
            //if (procEvent.lgQ>500*exp(-(microTime-100)/87))
            /*if (procEvent.lgQ<150)
            {
                continue;
            }*/


            if(procEvent.targetPos==0)        // discard events during target-changer movement
            {
                continue;
            }

            // Fill raw TOF before gates:

            TOFHistos[procEvent.targetPos-1]->Fill(microTime);
            energyHistos[procEvent.targetPos-1]->Fill(rKE);

            // Apply gates:
            if (timeDiff < MACRO_LENGTH && timeDiff > 0 // omit events outside beam-on period
                    //&& (procEvent.sgQ/(double)procEvent.lgQ < 0.495
                    //|| procEvent.sgQ/(double)procEvent.lgQ > 0.505)

                    /*&& procEvent.lgQ < 65000
                    && procEvent.sgQ < 32500

                    && procEvent.lgQ>50
                    && procEvent.sgQ>50*/
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

                timeDiffVEnergy1->Fill(eventTimeDiff,prevRKE);
                timeDiffVEnergy2->Fill(eventTimeDiff,rKE);
                energy1VEnergy2->Fill(prevRKE,rKE);
                time1Vtime2->Fill(prevMicroTime,microTime);
                lgQ1VlgQ2->Fill(prevlgQ,procEvent.lgQ);

                microNoH->Fill(microNo);

                /*****************************************************************/
                // Fill target-specific plots

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

        timeDiffHistoBlank->Write();
        timeDiffHistoTarget1->Write();
        timeDiffHistoTarget2->Write();
        timeDiffHistoTarget3->Write();
        timeDiffHistoTarget4->Write();
        timeDiffHistoTarget5->Write();

        timeDiffVEnergy1->Write();
        timeDiffVEnergy2->Write();
        energy1VEnergy2->Write();
        time1Vtime2->Write();
        lgQ1VlgQ2->Write();
        lgQ1VlgQ2_background->Write();

        for(auto histo : TOFHistos)
        {
            histo->Write();
        }

        for(auto histo : energyHistos)
        {
            histo->Write();
        }

        std::vector<long> macrosPerTarget;
        std::vector<long> microsPerTarget;

        histoFile->cd("/macroTime");

        for(unsigned int i=0; i<tarGates.size()-1; i++)
        {
            macrosPerTarget.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
            microsPerTarget.push_back(macrosPerTarget.back()*(MACRO_LENGTH/MICRO_LENGTH));
        }

        histoFile->cd("/");
        histoFile->cd(s.c_str());

        // perform iterative deadtime correction, until average deadtime changes
        // <0.1%
        for(unsigned int i=0; (unsigned int)i<microsPerTarget.size(); i++)
        {
            vector<double> deadtimeBins = generateDeadtimeCorrection(TOFHistos[i], microsPerTarget[i]);

            double averageDeadtimeDiff = applyDeadtimeCorrection(TOFHistos[i], deadtimeBins) - 0;

            /*while(averageDeadtimeDiff>0.001)
            {
                rawTOF = correctedTOF;

                vector<double> prevDeadtimeBins = deadtimeBins;
                deadtimeBins = generateDeadtimeCorrection(rawTOF, microsPerTarget[i]);

                for(unsigned int j=0; j<deadtimeBins.size(); j++)
                {
                    deadtimeBins[j] = deadtimeBins[j]-prevDeadtimeBins[j];
                }

                correctedTOF = ((TH1D*)plots[i]->getTOFHisto());
                deadtimeH = ((TH1D*)plots[i]->getDeadtimeHisto());
                averageDeadtimeDiff = applyDeadtimeCorrection(rawTOF, correctedTOF, deadtimeH, deadtimeBins) - averageDeadtimeDiff;
            }*/
        }

        // calculate likelihood of double peak in wavelet, for each TOF bin
        const int waveletBins = (TOF_BINS/TOF_RANGE)*WAVELET_PERIOD;

        vector<double> eventsPerMicroPerBin(TOF_BINS);
        vector<double> doublePeakOddsPerBin(TOF_BINS);

        //TH1D* tof = (TH1D*)plots[0]->getTOFHisto();
        /*for(int i=0; i<eventsPerMicroPerBin.size(); i++)
        {
            if(microsPerTarget[0]==0)
            {
                cerr << "Error: can't calculate double peak likelihood with 0 total micros" << endl;
                continue;
            }

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
        */
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

        TH1I* fineTimeHBlank = new TH1I("fineTimeHBlank","fineTimeBlank",6000,-10,10);
        fineTimeHBlank->GetXaxis()->SetTitle("fine time, blank target");

        TH1I* fineTimeHTarget1 = new TH1I("fineTimeHTarget1","fineTimeTarget1",6000,-10,10);
        fineTimeHTarget1->GetXaxis()->SetTitle("fine time, target 1");

        TH1I* fineTimeHTarget2 = new TH1I("fineTimeHTarget2","fineTimeTarget2",6000,-10,10);
        fineTimeHTarget2->GetXaxis()->SetTitle("fine time, target 2");

        TH1I* fineTimeHTarget3 = new TH1I("fineTimeHTarget3","fineTimeTarget3",6000,-10,10);
        fineTimeHTarget3->GetXaxis()->SetTitle("fine time, target 3");

        TH1I* fineTimeHTarget4 = new TH1I("fineTimeHTarget4","fineTimeTarget4",6000,-10,10);
        fineTimeHTarget4->GetXaxis()->SetTitle("fine time, target 4");

        TH1I* fineTimeHTarget5 = new TH1I("fineTimeHTarget5","fineTimeTarget5",6000,-10,10);
        fineTimeHTarget5->GetXaxis()->SetTitle("fine time, target 5");

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

        TH1I* fineTimeH = new TH1I("fineTimeH","fineTimeH",1000,-5,5);
        fineTimeH->GetXaxis()->SetTitle("fine time of event");

        TH1I* diffCompleteTimeH = new TH1I("diffCompleteTimeH","diffCompleteTime",pow(2,20),0,pow(2,26));
        diffCompleteTimeH->GetXaxis()->SetTitle("difference between complete time of consecutive events");

        TH1I* diffMacroCompleteTimesH = new TH1I("diffMacroCompleteTimesH","diffMacroCompleteTime",pow(2,20),0,pow(2,20));
        diffMacroCompleteTimesH->GetXaxis()->SetTitle("difference between complete time of event and its macrotime");

        // create a subdirectory for holding DPP-mode waveform data
        gDirectory->mkdir("waveformsDir","raw DPP waveforms");

        /*************************************************************************/
        // Loop through sorted trees to calculate advanced histogram variables

        string name = t->GetName();
        if(name=="macroTime")
        {
            setBranchesTC(t);
        }

        else
        {
            setBranchesHistos(t);
        }

        long totalEntries = t->GetEntries();

        cout << "Populating " << t->GetName() << " histograms..." << endl;

        TDirectory* waveformsDir = (TDirectory*)gDirectory->Get("waveformsDir");
        waveformsDir->cd();

        int prevMacroNo = 0;
        double prevCompleteTime = 0;

        double timeDiff;
        double microTime;

        // loop through the channel-specific tree and populate histos
        for(long j=0; j<totalEntries; j++)
        {
            t->GetEntry(j);

            timeDiff = procEvent.completeTime-procEvent.macroTime;
            microTime = fmod(timeDiff,MICRO_LENGTH);

            if (procEvent.targetPos == 0)         // discard events during target-changer movement
            {
                continue;
            }

            if(j%50000==0)
            {
                cout << "Processed " << j << " events through basic histos...\r";
                fflush(stdout);

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

            /*if (name!="targetChanger" &&
               ((procEvent.sgQ/(double)procEvent.lgQ > 0.498
                        && procEvent.sgQ/(double)procEvent.lgQ < 0.502)
               || procEvent.completeTime > procEvent.macroTime + MACRO_LENGTH
               || procEvent.completeTime < procEvent.macroTime))
            {
                    continue;
            }*/

            /*&& (prevsgQ/(double)prevlgQ > 0.498
              && prevsgQ/(double)prevlgQ < 0.502)
              && procEvent.lgQ < 65500
              && procEvent.sgQ < 32750
              && procEvent.lgQ > 100
              && procEvent.sgQ > 100
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
                    fineTimeHBlank->Fill(procEvent.fineTime);
                    macroNoHBlank->Fill(procEvent.macroNo);
                    break;
                case 2:
                    fineTimeHTarget1->Fill(procEvent.fineTime);
                    macroNoHTarget1->Fill(procEvent.macroNo);
                    break;
                case 3:
                    fineTimeHTarget2->Fill(procEvent.fineTime);
                    macroNoHTarget2->Fill(procEvent.macroNo);
                    break;
                case 4:
                    fineTimeHTarget3->Fill(procEvent.fineTime);
                    macroNoHTarget3->Fill(procEvent.macroNo);
                    break;
                case 5:
                    fineTimeHTarget4->Fill(procEvent.fineTime);
                    macroNoHTarget4->Fill(procEvent.macroNo);
                    break;
                case 6:
                    fineTimeHTarget5->Fill(procEvent.fineTime);
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

            fineTimeH->Fill(procEvent.fineTime);
            diffCompleteTimeH->Fill(procEvent.completeTime-prevCompleteTime);

            diffMacroCompleteTimesH->Fill(procEvent.completeTime-procEvent.macroTime);

            if(name=="macroTime" || name=="targetChanger")
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
