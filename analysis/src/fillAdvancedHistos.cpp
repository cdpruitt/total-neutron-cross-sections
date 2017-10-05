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
#include "../include/dataStructures.h"
#include "../include/branches.h"
#include "../include/plots.h"
#include "../include/fillAdvancedHistos.h"
#include "../include/waveform.h"
#include "../include/config.h"

using namespace std;

// experimentally-determined digitizer deadtime
const int DEADTIME_PERIOD = 150; // in ns
const int DEADTIME_TRANSITION_PERIOD = 15; // in ns

const unsigned int MACRO_EVENTS_LOW_THRESHOLD = 100;
const unsigned int MACRO_EVENTS_HIGH_THRESHOLD = 280;

extern Config config;

TH1D* convertTOFtoEnergy(TH1D* tof, string name)
{
    TH1D* energy = timeBinsToRKEBins(tof, name); 

    if(!tof)
    {
        cerr << "Error: cannot convert empty TOF histogram to energy units in convertTOFtoEnergy()" << endl;
        return energy;
    }

    unsigned int tofBins = tof->GetNbinsX()-2;

    TRandom3 *randomizeBin = new TRandom3();

    for(unsigned int j=1; j<tofBins+1; j++)
    {
        // convert time into neutron velocity based on flight path distance
        double velocity = pow(10.,7.)*(config.facilityConfig.FLIGHT_DISTANCE)
            /(tof->GetBinCenter(j)
                    /*+randomizeBin->Uniform(
                        -(config.plotConfig.TOF_RANGE)/(double)(2*tofBins),
                         (config.plotConfig.TOF_RANGE)/(double)(2*tofBins))*/
             ); // in meters/sec 

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
    vector<double> eventsPerMicroPerBin(config.plotConfig.TOF_BINS, 0);
    vector<double> deadtimesPerBin(config.plotConfig.TOF_BINS,0);

    string name = tof->GetName();
    string eventsPerMicroName = name + "EventsPerMicro";
    TH1D* eventsPerMicroPerBinH = new TH1D(eventsPerMicroName.c_str(), eventsPerMicroName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND);

    string deadtimeName = name + "Deadtime";
    TH1D* deadtimeHisto = new TH1D(deadtimeName.c_str(), deadtimeName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND);

    if(totalNumberOfMicros==0)
    {
        cerr << "Error: tried to create deadtimes, but totalNumberOfMicros = 0" << endl;
        return deadtimesPerBin;
    }

    const int deadtimeBins = ((double)config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*DEADTIME_PERIOD;
    const int deadtimeTransitionBins = ((double)config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*DEADTIME_TRANSITION_PERIOD;

    for(int j=0; j<eventsPerMicroPerBin.size(); j++)
    {
        eventsPerMicroPerBin[j] = tof->GetBinContent(j+1)/((double)totalNumberOfMicros);
        eventsPerMicroPerBinH->SetBinContent(j+1,eventsPerMicroPerBin[j]);
    }

    eventsPerMicroPerBinH->Write();

    // find the fraction of the time that the detector is dead for each bin in the micropulse
    // set deadtime fraction base case

    // use deadtime base case to calculate deadtime for remaining bins

    for(int j=0; j<config.plotConfig.TOF_BINS; j++)
    {
        for(int k=j-(deadtimeBins+deadtimeTransitionBins); k<j; k++)
        {
            if(k<(j-deadtimeBins))
            {
                if(k<0)
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k+config.plotConfig.TOF_BINS]*(double)(1-deadtimesPerBin[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }

                else
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimesPerBin[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if(k<0)
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k+config.plotConfig.TOF_BINS]*(double)(1-deadtimesPerBin[j]);
                }

                else
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimesPerBin[j]);
                }
            }
        }

        deadtimesPerBin[j] += (eventsPerMicroPerBin[j]/2)*(1-deadtimesPerBin[j]); // last bin contributes 1/2 its value

        // scale up events per micro on-the-fly
        //eventsPerMicroPerBin[j] *= (1-deadtimesPerBin[j]);

        deadtimeHisto->SetBinContent(j+1,deadtimesPerBin[j]);
    }

    deadtimeHisto->Write();

    return deadtimesPerBin;
}

int fillAdvancedHistos(string inputFileName, string treeName, string outputFileName)
{
    cout << "Filling advanced histograms for tree \"" << treeName << "\"..." << endl;

    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    TTree* tree = (TTree*)inputFile->Get(treeName.c_str());
    if(!tree)
    {
        cerr << "Error: tried to populate advanced histos, but failed to find " << treeName << " in " << inputFileName << endl;
        return 1;
    }

    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    TDirectory* directory = (TDirectory*)outputFile->GetDirectory(treeName.c_str());
    if(!directory)
    {
        directory = outputFile->mkdir(treeName.c_str(),treeName.c_str());
    }

    outputFile->cd(treeName.c_str());

    // create diagnostic histograms
    TH2I *triangle = new TH2I("triangle","TOF vs. lgQ",config.plotConfig.TOF_RANGE,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND/10,8192,0,65536);
    TH2I *triangleRKE = new TH2I("triangleRKE","relativistic KE vs. lgQ",config.plotConfig.NUMBER_ENERGY_BINS,config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND,2048,0,65536);
    TH2I *sgQlgQ = new TH2I("sgQlgQ","short gate Q vs. long gate Q",2048,0,65536,2048,0,65536);
    TH1I *QRatio = new TH1I("QRatio","short gate Q/long gate Q",1000,0,1);

    TH1I* timeDiffHistoBlank = new TH1I("time since last event, blank","time since last event",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH1I* timeDiffHistoTarget1 = new TH1I("time since last event, target 1","time since last event",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH1I* timeDiffHistoTarget2 = new TH1I("time since last event, target 2","time since last event",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH1I* timeDiffHistoTarget3 = new TH1I("time since last event, target 3","time since last event",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH1I* timeDiffHistoTarget4 = new TH1I("time since last event, target 4","time since last event",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH1I* timeDiffHistoTarget5 = new TH1I("time since last event, target 5","time since last event",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);

    TH1I* goodMacroNoH = new TH1I("good macros","good macros",200000,0,200000);
    TH1I* goodMacroNoHBlank = new TH1I("good macros, blank","good macros, blank",200000,0,200000);
    TH1I* goodMacroNoHTarget1 = new TH1I("good macros, target 1","good macros, target 1",200000,0,200000);
    TH1I* goodMacroNoHTarget2 = new TH1I("good macros, target 2","good macros, target 2",200000,0,200000);
    TH1I* goodMacroNoHTarget3 = new TH1I("good macros, target 3","good macros, target 3",200000,0,200000);
    TH1I* goodMacroNoHTarget4 = new TH1I("good macros, target 4","good macros, target 4",200000,0,200000);
    TH1I* goodMacroNoHTarget5 = new TH1I("good macros, target 5","good macros, target 5",200000,0,200000);

    TH2I *timeDiffVEnergy1 = new TH2I("time difference vs. energy of first","time difference vs. energy of first",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE,500,2,700);
    TH2I *timeDiffVEnergy2 = new TH2I("time difference vs. energy of second","time difference vs. energy of second",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE,config.plotConfig.NUMBER_ENERGY_BINS,config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND);
    TH2I *energy1VEnergy2 = new TH2I("energy of first vs. energy of second","energy of first vs. energy of second",10*config.plotConfig.NUMBER_ENERGY_BINS,config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND,10*config.plotConfig.NUMBER_ENERGY_BINS,config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND);
    TH2I *time1Vtime2 = new TH2I("time of first vs. time of second","time of first vs. time of second",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE,config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH2I *lgQ1VlgQ2 = new TH2I("lgQ of first vs. lgQ of second","lgQ of first vs. lgQ of second",1024,0,65536,1024,0,65536);
    TH2I *lgQ1VlgQ2_background = new TH2I("lgQ of first vs. lgQ of second_background","lgQ of first vs. lgQ of second_background",1024,0,65536,1024,0,65536);

    TH1I *doublePeakH = new TH1I("likelihood of double peak in wavelet","likelihood of double peak in wavelet",config.plotConfig.TOF_BINS,0,config.plotConfig.TOF_RANGE);

    TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
    microNoH->GetXaxis()->SetTitle("micropulse number of each event");

    TH1I *orderInMicroH = new TH1I("orderInMicroH","orderInMicro",8,0,8);
    orderInMicroH->GetXaxis()->SetTitle("order in micro");

    TH1D *gammaToGammaTimeH = new TH1D("gammaToGammaTimeH","time between consecutive gammas", 1000, -5 , 5);

    TH1I *gammaOffsetH = new TH1I("gammaOffsetH","gammaOffsetH",1000,-3,3);
    
    TH2I* gammaOffsetByGammaNumber = new TH2I("gammaOffsetByGammaNumber","gammaOffsetByGammaNumber",12,0,12,1000,-5,5);

    TH1I* numberOfGammasH = new TH1I("numberOfGammasH", "number of gammas in each macropulse", 20, 0, 20);

    TH1D* ratioBadMacrosH = new TH1D("ratio of bad macros","ratio of bad macros",7,0,7);

    vector<TH1D*> TOFHistos;
    vector<TH1D*> energyHistos;

    for(unsigned int i=0; i<config.targetConfig.TARGET_ORDER.size(); i++)
    {
        string TOFName = config.targetConfig.TARGET_ORDER[i] + "TOF";
        TOFHistos.push_back(new TH1D(TOFName.c_str(),TOFName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND));

        string energyName =  config.targetConfig.TARGET_ORDER[i] + "Energy";
        energyHistos.push_back(timeBinsToRKEBins(TOFHistos.back(),energyName));
    }

    // reattach to the channel-specific tree for reading out data
    ProcessedEvent procEvent;
    setBranchesHistos(tree, procEvent);

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

    vector<int> macrosByTarget;

    long totalEntries = tree->GetEntries();

    long prevMacroNo = 0;
    int prevTargetPos = 0;

    double timeOffset = 0;

    // Indicate the range of times considered to be gamma rays (for the purposes of
    // counting gamma rays)
    const double GAMMA_TIME = pow(10,7)*config.facilityConfig.FLIGHT_DISTANCE/C;
    cout << "Gamma time is " <<  GAMMA_TIME << "." << endl;

    tree->GetEntry(totalEntries-1);

    vector<double> gammaOffsets(procEvent.macroNo,0);

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if(procEvent.macroNo!=prevMacroNo)
        {
            numberOfGammasH->Fill(numberOfGammas);

            if(numberOfGammas>0)
            {
                // calculate average gamma offset for previous macro
                gammaOffsets[prevMacroNo] =
                    (gammaTimesSum/(double)numberOfGammas-GAMMA_TIME);
            }

            gammaOffsetByGammaNumber->Fill(numberOfGammas,gammaOffsets[prevMacroNo]);

            gammaOffsetH->Fill(gammaOffsets[prevMacroNo]);
            

            numberOfGammas=0;
            gammaTimesSum=0;
        }

        double timeDiff = procEvent.completeTime-procEvent.macroTime;
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        // test if gamma
        if(abs(microTime-GAMMA_TIME)<2.5)
        {
            numberOfGammas++;
            gammaTimesSum += microTime;
        }

        prevMacroNo = procEvent.macroNo;

        if(i%10000==0)
        {
            cout << "Processed " << i << " events through gamma offset correction...\r";
        }
    }

    gammaOffsetH->Write();
    gammaOffsetByGammaNumber->Write();

    numberOfGammasH->Write();

    int currentMacroNo = -1;
    unsigned int numberOfEventsInCurrentMacro = 0;
    unsigned int targetPosOfCurrentMacro = 0;
    vector<long> goodMacrosByTarget(7,0);
    vector<long> badMacrosByTarget(7,0);

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if(procEvent.targetPos==0)        // discard events during target-changer movement
        {
            continue;
        }

        /*if(currentMacroNo!=procEvent.macroNo)
          {
        // new macropulse; test if "good macropulse"
        currentMacroNo = procEvent.macroNo;
        long j = i;

        while(currentMacroNo==procEvent.macroNo)
        {
        numberOfEventsInCurrentMacro = procEvent.eventNo;
        targetPosOfCurrentMacro = procEvent.targetPos;

        j++;
        if(j>=totalEntries)
        {
        // reached end of tree
        break;
        }

        tree->GetEntry(j);
        }

        if(j>=totalEntries)
        {
        break;
        }

        //cout << "macroNo = " << currentMacroNo << ", number of events = " << numberOfEventsInCurrentMacro << endl;

        if(numberOfEventsInCurrentMacro > MACRO_EVENTS_LOW_THRESHOLD
        && numberOfEventsInCurrentMacro < MACRO_EVENTS_HIGH_THRESHOLD)
        {
        // this macropulse is OK; start again at the first event and populate histos
        goodMacrosByTarget[targetPosOfCurrentMacro]++;
        tree->GetEntry(i);
        }

        else
        {
        // this macropulse is not OK; start at the first event of the next macro
        badMacrosByTarget[targetPosOfCurrentMacro]++;
        i = j-1;
        continue;
        }
        }*/

        procEvent.completeTime -= gammaOffsets[procEvent.macroNo];

        double timeDiff = procEvent.completeTime-procEvent.macroTime;

        eventTimeDiff = procEvent.completeTime-prevCompleteTime;

        /*****************************************************************/
        // Calculate event properties

        // find which micropulse the event is in and the time since
        // the start of the micropulse

        // first, save previous event's micropulse number (we'll need this
        // to calculate event ordering in each micropulse)
        prevMicroNo = microNo;

        microNo = floor(timeDiff/config.facilityConfig.MICRO_LENGTH);
        microTime = fmod(timeDiff,config.facilityConfig.MICRO_LENGTH);

        // convert microTime into neutron velocity based on flight path distance
        double velocity = pow(10.,7.)*config.facilityConfig.FLIGHT_DISTANCE/microTime; // in meters/sec 

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


        // Fill raw TOF before gates:

        // Apply gates:
        if (timeDiff < config.facilityConfig.MACRO_LENGTH && timeDiff > 0 // omit events outside beam-on period
                //&& procEvent.sgQ >= procEvent.lgQ
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

            TOFHistos[procEvent.targetPos-1]->Fill(microTime);
            energyHistos[procEvent.targetPos-1]->Fill(rKE);

            goodMacroNoH->Fill(procEvent.macroNo);

            switch(procEvent.targetPos)
            {
                case 1:
                    timeDiffHistoBlank->Fill(eventTimeDiff);
                    goodMacroNoHBlank->Fill(procEvent.macroNo);
                    break;
                case 2:
                    timeDiffHistoTarget1->Fill(eventTimeDiff);
                    goodMacroNoHTarget1->Fill(procEvent.macroNo);
                    break;
                case 3:
                    timeDiffHistoTarget2->Fill(eventTimeDiff);
                    goodMacroNoHTarget2->Fill(procEvent.macroNo);
                    break;
                case 4:
                    timeDiffHistoTarget3->Fill(eventTimeDiff);
                    goodMacroNoHTarget3->Fill(procEvent.macroNo);
                    break;
                case 5:
                    timeDiffHistoTarget4->Fill(eventTimeDiff);
                    goodMacroNoHTarget4->Fill(procEvent.macroNo);
                    break;
                case 6:
                    timeDiffHistoTarget5->Fill(eventTimeDiff);
                    goodMacroNoHTarget5->Fill(procEvent.macroNo);
                    break;
            }

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
              temp << "macroNo " << procEvent.macroNo << ", eventNo " << procEvent.eventNo;
              TH1I* waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,procEvent.waveform->size());

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

    for(size_t i=1; i<goodMacrosByTarget.size(); i++)
    {
        if((badMacrosByTarget[i]==0) && goodMacrosByTarget[i]==(0))
        {
            continue;
        }

        ratioBadMacrosH->SetBinContent(i+1,
                (double)(badMacrosByTarget[i])/
                (double)(badMacrosByTarget[i]+goodMacrosByTarget[i]));
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

    goodMacroNoH->Write();
    goodMacroNoHBlank->Write();
    goodMacroNoHTarget1->Write();
    goodMacroNoHTarget2->Write();
    goodMacroNoHTarget3->Write();
    goodMacroNoHTarget4->Write();
    goodMacroNoHTarget5->Write();

    timeDiffVEnergy1->Write();
    timeDiffVEnergy2->Write();
    energy1VEnergy2->Write();
    time1Vtime2->Write();
    lgQ1VlgQ2->Write();
    lgQ1VlgQ2_background->Write();

    ratioBadMacrosH->Write();

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

    outputFile->cd("/macroTime");

    for(unsigned int i=0; i<config.targetConfig.TARGET_ORDER.size(); i++)
    {
        macrosPerTarget.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
        //double badMacroRatio = ((TH1D*)gDirectory->Get("ratioBadMacrosH"))->GetBinContent(i+2);
        microsPerTarget.push_back(macrosPerTarget.back()/*(1-badMacroRatio)*/
                *(config.facilityConfig.MACRO_LENGTH/config.facilityConfig.MICRO_LENGTH));
    }

    outputFile->cd();
    outputFile->cd("/summedDet");

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
    //const int waveletBins = (config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*WAVELET_PERIOD;

    vector<double> eventsPerMicroPerBin(config.plotConfig.TOF_BINS);
    vector<double> doublePeakOddsPerBin(config.plotConfig.TOF_BINS);

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

      for(int i=0; i<config.plotConfig.TOF_BINS; i++)
      {
      for(int j=i; j<i+waveletBins/2; j++)
      {
      if(j>config.plotConfig.TOF_BINS)
      {
      doublePeakOddsPerBin[i] += eventsPerMicroPerBin[j-config.plotConfig.TOF_BINS];
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

    outputFile->Close();

    return 0;
}


