#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectoryFile.h"
#include "TF1.h"

#include "../include/physicalConstants.h"
#include "../include/dataStructures.h"
#include "../include/branches.h"
#include "../include/plots.h"
#include "../include/fillCSHistos.h"
#include "../include/waveform.h"
#include "../include/config.h"
#include "../include/GammaCorrection.h"

using namespace std;

extern Config config;

const double Q_LOW_THRESHOLD = 0;
const double Q_HIGH_THRESHOLD = 65550;

int identifyGoodMacros(TTree* eventTree, DetectorEvent& event, vector<MacropulseEvent>& macropulseList, ofstream& logFile)
{
    cout << "Starting identification of good macropulses..." << endl;

    // tally the number of events in each macropulse
    eventTree->GetEntry(0);
    macropulseList.push_back(MacropulseEvent(event.macroNo,0,event.targetPos));

    unsigned int totalEntries = eventTree->GetEntries();
    for(long i=1; i<totalEntries; i++)
    {
        eventTree->GetEntry(i);

        if(macropulseList.back().macroNo<event.macroNo)
        {
            macropulseList.push_back(MacropulseEvent(event.macroNo, 0, event.targetPos));
        }

        macropulseList.back().numberOfEventsInMacro++;

        if(i%10000==0)
        {
            cout << "Found " << macropulseList.size() << " macropulses during good macropulse identification...\r";
            fflush(stdout);
        }
    }

    cout << endl;

    // calculate the average number of events in each macropulse, by target
    vector<double> averageEventsPerMacropulseByTarget(config.target.TARGET_ORDER.size(),0);
    vector<unsigned int> numberOfMacropulsesByTarget(config.target.TARGET_ORDER.size(),0);

    for(auto& macropulse : macropulseList)
    {
        averageEventsPerMacropulseByTarget[macropulse.targetPos]
            += macropulse.numberOfEventsInMacro;
        numberOfMacropulsesByTarget[macropulse.targetPos]++;
    }

    for(int i=0; i<averageEventsPerMacropulseByTarget.size(); i++)
    {
        averageEventsPerMacropulseByTarget[i] /= numberOfMacropulsesByTarget[i];
        logFile << "Target \"" << config.target.TARGET_ORDER[i]
            << "\" had, on average, " << averageEventsPerMacropulseByTarget[i]
            << " events per macropulse." << endl;
    }

    cout << "Finished calculating average events per macropulse."
        << endl;

    // using the average just calculate, identify macros with too few/too many events
    // per macro (indicating readout problem)
    for(auto& macropulse : macropulseList)
    {
        if(abs(averageEventsPerMacropulseByTarget[macropulse.targetPos]-macropulse.numberOfEventsInMacro)
                <3*sqrt(averageEventsPerMacropulseByTarget[macropulse.targetPos]))
        {
            // good macro
            macropulse.isGoodMacro = true;
        }
    }

    return 0;
}

int fillCSHistos(string inputFileName, string gammaCorrectionFileName, ofstream& logFile, string outputFileName)
{
    ifstream f(outputFileName);

    if(f.good())
    {
        cout << outputFileName << " already exists; skipping gated histogramming of events." << endl;
        logFile << outputFileName << " already exists; skipping gated histogramming of events." << endl;
        return 2;
    }

    logFile << endl << "*** Filling CS histos ***" << endl;

    // open input tree
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    // open gamma correction file
    TFile* gammaCorrectionFile = new TFile(gammaCorrectionFileName.c_str(),"READ");
    if(!gammaCorrectionFile->IsOpen())
    {
        cerr << "Error: failed to open " << gammaCorrectionFileName << "  to read gamma correction." << endl;
        return 1;
    }

    TDirectory* gammaDirectory = (TDirectory*)gammaCorrectionFile->Get("summedDet");
    if(!gammaDirectory)
    {
        cerr << "Error: failed to open summedDet directory in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;
        return 1;
    }

    gammaDirectory->cd();

    TH1D* gammaCorrectionHisto = (TH1D*)gammaDirectory->Get("gammaCorrection");
    if(!gammaCorrectionHisto)
    {
        cerr << "Error: failed to open gammaCorrections histo in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;
        return 1;
    }

    vector<double> gammaCorrectionList;

    unsigned int gammaCorrectionBins = gammaCorrectionHisto->GetNbinsX();
    for(unsigned int i=1; i<=gammaCorrectionBins; i++)
    {
        gammaCorrectionList.push_back(gammaCorrectionHisto->GetBinContent(i));
    }

    // create outputFile
    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    for(auto& channelName : config.cs.DETECTOR_NAMES)
    {
        cout << "Filling gated histograms for tree \"" << channelName << "\"..." << endl;

        TTree* tree = (TTree*)inputFile->Get(channelName.c_str());
        if(!tree)
        {
            cerr << "Error: tried to populate advanced histos, but failed to find " << channelName << " in " << inputFileName << endl;
            return 1;
        }

        // connect input tree to event data buffer
        DetectorEvent event;
        vector<int>* waveformPointer = 0;

        tree->SetBranchAddress("cycleNumber",&event.cycleNumber);
        tree->SetBranchAddress("macroNo",&event.macroNo);
        tree->SetBranchAddress("macroTime",&event.macroTime);
        tree->SetBranchAddress("fineTime",&event.fineTime);
        tree->SetBranchAddress("eventNo",&event.eventNo);
        tree->SetBranchAddress("completeTime",&event.completeTime);
        tree->SetBranchAddress("targetPos",&event.targetPos);
        tree->SetBranchAddress("sgQ",&event.sgQ);
        tree->SetBranchAddress("lgQ",&event.lgQ);
        tree->SetBranchAddress("vetoed",&event.vetoed);
        tree->SetBranchAddress("waveform",&waveformPointer);

        vector<MacropulseEvent> macropulseList;
        identifyGoodMacros(tree, event, macropulseList, logFile);

        TDirectory* directory = outputFile->mkdir(channelName.c_str(),channelName.c_str());
        directory->cd();

        vector<TH1D*> goodMacroHistos;
        for(string targetName : config.target.TARGET_ORDER)
        {
            string macroNumberName = targetName + "GoodMacros";
            goodMacroHistos.push_back(new TH1D(macroNumberName.c_str(),
                        macroNumberName.c_str(), 200000, 0, 200000));
        }

        unsigned long goodMacroNumber = 0;
        for(auto& macropulse : macropulseList)
        {
            if(macropulse.isGoodMacro)
            {
                goodMacroHistos[macropulse.targetPos]
                    ->SetBinContent(macropulse.macroNo+1, macropulse.numberOfEventsInMacro);
                goodMacroNumber++;
            }
        }

        logFile << "Fraction of all macropulses that were \"good\":"
            << 100*(double)goodMacroNumber/macropulseList.size() << "%." << endl;

        // create other diagnostic histograms used to examine run data
        TH1D* timeDiffHisto = new TH1D("time since last event","time since last event",
                config.plot.TOF_RANGE,0,config.plot.TOF_RANGE);
        TH2D* timeDiffVEnergy1 = new TH2D("time difference vs. energy of first",
                "time difference vs. energy of first",config.plot.TOF_RANGE,
                0,config.plot.TOF_RANGE,10*config.plot.NUMBER_ENERGY_BINS,2,700);

        TH2D* time1Vtime2 = new TH2D("time of first vs. time of second",
                "time of first vs. time of second",config.plot.TOF_RANGE,0,
                config.plot.TOF_RANGE,config.plot.TOF_RANGE,0,
                config.plot.TOF_RANGE);

        TH2D* energy1VEnergy2 = new TH2D("energy of first vs. energy of second",
                "energy of first vs. energy of second",
                10*config.plot.NUMBER_ENERGY_BINS, floor(config.plot.ENERGY_LOWER_BOUND), ceil(config.plot.ENERGY_UPPER_BOUND),
                10*config.plot.NUMBER_ENERGY_BINS, floor(config.plot.ENERGY_LOWER_BOUND), ceil(config.plot.ENERGY_UPPER_BOUND));

        TH1D* microNoH = new TH1D("microNoH","microNo",ceil(1.1*config.facility.MICROS_PER_MACRO)
                ,0,ceil(1.1*config.facility.MICROS_PER_MACRO));

        TH1D* vetoedTOF = new TH1D("vetoedTOF",
                "vetoedTOF",
                config.plot.TOF_BINS,
                config.plot.TOF_LOWER_BOUND,
                config.plot.TOF_UPPER_BOUND);

        TH2D* vetoedTriangle = new TH2D("vetoedTriangle",
                "vetoedTriangle",
                config.plot.TOF_RANGE,
                config.plot.TOF_LOWER_BOUND,
                config.plot.TOF_UPPER_BOUND/5,
                4096,0,65536);

        vector<TH1D*> TOFHistos;
        vector<TH2D*> triangleHistos;

        for(string targetName : config.target.TARGET_ORDER)
        {
            string TOFName = targetName + "TOF";
            TOFHistos.push_back(new TH1D(TOFName.c_str(),
                        TOFName.c_str(),
                        config.plot.TOF_BINS,
                        config.plot.TOF_LOWER_BOUND,
                        config.plot.TOF_UPPER_BOUND));

            string triangleName = targetName + "Triangle";
            triangleHistos.push_back(new TH2D(triangleName.c_str(),
                        triangleName.c_str(),
                        config.plot.TOF_RANGE,
                        config.plot.TOF_LOWER_BOUND,
                        config.plot.TOF_UPPER_BOUND/5,
                        4096,0,65536));
        }

        double prevCompleteTime = 0;
        double prevlgQ = 0;

        double microTime;
        double prevMicroTime = 0;
        int microNo;

        double timeDiff;
        double eventTimeDiff = 0;
        double velocity;
        double rKE;
        double prevRKE = 0;

        double prevAverageTime = 0;

        const double MACRO_LENGTH = config.facility.MICROS_PER_MACRO*config.facility.MICRO_LENGTH;

        unsigned long badMacroEvent = 0;
        unsigned long badChargeGateEvent = 0;
        unsigned long outsideMacro = 0;

        unsigned int totalEntries = tree->GetEntries();

        // fill advanced histos
        for(long i=0; i<totalEntries; i++)
        {
            tree->GetEntry(i);

            // throw away events during "bad" macropulses
            if(!macropulseList[event.macroNo].isGoodMacro)
            {
                badMacroEvent++;
                continue;
            }

            // charge gates:
            if(event.lgQ<Q_LOW_THRESHOLD || event.lgQ>Q_HIGH_THRESHOLD)
            {
                badChargeGateEvent++;
                continue;
            }

            /*****************************************************************/
            // Calculate event properties

            // find which micropulse the event is in and the time since the start of
            // the micropulse (the TOF)
            timeDiff = event.completeTime-event.macroTime;

            // correct times using average gamma time
            timeDiff -= gammaCorrectionList[event.macroNo];

            // timing gate
            if(timeDiff > MACRO_LENGTH)
            {
                outsideMacro++;
                continue;
            }

            eventTimeDiff = event.completeTime-prevCompleteTime;
            microNo = floor(timeDiff/config.facility.MICRO_LENGTH);
            microTime = fmod(timeDiff,config.facility.MICRO_LENGTH);

            // convert micropulse time into neutron velocity based on flight path distance
            velocity = (pow(10.,7.)*config.facility.FLIGHT_DISTANCE)/microTime; // in meters/sec 

            // convert velocity to relativistic kinetic energy
            rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            // veto gate
            if(event.vetoed)
            {
                vetoedTOF->Fill(microTime);
                vetoedTriangle->Fill(microTime, event.lgQ);

                continue;
            }

            TOFHistos[event.targetPos]->Fill(microTime);
            triangleHistos[event.targetPos]->Fill(microTime, event.lgQ);

            // fill detector histograms with event data
            timeDiffHisto->Fill(eventTimeDiff);
            timeDiffVEnergy1->Fill(eventTimeDiff,prevRKE);
            time1Vtime2->Fill(prevMicroTime,microTime);
            energy1VEnergy2->Fill(prevRKE,rKE);
            microNoH->Fill(microNo);

            prevlgQ = event.lgQ;
            prevMicroTime = microTime;
            prevCompleteTime = event.completeTime;
            prevRKE = rKE;

            if(i%10000==0)
            {
                cout << "Processed " << i << " " << channelName << " events into advanced CS histos...\r";
            }
        }

        cout << endl << "Finished populating \"" << channelName << "\" events into CS histos." << endl;
        cout << "Total events processed = " << totalEntries << endl;

        logFile << endl << "Fraction events filtered out by good macro gate: "
            << 100*(double)badMacroEvent/totalEntries << "%." << endl;

        logFile << "Fraction events filtered out by charge gate (" << Q_LOW_THRESHOLD
            << " < lgQ < " << Q_HIGH_THRESHOLD << "): "
            << 100*(double)badChargeGateEvent/totalEntries << "%." << endl;

        logFile << "Fraction events outside macropulse: "
            << 100*(double)outsideMacro/totalEntries << "%." << endl;

        // calculate width of gamma peak after gamma correction has been applied
        const double GAMMA_TIME = pow(10,7)*config.facility.FLIGHT_DISTANCE/C;
        const double GAMMA_WINDOW_WIDTH = config.timeOffsets.GAMMA_WINDOW_SIZE/2;

        TF1* gammaPeakFit = new TF1("gammaPeakFit","gaus",
                GAMMA_TIME-GAMMA_WINDOW_WIDTH, GAMMA_TIME+GAMMA_WINDOW_WIDTH);
        TOFHistos[1]->Fit("gammaPeakFit","Q0", "",
                GAMMA_TIME-GAMMA_WINDOW_WIDTH, GAMMA_TIME+GAMMA_WINDOW_WIDTH);
        logFile << endl << "Blank target gamma peak FWHM = "
            << 2.355*gammaPeakFit->GetParameter(2)
            << " ns (fit range: " << GAMMA_TIME-GAMMA_WINDOW_WIDTH
            << " - " << GAMMA_TIME+GAMMA_WINDOW_WIDTH << " ns)." << endl;

        for(auto& histo : TOFHistos)
        {
            histo->Write();
        }

        for(auto& histo : triangleHistos)
        {
            histo->Write();
        }

        vetoedTOF->Write();
        vetoedTriangle->Write();

        timeDiffHisto->Write();
        timeDiffVEnergy1->Write();
        time1Vtime2->Write();
        energy1VEnergy2->Write();
        microNoH->Write();

        for(auto& histo : goodMacroHistos)
        {
            histo->Write();
        }
    }

    outputFile->Close();
    inputFile->Close();

    logFile << endl << "*** Finished filling CS histos ***" << endl;

    return 0;
}
