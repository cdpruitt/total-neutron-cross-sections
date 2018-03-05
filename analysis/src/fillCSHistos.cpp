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

int fillCSHistos(string vetoedInputFileName, string nonVetoInputFileName, bool useVetoPaddle, string macropulseFileName, string gammaCorrectionFileName, ofstream& logFile, string outputFileName)
{
    ifstream f(outputFileName);

    if(f.good())
    {
        cout << outputFileName << " already exists; skipping gated histogramming of events." << endl;
        logFile << outputFileName << " already exists; skipping gated histogramming of events." << endl;
        return 2;
    }

    logFile << endl << "*** Filling CS histos ***" << endl;

    TFile* vetoedInputFile;
    if(useVetoPaddle)
    {
        // open vetoed input tree
        vetoedInputFile = new TFile(vetoedInputFileName.c_str(),"READ");
        if(!vetoedInputFile->IsOpen())
        {
            cerr << "Error: failed to open " << vetoedInputFileName << "  to fill histos." << endl;
            return 1;
        }
    }

    // open non-vetoed input tree
    TFile* nonVetoInputFile = new TFile(nonVetoInputFileName.c_str(),"READ");
    if(!nonVetoInputFile->IsOpen())
    {
        cerr << "Error: failed to open " << nonVetoInputFileName << "  to fill histos." << endl;
        return 1;
    }

    // open macropulse tree
    TFile* macropulseFile = new TFile(macropulseFileName.c_str(),"READ");
    if(!macropulseFile->IsOpen())
    {
        cerr << "Error: failed to open " << macropulseFileName << "  to fill histos." << endl;
        vetoedInputFile->Close();
        return 1;
    }

    TTree* macropulseTree = (TTree*)(macropulseFile->Get("macropulses"));
    if(!macropulseTree)
    {
        cerr << "Error: failed to open macropulses tree to gate histos." << endl;
        vetoedInputFile->Close();
        macropulseFile->Close();
        return 1;
    }

    MacropulseEvent me;

    macropulseTree->SetBranchAddress("cycleNumber",&me.cycleNumber);
    macropulseTree->SetBranchAddress("macroNo",&me.macroNo);
    macropulseTree->SetBranchAddress("macroTime",&me.macroTime);
    macropulseTree->SetBranchAddress("targetPos",&me.targetPos);
    macropulseTree->SetBranchAddress("numberOfEventsInMacro",&me.numberOfEventsInMacro);
    macropulseTree->SetBranchAddress("numberOfMonitorsInMacro",&me.numberOfMonitorsInMacro);
    macropulseTree->SetBranchAddress("isGoodMacro",&me.isGoodMacro);

    vector<MacropulseEvent> macropulseList;

    int numberOfEntries = macropulseTree->GetEntries();

    if(numberOfEntries==0)
    {
        cerr << "Error: no macropulses found in macropulseTree during fillCSHistos." << endl;

        if(vetoedInputFile)
        {
            vetoedInputFile->Close();
        }

        nonVetoInputFile->Close();
        macropulseFile->Close();
        return 1;
    }

    for(int i=0; i<numberOfEntries; i++)
    {
        macropulseTree->GetEntry(i);

        macropulseList.push_back(me);
    }

    // open gamma correction file
    TFile* gammaCorrectionFile = new TFile(gammaCorrectionFileName.c_str(),"READ");
    if(!gammaCorrectionFile->IsOpen())
    {
        cerr << "Error: failed to open " << gammaCorrectionFileName << "  to read gamma correction." << endl;

        if(vetoedInputFile)
        {
            vetoedInputFile->Close();
        }

        nonVetoInputFile->Close();
        macropulseFile->Close();
        return 1;
    }

    TDirectory* gammaDirectory = (TDirectory*)gammaCorrectionFile->Get(config.analysis.GAMMA_CORRECTION_TREE_NAME.c_str());
    if(!gammaDirectory)
    {
        cerr << "Error: failed to open summedDet directory in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;

        if(vetoedInputFile)
        {
            vetoedInputFile->Close();
        }

        nonVetoInputFile->Close();
        macropulseFile->Close();
        gammaCorrectionFile->Close();
        return 1;
    }

    gammaDirectory->cd();

    TH1D* gammaCorrectionHisto = (TH1D*)gammaDirectory->Get("gammaCorrection");
    if(!gammaCorrectionHisto)
    {
        cerr << "Error: failed to open gammaCorrections histo in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;

        if(vetoedInputFile)
        {
            vetoedInputFile->Close();
        }

        nonVetoInputFile->Close();
        macropulseFile->Close();
        gammaCorrectionFile->Close();

        return 1;
    }

    vector<double> gammaCorrectionList;

    int gammaCorrectionBins = gammaCorrectionHisto->GetNbinsX();
    for(int i=1; i<=gammaCorrectionBins; i++)
    {
        gammaCorrectionList.push_back(gammaCorrectionHisto->GetBinContent(i));
    }

    // define gamma times
    const double GAMMA_TIME = pow(10,7)*config.facility.FLIGHT_DISTANCE/C;
    const double GAMMA_WINDOW_WIDTH = config.time.GAMMA_WINDOW_SIZE/2;

    // create outputFile
    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    for(auto& channel : config.digitizer.CHANNEL_MAP)
    {
        if(channel.second == "-" || channel.second == "macroTime"
                || channel.second == "targetChanger")
        {
            continue;
        }

        bool isDetector = false;

        TTree* tree;

        for(auto& detName : config.cs.DETECTOR_NAMES)
        {
            if(channel.second == detName)
            {
                isDetector = true;
                break;
            }
        }

        cout << "Filling gated histograms for tree \"" << channel.second << "\"..." << endl;

        if(isDetector && useVetoPaddle)
        {
            tree = (TTree*)vetoedInputFile->Get(channel.second.c_str());
        }

        else
        {
            tree = (TTree*)nonVetoInputFile->Get(channel.second.c_str());
        }

        if(!tree)
        {
            cerr << "Error: tried to populate advanced histos, but failed to find " << channel.second << " in " << vetoedInputFileName << endl;
            vetoedInputFile->Close();
            macropulseFile->Close();
            gammaCorrectionFile->Close();

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
        tree->SetBranchAddress("waveform",&waveformPointer);

        if(isDetector)
        {
            tree->SetBranchAddress("vetoed",&event.vetoed);
        }

        TDirectory* directory = outputFile->mkdir(channel.second.c_str(),channel.second.c_str());
        directory->cd();

        vector<TH1D*> goodMacroHistos;
        for(string targetName : config.target.TARGET_ORDER)
        {
            string macroNumberName = targetName + "GoodMacros";
            goodMacroHistos.push_back(new TH1D(macroNumberName.c_str(),
                        macroNumberName.c_str(), 500000, 0, 500000));
        }

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

        TH1D* microNoH = new TH1D("microNoH","microNo",config.facility.MICROS_PER_MACRO+1
                ,0,config.facility.MICROS_PER_MACRO+1);

        vector<TH1D*> TOFHistos;
        vector<TH2D*> triangleHistos;
        vector<TH1D*> vetoTOFHistos;
        vector<TH2D*> vetoTriangleHistos;

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
                        config.plot.TOF_UPPER_BOUND,
                        pow(2,9),0,pow(2,15)));

            string vetoTOFName = "veto" + TOFName;
            vetoTOFHistos.push_back(new TH1D(vetoTOFName.c_str(),
                        vetoTOFName.c_str(),
                        config.plot.TOF_BINS,
                        config.plot.TOF_LOWER_BOUND,
                        config.plot.TOF_UPPER_BOUND));

            string vetoTriangleName = "veto" + triangleName;
            vetoTriangleHistos.push_back(new TH2D(vetoTriangleName.c_str(),
                        vetoTriangleName.c_str(),
                        config.plot.TOF_RANGE,
                        config.plot.TOF_LOWER_BOUND,
                        config.plot.TOF_UPPER_BOUND,
                        pow(2,9),0,pow(2,15)));
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

        long badMacroEvent = 0;
        long badChargeGateEvent = 0;
        long badChargeRatioEvent = 0;
        long outsideMacro = 0;

        int totalEntries = tree->GetEntries();

        int currentMacropulse = 0;
        bool endGatedHistoFill = false;

        int startOfCycleMacro = 0;
        int prevCycleNumber = 0;

        int facilityCounter = 0;

        vector<int> targetPositionPreviousMacro(7,-1);
        vector<int> targetPositionMacroCounter(7,0);

        // fill advanced histos
        for(long i=0; i<totalEntries; i++)
        {
            tree->GetEntry(i);

            if(event.cycleNumber > prevCycleNumber)
            {
                startOfCycleMacro = event.macroNo;
                prevCycleNumber = event.cycleNumber;
            }

            if(event.macroNo > macropulseList[currentMacropulse].macroNo)
            {
                currentMacropulse++;
                if(currentMacropulse>(macropulseList.size()-1))
                {
                    cout << "Reached end of gatedMacropulseList; ending fillCSHistos." << endl;
                    break;
                }

                facilityCounter++;

                if(macropulseList[currentMacropulse].macroTime - 8.4*pow(10,6) > macropulseList[currentMacropulse-1].macroTime)
                {
                    facilityCounter = 0;
                }

                continue;
            }

            // throw away events during "bad" macropulses
            if(!(macropulseList[currentMacropulse].isGoodMacro))
            {
                badMacroEvent++;
                continue;
            }

            if((int)event.macroNo > targetPositionPreviousMacro[event.targetPos])
            {
                targetPositionPreviousMacro[event.targetPos] = event.macroNo;
                targetPositionMacroCounter[event.targetPos]++;
            }

            // charge gates:
            if(isDetector)
            {
                if(event.lgQ<config.analysis.CHARGE_GATE_LOW_THRESHOLD
                        || event.lgQ>config.analysis.CHARGE_GATE_HIGH_THRESHOLD)
                {
                    badChargeGateEvent++;
                    continue;
                }

                /*if(event.sgQ/(double)event.lgQ < config.analysis.Q_RATIO_LOW_THRESHOLD
                        || event.sgQ/(double)event.lgQ > config.analysis.Q_RATIO_HIGH_THRESHOLD)
                {
                    badChargeRatioEvent++;
                    continue;
                }*/
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

            // micropulse gate:
            if(microNo < config.facility.FIRST_GOOD_MICRO
                    || microNo >= config.facility.LAST_GOOD_MICRO)
            {
                continue;
            }

            // veto gate: apply to neutron events only
            if(event.vetoed && microTime > GAMMA_TIME+GAMMA_WINDOW_WIDTH*2)
            {
                vetoTOFHistos[event.targetPos]->Fill(microTime);
                vetoTriangleHistos[event.targetPos]->Fill(microTime, event.lgQ);

                continue;
            }

            // convert micropulse time into neutron velocity based on flight path distance
            velocity = (pow(10.,7.)*config.facility.FLIGHT_DISTANCE)/microTime; // in meters/sec 

            // convert velocity to relativistic kinetic energy
            rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

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

            goodMacroHistos[event.targetPos]->Fill(event.macroNo+1);

            if(i%10000==0)
            {
                cout << "Processed " << i << " " << channel.second << " events into advanced CS histos...\r";
            }
        }

        cout << endl << "Finished populating \"" << channel.second << "\" events into CS histos." << endl;
        cout << "Total events processed = " << totalEntries << endl;

        logFile << endl << "Fraction events filtered out by good macro gate: "
            << 100*(double)badMacroEvent/totalEntries << "%." << endl;

        logFile << "Fraction events filtered out by charge gate (" << config.analysis.CHARGE_GATE_LOW_THRESHOLD
            << " < lgQ < " << config.analysis.CHARGE_GATE_HIGH_THRESHOLD << "): "
            << 100*(double)badChargeGateEvent/totalEntries << "%." << endl;

        logFile << "Fraction events filtered out by charge ratio gate (" << config.analysis.Q_RATIO_LOW_THRESHOLD
            << " < lgQ < " << config.analysis.Q_RATIO_HIGH_THRESHOLD << "): "
            << 100*(double)badChargeRatioEvent/totalEntries << "%." << endl;

        logFile << "Fraction events outside macropulse: "
            << 100*(double)outsideMacro/totalEntries << "%." << endl;

        for(auto& histo : TOFHistos)
        {
            histo->Write();
        }

        for(auto& histo : triangleHistos)
        {
            histo->Write();
        }

        for(auto& histo : vetoTOFHistos)
        {
            histo->Write();
        }

        for(auto& histo : vetoTriangleHistos)
        {
            histo->Write();
        }

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

    macropulseFile->Close();

    if(useVetoPaddle)
    {
        vetoedInputFile->Close();
    }

    nonVetoInputFile->Close();

    outputFile->Close();

    logFile << endl << "*** Finished filling CS histos ***" << endl;

    return 0;
}
