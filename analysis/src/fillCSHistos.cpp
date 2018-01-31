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

int fillCSHistos(string inputFileName, string macropulseFileName, string gammaCorrectionFileName, ofstream& logFile, string outputFileName)
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

    // open macropulse tree
    TFile* macropulseFile = new TFile(macropulseFileName.c_str(),"READ");
    if(!macropulseFile->IsOpen())
    {
        cerr << "Error: failed to open " << macropulseFileName << "  to fill histos." << endl;
        inputFile->Close();
        return 1;
    }

    TTree* macropulseTree = (TTree*)(macropulseFile->Get("macropulses"));
    if(!macropulseTree)
    {
        cerr << "Error: failed to open macropulses tree to gate histos." << endl;
        inputFile->Close();
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
        inputFile->Close();
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
        inputFile->Close();
        macropulseFile->Close();
        return 1;
    }

    TDirectory* gammaDirectory = (TDirectory*)gammaCorrectionFile->Get(config.analysis.GAMMA_CORRECTION_TREE_NAME.c_str());
    if(!gammaDirectory)
    {
        cerr << "Error: failed to open summedDet directory in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;
        inputFile->Close();
        macropulseFile->Close();
        gammaCorrectionFile->Close();
        return 1;
    }

    gammaDirectory->cd();

    TH1D* gammaCorrectionHisto = (TH1D*)gammaDirectory->Get("gammaCorrection");
    if(!gammaCorrectionHisto)
    {
        cerr << "Error: failed to open gammaCorrections histo in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;
        inputFile->Close();
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

    // create outputFile
    TFile* outputFile = new TFile(outputFileName.c_str(),"UPDATE");

    for(auto& channel : config.digitizer.CHANNEL_MAP)
    {
        if(channel.second == "-" || channel.second == "macroTime"
                || channel.second == "targetChanger")
        {
            continue;
        }

        cout << "Filling gated histograms for tree \"" << channel.second << "\"..." << endl;

        TTree* tree = (TTree*)inputFile->Get(channel.second.c_str());
        if(!tree)
        {
            cerr << "Error: tried to populate advanced histos, but failed to find " << channel.second << " in " << inputFileName << endl;
            inputFile->Close();
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
        tree->SetBranchAddress("vetoed",&event.vetoed);
        tree->SetBranchAddress("waveform",&waveformPointer);

        TDirectory* directory = outputFile->mkdir(channel.second.c_str(),channel.second.c_str());
        directory->cd();

        vector<TH1D*> goodMacroHistos;
        for(string targetName : config.target.TARGET_ORDER)
        {
            string macroNumberName = targetName + "GoodMacros";
            goodMacroHistos.push_back(new TH1D(macroNumberName.c_str(),
                        macroNumberName.c_str(), 500000, 0, 500000));
        }

        vector<TH1D*> macrosSinceCycleHistos;
        for(string targetName : config.target.TARGET_ORDER)
        {
            string macrosSinceCycleName = targetName + "MacrosSinceCycle";
            macrosSinceCycleHistos.push_back(new TH1D(macrosSinceCycleName.c_str(),
                        macrosSinceCycleName.c_str(), 1000, 0, 1000));
        }

        vector<TH1D*> macrosSinceFacilityGapHistos;
        for(string targetName : config.target.TARGET_ORDER)
        {
            string macrosSinceFacilityGapName = targetName + "MacrosSinceFacilityGap";
            macrosSinceFacilityGapHistos.push_back(new TH1D(macrosSinceFacilityGapName.c_str(),
                        macrosSinceFacilityGapName.c_str(), 100, 0, 100));
        }

        vector<TH1D*> postCycleHistos;
        for(string targetName : config.target.TARGET_ORDER)
        {
            string postCycleHistoName = targetName + "postCycle";
            postCycleHistos.push_back(new TH1D(postCycleHistoName.c_str(),
                        postCycleHistoName.c_str(), 500000, 0, 500000));
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
            /*if(event.lgQ<config.analysis.CHARGE_GATE_LOW_THRESHOLD
                    || event.lgQ>config.analysis.CHARGE_GATE_HIGH_THRESHOLD)
            {
                badChargeGateEvent++;
                continue;
            }

            if(event.sgQ/(double)event.lgQ < config.analysis.Q_RATIO_LOW_THRESHOLD
                    || event.sgQ/(double)event.lgQ > config.analysis.Q_RATIO_HIGH_THRESHOLD)
            {
                badChargeRatioEvent++;
                continue;
            }*/

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

            // veto gate
            if(event.vetoed)
            {
                vetoedTOF->Fill(microTime);
                vetoedTriangle->Fill(microTime, event.lgQ);

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

            macrosSinceCycleHistos[event.targetPos]
                ->Fill(event.macroNo-startOfCycleMacro);

            macrosSinceFacilityGapHistos[event.targetPos]
                ->Fill(facilityCounter);

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

        // calculate width of gamma peak after gamma correction has been applied
        const double GAMMA_TIME = pow(10,7)*config.facility.FLIGHT_DISTANCE/C;
        const double GAMMA_WINDOW_WIDTH = config.time.GAMMA_WINDOW_SIZE/2;

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

        for(auto& histo : macrosSinceCycleHistos)
        {
            histo->Write();
        }

        for(auto& histo : macrosSinceFacilityGapHistos)
        {
            histo->Write();
        }

        // calculate time autocorrelation of number of events in good macro histos
        /*vector<TH1D*> macroAutocorrelationHistos;
        for(string targetName : config.target.TARGET_ORDER)
        {
            string macroAutocorrelationName = targetName + "MacroAutocorrelation";
            macroAutocorrelationHistos.push_back(new TH1D(macroAutocorrelationName.c_str(),
                        macroAutocorrelationName.c_str(), 1000, 0, 1000));
        }

        // calculate average number of events per macro for each target
        vector<double> goodMacrosAverageByTarget(goodMacroHistos.size());
        vector<double> goodMacrosAverageVarianceByTarget(goodMacroHistos.size());

        for(int i=0; i<goodMacroHistos.size(); i++)
        {
            cout << "calculating macro average for " << config.target.TARGET_ORDER[i]
                << "..." << endl;

            TH1D* histo = goodMacroHistos[i];
            int numberOfBins = histo->GetNbinsX();
            int numberOfMacros = 0;

            for(int j=1; j<numberOfBins-1; j++)
            {
                double binContent = histo->GetBinContent(j);
                if(binContent>0)
                {
                    goodMacrosAverageByTarget[i] += binContent;
                    numberOfMacros++;
                }
            }

            goodMacrosAverageByTarget[i] /= numberOfMacros;

            for(int j=1; j<numberOfBins-1; j++)
            {
                double binContent = histo->GetBinContent(j);
                if(binContent>0)
                {
                    goodMacrosAverageVarianceByTarget[i] += pow(binContent-goodMacrosAverageByTarget[i],2);
                }
            }

            goodMacrosAverageVarianceByTarget[i] /= numberOfMacros;

            cout << "macro average for " << config.target.TARGET_ORDER[i]
                << " = " << goodMacrosAverageByTarget[i] << endl;

            cout << "macro average variance for " << config.target.TARGET_ORDER[i]
                << " = " << goodMacrosAverageVarianceByTarget[i] << endl;
        }

        for(int i=0; i<goodMacroHistos.size(); i++)
        {
            TH1D* macroHisto = goodMacroHistos[i];
            double average = goodMacrosAverageByTarget[i];
            double averageVariance = goodMacrosAverageVarianceByTarget[i];

            int numberOfBins = macroHisto->GetNbinsX();

            for(int delay = 1; delay<1000; delay+=1)
            {
                double correlation = 0;

                int macroCounter = 0;

                for(int j=0; j+delay<numberOfBins; j++)
                {
                    double binContent = macroHisto->GetBinContent(j);
                    double delayedBinContent = macroHisto->GetBinContent(j+delay);

                    if(binContent<=0 || delayedBinContent<=0)
                    {
                        continue;
                    }

                    correlation +=
                        (binContent-average)*
                        (delayedBinContent-average);

                    macroCounter++;

                    if(j%1000==0)
                    {
                        cout << "Calculated autocorrelation (delay = " << delay << ") through " << j << " good macropulses...\r";
                    }
                }

                correlation /= (macroCounter-1)*averageVariance;

                macroAutocorrelationHistos[i]->SetBinContent(delay, correlation);
            }
        }

        for(auto& histo : macroAutocorrelationHistos)
        {
            histo->Write();
        }*/
    }

    macropulseFile->Close();
    outputFile->Close();
    inputFile->Close();

    logFile << endl << "*** Finished filling CS histos ***" << endl;

    return 0;
}
