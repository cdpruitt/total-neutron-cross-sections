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
#include "../include/dataStructures.h"
#include "../include/branches.h"
#include "../include/plots.h"
#include "../include/fillCSHistos.h"
#include "../include/waveform.h"
#include "../include/config.h"
#include "../include/GammaCorrection.h"

using namespace std;

extern Config config;

int fillCSHistos(string inputFileName, string treeName, string outputFileName)
{
    vector<GammaCorrection> gammaCorrectionList;
    calculateGammaCorrection(inputFileName, treeName, gammaCorrectionList);

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

    TH1I* timeDiffHisto = new TH1I("time since last event","time since last event",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);

    TH1I* goodMacroNoH = new TH1I("good macros","good macros",200000,0,200000);

    TH2I *timeDiffVEnergy1 = new TH2I("time difference vs. energy of first","time difference vs. energy of first",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE,500,2,700);
    TH2I *timeDiffVEnergy2 = new TH2I("time difference vs. energy of second","time difference vs. energy of second",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE,config.plotConfig.NUMBER_ENERGY_BINS,config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND);
    TH2I *energy1VEnergy2 = new TH2I("energy of first vs. energy of second","energy of first vs. energy of second",10*config.plotConfig.NUMBER_ENERGY_BINS,config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND,10*config.plotConfig.NUMBER_ENERGY_BINS,config.plotConfig.ENERGY_LOWER_BOUND,config.plotConfig.ENERGY_UPPER_BOUND);
    TH2I *time1Vtime2 = new TH2I("time of first vs. time of second","time of first vs. time of second",config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE,config.plotConfig.TOF_RANGE,0,config.plotConfig.TOF_RANGE);
    TH2I *lgQ1VlgQ2 = new TH2I("lgQ of first vs. lgQ of second","lgQ of first vs. lgQ of second",1024,0,65536,1024,0,65536);
    TH2I *lgQ1VlgQ2_background = new TH2I("lgQ of first vs. lgQ of second_background","lgQ of first vs. lgQ of second_background",1024,0,65536,1024,0,65536);

    TH1I *microNoH = new TH1I("microNoH","microNo",360,0,360);
    TH1I *orderInMicroH = new TH1I("orderInMicroH","orderInMicro",8,0,8);
    TH1I *gammaOffsetH = new TH1I("gammaOffsetH","gammaOffsetH",1000,-3,3);
    TH2I* gammaOffsetByGammaNumber = new TH2I("gammaOffsetByGammaNumber","gammaOffsetByGammaNumber",12,0,12,1000,-5,5);
    TH1I* numberOfGammasH = new TH1I("numberOfGammasH", "number of gammas in each macropulse", 20, 0, 20);

    vector<TH1D*> TOFHistos;

    for(unsigned int i=0; i<config.targetConfig.TARGET_ORDER.size(); i++)
    {
        string TOFName = config.targetConfig.TARGET_ORDER[i] + "TOF";
        TOFHistos.push_back(new TH1D(TOFName.c_str(),TOFName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND));
    }

    // re-attach to the channel-specific tree for reading out data
    DetectorEvent event;

    tree->SetBranchAddress("cycleNumber",&event.cycleNumber);
    tree->SetBranchAddress("macroNo",&event.macroNo);
    tree->SetBranchAddress("macroTime",&event.macroTime);
    tree->SetBranchAddress("fineTime",&event.fineTime);
    tree->SetBranchAddress("eventNo",&event.eventNo);
    tree->SetBranchAddress("completeTime",&event.completeTime);
    tree->SetBranchAddress("targetPos",&event.targetPos);
    tree->SetBranchAddress("sgQ",&event.sgQ);
    tree->SetBranchAddress("lgQ",&event.lgQ);
    tree->SetBranchAddress("waveform",&event.waveform);

    double microTime;
    int microNo, prevMicroNo;
    long runningMicroNo = 0;

    // create variables for DESCRIBING MICROPULSE TYPE (i.e., whether a
    // micropulse has a gamma at the start and keep track of the influence of
    // early micropulse events on later ones)
    int orderInMicro = 0;
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

    double rKE = 0;


    gammaOffsetH->Write();
    gammaOffsetByGammaNumber->Write();

    numberOfGammasH->Write();

    unsigned int totalEntries = tree->GetEntries();

    for(long i=0; i<totalEntries; i++)
    {
        tree->GetEntry(i);

        if(event.targetPos==0)        // discard events during target-changer movement
        {
            continue;
        }

        event.completeTime -= gammaCorrectionList[event.macroNo].averageGammaOffset;

        double timeDiff = event.completeTime-event.macroTime;

        eventTimeDiff = event.completeTime-prevCompleteTime;

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
        //if (event.lgQ>500*exp(-(microTime-100)/87))
        /*if (event.lgQ<150)
          {
          continue;
          }*/


        // Fill raw TOF before gates:

        // Apply gates:
        if (timeDiff < config.facilityConfig.MACRO_LENGTH && timeDiff > 0 // omit events outside beam-on period
                //&& event.sgQ >= event.lgQ
                //&& (event.sgQ/(double)event.lgQ < 0.495
                //|| event.sgQ/(double)event.lgQ > 0.505)

                /*&& event.lgQ < 65000
                  && event.sgQ < 32500

                  && event.lgQ>50
                  && event.sgQ>50*/
           )
        {
            /*****************************************************************/
            // Fill troubleshooting plots with event variables (rKE, microtime, etc.)

            TOFHistos[event.targetPos-1]->Fill(microTime);

            triangle->Fill(microTime,event.lgQ);
            sgQlgQ->Fill(event.sgQ,event.lgQ);
            QRatio->Fill(event.sgQ/(double)event.lgQ);
            triangleRKE->Fill(rKE,event.lgQ);
            orderInMicroH->Fill(orderInMicro);
            timeDiffHisto->Fill(eventTimeDiff);

            timeDiffVEnergy1->Fill(eventTimeDiff,prevRKE);
            timeDiffVEnergy2->Fill(eventTimeDiff,rKE);
            energy1VEnergy2->Fill(prevRKE,rKE);
            time1Vtime2->Fill(prevMicroTime,microTime);
            lgQ1VlgQ2->Fill(prevlgQ,event.lgQ);

            microNoH->Fill(microNo);

        }

        /*****************************************************************/

        prevlgQ = event.lgQ;
        prevsgQ = event.sgQ;
        prevMicroTime = microTime;
        prevCompleteTime = event.completeTime;
        prevRKE = rKE;

        if(i%10000==0)
        {
            cout << "Processed " << i << " events...\r";
        }
    }

    cout << endl << "Finished populating \"" << treeName << "\" events into advanced histos." << endl;
    cout << "Total events processed = " << totalEntries << endl;

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

    for(auto histo : TOFHistos)
    {
        histo->Write();
    }

    outputFile->Close();

    return 0;
}
