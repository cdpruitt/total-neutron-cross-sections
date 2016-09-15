// project-specific classes
#include "../include/target.h"
#include "../include/driver.h"
#include "../include/analysisConstants.h"
#include "../include/resort.h"
#include "../include/histos.h"
#include "../include/plots.h"
#include "../include/waveform.h"

// ROOT library classes
#include "TFile.h"
#include "TTree.h"

// STL classes
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
    /*************************************************************************/
    /* Parse input parameters */

    // location of digitizer-produced data file
    string inputFileName = argv[1];

    // location of directory where all output will be stored
    string outputDirectoryName = argv[2];

    // run number of data - to figure out which targets were in which positions
    // during this run
    int runNumber = stoi(argv[3]);

    // flag indicating whether ./waveform should be run
    bool runWaveform = argv[4];

    // flag indicating whether ./DPPwaveform should be run
    //bool runDPPFitting = argv[5];


    /*************************************************************************/
    /* Retrieve target data */

    vector<string> targetOrder;

    if(runNumber>5 && runNumber<=151)
    {
        targetOrder = {"blank","shortCarbon","longCarbon","Sn112","NatSn","Sn124"};
    }

    else if(runNumber==152)
    {
        targetOrder = {"blank","Sn112","NatSn","Sn124","shortCarbon","longCarbon"};
    }
        
    else if(runNumber>=153 && runNumber<=168)
    {
        targetOrder = {"blank","Sn112","NatSn","Sn124"};
    }

    else if(runNumber>=169 && runNumber<=180)
    {
        targetOrder = {"blank","Sn112","NatSn","Sn124","shortCarbon"};
    }

    else
    {
        cout << "Error - bad run number given to ./driver" << endl;
    }

    vector<Target*> targets;
    for(string s : targetOrder)
    {
        targets.push_back(
                new Target(TARGET_DATA_FILE_PATH + s + TARGET_DATA_FILE_EXTENSION));
    }


    /*************************************************************************/
    /* Analyze data */

    // Set filenames for output files
    string analysisDirectory = argv[2];
    string rawFileName = analysisDirectory + "raw.root";
    string sortedFileName = analysisDirectory + "sorted.root";
    string waveformFileName = analysisDirectory + "waveform.root";
    string DPPwaveformFileName = analysisDirectory + "DPPwaveform.root";
    string histoFileName = analysisDirectory + "histos.root";
    string CSFileName = analysisDirectory + "cross-sections.root";

    string errorFileName = analysisDirectory + "error.txt";
    string tempFileName = analysisDirectory + "temp.root";

    /*************************************************************************/
    /* extract raw data from binary file */

    // Check to see if raw trees already exist
    TFile *rawFile;
    rawFile = new TFile(rawFileName.c_str(),"UPDATE");
    if(!rawFile->Get("tree"))
    {
        // we need to (re)sort the raw data into trees
        extractRawData(inputFileName,rawFileName);
    }

    else
    {
        cout << "Found previously existing raw sort " << rawFileName << ". Skipping raw sort." << endl;
    }

    rawFile->Close();

    /*************************************************************************/
    /* separate data by channel and assign to macropulse */

    // Check to see if sorted trees already exist
    TFile* sortedFile = new TFile(sortedFileName.c_str(),"READ");
    if(!sortedFile->IsOpen())
    {
        vector<TTree*> orchardRaw;       // channel-specific DPP events NOT assigned to macropulses
        vector<TTree*> orchardRawW;      // channel-specific waveform events NOT assigned to macropulses

        // separate all data by channel
        separateByChannel(rawFileName, tempFileName, orchardRaw, orchardRawW);

        // Create an error log where sorting errors can be recorded for review
        ofstream errorFile;
        errorFile.open(errorFileName);
        errorFile.precision(13);

        vector<TTree*> orchardProcessed; // channel-specific DPP events assigned to macropulses
        vector<TTree*> orchardProcessedW;// channel-specific waveform events assigned to macropulses

        sortedFile = new TFile(sortedFileName.c_str(),"CREATE");

        // Next, extract target changer events from the input tree and add to the
        // target changer trees (DPP and waveform), assigning a macropulse to each
        // target changer event. 

        processTargetChanger(rawFileName, sortedFile, errorFile);

        // Last, now that the macropulse structure is assigned by the target changer
        // events, we can assign detector events to the correct macropulse.
        processDPPEvents(sortedFile, orchardRaw, orchardProcessed, errorFile);
        processWaveformEvents(sortedFile, orchardRawW, orchardProcessedW, errorFile);

        //cout << "Total number of ch0 waveform-mode events processed = " << numberOfCh0Waveforms << endl;
       // cout << "Total number of ch2 waveform-mode events processed = " << numberOfCh2Waveforms << endl;
       //   cout << "Total number of ch4 waveform-mode events processed = " << numberOfCh4Waveforms << endl;

        sortedFile->Write();
        sortedFile->Close();
    }

    else
    {
        cout << "Found previously existing file " << sortedFileName << ". Skipping resort..." << endl;
    }

    /*************************************************************************/
    /* perform fits on DPP wavelets and waveform-mode data */

    TFile* histoFile = new TFile(histoFileName.c_str(),"READ");
    TFile* CSFile = new TFile(CSFileName.c_str(),"READ");
    TFile* waveformFile = new TFile(waveformFileName.c_str(),"READ");

    vector<Plots*> plots;
    for(string s : targetOrder)
    {
        plots.push_back(new Plots(s, histoFile, waveformFile));
    }

    // analyze the waveform-mode data, including peak-fitting and deadtime extraction
    if(runWaveform)
    {
        waveform(sortedFileName, waveformFile, targets, plots);
    }

    // perform peak-fitting on the DPP-mode wavelets
    //if(DPPwaveform)
    //{
    //    DPPwaveform(sortedFileName,DPPWaveformFileName);
    //}

    // fill histograms with DPP data
    histos(sortedFileName, histoFile, CSFileName, plots);

    // calculate cross sections
    //calculateCS(sortedFileName,CSFileName);

    return 0;
}
