// project-specific classes
#include "../include/analysisConstants.h"
#include "../include/plottingConstants.h"

#include "../include/raw.h"
#include "../include/separate.h"
#include "../include/histos.h"
#include "../include/plots.h"
#include "../include/waveform.h"
#include "../include/crossSection.h"
#include "../include/vetoEvents.h"

// ROOT library classes
#include "TFile.h"
#include "TTree.h"

// STL classes
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

int main(int, char* argv[])
{
    /*************************************************************************/
    /* Set up input/output filenames */
    /*************************************************************************/

    // location of digitizer-produced data file
    string rawDataFileName = argv[1];

    // location of directory where all analysis output will be stored
    string outputDirectoryName = argv[2];

    // name of experimental directory where target information is located
    string expName = argv[3];

    // run number of data - to figure out which targets were in which positions
    // during this run
    int runNumber = stoi(argv[4]);

    string analysisDirectory = argv[2];
    string rawTreeFileName = analysisDirectory + "raw.root";
    string processedFileName = analysisDirectory + "sorted.root";
    string waveformFileName = analysisDirectory + "waveform.root";
    string DPPwaveformFileName = analysisDirectory + "DPPwaveform.root";
    string histoFileName = analysisDirectory + "histos.root";
    string CSFileName = analysisDirectory + "cross-sections.root";
    string CSFileNameLowThresh = analysisDirectory + "cross-sections_low.root";

    string errorFileName = analysisDirectory + "error.txt";
    string tempFileName = analysisDirectory + "temp.root";

    /*************************************************************************/
    /* Start analysis */
    /*************************************************************************/

    /*************************************************************************/
    /* Populate raw event data into a tree */
    /*************************************************************************/
    TFile* rawTreeFile = new TFile(rawTreeFileName.c_str(),"READ");
    if(!rawTreeFile->IsOpen())
    {
        // create a raw data tree
        readRawData(rawDataFileName,rawTreeFileName);
    }

    else
    {
        cout << "Found previously existing raw data tree. Skipping import of raw data file..." << endl;
        rawTreeFile->Close();
    }

    /*************************************************************************/
    /* "Process" raw events by assigning time and target data */
    /*************************************************************************/
    TFile* processedFile = new TFile(processedFileName.c_str(),"READ");
    if(!processedFile->IsOpen())
    {
        vector<TTree*> orchardRaw;
        vector<TTree*> orchardRawW;

        // separate all data by channel and event type
        separateByChannel(rawTreeFileName, tempFileName, orchardRaw, orchardRawW);

        vector<TTree*> orchardProcessed; // channel-specific DPP events assigned to macropulses
        vector<TTree*> orchardProcessedW;// channel-specific waveform events assigned to macropulses

        processedFile = new TFile(processedFileName.c_str(),"CREATE");

        // Next, extract target changer events from the input tree and add to the
        // target changer trees (DPP and waveform), assigning a macropulse to each
        // target changer event. 

        processTargetChanger(rawTreeFileName, processedFile);

        // Last, now that the macropulse structure is assigned by the target changer
        // events, we can assign detector events to the correct macropulse.
        processDPPEvents(processedFile, orchardRaw, orchardProcessed);
        processWaveformEvents(processedFile, orchardRawW, orchardProcessedW);

        vetoEvents(orchardProcessed[2],orchardProcessed[3]);

        //cout << "Total number of ch0 waveform-mode events processed = " << numberOfCh0Waveforms << endl;
       // cout << "Total number of ch2 waveform-mode events processed = " << numberOfCh2Waveforms << endl;
       //   cout << "Total number of ch4 waveform-mode events processed = " << numberOfCh4Waveforms << endl;

        processedFile->Write();
        processedFile->Close();
    }

    else
    {
        cout << "Found previously existing processed tree file. Skipping tree processing..." << endl;
    }

    /*************************************************************************/
    /* perform fits on DPP wavelets and waveform-mode data */

    // analyze the waveform-mode data, including peak-fitting and deadtime extraction
    //if(runWaveform)

    //waveform(processedFileName, waveformFileName);

    histos(processedFileName, histoFileName);

    // Apply deadtime correction to DPP-mode data
    correctForDeadtime(histoFileName, histoFileName, dirs[2]);
    correctForDeadtime(histoFileName, histoFileName, dirs[3]);

    // calculate cross sections
    calculateCS(histoFileName,dirs[2],CSFileName,expName,runNumber);
    calculateCS(histoFileName,dirs[3],CSFileNameLowThresh,expName,runNumber);

    return 0;
}
