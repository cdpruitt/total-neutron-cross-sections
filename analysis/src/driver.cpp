// project-specific classes
#include "../include/analysisConstants.h"
#include "../include/raw.h"
#include "../include/separate.h"
#include "../include/histos.h"
#include "../include/plots.h"
#include "../include/waveform.h"
#include "../include/crossSection.h"

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

    // run number of data - to figure out which targets were in which positions
    // during this run
    int runNumber = stoi(argv[3]);

    string analysisDirectory = argv[2];
    string rawTreeFileName = analysisDirectory + "raw.root";
    string processedTreeFileName = analysisDirectory + "sorted.root";
    string waveformFileName = analysisDirectory + "waveform.root";
    string DPPwaveformFileName = analysisDirectory + "DPPwaveform.root";
    string histoFileName = analysisDirectory + "histos.root";
    string CSFileName = analysisDirectory + "cross-sections.root";

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
    TFile* processedTreeFile = new TFile(processedTreeFileName.c_str(),"READ");
    if(!processedTreeFile->IsOpen())
    {
        vector<TTree*> orchardRaw;
        vector<TTree*> orchardRawW;

        // separate all data by channel and event type
        separateByChannel(rawTreeFileName, tempFileName, orchardRaw, orchardRawW);

        vector<TTree*> orchardProcessed; // channel-specific DPP events assigned to macropulses
        vector<TTree*> orchardProcessedW;// channel-specific waveform events assigned to macropulses

        processedTreeFile = new TFile(processedTreeFileName.c_str(),"CREATE");

        // Next, extract target changer events from the input tree and add to the
        // target changer trees (DPP and waveform), assigning a macropulse to each
        // target changer event. 

        processTargetChanger(rawTreeFileName, processedTreeFile);

        // Last, now that the macropulse structure is assigned by the target changer
        // events, we can assign detector events to the correct macropulse.
        processDPPEvents(processedTreeFile, orchardRaw, orchardProcessed);
        processWaveformEvents(processedTreeFile, orchardRawW, orchardProcessedW);

        //cout << "Total number of ch0 waveform-mode events processed = " << numberOfCh0Waveforms << endl;
       // cout << "Total number of ch2 waveform-mode events processed = " << numberOfCh2Waveforms << endl;
       //   cout << "Total number of ch4 waveform-mode events processed = " << numberOfCh4Waveforms << endl;

        processedTreeFile->Write();
        processedTreeFile->Close();
    }

    else
    {
        cout << "Found previously existing processed tree file. Skipping tree processing..." << endl;
    }

    /*************************************************************************/
    /* perform fits on DPP wavelets and waveform-mode data */

    // analyze the waveform-mode data, including peak-fitting and deadtime extraction
    //if(runWaveform)

    waveform(processedTreeFileName, waveformFileName);

    histos(processedTreeFileName, histoFileName);

    // Apply deadtime correction to DPP-mode data
    correctForDeadtime(histoFileName, waveformFileName);

    // calculate cross sections
    calculateCS(histoFileName,CSFileName,runNumber);

    return 0;
}
