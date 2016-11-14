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

    string analysisDirectory = argv[2];
    string rawTreeFileName = analysisDirectory + "raw.root";
    string sortedFileName = analysisDirectory + "sorted.root";
    string waveformFileName = analysisDirectory + "waveform.root";
    string DPPwaveformFileName = analysisDirectory + "DPPwaveform.root";
    string histoFileName = analysisDirectory + "histos.root";
    string CSFileNameHighThresh = analysisDirectory + "cross-sections_highT.root";
    string CSFileNameLowThresh = analysisDirectory + "cross-sections_lowT.root";

    string errorFileName = analysisDirectory + "error.txt";
    string tempFileName = analysisDirectory + "temp.root";

    /*************************************************************************/
    /* Start analysis */
    /*************************************************************************/

    /*************************************************************************/
    /* Populate raw event data into a tree */
    /*************************************************************************/
    ifstream f(rawTreeFileName);
    if(!f.good())
    {
        // create a raw data tree for this subrun
        readRawData(rawDataFileName,rawTreeFileName);
    }

    else
    {
        cout << "Raw data tree already exists." << endl;
        f.close();
    }

    /*************************************************************************/
    /* "Process" raw events by assigning time and target data */
    /*************************************************************************/
    ifstream p(sortedFileName);
    if(!p.good())
    {
        // convert the raw data tree into a processed data tree
        TFile* sortedFile = new TFile(sortedFileName.c_str(),"CREATE");

        vector<TTree*> orchardProcessed; // channel-specific DPP events assigned to macropulses
        vector<TTree*> orchardProcessedW;// channel-specific waveform events assigned to macropulses

        // separate all data by channel and event type
        separateByChannel(rawTreeFileName, sortedFile, orchardProcessed, orchardProcessedW);

        // uncomment for NEVT_AGGR = 10 behavior
        /*

        // Next, extract target changer events from the input tree and add to the
        // target changer trees (DPP and waveform), assigning a macropulse to each
        // target changer event. 

        processTargetChanger(rawTreeFileName, sortedFile);

        // Last, now that the macropulse structure is assigned by the target changer
        // events, we can assign detector events to the correct macropulse.
        processDPPEvents(sortedFile, orchardRaw, orchardProcessed);
        processWaveformEvents(sortedFile, orchardRawW, orchardProcessedW);
        */

        vetoEvents(orchardProcessed[4],orchardProcessed[6], "highThreshold");
        vetoEvents(orchardProcessed[5],orchardProcessed[6], "lowThreshold");

        //cout << "Total number of ch0 waveform-mode events processed = " << numberOfCh0Waveforms << endl;
        //cout << "Total number of ch2 waveform-mode events processed = " << numberOfCh2Waveforms << endl;
        //cout << "Total number of ch4 waveform-mode events processed = " << numberOfCh4Waveforms << endl;

        sortedFile->Write();
        sortedFile->Close();
    }

    else
    {
        cout << "Found previously existing processed tree file. Skipping tree processing..." << endl;
    }

    /*************************************************************************/
    /* perform fits on DPP wavelets and waveform-mode data */

    // analyze the waveform-mode data, including peak-fitting and deadtime extraction
    //if(runWaveform)

    //waveform(sortedFileName, waveformFileName);

    histos(sortedFileName, histoFileName);

    // Apply deadtime correction to DPP-mode data
    correctForDeadtime(histoFileName, histoFileName, get<1>(channelMap[4]));
    correctForDeadtime(histoFileName, histoFileName, get<1>(channelMap[5]));
    return 0;
}
