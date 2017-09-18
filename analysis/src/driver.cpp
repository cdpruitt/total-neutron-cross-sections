// project-specific classes
#include "../include/analysisConstants.h"
#include "../include/plottingConstants.h"

#include "../include/raw.h"
#include "../include/assignMacropulses.h"
#include "../include/histos.h"
#include "../include/plots.h"
#include "../include/waveform.h"
#include "../include/veto.h"
#include "../include/experiment.h"

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
    /* Figure out where analysis input/output should be */
    /*************************************************************************/

    // location of digitizer-produced data file
    string rawDataFileName = argv[1];

    // location of directory where all analysis output will be stored
    string analysisDirectory = argv[2];

    // location of experimental configuration files
    string experimentName = argv[3];


    /*************************************************************************/
    /* Configure how this run should be analyzed */
    /*************************************************************************/

    // read the mapping from detectors to digitizer channels
    // (see <experiment-name>/channelMap.txt).
    int runNumber = atoi(argv[4]);
    vector<string> channelMap = getChannelMap(experimentName, runNumber);

    // toggle use of the charged-particle veto paddle to reject DPP events
    string useVetoPaddleString = argv[5];
    bool useVetoPaddle = false;
    if(useVetoPaddleString=="true")
    {
        useVetoPaddle = true;
    }

    // toggle use of waveform event data during analysis
    string processWaveformEventsString = argv[6];
    bool processWaveformEvents = false;
    if(processWaveformEventsString=="true")
    {
        processWaveformEvents = true;
    }

    // toggle use of DPP wavelet data during analysis
    string processDPPWaveformEventsString = argv[7];
    bool processDPPWaveformEvents = false;
    if(processDPPWaveformEventsString=="true")
    {
        processDPPWaveformEvents = true;
    }

    /*************************************************************************/
    /* Start analysis:
     * Separate raw event data by channel and event type and store the results
     * into ROOT trees */
    /*************************************************************************/
    string rawTreeFileName = analysisDirectory + "raw.root";
    ifstream f(rawTreeFileName);
    if(!f.good())
    {
        // create a raw data tree for this subrun
        readRawData(rawDataFileName,rawTreeFileName,channelMap);
    }

    else
    {
        cout << "Raw data tree already exists." << endl;
        f.close();
    }

    /*************************************************************************/
    /* Assign macropulse number to each event */
    /*************************************************************************/
    string sortedFileName = analysisDirectory + "sorted.root";
    ifstream p(sortedFileName);
    if(!p.good())
    {
        // separate all data by channel and event type
        assignMacropulses(rawTreeFileName, sortedFileName, channelMap);

        cout << "finished assignment" << endl;
    }

    else
    {
        cout << "Sorted data tree already exists." << endl;
        p.close();
    }

    /*************************************************************************/
    /* Veto detector events using the charged-particle paddle */
    /*************************************************************************/
    string vetoedFileName = analysisDirectory + "vetoed.root";
    if(useVetoPaddle)
    {
        ifstream v(vetoedFileName);
        if(!v.good())
        {
            vetoEvents(sortedFileName, vetoedFileName, detectorNames, "veto");
        }

        else
        {
            cout << "Vetoed data tree already exists." << endl;
            v.close();
        }
    }

    /*************************************************************************/
    /* Process waveform events */
    /*************************************************************************/
    if(processWaveformEvents)
    {
        string waveformFileName = analysisDirectory + "waveform.root";
        ifstream w(waveformFileName);
        if(!w.good())
        {
            waveform(sortedFileName, waveformFileName, channelMap, "waveform");
        }

        else
        {
            cout << "Waveform events already processed." << endl;
            w.close();
        }
    }

    /*************************************************************************/
    /* Process DPP waveform events */
    /*************************************************************************/
    if(processDPPWaveformEvents)
    {
        string DPPWaveformFileName = analysisDirectory + "DPPwaveform.root";
        ifstream dppW(DPPWaveformFileName);
        if(!dppW.good())
        {
            waveform(sortedFileName, DPPWaveformFileName, channelMap, "DPP");
        }

        else
        {
            cout << "DPP waveform events already processed." << endl;
            dppW.close();
        }
    }

    /*************************************************************************/
    /* Populate each event into TOF histograms in preparation for cross section
     * calculation */
    /*************************************************************************/
    string histoFileName = analysisDirectory + "histos.root";
    ifstream h(histoFileName);
    if(!h.good())
    {
        //Uncomment to use vetoed trees
        //histos(sortedFileName, vetoedFileName, histoFileName, channelMap);
        //Uncomment to use unvetoed trees
        histos(sortedFileName, sortedFileName, histoFileName, channelMap);

        /*************************************************************************/
        /* Process events into histograms in preparation for cross section
         * calculation */
        /*************************************************************************/

        TFile* histoFile = new TFile(histoFileName.c_str(),"UPDATE");

        for(string name : detectorNames)
        {
            histoFile->cd("/");
            histoFile->cd(name.c_str());

            for(string positionName : positionNames)
            {
                string histoName = positionName + "TOFCorrected";
                TH1D* tof = (TH1D*)gDirectory->Get(histoName.c_str());
                convertTOFtoEnergy(tof, positionName + "CorrectedEnergy");
            }
        }

        histoFile->Write();
        histoFile->Close();
    }

    else
    {
        cout << "Histos already exist." << endl;
        h.close();
    }

    return 0;
}
