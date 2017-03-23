// project-specific classes
#include "../include/analysisConstants.h"
#include "../include/plottingConstants.h"

#include "../include/raw.h"
#include "../include/separate.h"
#include "../include/histos.h"
#include "../include/plots.h"
#include "../include/waveform.h"
#include "../include/crossSection.h"
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
    /* Set up input/output filenames */
    /*************************************************************************/

    // location of digitizer-produced data file
    string rawDataFileName = argv[1];

    // location of directory where all analysis output will be stored
    string analysisDirectory = argv[2];
 
    // location of experimental configuration files
    string experimentName = argv[3];

    int runNumber = atoi(argv[4]);

    // create names of all output files
    string rawTreeFileName = analysisDirectory + "raw.root";
    string sortedFileName = analysisDirectory + "sorted.root";
    string vetoedFileName = analysisDirectory + "vetoed.root";
    string histoFileName = analysisDirectory + "histos.root";
    string waveformFileName = analysisDirectory + "waveform.root";
    string DPPWaveformFileName = analysisDirectory + "DPPwaveform.root";
    string errorFileName = analysisDirectory + "error.txt";

    vector<string> channelMap = getChannelMap(experimentName, runNumber);

    /*************************************************************************/
    /* Start analysis:
     * Separate raw event data by channel and event type into ROOT trees */
    /*************************************************************************/
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
    /* Process raw events: assign time and target data */
    /*************************************************************************/
    ifstream p(sortedFileName);
    if(!p.good())
    {
        // separate all data by channel and event type
        assignMacropulses(rawTreeFileName, sortedFileName, channelMap);
    }

    else
    {
        cout << "Sorted data tree already exists." << endl;
        p.close();
    }

    /*************************************************************************/
    /* Eliminate detector events using the veto paddle */
    /*************************************************************************/
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

    /*************************************************************************/
    /* Process waveform data */
    /*************************************************************************/

    ifstream w(waveformFileName);
    if(!w.good())
    {
        waveform(sortedFileName, waveformFileName, channelMap, "waveform");
    }

    else
    {
        cout << "Waveform data already processed." << endl;
        w.close();
    }

    /*************************************************************************/
    /* Process DPP waveform data */
    /*************************************************************************/

    ifstream dppW(DPPWaveformFileName);

    if(!dppW.good())
    {
        waveform(sortedFileName, DPPWaveformFileName, channelMap, "DPP");
    }

    else
    {
        cout << "DPP waveform data already processed." << endl;
        dppW.close();
    }

    /*************************************************************************/
    /* Process events into histograms in preparation for cross section
     * calculation */
    /*************************************************************************/
    ifstream h(histoFileName);
    if(!h.good())
    {
        //Uncomment to use vetoed trees
        histos(sortedFileName, vetoedFileName, histoFileName, channelMap);
        //Uncomment to use unvetoed trees
        //histos(sortedFileName, sortedFileName, histoFileName, channelMap);
    }

    else
    {
        cout << "Histos already exist." << endl;
        h.close();
    }

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
            string histoName = positionName + "TOF";
            TH1I* tof = (TH1I*)gDirectory->Get(histoName.c_str());
            convertTOFtoEn(tof, positionName + "CorrectedEnergy");
        }
    }
    histoFile->Write();
    histoFile->Close();

    return 0;
}
