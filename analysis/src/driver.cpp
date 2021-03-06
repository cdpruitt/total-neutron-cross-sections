// project-specific classes
#include "../include/raw.h"
#include "../include/identifyMacropulses.h"
#include "../include/assignEventsToMacropulses.h"
#include "../include/fillBasicHistos.h"
#include "../include/fillCSHistos.h"
#include "../include/correctForDeadtime.h"
#include "../include/produceEnergyHistos.h"
#include "../include/plots.h"
#include "../include/waveform.h"
#include "../include/veto.h"
#include "../include/experiment.h"
#include "../include/config.h"
#include "../include/GammaCorrection.h"
#include "../include/identifyGoodMacros.h"

// ROOT library classes
#include "TFile.h"
#include "TTree.h"

// STL classes
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

Config config;

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

    // toggle use of the charged-particle veto paddle to reject DPP events
    string useVetoPaddleString = argv[5];
    bool useVetoPaddle = false;
    if(useVetoPaddleString=="true")
    {
        useVetoPaddle = true;
    }

    config = Config(experimentName, runNumber);

    string logFileName = analysisDirectory + "log.txt";
    ofstream log(logFileName);

    /*************************************************************************/
    /* Start analysis:
     * Separate raw event data by channel and event type and store the results
     * into ROOT trees */
    /*************************************************************************/
    cout << endl << "Start processing event data into raw data tree..." << endl;

    string rawTreeFileName = analysisDirectory + config.analysis.RAW_TREE_FILE_NAME;
    if(readRawData(rawDataFileName, rawTreeFileName, log))
    {
        return 1;
    }

    /*************************************************************************/
    /* Assign each event to its macropulse */
    /*************************************************************************/
    cout << endl << "Start macropulse identification..." << endl;

    string sortedFileName = analysisDirectory + config.analysis.MACROPULSE_ASSIGNED_FILE_NAME;
    string macropulseFileName = analysisDirectory + config.analysis.MACROPULSES_FILE_NAME;

    vector<MacropulseEvent> macropulseList;

    switch(identifyMacropulses(rawTreeFileName, sortedFileName, log, macropulseList))
    {
        case 0:
            // recreate sorted.root file

            assignEventsToMacropulses(
                    rawTreeFileName,
                    sortedFileName,
                    log,
                    macropulseList
                    );

            /******************************************************************/
            /* Identify "good" macropulses */
            /******************************************************************/
            identifyGoodMacros(macropulseFileName, macropulseList, log);

            break;

        case 1:
            // error state - end analysis
            return 1;

        case 2:
            // sorted.root already exists; skip to next analysis step
            break;
    }

    
    /*************************************************************************/
    /* Veto detector events using the charged-particle paddle */
    /*************************************************************************/
    string vetoedFileName = analysisDirectory + config.analysis.PASSED_VETO_FILE_NAME;
    if(useVetoPaddle)
    {
        cout << endl << "\"Veto Events\" flag enabled; start processing detector events through veto..." << endl;
        vetoEvents(sortedFileName, vetoedFileName, log, "veto");
    }

    /******************************************************************/
    /* Populate events into basic histograms */
    /******************************************************************/
    string histoFileName = analysisDirectory + config.analysis.HISTOGRAM_FILE_NAME;
    fillBasicHistos(sortedFileName, log, histoFileName);

    /*****************************************************/
    /* Use raw TOF histos to create deadtime correction  */
    /*****************************************************/
    string deadtimeFileName = analysisDirectory + "deadtime.root";
    generateDeadtimeCorrection(histoFileName, log, deadtimeFileName);

    /*****************************************************/
    /* Calculate macropulse time correction using gammas */
    /*****************************************************/
    string gammaCorrectionFileName = analysisDirectory + "gammaCorrection.root";

    calculateGammaCorrection(
            sortedFileName,
            log,
            config.analysis.GAMMA_CORRECTION_TREE_NAME,
            gammaCorrectionFileName);

    /******************************************************************/
    /* Populate events into gated histograms, using time correction   */
    /******************************************************************/
    string gatedHistoFileName = analysisDirectory + "gatedHistos.root";
    fillCSHistos(vetoedFileName, sortedFileName, useVetoPaddle, macropulseFileName, gammaCorrectionFileName, log, gatedHistoFileName);

    /*****************************************************/
    /* Apply deadtime correction to gated histograms     */
    /*****************************************************/
    string correctedHistoFileName = analysisDirectory + "correctedHistos.root";
    applyDeadtimeCorrection(gatedHistoFileName, deadtimeFileName, histoFileName, gammaCorrectionFileName, log, correctedHistoFileName);

    /*************************************************************************/
    /* Convert TOF histograms into energy in preparation for cross section
     * calculation */
    /*************************************************************************/
    //string energyFileName = analysisDirectory + config.analysis.ENERGY_PLOTS_FILE_NAME;
    //produceEnergyHistos(correctedHistoFileName, log, energyFileName);

    return 0;
}
