#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TMath.h"

#include "../include/physicalConstants.h"
#include "../include/analysisConstants.h"
#include "../include/plottingConstants.h"
#include "../include/dataStructures.h"
#include "../include/waveformFitting.h"
#include "../include/waveform.h"
#include "../include/plots.h"
#include "../include/branches.h"
#include "../include/experiment.h"

using namespace std;

TH1I *relativeTriggerSampleHisto;
TH2I *triggerWalk;
TH1I *peakHisto;
TF1 *fittingFunc;

extern struct ProcessedEvent procEvent;

// Declare variables to be used for calculating neutron TOFs
int microNo;
double microTime, velocity, rKE;

/*****************************************************************************/

// Keep track of sample number where triggers are found
vector<double> triggerValues;

// for creating a histogram of each waveform
TH1* waveformH;

// for creating a 'wrapped' histogram of each waveform that plots all
// micropulses from the waveform on top of each other
TMultiGraph* waveformWrap;

// To show the trigger locations found on each waveform, we'd like a generic 
// 'holder' that will point to whatever waveform we want to process
TH1 *triggerH;

struct fitData
{
    double trigger1Time = 0;
    double trigger2Time = 0;
    double peak1Amplitude = 0;
    double peak2Amplitude = 0;
    double peak1Derivative = 0;
    double peak2Derivative = 0;
    double chiSquare = 0;
    bool goodFit = false;

    void clear()
    {
        trigger1Time = 0;
        trigger2Time = 0;
        peak1Amplitude = 0;
        peak2Amplitude = 0;
        peak1Derivative = 0;
        peak2Derivative = 0;
        chiSquare = 0;
        goodFit = false;
    }
} data;
