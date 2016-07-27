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

using namespace std;

/*****************************************************************************/
// General waveform visualization and fitting parameters

const int SAMPLE_PERIOD = 2;   // digitizer sample length, in ns
//const int ADC_RANGE  = 16383;  // Range of ADC voltage steps (internal units)
const int THRESHOLD = 25;      // displacement from baseline needed to trigger
                               // software threshold
const int DERIVATIVE_THRESHOLD = -10; // derivative needed to trigger software
                                      // threshold
const int PEAKFIT_WINDOW = 10; // set the size of the window (in number of
                               // samples) where peak-fitting is done on raw
                               // waveforms
const int PEAKFIT_OFFSET = -8; // set the offset of the peak-fitting window,
                               // relative to the raw trigger time (in samples)
// Define time offset of waveforms relative to macropulse clock time
const int WAVEFORM_OFFSET = -960; // in ns

// After a trigger, prevent re-triggering on TRIGGER_HOLDOFF samples
const int TRIGGER_HOLDOFF = 10; // in samples

// If chi-square of a single-peak fit is worse than this, try fitting as a double-peak
const float ERROR_LIMIT = 5.0;

// If chi-square of the double peak fit is worse than this, throw trigger away
// and generate error
const float ERROR_LIMIT_2 = 50.0;

/*****************************************************************************/
// Variables for calculating the baseline:

double BASELINE = 14830;
const int BASELINE_SAMPLES = 10; // number of samples averaged to calculate baseline
const int BASELINE_THRESHOLD = 25; // used to reject samples from being used to calculate the baseline average 
const int BASELINE_LIMIT = 128; // used to abort waveform fitting if baseline can't be established within this # of samples

/*****************************************************************************/
// Fitting function-specific variables

const int nParamsOnePeak = 6;
const int nParamsOnePeakExpBack = 8;
const int nParamsTwoPeaks = 8;

// Initial fitting function parameters
float A_init  = -300;     // amplitude of peak
float B_init  = -300;     // amplitude of peak 2
float trig1_init = 0;     // first peak trigger time
float trig2_init = 40;    // second peak trigger time
float n_init  = 8;        // exponent of peak monomial
float d_init  = 1.4;      // decay constant of peak exponential
float C_init  = BASELINE; // background offset
float m_init = 0;         // flat baseline
float E_init  = 0.0;      // previous-peak background: amplitude
float tE_init = 0; // previous-peak background: time offset

// Fitting function parameter limits:
// One-peak fitting
const float onePeakParamMin[nParamsOnePeak] = {-30000., -20.,  5, 0.2, 14800., -0.2};
                                             //       A    t1   n    d       C     m
const float onePeakParamMax[nParamsOnePeak] = {    -30.,  0., 20,   3, 14860.,  0.2};

// Two-peak fitting
const float twoPeakParamMin[nParamsTwoPeaks] = {-30000., -30000.,   -20.,   0,  5, 0.2, 14800., -0.2};
                                             //        A         B    t1   t2   n    d       C     m
const float twoPeakParamMax[nParamsTwoPeaks] = {    -30.,     -30.,   0,  120, 20,   3, 14860.,  0.2};

// Variables for manual CFD waveform fitting (diagnostic mode)
//const double CFD_SCALEDOWN = 2;
//const double CFD_DELAY = 3;

// For functions that use a Fermi function, define the x-axis offset
const double FERMI_OFFSET = 11.5; // in ns

// Indicate the range of times considered to be gamma rays (for the purposes of
// counting gamma rays)
const double GAMMA_WINDOW[2] = {80,100};

int numberGoodFits = 0;
int numberBadFits = 0;
int numberOnePeakFits = 0;
int numberOnePeakExpBackFits = 0; // Successfully fit as one peak riding on
// an exponential tail
int numberTwoPeakFits = 0;        // Successfully fit as two peaks

/*****************************************************************************/
// Cross-section calculation variables and parameters

/* Experimental constants */

const double FLIGHT_DISTANCE = 2672; // detector distance from source, in cm
// 2080 for Neon  
const double MICRO_PERIOD = 1788.814; // in ns
const double C = 299792458; // speed of light in m/s
const double NEUTRON_MASS = 939.56536; // in MeV/c^2
double avo = 6.022*pow(10.,23.); // Avogadro's number, in atoms/mol

/* Target data */

const int NUMBER_OF_TARGETS = 6;

const vector<string> targetNamesWaveform = {"blankWaveform", "shortCarbonWaveform", "longCarbonWaveform", "Sn112Waveform", "NatSnWaveform", "Sn124Waveform"}; 

// physical target data, listed in order:
// {blank, Sn112, Natural Sn, Sn124, short carbon, long carbon} 

double targetlength[6] = {0,1.37,2.74,1.365,1.370,1.370}; //cm
double targetMolMass[6] = {0,12.01,12.01,112,118.7,124}; //g/mol
double targetdensity[6] = {0,2.2,2.2,6.89,7.31,7.63}; //g/cm^3

/* Plotting variables */
const double ENERGY_LOWER_BOUND = 1; // cross-section plots' lower bound, in MeV
const double ENERGY_UPPER_BOUND = 500; // cross-section plots' upper bound, in MeV

// keep track of which order the targets are in, based on which run number we're
// sorting
vector<int> order;

const int TOF_RANGE = 1800; // in ns
const int TOF_BINS = 18000;

struct Plots
{
    vector<TH1I*> TOFHistos;
    vector<TH1I*> energyHistos;
    vector<TH1I*> correctedEnergyHistos;
    vector<TGraph*> CSGraphs;
} plots;

TH1I *relativeTriggerSampleHisto;

TH2I *triggerWalk;

TH1I *peakHisto;

TF1 *fittingFunc;
TF1Convolution *convolvedPeakFunc;

// Set number of bins for energy histograms
const int NUMBER_ENERGY_BINS = 50;

// declare vectors to hold the scaled cross-sections of each target; i = target #, j = bin
vector<vector<double>*> sigma;

vector<double> sigmaXAxis;

// Declare variables to be used for calculating neutron TOFs
int microNo;
double microTime, velocity, rKE;

/*****************************************************************************/

// Keep track of sample number where triggers are found
vector<double> triggerList;

float chisqMax = 1.5;
float chisqThresh = 0.4;

/*****************************************************************************/
/* Define the peak-fitting functional form here:
 *
 *  onePeakForm = A * (t-t0)^n * exp[-((t-t0)^1)/w] + [C + m*(t-t0)]
 *       + B * (t-t1)^n * exp[-((t-t1)^1)/w]
 *
 *  par[0] = A = amplitude of first peak
 *  par[1] = B = amplitude of second peak
 *  par[2] = t1 = starting time of first peak
 *  par[3] = t2 = starting time of second peak
 *  par[4] = n = polynomial order of peak (rise rate)
 *  par[5] = w = exponential decay constant of peak (fall rate)
 *  par[6] = C = baseline offset
 *  par[7] = m = slope of baseline
 *
 *  During the first attempt at fitting, set B = 0 and allow A, t1, C, and m
 *  to vary. If this produces a satisfactory fit (chi squared < ERROR_LIMIT),
 *  then we end fitting and accept the fit's trigger time, amplitude, etc for
 *  this peak.
 *  If this failes to produce a satisfactory fit, then we reset the parameters
 *  to their initial states and allow B and t2 to vary as well (i.e., allow
 *  two peaks in the fit). If this fitting also fails, we skip the peak region
 *  and indicate an error.

Double_t onePeakForm(Double_t *x, Double_t *par)
{
    Double_t fitval;    // fitted value of function
    Double_t arg1 = 0;   // argument of exponential 1
    Double_t arg2 = 0;   // argument of exponential 2

    if (par[5]!=0)
    {
        arg1 = pow((x[0]-par[2]),2)/par[5];
        arg2 = pow((x[0]-par[3]),2)/par[5];
    }

    fitval = par[6] + par[7]*(x[0]-par[2]);
    fitval += par[0]*pow((x[0]-par[2]),par[4])*exp(-arg1);
    fitval += par[1]*pow((x[0]-par[3]),par[4])*exp(-arg2);

    // If before first peak start, set to background + peak's gaussian
    if (x[0]<par[2])
    {
        fitval = par[6] + par[7]*(x[0]-par[2]);
    }

    // If before second peak start, set to background + first peak
    if (x[0]<par[3] && x[0]>=par[2])
    {
        fitval = par[6] + par[7]*(x[0]-par[2]);
        fitval += par[0]*pow((x[0]-par[2]),par[4])*exp(-arg1);
    }

    return fitval;
}
*/

/*****************************************************************************/
/* The peak-fitting functional form is defined as a linear background plus
 * (up to) two peaks, with each peak a convolution of a Maxwell-Boltzmann and a
 * Gaussian.
 *
 * Thus, convolutedPeakFunc =
 *   Integral (A * (i-t1)^n * exp[-((i-t1)^2)/d + (t-i-t1)^2/(2w^2)]) di
 * + Integral (B * (i-t2)^n * exp[-((i-t2)^2)/d + (t-i-t2)^2/(2w^2)]) di
 * + C + m*(t-t0)
 *
 * This form is realized as a TF1Convolution of a user-defined function,
 * fittingFunc, and a generic Gaussian provided by ROOT.
 *
 * fittingFunc = Maxwell-Boltzmann-like distribution + background
 *          = A * (t-t1)^n * exp[-((t-t1)^2)/d] + C + m*(t-t1)
 *
 * Two fittingFuncs are present in the final convolutedPeakFunc expression,
 * and they each have independent amplitudes and time zeroes.
 *
 *  par[0] = A = amplitude of first peak
 *  par[1] = t1 = starting time of first peak
 *  par[2] = n = polynomial order of peak (rise rate)
 *  par[3] = d = exponential decay constant of peak (fall rate)
 *  par[4] = C = baseline offset
 *  par[5] = m = slope of baseline
 *
 *  During the first attempt at fitting, set B = 0 and allow A, t1, C, and m
 *  to vary. If this produces a satisfactory fit (chi squared < ERROR_LIMIT),
 *  then we end fitting and accept the fit's trigger time, amplitude, etc for
 *  this peak.
 *  If this fails to produce a satisfactory fit, allow B and t2 to vary as well
 *  (i.e., allow two peaks in the fit). If this fitting also fails, we ignore
 *  this trigger.
 */

Double_t onePeakForm(Double_t *x, Double_t *par)
{
    Double_t fitval;    // fitted value of function

    if (par[2]==0 || par[3]==0)
    {
        cout << "Error: divide by 0 in onePeakForm" << endl;
        exit(1);
    }

    // Define exponential*fermi function as peak shape
    Double_t arg1 = pow((x[0]-par[1]),1)/par[2];
    Double_t arg2 = pow((x[0]-(par[1]+FERMI_OFFSET)),1)/par[3];

    // Start with linear background
    fitval = par[4] + par[5]*(x[0]-par[1]);

    // add peak
    fitval += par[0]*exp(-arg1)/(1+exp(-arg2));

    return fitval;
}

Double_t onePeakExpBackForm(Double_t *x, Double_t *par)
{
    Double_t fitval;    // fitted value of function

    if (par[2]==0 || par[3]==0)
    {
        cout << "Error: divide by 0 in onePeakOnExpForm" << endl;
        exit(1);
    }

    // Define exponential*fermi function as peak shape
    Double_t arg1 = pow((x[0]-par[1]),1)/par[2];
    Double_t arg2 = pow((x[0]-(par[1]+FERMI_OFFSET)),1)/par[3];

    // Define exponential decay as previous peak background
    Double_t arg3 = pow((x[0]-(par[7])),1)/par[3];

    // Start with linear background
    fitval = par[4] + par[5]*(x[0]-par[1]);

    // Add exponential background of previous peak
    fitval += par[6]*exp(-arg3);

    // add peak
    fitval += par[0]*exp(-arg1)/(1+exp(-arg2));

    return fitval;
}

Double_t twoPeakForm(Double_t *x, Double_t *par)
{
    Double_t fitval;    // fitted value of function

    if (par[4]==0 || par[5]==0)
    {
        cout << "Error: divide by 0 in twoPeakForm" << endl;
        exit(1);
    }

    // Define exponential*fermi function as peak shape
    Double_t arg1 = pow((x[0]-par[2]),1)/par[4];
    Double_t arg2 = pow((x[0]-(par[2]+FERMI_OFFSET)),1)/par[5];

    Double_t arg3 = pow((x[0]-par[3]),1)/par[4];
    Double_t arg4 = pow((x[0]-(par[3]+FERMI_OFFSET)),1)/par[5];


    // Start with linear background
    fitval = par[6] + par[7]*(x[0]-par[2]);

    // add peaks
    fitval += par[0]*exp(-arg1)/(1+exp(-arg2));
    fitval += par[1]*exp(-arg3)/(1+exp(-arg4));

    return fitval;
}

/*Double_t detForm(Double_t *x, Double_t *par)
{
    Double_t fitval;    // fitted value of function
    Double_t arg = 0;   // argument of exponential

    // Uncomment for 'real' detForm
    if (par[2]==0)
    {
        cout << "Error: divide by 0 in detForm" << endl;
        exit(1);
    }

    fitval = par[0]*exp(-pow((x[0]-par[1]),2)/(2*pow(par[2],2)));

    return fitval;
}*/

/*Double_t cfdForm(Double_t *x, Double_t *par)
{
    Double_t fitval;    // fitted value of function

    // uncomment for gaussian peakform definition
    //fitval = par[0]*exp(-pow((x[0]-par[1]),2)/(2*pow(par[2],2)));

    // uncomment for 'correct' peakform definition
    Double_t arg1 = 0;   // argument of peak 1 exponential
    Double_t arg2 = 0;   // argument of peak 2 exponential

    if (par[5]==0)
    {
        cout << "Error: divide by 0 in onePeakForm" << endl;
        exit(1);
    }

    arg1 = pow((x[0]-par[2]),1)/par[4];
    arg2 = pow((x[0]-(par[2]+11.5)),1)/par[5];

    // define background
    fitval = par[6] + par[7]*(x[0]-par[2]);

    // add scaled-down part of CFD form
    fitval += (par[0]/CFD_SCALEDOWN)*exp(-arg1)/(1+par[1]*exp(-arg2));

    // add delayed part of CFD form
    arg1 = pow((x[0]-(par[2]+CFD_DELAY)),1)/par[4];
    arg2 = pow((x[0]-(par[2]+11.5+CFD_DELAY)),1)/par[5];
    fitval += -par[0]*exp(-arg1)/(1+par[1]*exp(-arg2));

    return fitval;
}*/


/*****************************************************************************/
/* event variables */

// variables for holding tree event data to be histogrammed
unsigned int evtNo, macroNo, targetPos;
unsigned int sgQ, lgQ;
double completeTime, macroTime;

vector<int> *waveform; // for holding one event's waveform data

// for creating a histogram of each waveform
TH1* waveformH;

// for creating a 'wrapped' histogram of each waveform that plots all
// micropulses from the waveform on top of each other
TMultiGraph* waveformWrap;

// To show the trigger locations found on each waveform, we'd like a generic 
// 'holder' that will point to whatever waveform we want to process
TH1 *triggerH;


/*****************************************************************************/
/* Tree variables */

vector<TTree*> orchard; // holds waveform-mode channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree


/*****************************************************************************/
// Re-link to an already-existing tree's data so we can read the tree
void setBranchesW(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("macroTime",&macroTime);
    tree->SetBranchAddress("targetPos",&targetPos);
    tree->SetBranchAddress("waveform",&waveform);
}

/*****************************************************************************/
// Check a set of waveform samples to see if they are suitable for calculating
// the baseline. If so, return the average of the samples and treat as the
// baseline.

// Pass in a set of BASELINE_SAMPLES consecutive samples
float baselineWithinThreshold(vector<int> bw)
{
    // find the average of these consecutive samples
    float baselineAvg=0;

    for(int i=0; (size_t)i<bw.size(); i++)
    {
        baselineAvg += bw.at(i);
    }

    baselineAvg /= BASELINE_SAMPLES;

    // test each sample to see if it deviates too far from the baseline
    // average (i.e., more than +/- BASELINE_THRESHOLD).
    for(int i=0; (size_t)i<bw.size(); i++)
    {
        if(bw.at(i) > baselineAvg+BASELINE_THRESHOLD || bw.at(i) < baselineAvg-BASELINE_THRESHOLD)
        {
            // if one of the samples deviates too far from the baseline average,
            // then the set of samples being used to calculate the baseline is
            // NOT representative of the baseline. Thus, return false, and a new
            // set of baseline samples will be used.
            return 0;
        }
    }
    // All the samples in the baseline window are within BASELINE_THRESHOLD of
    // each other, so they should be appropriate for calculating the baseline
    // for this waveform.
    return baselineAvg;
}
/*****************************************************************************/

/*****************************************************************************/
// Starting at the beginning of a waveform, try to find a window of samples
// BASELINE_LIMIT long that are all within +/- BASELINE_THRESHOLD ADC units)
// of each other. Average these points to find the baseline.
float calculateBaseline()
{
    vector<int> baselineWindow;

    // Load the first BASELINE_SAMPLES waveform samples into the baseline window
    for(int i=0; i<BASELINE_SAMPLES; i++)
    {
        baselineWindow.push_back(waveform->at(i));
    }

    // Loop through the waveform to try to find a BASELINE_SAMPLES-long stretch
    // that can be used to calculate the baseline. Truncate the search at
    // BASELINE_LIMIT samples.
    for(int i=BASELINE_SAMPLES; i<waveform->size(); i++)
    {
        BASELINE = baselineWithinThreshold(baselineWindow);

        // test to see if the baseline returned by baselineWithinThreshold is
        // non-zero.

        if(BASELINE==0)
        {
            // Baseline calculation failed, so move the baseline window forward
            // one step and try to calculate the baseline again.
            baselineWindow.erase(baselineWindow.begin());
            baselineWindow.push_back(waveform->at(i));
        }

        else
        {
            // Successful baseline calculation. End looping to find a baseline
            // window.
            break;
        }
    }

    /*if(BASELINE < paramMin[6] || BASELINE > paramMax[6])
    {
        cout << "Baseline outside of bounds: " << BASELINE << endl;
    }*/

    return BASELINE;
}
/*****************************************************************************/

/*****************************************************************************/
// Using the waveform value at time i, check for a peak, and return true if
// there is
bool isTrigger(int i)
{
    // Check for triggers using:
    //      - a raw threshold above the baseline
    //      - a derivative threshold above the value DERIVATIVE_THRESHOLD
    //      - by rejecting triggers if the PREVIOUS point was already above
    //      these thresholds

    if((waveform->at(i) <= BASELINE-THRESHOLD
       && (waveform->at(i)-waveform->at(i-1))/(double)SAMPLE_PERIOD <= DERIVATIVE_THRESHOLD)

       && (waveform->at(i-1) > BASELINE-THRESHOLD
       || (waveform->at(i-1)-waveform->at(i-2))/(double)SAMPLE_PERIOD > DERIVATIVE_THRESHOLD))
    {
        return true;
    }

    return false;
}

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

/*double cfd(TF1* fit, double triggerTime)
{
    // create scaled-down part of CFD
    TF1* cfdFunc = new TF1("cfdFunc",cfdForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW),8);
    cfdFunc->SetParameters(fit->GetParameters());

    return cfdFunc->GetX(BASELINE,triggerTime-5,triggerTime+5);
}*/

fitData fitTrigger(int waveformNo, vector<double> rawTriggerList)
{
    // A trigger has been detected on the current waveform; fitTrigger attempts
    // to fit the region around this trigger using a series of progressively
    // more complicated fitting functions.
    // If any of these functions fits the waveform to within an error boundary,
    // fitTrigger accepts it as a good fit and exits.

    // reset fit data
    data.clear();

    // extract the waveform chunk we'd like to fit
    stringstream temp;
    temp << "waveform" << waveformNo << "_peak";

    peakHisto = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size(),SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+waveform->size()));

    for (int i=0; i<waveform->size(); i++)
    {
        peakHisto->SetBinContent(i,waveform->at(i));
    }

    /*************************************************************************/
    /*if(rawTriggerList.size() == 1)
    {
        fittingFunc = new TF1("fittingFunc",onePeakForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+waveform->size()),nParamsOnePeak);
        fittingFunc->SetParameters(A_init,trig1_init,n_init,d_init,C_init,m_init);

        // Lock in the shape of the peak (manually chosen to match real peak shape)
        fittingFunc->FixParameter(2,n_init);
        fittingFunc->FixParameter(3,d_init);

        // Lock in a constant background
        fittingFunc->FixParameter(5,m_init);

        // Set parameter boundaries
        for (int i=0; i<fittingFunc->GetNpar(); i++)
        {
            fittingFunc->SetParLimits(i,onePeakParamMin[i],onePeakParamMax[i]);
        }

        // fit peak 
        peakHisto->Fit("fittingFunc","RQ");

        if(fittingFunc->GetChisquare() < ERROR_LIMIT)
        {
            // Success - we've achieved a good fit with just one peak
            // Output fit data

            data.peak1Amplitude = fittingFunc->GetMinimum(fittingFunc->GetParameter(1),fittingFunc->GetParameter(1)+10);

            double relativeTriggerSample = fittingFunc->GetX((fittingFunc->GetParameter(4)+data.peak1Amplitude)/2,fittingFunc->GetParameter(1)-10,fittingFunc->GetParameter(1)+20);

            data.trigger1Time = SAMPLE_PERIOD*(triggerSample+relativeTriggerSample);
            data.peak1Derivative = fittingFunc->Derivative(relativeTriggerSample);
            data.chiSquare = fittingFunc->GetChisquare();
            data.goodFit = true;

            relativeTriggerSampleHisto->Fill(relativeTriggerSample);

            //cout << "getX from peakFit = " << relativeTriggerSample << endl;
            //cout << "derivative = " << data.peak1Derivative << endl;

            //cout << "monomial order = " << fittingFunc->GetParameter(2) << endl;

            numberOnePeakFits++;
        }
    }

    else if(rawTriggerList.size() == 2)
    {
        fittingFunc = new TF1("fittingFunc",twoPeakForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+waveform->size()),nParamsTwoPeaks);
        fittingFunc->SetParameters(A_init,B_init,trig1_init,trig2_init,n_init,d_init,C_init,m_init);

        // Lock in the shape of the peak (manually chosen to match real peak shape)
        fittingFunc->FixParameter(4,n_init);
        fittingFunc->FixParameter(5,d_init);

        // Lock in a constant background
        fittingFunc->FixParameter(7,m_init);

        // Set parameter boundaries
        for (int i=0; i<fittingFunc->GetNpar(); i++)
        {
            fittingFunc->SetParLimits(i,twoPeakParamMin[i],twoPeakParamMax[i]);
        }

        // fit peak 
        peakHisto->Fit("fittingFunc","RQ");

        if(fittingFunc->GetChisquare() < ERROR_LIMIT_2)
        {
            data.peak1Amplitude = fittingFunc->GetMinimum(fittingFunc->GetParameter(2),fittingFunc->GetParameter(2)+10);
            data.peak2Amplitude = fittingFunc->GetMinimum(fittingFunc->GetParameter(3),fittingFunc->GetParameter(3)+10);

            double relativeTriggerSample = fittingFunc->GetX((fittingFunc->GetParameter(6)+data.peak1Amplitude)/2,fittingFunc->GetParameter(2)-10,fittingFunc->GetParameter(2)+10);
            data.trigger1Time = SAMPLE_PERIOD*(triggerSample+relativeTriggerSample);
            data.peak1Derivative = fittingFunc->Derivative(relativeTriggerSample);
            relativeTriggerSampleHisto->Fill(relativeTriggerSample);

            relativeTriggerSample = fittingFunc->GetX((fittingFunc->GetParameter(6)+data.peak2Amplitude)/2,fittingFunc->GetParameter(3)-5,fittingFunc->GetParameter(2)+5);
            data.trigger2Time = SAMPLE_PERIOD*(triggerSample+relativeTriggerSample);
            data.peak2Derivative = fittingFunc->Derivative(relativeTriggerSample);
            relativeTriggerSampleHisto->Fill(relativeTriggerSample);

            data.chiSquare = fittingFunc->GetChisquare();
            data.goodFit = true;

            numberTwoPeakFits++;
        }
    }
*/

    
    /*************************************************************************/

            /*
        else
        {
            // failed to fit w/ two peaks; try allowing three peaks
            fittingFunc = new TF1("fittingFunc",threePeakForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW),nParamsthreePeaks);
            fittingFunc->SetParameters(A_init,B_init,trig1_init,trig2_init,n_init,d_init,C_init,m_init);

            // Lock in the shape of the peak (manually chosen to match real peak shape)
            fittingFunc->FixParameter(4,n_init);
            fittingFunc->FixParameter(5,d_init);

            // Lock in a constant background
            fittingFunc->FixParameter(7,m_init);

            // Set parameter boundaries
            for (int i=0; i<fittingFunc->GetNpar(); i++)
            {
                fittingFunc->SetParLimits(i,threePeakParamMin[i],threePeakParamMax[i]);
            }

            // fit peak 
            peakHisto->Fit("fittingFunc","RQ");

            // see if peak fitting was successful
            if(fittingFunc->GetChisquare() < ERROR_LIMIT_2)
            {
                data.peak1Amplitude = fittingFunc->GetMinimum(fittingFunc->GetParameter(2),fittingFunc->GetParameter(2)+10);
                data.peak2Amplitude = fittingFunc->GetMinimum(fittingFunc->GetParameter(3),fittingFunc->GetParameter(3)+10);

                double relativeTriggerSample = fittingFunc->GetX((fittingFunc->GetParameter(6)+data.peak1Amplitude)/2,fittingFunc->GetParameter(2)-10,fittingFunc->GetParameter(2)+10);
                data.trigger1Time = SAMPLE_PERIOD*(triggerSample+relativeTriggerSample);
                data.peak1Derivative = fittingFunc->Derivative(relativeTriggerSample);
                relativeTriggerSampleHisto->Fill(relativeTriggerSample);

                relativeTriggerSample = fittingFunc->GetX((fittingFunc->GetParameter(6)+data.peak2Amplitude)/2,fittingFunc->GetParameter(3)-5,fittingFunc->GetParameter(2)+5);
                data.trigger2Time = SAMPLE_PERIOD*(triggerSample+relativeTriggerSample);
                data.peak2Derivative = fittingFunc->Derivative(relativeTriggerSample);
                relativeTriggerSampleHisto->Fill(relativeTriggerSample);

                data.chiSquare = fittingFunc->GetChisquare();
                data.goodFit = true;

                numberthreePeakFits++;
            }
            */

    //delete peakHisto;
    //delete fittingFunc;

    /*if(data.chiSquare<5)
      {
      delete peakHisto;
      }*/

    // failed to fit
    return data;
}

void fillTriggerHistos(double triggerTime, int waveformNo)
{
    // Calculate time of flight from trigger time
    microTime = fmod((triggerTime),MICRO_PERIOD);

    microNo = floor((triggerTime)/MICRO_PERIOD);

    // convert microTime into neutron velocity based on flight path distance
    velocity = pow(10.,7.)*FLIGHT_DISTANCE/microTime; // in meters/sec 

    // convert velocity to relativistic kinetic energy
    rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

    triggerWalk->Fill(microTime,waveformNo);

    if (targetPos>0 && targetPos<=NUMBER_OF_TARGETS)
    {
        plots.TOFHistos[targetPos-1]->Fill(microTime);
        plots.energyHistos[targetPos-1]->Fill(rKE);
    }
}

/*****************************************************************************/
void processTriggers(int waveformNo, vector<double> rawTriggerList)
{
    // Uncomment to use raw trigger sample as trigger time
    for(int i=0; i<rawTriggerList.size(); i++)
    {
        fillTriggerHistos(rawTriggerList[i]+completeTime-macroTime, waveformNo);
    }

    // Uncomment to use fitted peak threshold-intercept as trigger time
    /*if(fitTrigger(waveformNo, rawTriggerList).goodFit)
    {
        triggerList.push_back(data.trigger1Time+completeTime);
        numberGoodFits++;
    }

    else
    {
        numberBadFits++;
    }
    */

    if(waveformNo%1000==0)
    {
        cout << "Processing triggers on waveform " << waveformNo << "\r";
        fflush(stdout);
    }
}


/*****************************************************************************/
double tofToRKE(double TOF)
{
    double velocity = pow(10.,7.)*FLIGHT_DISTANCE/TOF; // in meters/sec 
    
    if (velocity>C)
    {
        return -1;
    }

    // convert velocity to relativistic kinetic energy
    double RKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV
    if(RKE<0)
    {
        return -1;
    }
    return RKE;
}

vector<double> scaleBins(vector<double> inputBins, int nInputBins, int scaledown)
{
    vector<double> outputBins((int)floor(nInputBins/scaledown));
    for(int i=0; i<nInputBins; i++)
    {
        if(i>floor(nInputBins/scaledown)*scaledown)
        {
            break;
        }
        outputBins[(int)floor(i/scaledown)] += inputBins[i];
    }

    for(int i=0; i<nInputBins/scaledown; i++)
    {
        outputBins[i] /= scaledown;
    }

    return outputBins;
}

// Map an input histogram with bins in the time domain to equivalent bins in the relativistic
// kinetic energy domain
TH1* timeBinsToRKEBins(TH1 *inputHisto)
{
    string newName;
    newName = inputHisto->GetName();
    newName += "RKE";

    TAxis* oldAxis = inputHisto->GetXaxis();
    int nOldBins = oldAxis->GetNbins();

    double minimumTime = (((TAxis*)inputHisto->GetXaxis())->GetXmin());
    int minimumBin = 0;

    for(int i=0; i<nOldBins; i++)
    {
        if(tofToRKE(minimumTime)>0 && tofToRKE(minimumTime)<ENERGY_UPPER_BOUND)
        {
            break;
        }

        minimumTime = inputHisto->GetBinCenter(i);
        minimumBin = i;
    }

    int nUnscaledEnergyBins = nOldBins-minimumBin;
    double tentativeEnergy = tofToRKE(minimumTime);
    if(tentativeEnergy==-1)
    {
        cout << "Error: energy of old min time " << minimumTime << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    double newXMax = tentativeEnergy;

    double oldXMax = (((TAxis*)inputHisto->GetXaxis())->GetXmax());

    tentativeEnergy = tofToRKE(oldXMax);
    if(tentativeEnergy==-1)
    {
        cout << "Error: energy of old max time " << oldXMax << " was not finite: " << tentativeEnergy << " (MeV)" << endl;
        exit(1);
    }

    double newXMin = tentativeEnergy;

    // Remap bins from old histo to new histo
    vector<double> unscaledEnergyBins(nUnscaledEnergyBins);

    // Reorder bins to go from lowest energy (shortest time) to highest energy (longest time)
    for(int i=0; i<unscaledEnergyBins.size(); i++)
    {
        unscaledEnergyBins[i] = tofToRKE(oldAxis->GetBinCenter(nOldBins-i-1));
    }

    // Downscale bins to desired granularity
    vector<double> scaledEnergyBins = scaleBins(unscaledEnergyBins, nUnscaledEnergyBins, TOF_BINS/NUMBER_ENERGY_BINS);
    // n bins are defined n+1 points (like fence sections and fence posts)
    scaledEnergyBins.push_back(newXMax);

    int nScaledEnergyBins = floor(nUnscaledEnergyBins/(TOF_BINS/NUMBER_ENERGY_BINS));

    TH1* outputHisto = new TH1D(newName.c_str(),
                                newName.c_str(),
                                nScaledEnergyBins,
                                newXMin,
                                newXMax);

    // Assign the remapped bins to the new histo
    ((TAxis*)outputHisto->GetXaxis())->Set(nScaledEnergyBins,
&scaledEnergyBins[0]);

    return outputHisto;
}

void processWaveforms()
{
    TH1I *blankTOF = new TH1I("blankTOF","blank TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target1TOF = new TH1I("target1TOF","target 1 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target2TOF = new TH1I("target2TOF","target 2 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target3TOF = new TH1I("target3TOF","target 3 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target4TOF = new TH1I("target4TOF","target 4 TOF",TOF_BINS,0,TOF_RANGE);
    TH1I *target5TOF = new TH1I("target5TOF","target 5 TOF",TOF_BINS,0,TOF_RANGE);

    plots.TOFHistos.push_back(blankTOF);
    plots.TOFHistos.push_back(target1TOF);
    plots.TOFHistos.push_back(target2TOF);
    plots.TOFHistos.push_back(target3TOF);
    plots.TOFHistos.push_back(target4TOF);
    plots.TOFHistos.push_back(target5TOF);

    triggerWalk = new TH2I("triggerWalk","trigger time vs. waveform chunk #",200,0,200,1000,0,1000);

    relativeTriggerSampleHisto = new TH1I("relativeTriggerSampleHisto","relative trigger time, from start of fitted wavelet",100,PEAKFIT_OFFSET*SAMPLE_PERIOD,(PEAKFIT_OFFSET+PEAKFIT_WINDOW)*SAMPLE_PERIOD);

    // create raw (unnormalized) neutron energy plots
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        plots.energyHistos.push_back((TH1I*)timeBinsToRKEBins(plots.TOFHistos[i])); 
    }

    // we'll need to calculate the baseline of each waveform to know when to
    // trigger peaks
    float BASELINE;

    // Loop through all channel-specific trees
    for(int i=0; (size_t)i<orchard.size(); i++)
    {
        // point event variables at the correct tree in preparation for reading
        // data
        setBranchesW(orchard[i]);

        int totalEntries = orchard[i]->GetEntries();
        cout << "Total DPP events on ch. " << 2*i << ": " << totalEntries << endl;

        TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan");

        if(!mycan)
        {
            mycan = new TCanvas("mycan","mycan");
        }

        // EVENT LOOP for sorting through channel-specific waveforms
        for(int j=0; j<totalEntries; j++)
        {
            triggerList.clear();

            // pull individual waveform event
            orchard[i]->GetEntry(j);

            //cout << "waveform chunk time = " << completeTime << endl;

            // calculate the baseline for this waveform
            BASELINE = calculateBaseline();

            vector<double> rawTriggerList;

            // Loop through all points in the waveform and fit peaks
            for(int k=2; (size_t)k<waveform->size(); k++)
            {
                // Check to see if this point creates a new trigger
                if(isTrigger(k))
                {
                    // trigger found - plot/fit/extract time
                    rawTriggerList.push_back(k);

                    // shift waveform index ahead by TRIGGER_HOLDOFF to prevent
                    // retriggering
                    //k += TRIGGER_HOLDOFF;
                }
                /*if(triggerList.size()>10)
                {
                    for(int m=0; m<triggerList.size(); m++)
                    {
                        cout << "triggerList[" << m << "] = " << triggerList[m] << endl;
                    }
                    break;
                }*/
            }

            processTriggers(j, rawTriggerList);

            /*stringstream temp;
            temp << "waveform " << j;
            waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size()+10,0,2*(waveform->size()+10));

            for(int k=0; k<waveform->size(); k++)
            {
                waveformH->SetBinContent(k,waveform->at(k));
            }
            */

            /*temp.str("");
            temp << "waveformWrap" << j;

            waveformWrap = new TMultiGraph(temp.str().c_str(), temp.str().c_str());

            vector<TGraph*> microGraphs;

            // Create a new graph for each micropulse period to be plotted
            for (int m = 0; m<floor(2*waveform->size()/MICRO_PERIOD); m++)
            {
                microGraphs.push_back(new TGraph());
            }

            // Fill each micropulse graph with waveform samples
            for (int l = 0; l<waveform->size(); l++)
            {
                microGraphs[(int)floor(l/(double)MICRO_PERIOD)]->SetPoint(microGraphs[(int)floor(l/(double)MICRO_PERIOD)]->GetN(),fmod(2*l+WAVEFORM_OFFSET,MICRO_PERIOD),waveform->at(l));
                //cout << "Adding value " << waveform->at(l) << " to position " << fmod(l,MICRO_PERIOD) << " in microGraph " << floor(l/(double)MICRO_PERIOD) << endl;
            }

            // Add each graph to the MultiGraph
            for (int m = 0; m<microGraphs.size(); m++)
            {
                //cout << "adding graph " << m << " to multigraph" << endl;
                microGraphs[m]->Draw();
                waveformWrap->Add(microGraphs[m],"l");
            }

            waveformWrap->Write();
            */

            //gPad->Modified();
            //mycan->Update();

            // Fill trigger histogram
            /*for (int l = 0; l<triggerList.size(); l++)
              {
              cout << "trigger " << l << " = " << triggerList[l] << ", " << triggerValues[l] << endl;
              }*/

            /*
            temp << "triggers";
            triggerH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size()+10,0,2*(waveform->size()+10));

            for(int k=0; k<triggerList.size(); k++)
            {
                triggerH->SetBinContent(triggerList[k],triggerValues[k]);
            }

            TCanvas *c1 = new TCanvas;
            c1->DrawFrame(0,0,waveform->size()+10,16383);

            //waveformH->Draw("same");
            triggerH->SetOption("P");

            triggerH->SetMarkerStyle(29);
            triggerH->SetMarkerSize(2);
            triggerH->SetMarkerColor(2);
            */
            //triggerH->Draw();

            //triggerH->Write();
            //cout << "Finished processing waveform " << j << endl << endl;

            /*if(j>1)
            {
                break;
            }*/
        }
    }

    cout << endl;
}

void calculateCS(string monitorFileName)
{
    TFile* monitorFile = new TFile(monitorFileName.c_str(),"READ");
    if(!monitorFile->IsOpen())
    {
        cout << "Error: failed to open histos.root for monitor counts" << endl;
        exit(1);
    }

    // switch to the monitor directory
    gDirectory->cd("/");
    gDirectory->GetDirectory("monitor")->cd();

    // to normalize flux between all channels, we'll need to keep track of the
    // number of counts that came in through the monitor during that target's
    // beam time
    vector<long> monCounts;

    for(int i=2; i<8; i++)
    {
        monCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i));
        cout << "target position " << i << " monitor counts = " << monCounts[i-2] << endl;
    }
    monitorFile->Close();

    // Loop through the relativistic kinetic energy histograms and use them
    // to populate cross-section and other histograms for each target

    int numberOfBins = ((TAxis*)plots.energyHistos[0]->GetXaxis())->GetNbins();

    for(int i=0; (size_t)i<order.size(); i++)
    {
        sigma.push_back(new vector<double>);
    }

    for(int i=0; (size_t)i<order.size(); i++)
    {
        for(int j=0; j<numberOfBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(plots.energyHistos[0]->GetBinContent(j) <= 0 || plots.energyHistos[i]->GetBinContent(j) <= 0)
            {
                sigma[order[i]]->push_back(0);
            }

            else
            {
                // calculate the cross section and
                sigma[order[i]]->push_back(-log((plots.energyHistos[i]->GetBinContent(j)/(double)plots.energyHistos[0]->GetBinContent(j))*(monCounts[0]/(double)monCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)avo*pow(10.,-24)/(double)targetMolMass[order[i]])); // in barns
            }

            if(i==0)
            {
                sigmaXAxis.push_back(plots.energyHistos[0]->GetBinCenter(j));
            }
        }
    }
}

void fillCSGraphs(TFile* outFile)
{
    for(int i=0; i<order.size(); i++)
    {
        outFile->cd();
        plots.CSGraphs.push_back(new TGraphErrors(sigmaXAxis.size(),&sigmaXAxis[0],&(sigma[order[i]]->at(0))));
        plots.CSGraphs[i]->SetNameTitle(targetNamesWaveform[order[i]].c_str(),targetNamesWaveform[order[i]].c_str());
        plots.CSGraphs[i]->Write();
    }
}

int main(int argc, char* argv[])
{
    string inFileName = argv[1];
    inFileName += "resort.root";

    string outFileName = argv[1];
    outFileName += "DPPwaveform.root";

    string monitorFileName = argv[1];
    monitorFileName += "histos.root";

    string runNumber = argv[2];

    TFile* inFile = new TFile(inFileName.c_str(),"READ");
    if(!inFile->IsOpen())
    {
        cout << "Error: failed to open resort.root" << endl;
        exit(1);
    }

    if(inFile->Get("ch4ProcessedTree"))
    {
        cout << "Located DPP tree in " << inFileName << "." << endl;
    }

    //TTree* ch0Tree = (TTree*)file->Get("targetChangerTree");
    //TTree* ch2Tree = (TTree*)file->Get("ch2ProcessedTreeW");
    TTree* ch4Tree = (TTree*)inFile->Get("ch4ProcessedTree");

    //orchard.push_back(ch0TreeW);
    //orchard.push_back(ch2TreeW);
    orchard.push_back(ch4Tree);

    // Target order changes between runs. So use the run number to map the
    // correct target to where it was during that run in the target changer

    // Targets are labeled by number as follows:
    // blank = 0, sc = 1, lc = 2, Sn112 = 3, NatSn = 4, Sn124 = 5

    if(stoi(runNumber)<=5)
    {
        cout << "Neon run - stopping sort." << endl;
        exit(0);
    }

    if(stoi(runNumber)<=151)
    {
        // blank, short carbon, long carbon, Sn112, NatSn, Sn124
        order.push_back(0);
        order.push_back(1);
        order.push_back(2);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
    }

    else if(stoi(runNumber)==152)
    {
        // blank, Sn112, NatSn, Sn124, short carbon, long carbon
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
        order.push_back(1);
        order.push_back(2);
    }

    else if(stoi(runNumber)>=153 && stoi(runNumber)<=168)
    {
        // blank, Sn112, NatSn, Sn124
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
    }

    else if(stoi(runNumber)>=169 && stoi(runNumber)<=180)
    {
        // blank, Sn112, NatSn, Sn124, short carbon
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
        order.push_back(1);
    }

    // open output file to contain waveform histos
    TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

    // Extract triggers from waveforms
    processWaveforms();

    cout << "Number of good fits: " << numberGoodFits << endl;
    cout << "onePeak = " << numberOnePeakFits << endl; 
    cout << "onePeakExpBack = " << numberOnePeakExpBackFits << endl; 
    cout << "twoPeaks = " << numberTwoPeakFits << endl << endl; 
    cout << "Number of bad fits: " << numberBadFits << endl;

    // Calculate cross-sections from waveforms' trigger time data
    calculateCS(monitorFileName);
 
    // Plot cross-sections
    fillCSGraphs(outFile);

    outFile->Write();
    inFile->Close();
}
