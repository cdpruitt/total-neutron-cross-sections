#ifndef WAVEFORM_FITTING_H
#define WAVEFORM_FITTING_H

/*****************************************************************************/
// Waveform-fitting parameters

const int ADC_RANGE  = 16383;  // Range of ADC voltage steps (internal units)

const int THRESHOLD = 25;      // displacement from baseline needed to trigger
                               // software threshold
const int DERIVATIVE_THRESHOLD = -10; // derivative needed to trigger software
                                      // threshold
const int PEAKFIT_WINDOW = 40; // set the size of the window (in number of
                               // samples) where peak-fitting is done on raw
                               // waveforms
const int PEAKFIT_OFFSET = -8; // set the offset of the peak-fitting window,
                               // relative to the raw trigger time (in samples)
// Define time offset of waveforms relative to macropulse clock time
const int WAVEFORM_OFFSET = -960; // in ns

// After a trigger, prevent re-triggering on TRIGGER_HOLDOFF samples
const int TRIGGER_HOLDOFF = 10; // in samples

// If chi-square of a single-peak fit is worse than this, try fitting as a double-peak
const float ERROR_LIMIT = 20.0;

// If chi-square of the double peak fit is worse than this, throw trigger away
// and generate error
const float ERROR_LIMIT_2 = 50.0;

/*****************************************************************************/
// Variables for calculating the baseline:

double BASELINE = 14830;
const int BASELINE_SAMPLES = 10; // number of samples averaged to calculate baseline
const int BASELINE_THRESHOLD = 25; // used to reject samples from being used to calculate the baseline average 
const int BASELINE_LIMIT = 1000; // used to abort waveform fitting if baseline can't be established within this # of samples

/*****************************************************************************/
// Waveform parameters

int numberGoodFits = 0;
int numberBadFits = 0;
int numberOnePeakFits = 0;
int numberOnePeakExpBackFits = 0; // Successfully fit as one peak riding on
                                  // an exponential tail
int numberTwoPeakFits = 0;        // Successfully fit as two peaks

// total number of micropulses processed per target (for performing dead time
// calculation)
std::vector<long> microsPerTarget(6,0);



double onePeakForm(double *x, double *par);
double onePeakExpBackForm(double *x, double *par);
double twoPeakForm(double *x, double *par);

/*****************************************************************************/
// Fitting-function-specific variables

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
const float onePeakParamMin[nParamsOnePeak] = {-300000., -10.,  5, 0.2, 14800., -0.2};
                                             //       A    t1   n    d       C     m
const float onePeakParamMax[nParamsOnePeak] = {    -30.,   2., 20,   3, 14860.,  0.2};

// One-peak-plus-exponential-tail fitting
const float onePeakExpBackParamMin[nParamsOnePeakExpBack] = {-300000., -10.,  5, 0.2, 14800., -0.2, -30000, -5};
                                                           //       A    t1   n    d       C     m       E  tE
const float onePeakExpBackParamMax[nParamsOnePeakExpBack] = {    -30.,   6., 20,   3, 14860.,  0.2,  30000, 20};

// Two-peak fitting
const float twoPeakParamMin[nParamsTwoPeaks] = {-300000., -300000., -10., 10,  5, 0.2, 14800., -0.2};
                                             //        A         B    t1  t2   n    d       C     m
const float twoPeakParamMax[nParamsTwoPeaks] = {    -30.,     -30.,    6, 50, 20,   3, 14860.,  0.2};

// Variables for manual CFD waveform fitting (diagnostic mode)
//const double CFD_SCALEDOWN = 2;
//const double CFD_DELAY = 3;

#endif