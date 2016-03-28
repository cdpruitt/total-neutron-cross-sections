#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
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
// Waveform fitting functions and parameters

const float ch2ns   = 2.;     // time spacing between data points

const int ADC_RANGE  = 16383;  // Number of ADC voltage steps (internal units)
const int THRESHOLD = 25; // displacement from baseline needed to trigger
// slope needed to trigger
const int D_THRESHOLD = -10;

const int npts = 25;    // number of points per pulse

// set the size of the window (in number of samples) where peak-fitting is done
// on raw waveforms
const int PEAKFIT_WINDOW = 20; // in samples

// Set the offset of the processing window relative to the trigger sample
// (i.e., -4 means peak fitting includes 4 samples before the trigger when
// filling the fitting histogram)
const int PEAKFIT_START = -3; // in samples

// time spacing between samples
const int SAMPLE_PERIOD = 2; // in ns

const int BASELINE_SAMPLES = 10;

const int BASELINE_THRESHOLD = 25;

const int BASELINE_LIMIT = 1000;

const int BASELINE_WINDOW = 40;

const int WAVEFORM_OFFSET = -960; // in ns

const int TRIGGER_HOLDOFF = 10; // in samples

// If chi-square of the peak fit is worse than this, try fitting as a double-peak
const float ERROR_LIMIT = 10.0;

// Declare a float to hold the calculated baseline value
float baseline = 14800;

const int nParams = 7;   // number of parameters for fit

// Initial single-peak fitting function parameters
float A_init  = -10;   // amplitude of peak
float trig_init = 0; // trigger time
float n_init  = 1.5;   // exponent of monomial
float d_init  = 3.7;   // decay constant of exponential (in samples)
float C_init  = baseline;  // background is zero
float m_init  = 0;   // background is zero

/*****************************************************************************/
// Cross-section calculation variables and parameters

/* Experimental constants */

const double FLIGHT_DISTANCE = 2672; // detector distance from source, in cm
// 2080 for Neon  

// Period of micropulses
const double MICRO_PERIOD = 1788.820; // in ns

// Physical constants
const double C = 299792458; // speed of light in m/s

const double NEUTRON_MASS = 939.56536; // in MeV/c^2

double avo = 6.022*pow(10.,23.); // Avogadro's number, in atoms/mol


/* Target data */

// State the number of targets currently being cycled over to indicate how many
// histograms should be populated
const int noTargets = 5;

// physical target data, listed in order:
// {blank, Sn112, Natural Sn, Sn124, short carbon, long carbon} 

// lengths of each target:
double targetlength[6] = {0,1.37,2.74,1.365,1.370,1.370}; //cm

// molar mass of each target:
double targetMolMass[6] = {0,12.01,12.01,112,118.7,124}; //g/mol

// density of each target:
double targetdensity[6] = {0,2.2,2.2,6.89,7.31,7.63}; //g/cm^3

// keep track of which order the targets are in, based on which run number we're
// sorting
vector<int> order;


// Declare histograms to hold raw neutron energies, calculated from trigger
// times of waveform peaks
TH1I *blankRaw;
TH1I *target1Raw;
TH1I *target2Raw;
TH1I *target3Raw;
TH1I *target4Raw;
TH1I *target5Raw;

TH1I *blankRawLog;
TH1I *target1RawLog;
TH1I *target2RawLog;
TH1I *target3RawLog;
TH1I *target4RawLog;
TH1I *target5RawLog;

TH1I *TOF;

TH2I *triggerWalk;

TH1I *peakHisto;

TF1 *onePeakFunc;

// Set number of bins for energy histograms
const int noBins = 30;

// declare arrays to hold the scaled cross-sections of each target; i = target #, j = bin
double sigma[6][noBins] = {0};
double sigmaLog[6][noBins] = {0};

// Declare variables to be used for calculating neutron TOFs
int microNo;
double microTime, velocity, rKE;

// Keep track of the number of waveforms collected during each target's period
// in the beam
int targetCounts[noTargets] = {0};

// minimum bound on permissible peak-fit parameter values
float paramMin[nParams] = {0., -6., 0., 2., 14700., -0.5};

// maximum bound on permissible peak-fit parameter values
float paramMax[nParams] = {-2000., 5., 4., 5., 14950., 0.5}; 
// (from Bec) digitizer bits : 0 volts = 128

// Keep track of sample number where triggers are found
vector<int> triggerList;
vector<int> triggerValues;

float chisqMax = 1.5;
float chisqThresh = 0.4;


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

vector<TTree*> orchardW; // holds waveform-mode channel-specific trees
// so that fillHistos can loop through all the trees and make the same histos
// for each tree


/*****************************************************************************/
// Re-link to an already-existing tree's data so we can read the tree
void setBranchesW(TTree* tree)
{
    tree->SetBranchAddress("macroNo",&macroNo);
    tree->SetBranchAddress("evtNo",&evtNo);
    tree->SetBranchAddress("completeTime",&completeTime);
    tree->SetBranchAddress("targetPos",&targetPos);
    tree->SetBranchAddress("waveform",&waveform);
}
/*****************************************************************************/


/*
// Main fitting function
void fitDataPulse()
{

  gROOT->Reset();
  gStyle->SetOptFit(1111);

  int npulses = 200; // number of pulses to analyze
  int peakTime[200][2];
  int peakValu[200][2];
  int nPeaks;
  int junk;

  //float parameters[nparam][npulses];
  float parameters[6][240];
  float chisq[240];

  //float realTime[npts];  // holds true tof
  //float time[npts];      // holds time values for plotting
  //float wave[npts];      // holds voltage values
  float realTime[15];  // holds true tof
  float time[15];      // holds time values for plotting
  float wave[15];      // holds voltage values
  float peak;

  float tmin = -20;
  float tmax =  tmin + ch2ns*npts; // max time on plot default set to one pulse

  //float pulseStart[npulses];
  //float pulseEnd[npulses];

  ifstream ifFile("pulseNeut.out");
  ostringstream titlestring;
  string title;

  TH1F *histChisq = new TH1F("chisq","",50,0,chisqMax);
  TH1F *histParam[6];
  for (int j=0; j<nparam; j++) {
    titlestring.str("");
    if (j==0)  titlestring << "A";
    if (j==1)  titlestring << "t0";
    if (j==2)  titlestring << "n";
    if (j==3)  titlestring << "w";
    if (j==4)  titlestring << "C";
    if (j==5)  titlestring << "m";
    title = titlestring.str();
    histParam[j]= new TH1F(title.c_str(),"",50,paramMin[j],paramMax[j]);
    titlestring << " value";
    title = titlestring.str();
    histParam[j]->GetXaxis()->SetTitle(title.c_str());
    histParam[j]->GetYaxis()->SetTitle("Counts");
  }


  TCanvas *pulsepad = new TCanvas("pulsepad");
  pulsepad->Divide(15,npulses/15);

  TH1F *plotData;

  for (int iP=0;iP<npulses;iP++) {

    // read number of peaks and positions
    ifFile >> nPeaks;
    for (int nP=0; nP<2; nP++) ifFile >> peakTime[iP][nP] >> peakValu[iP][nP];


    for (int i=0; i<npts; i++) {
      ifFile >> realTime[i] >> wave[i];
      if(ifFile.eof()) break;
      if(ifFile.bad()) break;

    // read pulse shape
      if (i==0) peak = wave[i];
      else if (wave[i]>peak) peak = wave[i];       
      time[i] = ch2ns*i;
    } // end loop over points in pulse

    pulsepad->cd(iP+1);  
    plotData = new TH1F("plotData","",npts-1,tmin,tmax);
    plotData->GetXaxis()->SetTitle("time (ns)");
    plotData->GetYaxis()->SetTitle("pulse (V)");
    plotData->SetStats(kFALSE);
    plotData->SetLineColor(2);
    plotData->SetMarkerStyle(7);

    for (int i=0; i<npts; i++) plotData->SetBinContent(i,wave[i]);



//------------------------ Fit pulse -------------------------//

    chisq[iP] = func->GetChisquare();

    for (int j=0; j<nparam; j++) {
      parameters[j][iP] = func->GetParameter(j);
      //if (chisq[iP]<chisqThresh) histParam[j]->Fill(parameters[j][iP]);
      histParam[j]->Fill(parameters[j][iP]);
      //cout << parameters[j][iP]<<endl;
    }

    ostringstream drawParam;
    string drawPar;

    TLatex params;
    params.SetTextSize(0.02);
    for (int j=0; j<nparam; j++) {
      drawParam.str("");
      drawParam << func->GetParName(j) << " ";
      drawParam << parameters[j][iP];
      drawPar = drawParam.str();
      params.DrawLatex( 30,250-7*j,drawPar.c_str() );
    }


    histChisq->Fill(-chisq[iP]/parameters[0][iP]);
    if (-chisq[iP]/parameters[0][iP]>chisqThresh) plotData->Draw("same");
    
    //histChisq->Fill(chisq[iP]);
    //if (chisq[iP]>chisqThresh) plotData->Draw("same");

    TGraph *plotPeaks = new TGraph(nPeaks,peakTime[iP],peakValu[iP]);  
    plotPeaks->SetMarkerStyle(8);  // small dots
    plotPeaks->SetMarkerColor(4);
    plotPeaks->SetLineColor(1);
    plotPeaks->Draw("p");


  } // end loop over pulses

  TCanvas *paramPlot = new TCanvas("paramPlot");

  paramPlot->Divide(2,3);
  for (int j=0; j<nparam; j++) {
    paramPlot->cd(j+1);
    histParam[j]->Draw();
  }

  TLatex fitfunction;
  fitfunction.SetTextSize(0.08);
  fitfunction.SetTextColor(2);
  fitfunction.DrawLatex(0.1,30,
			"fitf = A*(t-t0)^n*exp[-(t-t0)/w] + [C+m*(t-t0)]");


  TCanvas *chisqPlot = new TCanvas("chisqPlot");
  histChisq->Draw();

}*/

/*****************************************************************************/
// Check a set of waveform samples to see if they are suitable for calculating
// the baseline. If so, return the average of the samples and treat as the
// baseline.

// Pass in a set of BASELINE_SAMPLES consecutive samples
float baselineWithinThreshold(vector<int> bw)
{
    // find the average of these consecutive samples
    float baselineAvg=0;

    for(int i=0; i<bw.size(); i++)
    {
        baselineAvg += bw.at(i);
    }

    baselineAvg /= BASELINE_SAMPLES;

    // test each sample to see if it deviates too far from the baseline
    // average (i.e., more than +/- BASELINE_THRESHOLD).
    for(int i=0; i<bw.size(); i++)
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
    for(int i=BASELINE_SAMPLES; i<BASELINE_LIMIT; i++)
    {
        baseline = baselineWithinThreshold(baselineWindow);

        // test to see if the baseline returned by baselineWithinThreshold is
        // non-zero.

        if(baseline==0)
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

    if(baseline < paramMin[4] || baseline > paramMax[4])
    {
        cout << "Baseline outside of bounds: " << baseline << endl;
    }

    return baseline;
}
/*****************************************************************************/

/*****************************************************************************/
// Using the waveform value at time i, check for a peak, and return true if
// there is
bool trigger(int i)
{
    // Check for triggers using:
    //      - a raw threshold above the baseline
    //      - a derivative threshold above the value D_THRESHOLD
    //      - by rejecting triggers if the PREVIOUS point was already above
    //      these thresholds

    if(
            waveform->at(i) <= baseline-THRESHOLD
            && (waveform->at(i)-waveform->at(i-1))/(double)SAMPLE_PERIOD <= D_THRESHOLD

            && (waveform->at(i-1) > baseline-THRESHOLD
            || (waveform->at(i-1)-waveform->at(i-2))/(double)SAMPLE_PERIOD > D_THRESHOLD))
    {
        return true;
    }

    else
    {
        return false;
    }
}


/*****************************************************************************/
void processTrigger(int waveformNo, float triggerSample)
{
    //float triggerTime = triggerSample;

    // Uncomment to use peak fitting to extract trigger times
    if(triggerSample+PEAKFIT_WINDOW < waveform->size())
    {
        stringstream temp;
        temp << "waveform" << waveformNo << "_peak" << triggerList.size();

        // fill a histogram with peak and surrounding environment
        peakHisto = new TH1I(temp.str().c_str(),temp.str().c_str(),PEAKFIT_WINDOW,SAMPLE_PERIOD*PEAKFIT_START,SAMPLE_PERIOD*(PEAKFIT_START+PEAKFIT_WINDOW));

        for (int i=0; i<PEAKFIT_WINDOW; i++)
        {
            peakHisto->SetBinContent(i,waveform->at(triggerSample+PEAKFIT_START+i));
        }

        // reset fitting function to initial parameters in preparation for fitting
        onePeakFunc->SetParameters(A_init,trig_init,n_init,d_init,C_init,m_init);

        onePeakFunc->FixParameter(2,n_init);
        onePeakFunc->FixParameter(3,d_init);

        // fit peak 

        peakHisto->Fit("onePeakFunc","RQ");

        // get fit information
        TF1* onePeakFunc = peakHisto->GetFunction("onePeakFunc");

        //cout << onePeakFunc->GetChisquare() << endl;

        if(onePeakFunc->GetChisquare() < ERROR_LIMIT)
        {
            // we've achieved a good fit with just one peak
            // Extract trigger time from fit
            float triggerTime = triggerSample+onePeakFunc->GetX(baseline-THRESHOLD,onePeakFunc->GetParameter(1),onePeakFunc->GetParameter(1)+10);
            //cout << "Peak " << triggerList.size() << endl;
            //cout << "Trigger time = " << triggerTime << ", parameter(1) = " << onePeakFunc->GetParameter(1) << endl << endl;

            // Extract peak amplitude from peak (in ADC units relative to baseline)
            float peakHeight = onePeakFunc->GetMinimum(PEAKFIT_START,onePeakFunc->GetParameter(1)+10);

            // Calculate derivative of fitted peak at trigger time
            float derivative = onePeakFunc->Derivative(triggerTime);

            // Calculate time of flight from trigger time
            microNo = floor((2*triggerTime+WAVEFORM_OFFSET)/MICRO_PERIOD);
            microTime = fmod((2*triggerTime+WAVEFORM_OFFSET),MICRO_PERIOD);

            // convert microTime into neutron velocity based on flight path distance
            velocity = pow(10.,7.)*FLIGHT_DISTANCE/microTime; // in meters/sec 

            // convert velocity to relativistic kinetic energy
            rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

            TOF->Fill(microTime);

            triggerWalk->Fill(triggerTime,derivative);

            // Fill target-specific plots

            switch (targetPos)
            {
                case 1:
                    // BLANK
                    blankRaw->Fill(rKE);
                    blankRawLog->Fill(TMath::Log10(rKE));
                    break;

                case 2:
                    // TARGET 1
                    target1Raw->Fill(rKE);
                    target1RawLog->Fill(TMath::Log10(rKE));
                    break;

                case 3:
                    // TARGET 2
                    target2Raw->Fill(rKE);
                    target2RawLog->Fill(TMath::Log10(rKE));
                    break;

                case 4:
                    // TARGET 3
                    target3Raw->Fill(rKE);
                    target3RawLog->Fill(TMath::Log10(rKE));
                    break;

                case 5:
                    // TARGET 4
                    target4Raw->Fill(rKE);
                    target4RawLog->Fill(TMath::Log10(rKE));
                    break;

                case 6:
                    // TARGET 5
                    target5Raw->Fill(rKE);
                    target5RawLog->Fill(TMath::Log10(rKE));
                    break;

                default:
                    break;
            }
        }
        
        delete peakHisto;
    }

    triggerList.push_back(triggerSample);
    triggerValues.push_back(waveform->at(triggerSample));

    cout << "Processing trigger " << triggerList.size() << " on waveform " << waveformNo << "\r";
    fflush(stdout);
}
/*****************************************************************************/
/*  Functional form of fit (six parameters)
 *
 *  fitf = -A * (t-t0)^n * exp[-((t-t0)^1)/w] + [C + m*(t-t0)]
 *
 *  par[0] = A
 *  par[1] = t0
 *  par[2] = n
 *  par[3] = w
 *  par[4] = C
 *  par[5] = m
 *
 *  Hope to keep n, w, C? fixed
 *
 */

// Returns the value of the fitted function at a given time
Double_t fitf(Double_t *x, Double_t *par)
{
    Double_t fitval;    // fitted value of function
    Double_t arg = 0;   // argument of exponential

    if (par[3]!=0) arg = pow((x[0]-par[1]),1)/par[3];
    fitval  = exp(-arg); 
    fitval *= par[0] * pow((x[0]-par[1]),par[2]);
    fitval += par[4] + par[5]*(x[0]-par[1]);
    if (x[0]<par[1]) fitval = par[4] + par[5]*(x[0]-par[1]);
;
    return fitval;
}

/*****************************************************************************/
void processWaveforms()
{
    // create raw (unnormalized) neutron energy plots
    blankRaw = new TH1I("blank","blank",noBins,0,700);
    target1Raw = new TH1I("target1","target1",noBins,0,700);
    target2Raw = new TH1I("target2","target2",noBins,0,700);
    target3Raw = new TH1I("target3","target3",noBins,0,700);
    target4Raw = new TH1I("target4","target4",noBins,0,700);
    target5Raw = new TH1I("target5","target5",noBins,0,700);

    // create raw log-scaled neutron energy plots
    blankRawLog = new TH1I("blankLog","blank",noBins,0,TMath::Log10(700));
    target1RawLog = new TH1I("target1Log","target1",noBins,0,TMath::Log10(700));
    target2RawLog = new TH1I("target2Log","target2",noBins,0,TMath::Log10(700));
    target3RawLog = new TH1I("target3Log","target3",noBins,0,TMath::Log10(700));
    target4RawLog = new TH1I("target4Log","target4",noBins,0,TMath::Log10(700));
    target5RawLog = new TH1I("target5Log","target5",noBins,0,TMath::Log10(700));

    TOF = new TH1I("TOF","Summed-detector time of flight",1800,0,MICRO_PERIOD*1.05);

    triggerWalk = new TH2I("triggerWalk","Derivative at trigger time vs. peak amplitude",floor(MICRO_PERIOD*1.05),0,MICRO_PERIOD*1.05,200,-200,0);

    // Define peak-fitting function
    onePeakFunc = new TF1("onePeakFunc",fitf,SAMPLE_PERIOD*PEAKFIT_START,SAMPLE_PERIOD*(PEAKFIT_START+PEAKFIT_WINDOW),6);
    onePeakFunc->SetParNames("A","trig","n","d","C","m");

    // Set limits on parameter values (defined above)
    for (int j=0; j<nParams; j++) {
        onePeakFunc->SetParLimits(j,paramMin[j],paramMax[j]);
    }

    // we'll need to calculate the baseline of each waveform to know when to
    // trigger peaks
    float baseline;

    // Loop through all channel-specific trees
    for(int i = 0; i<orchardW.size(); i++)
    {
        // point event variables at the correct tree in preparation for reading
        // data
        setBranchesW(orchardW[i]);

        // Find total number of events to loop over
        int totalEntries = orchardW[i]->GetEntries();
        cout << "Processing ch. " << 2*i << " waveforms" << endl;

        // EVENT LOOP for sorting through channel-specific waveforms
        for(int j=0; j<totalEntries; j++)
        {
            // reset trigger list
            triggerList.clear();
            triggerValues.clear();

            // pull individual waveform event
            orchardW[i]->GetEntry(j);

            // calculate the baseline for this waveform
            baseline = calculateBaseline();

            // add event to target-position tracker
            if (targetPos > 0)
            {
                targetCounts[targetPos-1]++;
            }

            // Loop through all points in the waveform and fit peaks
            for(int k=BASELINE_WINDOW; k<waveform->size(); k++)
            {
                // Check to see if this point creates a new trigger
                if(trigger(k))
                {
                    // trigger found - plot/fit/extract time

                    processTrigger(j, k);

                    // shift waveform index ahead by TRIGGER_HOLDOFF to prevent
                    // retriggering
                    k += TRIGGER_HOLDOFF;
                }
            }

            stringstream temp;
            temp << "waveform " << j;
            waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform->size()+10,0,2*(waveform->size()+10));

            for(int k=0; k<waveform->size(); k++)
            {
                waveformH->SetBinContent(k,waveform->at(k));
            }
            
            /*
            temp.str("");
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
            */

            TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan");

            if(!mycan)
            {
                mycan = new TCanvas("mycan","mycan");
            }

            //waveformWrap->Write();

            //gPad->Modified();
            //mycan->Update();
            
            /*// Fill trigger histogram
            for (int l = 0; l<triggerList.size(); l++)
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
            triggerH->Draw();
            */

            //triggerH->Write();
            //break;

            //cout << "Finished processing waveform " << j << endl << endl;
        }
    }
}

void calculateCS()
{

    // holds the raw target-specific energy histograms in preparation for populating
    // cross-sections
    vector<TH1I*> rawHistos;
    vector<TH1I*> rawLogHistos;

    rawHistos.push_back(blankRaw);
    rawHistos.push_back(target1Raw);
    rawHistos.push_back(target2Raw);
    rawHistos.push_back(target3Raw);
    rawHistos.push_back(target4Raw);
    rawHistos.push_back(target5Raw);

    rawLogHistos.push_back(blankRawLog);
    rawLogHistos.push_back(target1RawLog);
    rawLogHistos.push_back(target2RawLog);
    rawLogHistos.push_back(target3RawLog);
    rawLogHistos.push_back(target4RawLog);
    rawLogHistos.push_back(target5RawLog);

    for(int k=0; k<noTargets; k++)
    {
        cout << "target position " << k+1 << " counts = " << targetCounts[k] << endl;
    }

    // Loop through the relativistic kinetic energy histograms and use them
    // to populate cross-section and other histograms for each target
    for(int i=0; i<noTargets; i++)
    {
        // Calculate the cross-section for each bin of the energy plots
        // (number of bins set at top of this file)
        for(int j=0; j<noBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(rawHistos[0]->GetBinContent(j) <= 0 || rawHistos[i]->GetBinContent(j) <= 0)
            {
                sigma[i][j] = 0;
            }

            else
            {
                // we must have found positive-definite values found for the raw
                // histogram bins in questions
                // calculate the cross-section and
                // fill the relevant csHisto
                sigma[i][j] = -log((rawHistos[i]->GetBinContent(j)/(double)rawHistos[0]->GetBinContent(j))*(targetCounts[0]/(double)targetCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)avo*pow(10.,-24)/(double)targetMolMass[order[i]]); // in barns
                //cout << "sigma = " << i << ", bin content at " << j << " = " << sigma[i][j] << endl;

                // first, test to make sure we're not about to take log of 0 or
                // divide by 0
                if(rawLogHistos[0]->GetBinContent(j) <= 0 || rawLogHistos[i]->GetBinContent(j) <= 0)
                {
                    sigmaLog[i][j] = 0;
                }

                else
                {
                    // we must have found positive-definite values found for the rawLog
                    // histogram bins in questions; calculate the cross-section and
                    // fill the relevant csHisto

                    sigmaLog[i][j] = -log((rawLogHistos[i]->GetBinContent(j))/((double)rawLogHistos[0]->GetBinContent(j))*(targetCounts[0]/(double)targetCounts[i]))/((double)targetlength[order[i]]*(double)targetdensity[order[i]]*(double)avo*pow(10.,-24)/(double)targetMolMass[order[i]]); // in barns
                }
                //cout << "sigmaLog = " << i << ", bin content at " << j << " = " << sigmaLog[i][j] << endl;

            }
        }
    }
}

void fillCShistos()
{
    // declare the cross-section histograms to be filled
    TH1D *blankcs = new TH1D("blankcs","blank cross-section",noBins,0,700);
    TH1D *target1cs = new TH1D("target1cs","target 1 cross-section",noBins,0,700);
    TH1D *target2cs = new TH1D("target2cs","target 2 cross-section",noBins,0,700);
    TH1D *target3cs = new TH1D("target3cs","target 3 cross-section",noBins,0,700);
    TH1D *target4cs = new TH1D("target4cs","target 4 cross-section",noBins,0,700);
    TH1D *target5cs = new TH1D("target5cs","target 5 cross-section",noBins,0,700);

    TH1D *blankcsLog = new TH1D("blankcsLog","blank cross-section",noBins,0,TMath::Log10(700));
    TH1D *target1csLog = new TH1D("target1csLog","target 1 cross-section",noBins,0,TMath::Log10(700));
    TH1D *target2csLog = new TH1D("target2csLog","target 2 cross-section",noBins,0,TMath::Log10(700));
    TH1D *target3csLog = new TH1D("target3csLog","target 3 cross-section",noBins,0,TMath::Log10(700));
    TH1D *target4csLog = new TH1D("target4csLog","target 4 cross-section",noBins,0,TMath::Log10(700));
    TH1D *target5csLog = new TH1D("target5csLog","target 5 cross-section",noBins,0,TMath::Log10(700));

    // use holder for cross-section histograms to make looping through
    // histograms easier when we calculate cross-sections below
    vector<TH1D*> csHistos, csLogHistos;

    csHistos.push_back(blankcs);
    csHistos.push_back(target1cs);
    csHistos.push_back(target2cs);
    csHistos.push_back(target3cs);
    csHistos.push_back(target4cs);
    csHistos.push_back(target5cs);

    csLogHistos.push_back(blankcsLog);
    csLogHistos.push_back(target1csLog);
    csLogHistos.push_back(target2csLog);
    csLogHistos.push_back(target3csLog);
    csLogHistos.push_back(target4csLog);
    csLogHistos.push_back(target5csLog);

    for(int i=0; i<noTargets; i++)
    {
        for(int j=1; j<noBins; j++)
        {
            // first, test to make sure we're not about to take log of 0 or
            // divide by 0
            if(sigma[i][j] == 0)
            {
                continue;
            }

            else
            {
                //cout << "i = " << i << ", j = " << j << ", sigma[i][j] = " << sigma[i][j] << endl;
                csHistos[i]->SetBinContent(j,sigma[i][j]);
                csLogHistos[i]->SetBinContent(j,sigmaLog[i][j]);
            }
        }
    }
}

int main(int argc, char* argv[])
{
    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    string runDir = argv[1];
    string runNo = argv[2];
    string outpath = argv[3];

    stringstream treeName;
    stringstream scavengerEventsName;
    stringstream summedDetEventsName;

    treeName << "run" << runDir << "-" << runNo; 

    // fileIn/fileOut names to be accessed to open files
    stringstream fileInName, fileOutName, fileCSName;

    fileInName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_sorted.root";
    fileOutName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_waveform.root";

    TFile* file = new TFile(fileInName.str().c_str(),"READ");

    if(file->Get("ch0TreeW"))
    {
        cout << "Located waveform trees in " << fileInName << "." << endl;
    }

    TTree* ch0TreeW = (TTree*)file->Get("ch0TreeW");
    TTree* ch2TreeW = (TTree*)file->Get("ch2TreeW");
    TTree* ch4TreeW = (TTree*)file->Get("ch4TreeW");

    //orchardW.push_back(ch0TreeW);
    //orchardW.push_back(ch2TreeW);
    orchardW.push_back(ch4TreeW);

    // Target order changes between runs. So use the run number to map the
    // correct target to where it was during that run in the target changer

    // Targets are labeled by number as follows:
    // blank = 0, sc = 1, lc = 2, Sn112 = 3, NatSn = 4, Sn124 = 5

    if(stoi(runDir)<=5)
    {
        cout << "Neon run - stopping sort." << endl;
        exit(0);
    }

    else if(stoi(runDir)<=151)
    {
        // blank, short carbon, long carbon, Sn112, NatSn, Sn124
        order.push_back(0);
        order.push_back(1);
        order.push_back(2);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
    }

    else if(stoi(runDir)==152)
    {
        // blank, Sn112, NatSn, Sn124, short carbon, long carbon
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
        order.push_back(1);
        order.push_back(2);
    }

    else if(stoi(runDir)>=153 && stoi(runDir)<=172)
    {
        // blank, Sn112, NatSn, Sn124
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
    }

    else if(stoi(runDir)>=173 && stoi(runDir)<=180)
    {
        // blank, Sn112, NatSn, Sn124, short carbon
        order.push_back(0);
        order.push_back(3);
        order.push_back(4);
        order.push_back(5);
        order.push_back(1);
    }

    // open output file to contain waveform histos
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"RECREATE");

    // Extract triggers from waveforms
    processWaveforms();

    // Calculate cross-sections from waveforms' trigger time data
    calculateCS();
 
    // Plot cross-sections
    fillCShistos();

    fileOut->Write();
    file->Close();
}
