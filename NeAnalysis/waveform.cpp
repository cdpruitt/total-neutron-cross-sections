#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TDirectoryFile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TMath.h"

using namespace std;

/*****************************************************************************/
// Waveform fitting functions and parameters

const float ch2ns   = 2.;     // time spacing between data points

const int ADC_RANGE  = 65535;  // Number of ADC voltage steps (internal units)
const int THRESHOLD = 30; // displacement from baseline needed to trigger
// slope needed to trigger
const int D_THRESHOLD = -0.2;

const int npts = 25;    // number of points per pulse

// set the size of the window (in number of samples) where peak-fitting is done
// on raw waveforms
const int PROCESSING_WINDOW;

const int nParam = 6;   // number of parameters for fit

// minimum bound on permissible parameter values
float paramMin[nParam] = {0., tMin, 0., 0., 119., 0.};

// maximum bound on permissible parameter values
float paramMax[nParam] = {-200., tMax, 4., (tMax-tMin), 126., 0.1}; 
// (from Bec) digitizer bits : 0 volts = 128

// Keep track of sample number where triggers are found
vector<int> triggerList;
vector<int> triggerValues;

float chisqMax = 1.5;
float chisqThresh = 0.4;


/*****************************************************************************/

/*****************************************************************************/
/* event variables */

// variables for holding tree event data to be histogrammed
unsigned int evtNo, macroNo, targetPos;
unsigned int sgQ, lgQ;
double completeTime, macroTime;

vector<int> *waveform; // for holding one event's waveform data

// for creating a histogram of each waveform
TH1* waveformH;
/*****************************************************************************/

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

/*  Functional form of fit (six parameters)
 *
 *  fitf = A * (t-t0)^n * exp[-(t-t0)/w] + [C + m*(t-t0)]
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
    if (x[0]<par[1]) fitval = par[4];
    return fitval;
}

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

    // create fitting object
    TF1 *func  = new TF1("func",fitf,tmin,tmax,6);
    //TF1 *func  = new TF1("func",gaus,tmin,tmax,3);

    // Set initial parameters:  a, t0, n, w, C
    float A_init  = 50;                // half the range
    float t0_init = tmin+2*tmax/npts;  // a couple points into pulse
    float n_init  = 2;                 // squared
    float w_init  = (tmax-tmin)/3;     // a third of whole pulse
    float C_init  = 122;               // background is zero
    float m_init  = 0.5;                 // background is zero
    func->SetParameters(A_init,t0_init,n_init,w_init,C_init,m_init);
    func->SetParNames("A","t0","n","w","C","m");

    // Set limits on parameter values (defined above)
    for (int j=0; j<1; j++) {
       func->SetParLimits(j,paramMin[j],paramMax[j]);
    }

    // Fix certain parameters
    func->FixParameter(2,0.5);         // for new runs
    //func->FixParameter(2,2);         // for new runs
    //func->FixParameter(2,1.457);
    //func->FixParameter(2,1.457);
    //func->FixParameter(2,1.12);
    func->FixParameter(3,10);          // for new runs, exp = 1
    //func->FixParameter(3,13.55); // for exp power = 1.
    //func->FixParameter(3,2.6);   // for exp power = 0.7
    //func->FixParameter(3,1);     // for exp power = 0.5
    //func->FixParameter(4,121);
    func->FixParameter(4,121);
    func->FixParameter(5,0);
    //func->SetParameters(A_init,t0_init,n_init,w_init,C_init);


    plotData->Fit("func","R Q ","p");


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

    for(int i=0; i<bw->size(); i++)
    {
        baselineAvg += bw->at(i);
    }

    baselineAvg /= BASELINE_SAMPLES;

    // test each sample to see if it deviates too far from the baseline
    // average (i.e., more than +/- BASELINE_THREHOLD).
    for(int i=0; i<bw->size(); i++)
    {
        if(bw->at(i) > baselineAvg+BASELINE_THRESHOLD || bw->at(i) < baselineAvg-BASELINE_THRESHOLD)
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
// BASELINE_LIMITS long that are all within +/- BASELINE_THRESHOLD ADC units)
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
        if(!baseline)
        {
            // Baseline calculation failed, so move the baseline window forward
            // one step and try to calculate the baseline again.
            baselineWindow.erase(baselineWindow.begin);
            baselineWindow.push_back(waveform->at(i));
        }

        else
        {
            // Successful baseline calculation. End looping to find a baseline
            // window.
            break;
        }
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
            waveform->at(i) >= baseline+THRESHOLD &&
            (waveform->at(i)-waveform->at(i-1))/(double)SAMPLE_PERIOD=<D_THRESHOLD &&
            waveform->at(i-1) < baseline+THRESHOLD &&
            (waveform->at(i-1)-waveform->at(i-2))/(double)SAMPLE_PERIOD>D_THRESHOLD)
    {
        return true;
    }

    else
    {
        return false;
    }
}
/*****************************************************************************/


/*****************************************************************************/
void processTrigger(k)
{
    triggerList.push_back(k);
    triggerValues.push_back(waveform->at(k));
}
/*****************************************************************************/

/*****************************************************************************/
void processWaveforms()
{
    // we'll need to calculate the baseline of each waveform to know when to
    // trigger peaks
    float baseline;

    // Loop through all channel-specific trees
    for(int i = 0; i<orchardW.size(); i++)
    {
        // point event variables at the correct tree in preparation for reading
        // data
        setBranchesW(orchardW[i]);

        // Find event loop maximum
        int totalEntries = orchardW[i]->GetEntries();
        cout << "Processing ch. " << 2*i << " waveforms" << endl;

        // EVENT LOOP for sorting through channel-specific waveforms
        for(int j=0; j<totalEntries; j++)
        {
            // pull individual waveform event
            orchard[i]->GetEntry(j);

            // calculate the baseline for this waveform
            baseline = calculateBaseline();

            // Loop through all points in the waveform and fit peaks
            for(int k=PROCESSING_WINDOW; k<waveform->size(); k++)
            {
                // Check to see if this point creates a new trigger
                if(trigger(waveform->at(k)))
                {
                    // trigger found - plot/fit/extract time
                    processTrigger(k);
                }
            }


            stringstream temp;
            temp << "waveform " << j;
            waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),32000,0,64000);

            for(int k=0; k<waveform->size(); k++)
            {
                waveformH->SetBinContent(k,waveform->at(k));
            }

            TGraph* triggerGraph = new TGraph(waveform->size(),&triggerList[0],&triggerValues[0])
        }
    }
}


int main(int argc, char* argv[])
{
    // Open the raw tree from the initial sort. If it doesn't exist, exit.
    stringstream treeName;
    stringstream scavengerEventsName;
    stringstream summedDetEventsName;

    treeName << "run" << runDir << "-" << runNo; 

    fileInName << analysispath <<"analysis/run" << runDir << "/" << treeName.str() << "_histos.root";
    fileOutName << analysispath <<"analysis/run" << runDir << "/" << treeName.str() << "_waveform.root";

    TFile* file = new TFile(fileInName.str().c_str(),"READ");

    if(file->Get("ch0TreeW"))
    {
        cout << "Located waveform trees in " << fileInName << "." << endl;
    }

    TTree* ch0TreeW = (TTree*)file->Get("ch0TreeW");
    TTree* ch2TreeW = (TTree*)file->Get("ch2TreeW");
    TTree* ch4TreeW = (TTree*)file->Get("ch4TreeW");

    orchardW.push_back(ch0TreeW);
    orchardW.push_back(ch2TreeW);
    orchardW.push_back(ch4TreeW);

    // open output file to contain histos
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"RECREATE");
    if(!fileOut->IsOpen())
    {
        // No histogram file - need to create it and fill it before moving on to
        // cross-section histos
        TFile* fileOut = new TFile(fileOutName.str().c_str(),"CREATE");

        processWaveforms();
        fileOut->Write();
        file->Close();
    }
}
