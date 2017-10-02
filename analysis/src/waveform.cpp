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
#include "TVirtualFFT.h"
#include "../include/physicalConstants.h"
#include "../include/dataStructures.h"
#include "../include/target.h"
#include "../include/waveformFitting.h"
#include "../include/waveform.h"
#include "../include/plots.h"
#include "../include/branches.h"
#include "../include/experiment.h"

using namespace std;

TH2I *triggerWalk;
TH1I *peakHisto;
TF1 *fittingFunc;

unsigned int DPP_WAVEFORM_SAMPLES = 60;
unsigned int DPP_PEAKFIT_START = 2;
double DPP_PEAKFIT_OFFSET = -48.3;

extern struct ProcessedEvent procEvent;

// Declare variables to be used for calculating neutron TOFs
int microNo;
double microTime, velocity, rKE;

double prevTriggerTime = 0;
double eventTimeDiff = 0;

/*****************************************************************************/

// Keep track of sample number where triggers are found
vector<double> triggerValues;

/*****************************************************************************/

// for creating a histogram of each waveform
TH1* waveformH;

// for creating a 'wrapped' histogram of each waveform that plots all
// micropulses from the waveform on top of each other
TMultiGraph* waveformWrap;

// To show the trigger locations found on each waveform, we'd like a generic 
// 'holder' that will point to whatever waveform we want to process
TH1 *triggerH;


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
float calculateBaseline(const vector<int>& waveform)
{
    if(BASELINE_SAMPLES > waveform.size() ||
       BASELINE_LIMIT > waveform.size())
    {
        return BASELINE;
    }

    vector<int> baselineWindow;

    // Load the first BASELINE_SAMPLES waveform samples into the baseline window
    for(int i=0; i<BASELINE_SAMPLES; i++)
    {
        baselineWindow.push_back(waveform[i]);
    }

    // Loop through the waveform to try to find a BASELINE_SAMPLES-long stretch
    // that can be used to calculate the baseline. Truncate the search at
    // BASELINE_LIMIT samples.

    for(int i=BASELINE_SAMPLES; i<BASELINE_LIMIT; i++)
    {
        //BASELINE = baselineWithinThreshold(baselineWindow);

        // test to see if the baseline returned by baselineWithinThreshold is
        // non-zero.

        if(BASELINE==0)
        {
            // Baseline calculation failed, so move the baseline window forward
            // one step and try to calculate the baseline again.
            baselineWindow.erase(baselineWindow.begin());
            baselineWindow.push_back(waveform[i]);
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
bool isTrigger(int i, const vector<int>& waveform)
{
    // Check for triggers using:
    //      - a raw threshold above the baseline
    //      - a derivative threshold above the value DERIVATIVE_THRESHOLD
    //      - by rejecting triggers if the PREVIOUS point was already above
    //      these thresholds

    if((waveform[i] <= BASELINE-THRESHOLD
       /*&& (waveform[i]-waveform[i-1])/(double)SAMPLE_PERIOD <= DERIVATIVE_THRESHOLD*/)

       && (waveform[i-1] > BASELINE-THRESHOLD
       /*|| (waveform[i-1]-waveform[i-2])/(double)SAMPLE_PERIOD > DERIVATIVE_THRESHOLD*/))
    {
        return true;
    }

    else
    {
        return false;
    }
}

struct fitData
{
    double trigger1Time = 0;
    double trigger2Time = 0;
    double peak1Amplitude = 0;
    double peak2Amplitude = 0;
    double peak1TriggerAmplitude = 0;
    double peak2TriggerAmplitude = 0;
    double peak1Derivative = 0;
    double peak2Derivative = 0;
    double chiSquare = 10000;
    bool goodFit = false;

    void clear()
    {
        trigger1Time = 0;
        trigger2Time = 0;
        peak1Amplitude = 0;
        peak2Amplitude = 0;
        peak1TriggerAmplitude = 0;
        peak2TriggerAmplitude = 0;
        peak1Derivative = 0;
        peak2Derivative = 0;
        chiSquare = 10000;
        goodFit = false;
    }
} data;

bool testForEcho(vector<int>* waveform, int triggerSample)
{
    for(int i=triggerSample; i<triggerSample+PEAKFIT_WINDOW; i++)
    {
        if(waveform->at(i)>BASELINE)
        {
            return true;
        }
    }

    return false;
}

fitData fitTrigger(int waveformNo, double triggerSample, const vector<int>& waveform)
{
    // A trigger has been detected on the current waveform; fitTrigger attempts
    // to fit the region around this trigger using a series of progressively
    // more complicated fitting functions.
    // If any of these functions fits the waveform to within an error boundary,
    // fitTrigger accepts it as a good fit and exits.

    // reset fit data
    data.clear();

    // check to make sure we don't run off the end of the waveform
    if(triggerSample+PEAKFIT_WINDOW+PEAKFIT_OFFSET >= waveform.size() || triggerSample+PEAKFIT_OFFSET<0)
    {
        return data;
    }

    // extract the waveform chunk we'd like to fit
    stringstream temp;
    temp << "waveform" << waveformNo << "_peak" << triggerSample;

    //peakHisto = new TH1I(temp.str().c_str(),temp.str().c_str(),PEAKFIT_WINDOW,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW));
    peakHisto = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform.size(),0,waveform.size()*SAMPLE_PERIOD);

    for (int i=0; i<waveform.size(); i++)
    {
        peakHisto->SetBinContent(i,waveform[i]);
    }

    /*************************************************************************/
    // first, try fitting with one peak
    fittingFunc = new TF1("fittingFunc",onePeakForm,0,waveform.size()*SAMPLE_PERIOD,nParamsOnePeak);
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

    /*if(waveformNo<5000)
    {
        peakHisto->Write(); // uncomment to produce a fitted histo for each peak
    }*/

    if(fittingFunc->GetChisquare() < ERROR_LIMIT)
    {
        // Success - we've achieved a good fit with just one peak
        // Output fit data

        data.peak1Amplitude = fittingFunc->GetMinimum(fittingFunc->GetParameter(1),fittingFunc->GetParameter(1)+20);
        
        /*TF1* CFD = new TF1("CFD",CFDForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW),nParamsOnePeak);

        for(int i=0; i<fittingFunc->GetNpar(); i++)
        {
            CFD->SetParameter(i,fittingFunc->GetParameter(i));
        }*/

        //CFD->Write();

        //double relativeTriggerTime = CFD->GetX(fittingFunc->GetParameter(4),fittingFunc->GetParameter(1),fittingFunc->GetParameter(1)+10);

        //double relativeTriggerTime = fittingFunc->GetX(data.peak1Amplitude,fittingFunc->GetParameter(1),fittingFunc->GetParameter(1)+30);
        double relativeTriggerTime = fittingFunc->GetX((fittingFunc->GetParameter(4)+data.peak1Amplitude)/2,fittingFunc->GetParameter(1),fittingFunc->GetParameter(1)+20);
        //double relativeTriggerTime = fittingFunc->GetParameter(1);
        //double relativeTriggerTime = triggerSample*(double)SAMPLE_PERIOD;

        data.peak1TriggerAmplitude = fittingFunc->Eval(relativeTriggerTime);

        data.trigger1Time = relativeTriggerTime;

        data.peak1Derivative = fittingFunc->Derivative(relativeTriggerTime);
        data.chiSquare = fittingFunc->GetChisquare();
        data.goodFit = true;

        numberOnePeakFits++;
    }

    /*if(data.trigger1Time+DPP_PEAKFIT_OFFSET>4)
    {
        peakHisto->Write();
    }*/

    /*************************************************************************/

    /*else
    {
        fittingFunc = new TF1("fittingFunc",onePeakExpBackForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW),nParamsOnePeakExpBack);
        fittingFunc->SetParameters(A_init,trig1_init,n_init,d_init,C_init,m_init,E_init,tE_init);

        // Lock in the shape of the peak (manually chosen to match real peak shape)
        fittingFunc->FixParameter(2,n_init);
        fittingFunc->FixParameter(3,d_init);

        // Lock in a constant background
        fittingFunc->FixParameter(5,m_init);

        // Set parameter boundaries
        for (int i=0; i<fittingFunc->GetNpar(); i++)
        {
            fittingFunc->SetParLimits(i,onePeakExpBackParamMin[i],onePeakExpBackParamMax[i]);
        }

        // fit peak 
        peakHisto->Fit("fittingFunc","RQ");

        // see if peak fitting was successful
        if(fittingFunc->GetChisquare() < ERROR_LIMIT_2)
        {
            data.peak1Amplitude = fittingFunc->GetMinimum(fittingFunc->GetParameter(2),fittingFunc->GetParameter(2)+10);

            double relativeTriggerSample = fittingFunc->GetX((fittingFunc->GetParameter(6)+data.peak1Amplitude)/2,fittingFunc->GetParameter(2)-10,fittingFunc->GetParameter(2)+10);
            data.trigger1Time = SAMPLE_PERIOD*(triggerSample+relativeTriggerSample);
            data.peak1Derivative = fittingFunc->Derivative(relativeTriggerSample);
            relativeTriggerSampleHisto->Fill(relativeTriggerSample);

            data.chiSquare = fittingFunc->GetChisquare();
            data.goodFit = true;

            numberOnePeakExpBackFits++;
        }

        else
        {
            // failed to fit w/ one peak; try allowing two peaks
            fittingFunc = new TF1("fittingFunc",twoPeakForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW),nParamsTwoPeaks);
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

                numberTwoPeakFits++;
            }
        }
    }
*/

    delete peakHisto;
    delete fittingFunc;

    /*if(data.chiSquare<5)
      {
      delete peakHisto;
      }*/

    // failed to fit
    return data;
}

double calculateGammaOffset(const vector<double>& triggerList)
{
    double gammaOffset = 0;
    unsigned int numberOfGammas = 0;

    for(int i=0; (size_t)i<triggerList.size(); i++)
    {
        // Calculate time of flight from trigger time
        double rawMicroTime = fmod((triggerList[i]+WAVEFORM_OFFSET),MICRO_LENGTH);

        if(rawMicroTime>GAMMA_WINDOW[0] && rawMicroTime<GAMMA_WINDOW[1])
        {
            gammaOffset+=rawMicroTime;
            numberOfGammas++;
        }
    }

    if(numberOfGammas > 0)
    {
        gammaOffset/=numberOfGammas;
        gammaOffset -= pow(10,7)*FLIGHT_DISTANCE/C;
    }

    return gammaOffset;
}

//, int waveformNo
//triggerWalk->Fill(microTime,waveformNo);

void fillTriggerHistos(double triggerTime, vector<Plots>& plots)
{
    // Calculate time of flight from trigger time
    microTime = fmod((triggerTime),MICRO_LENGTH);

    microNo = floor((triggerTime)/MICRO_LENGTH);

    // convert microTime into neutron velocity based on flight path distance
    velocity = pow(10.,7.)*FLIGHT_DISTANCE/microTime; // in meters/sec 

    // convert velocity to relativistic kinetic energy
    rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

    if (procEvent.targetPos>0 && procEvent.targetPos<=tarGates.size()-1)
    {
        TH1D* tof = plots[procEvent.targetPos-1].getTOFHisto();
        TH1D* en = plots[procEvent.targetPos-1].getEnergyHisto();

        tof->Fill(microTime);
        en->Fill(rKE);
    }

    eventTimeDiff = triggerTime-prevTriggerTime;

    prevTriggerTime = triggerTime;
}

/*****************************************************************************/
void processTrigger(int waveformNo, int triggerSample, vector<double>& triggerList, const vector<int>& waveform)
{
    // Uncomment to use raw trigger sample as trigger time
    //double triggerTime = triggerSample*2;
    //triggerList.push_back(triggerTime);

    double timeOffset = 0;

    // Uncomment to use fitted peak threshold-intercept as trigger time
    if(fitTrigger(waveformNo, triggerSample, waveform).goodFit)
    {
        triggerList.push_back(data.trigger1Time);
        triggerValues.push_back(data.peak1TriggerAmplitude);

        /*if(data.peak2Amplitude && data.peak1Amplitude > 13000)
        {
            triggerList.push_back(data.trigger2Time);
        }*/

        numberGoodFits++;
    }

    else
    {
        //triggerList.push_back(triggerSample*SAMPLE_PERIOD-timeOffset);
        numberBadFits++;
    }

    //cout << data.chiSquare << endl;
}

void produceTriggerOverlay(int j, vector<double>& triggerList, vector<int>& waveform)
{
    stringstream temp;
    temp << "waveform " << j;
    waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform.size(),0,SAMPLE_PERIOD*waveform.size());

    for(int k=0; (size_t)k<waveform.size(); k++)
    {
        waveformH->SetBinContent(k,waveform.at(k));
    }

    waveformH->Write();

    temp << "triggers";
    triggerH = new TH1I(temp.str().c_str(),temp.str().c_str(),waveform.size(),0,SAMPLE_PERIOD*(waveform.size()));

    for(int k=0; (size_t)k<triggerList.size(); k++)
    {
        triggerH->SetBinContent((triggerList[k]/*+DPP_PEAKFIT_OFFSET*/)/(double)SAMPLE_PERIOD,triggerValues[k]);
    }

    TCanvas *c1 = new TCanvas;
    c1->DrawFrame(0,0,waveform.size()+10,16383);

    triggerH->SetOption("P");

    triggerH->SetMarkerStyle(29);
    triggerH->SetMarkerSize(3);
    triggerH->SetMarkerColor(2);

    triggerH->Write();
}

void processWaveforms(TTree* treeToSort, vector<Plots>& targetPlots, TFile* outFile, string mode)
{
    TH1D* fittedTimeHisto = new TH1D("raw fitted times", "raw fitted times", TOF_BINS,TOF_LOWER_BOUND,TOF_RANGE);

    TH2D* deltaTVsPulseHeight = new TH2D("delta T vs. pulse height","delta T vs. pulse height", 200, -10, 10, 1580, 0, 15800);

    TH2D* deltaTVsPulseIntegral0 = new TH2D("delta T vs. pulse integral, target 0","delta T vs. pulse integral, target 0",100, -10, 10, pow(2,15), 0, pow(2,15));
    TH2D* deltaTVsPulseIntegral1 = new TH2D("delta T vs. pulse integral, target 1","delta T vs. pulse integral, target 1",100, -10, 10, pow(2,15), 0, pow(2,15));
    TH2D* deltaTVsPulseIntegral2 = new TH2D("delta T vs. pulse integral, target 2","delta T vs. pulse integral, target 2",100, -10, 10, pow(2,15), 0, pow(2,15));
    TH2D* deltaTVsPulseIntegral3 = new TH2D("delta T vs. pulse integral, target 3","delta T vs. pulse integral, target 3",100, -10, 10, pow(2,15), 0, pow(2,15));
    TH2D* deltaTVsPulseIntegral4 = new TH2D("delta T vs. pulse integral, target 4","delta T vs. pulse integral, target 4",100, -10, 10, pow(2,15), 0, pow(2,15));
    TH2D* deltaTVsPulseIntegral5 = new TH2D("delta T vs. pulse integral, target 5","delta T vs. pulse integral, target 5",100, -10, 10, pow(2,15), 0, pow(2,15));

    TH1D* triggerAmplitudeHisto = new TH1D("triggerAmplitudeHisto","triggerAmplitudeHisto",pow(2,14),0,pow(2,14));
    TH1D* relativeTriggerTimeHisto = new TH1D("relativeTriggerTimeHisto","relative trigger time, from start of fitted wavelet",200,-5,5);
    TH2D* relativeTriggerTimeVsAmplitude = new TH2D("relativeTriggerTimeVSAmplitude","relative trigger time vs amplitude", 200, -5, 5, pow(2,14), 0, pow(2,14));

    TH1D* gammaToGammaTimeH = new TH1D("gammaToGammaTimeH","time between consecutive gammas", 1000, -5, 5);

    if(mode=="DPP")
    {
        setBranchesHistos(treeToSort);

        int totalEntries = treeToSort->GetEntries();
        cout << "Total waveforms = " << totalEntries << endl;

        vector<double> triggerList;

        waveformWrap = new TMultiGraph("DPP waveforms", "DPP waveforms");
        vector<TGraph*> waveletGraphs;
        vector<TGraph*> triggerGraphs;

        int gammaGate[2] = {80,90};

        double prevGammaTime = 0;

        for(int j=1; j<totalEntries; j++)
        {
            if(j%1000==0)
            {
                cout << "Processing triggers on waveform " << j << "\r";
                fflush(stdout);
            }

            /*if(j>500)
            {
                break;
            }*/

            triggerList.clear();
            triggerValues.clear();

            // pull individual waveform event
            treeToSort->GetEntry(j);

            // calculate the baseline for this waveform
            BASELINE = calculateBaseline(*procEvent.waveform);

            // Loop through all points in the waveform and fit peaks
            for(int k=DPP_PEAKFIT_START; (size_t)k<procEvent.waveform->size(); k++)
            {
                // Check to see if this point creates a new trigger
                if(isTrigger(k, *procEvent.waveform))
                {
                    // trigger found - plot/fit/extract time

                    double timeDiff = procEvent.completeTime-procEvent.macroTime;
                    double microTime = fmod(timeDiff,MICRO_LENGTH);

                    if(microTime > gammaGate[0] && microTime < gammaGate[1])
                    {
                        processTrigger(j, k, triggerList, *procEvent.waveform);

                        double fullTime = procEvent.completeTime-procEvent.macroTime+data.trigger1Time+DPP_PEAKFIT_OFFSET;
                        gammaToGammaTimeH->Fill(fmod(fullTime,MICRO_LENGTH)-prevGammaTime);
                        prevGammaTime=fmod(fullTime,MICRO_LENGTH);

                        fillTriggerHistos(fullTime, targetPlots);
                        fittedTimeHisto->Fill(data.trigger1Time);
                    }

                    //processTrigger(j, k, triggerList, *procEvent.waveform);

                    break;
                }
            }

            //produceTriggerOverlay(j, triggerList, *procEvent.waveform);

            TH2D* deltaTVsPulseIntegralHisto;

            switch(procEvent.targetPos-1)
            {
                case 0:
                    deltaTVsPulseIntegralHisto = deltaTVsPulseIntegral0;
                    break;

                case 1:
                    deltaTVsPulseIntegralHisto = deltaTVsPulseIntegral1;
                    break;

                case 2:
                    deltaTVsPulseIntegralHisto = deltaTVsPulseIntegral2;
                    break;

                case 3:
                    deltaTVsPulseIntegralHisto = deltaTVsPulseIntegral3;
                    break;

                case 4:
                    deltaTVsPulseIntegralHisto = deltaTVsPulseIntegral4;
                    break;

                case 5:
                    deltaTVsPulseIntegralHisto = deltaTVsPulseIntegral5;
                    break;
            }

            triggerAmplitudeHisto->Fill(data.peak1Amplitude);
            relativeTriggerTimeHisto->Fill(data.trigger1Time+DPP_PEAKFIT_OFFSET);
            relativeTriggerTimeVsAmplitude->Fill(data.trigger1Time+DPP_PEAKFIT_OFFSET,data.peak1Amplitude);

            deltaTVsPulseIntegralHisto->Fill(data.trigger1Time+DPP_PEAKFIT_OFFSET, procEvent.lgQ);
            deltaTVsPulseHeight->Fill(data.trigger1Time+DPP_PEAKFIT_OFFSET, data.peak1Amplitude);

            // Create a new graph for each wavelet
            //waveletGraphs.push_back(new TGraph());

            // Fill each micropulse graph with waveform samples
            /*for (int l=0; l<procEvent.waveform->size(); l++)
              {
              waveletGraphs.back()->SetPoint(l,l*SAMPLE_PERIOD,procEvent.waveform->at(l));
              }*/

            // Create a new graph for each wavelet
            //triggerGraphs.push_back(new TGraph());

            // Fill each micropulse graph with waveform samples
            //triggerGraphs.back()->SetPoint(0,triggerList[0],triggerValues[0]);

            fill(procEvent.waveform->begin(),procEvent.waveform->end(),BASELINE);
        }

        TGraph* exponentialFit = new TGraph();
        exponentialFit->SetPoint(0,-0.78,140);
        exponentialFit->SetPoint(1,0.16,206);
        exponentialFit->SetPoint(2,1.79,305);
        exponentialFit->SetPoint(3,2.82,417);
        exponentialFit->SetPoint(4,3.59,566);
        exponentialFit->SetPoint(5,4.4,929);
        exponentialFit->SetPoint(6,5.2,1482);
        exponentialFit->SetPoint(7,5.7,2149);
        exponentialFit->SetPoint(8,6.9,5319);
        exponentialFit->SetPoint(9,7.3,7808);
        exponentialFit->SetPoint(10,7.7,11395);
        exponentialFit->SetPoint(11,8.0,16200);
        exponentialFit->Write();

        // Add each wavelet graph to the MultiGraph
        for (int m=0; m<waveletGraphs.size(); m++)
        {
            //cout << "adding graph " << m << " to multigraph" << endl;
            waveletGraphs[m]->Draw();
            waveformWrap->Add(waveletGraphs[m],"l");
        }

        // Add each trigger graph to the MultiGraph
        for (int m=0; m<triggerGraphs.size(); m++)
        {
            //cout << "adding graph " << m << " to multigraph" << endl;
            triggerGraphs[m]->SetMarkerSize(2);
            triggerGraphs[m]->SetMarkerColor(2);
            triggerGraphs[m]->Draw();
            waveformWrap->Add(triggerGraphs[m],"*");
        }

        waveformWrap->Write();

        for(Plots p : targetPlots)
        {
            p.getTOFHisto()->Write();
            p.getEnergyHisto()->Write();
        }

        fittedTimeHisto->Write();

        deltaTVsPulseIntegral0->Write();
        deltaTVsPulseIntegral1->Write();
        deltaTVsPulseIntegral2->Write();
        deltaTVsPulseIntegral3->Write();
        deltaTVsPulseIntegral4->Write();
        deltaTVsPulseIntegral5->Write();

        deltaTVsPulseHeight->Write();

        relativeTriggerTimeHisto->Write();
        triggerAmplitudeHisto->Write();
        relativeTriggerTimeVsAmplitude->Write();

        gammaToGammaTimeH->Write();
    }

    else if(mode=="waveform")
    {
        setBranchesHistosW(treeToSort);

        triggerWalk = new TH2I("triggerWalk","trigger time vs. waveform chunk #",200,0,200,1000,0,1000);

        TH1I* monitorHisto = new TH1I("targetPosH", "targetPos", 7, 0, 7);
        monitorHisto->GetXaxis()->SetTitle("target position of each waveform");

        long triggersWithGamma = 0;
        vector<int> fullWaveform(MACRO_LENGTH/2, BASELINE);
        vector<double> triggerList;

        int prevTargetPos = 0;
        double firstTimetagInSeries = 0;

        int totalEntries = treeToSort->GetEntries();
        cout << "Total waveforms = " << totalEntries << endl;

        /*TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan");

          if(!mycan)
          {
          mycan = new TCanvas("mycan","mycan");
          }*/

        // EVENT LOOP for sorting through channel-specific waveforms
        for(int j=1; j<totalEntries; j++)
        {
            cout << "Processing triggers on waveform " << j << "\r";
            fflush(stdout);

            if(j>30)
            {
                break;
            }

            triggerList.clear();
            triggerValues.clear();

            prevTargetPos = procEvent.targetPos;

            // pull individual waveform event
            treeToSort->GetEntry(j);

            if(procEvent.evtNo==1)
            {
                // new macropulse; process the previous macropulse's waveform

                // calculate the baseline for this waveform
                BASELINE = calculateBaseline(fullWaveform);

                // Loop through all points in the waveform and fit peaks
                for(int k=PEAKFIT_WINDOW; (size_t)k<fullWaveform.size(); k++)
                {
                    // Check to see if this point creates a new trigger
                    if(isTrigger(k, fullWaveform))
                    {
                        // trigger found - plot/fit/extract time
                        processTrigger(j, k, triggerList, fullWaveform);

                        // shift waveform index past the end of this fitting window
                        // so that we don't refit the same data
                        k += PEAKFIT_WINDOW;
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

                // Use the gamma peaks to find the time offset for this waveform, and
                // adjust the microTime with this offset

                double gammaOffset = calculateGammaOffset(triggerList);

                /*for(int m=0; (size_t)m<triggerList.size(); m++)
                  {
                  fillTriggerHistos(triggerList[m]-gammaOffset, targetPlots);
                  }*/

                if(gammaOffset!=0)
                {
                    for(int m=0; (size_t)m<triggerList.size(); m++)
                    {
                        //fillTriggerHistos(triggerList[m]-gammaOffset+WAVEFORM_OFFSET, targetPlots);
                    }

                    triggersWithGamma++;
                }

                /*if(j<30)
                  {
                  produceTriggerOverlay(j, triggerList, fullWaveform);
                  }*/

                /*temp.str("");
                  temp << "waveformWrap" << j;

                  waveformWrap = new TMultiGraph(temp.str().c_str(), temp.str().c_str());

                  vector<TGraph*> microGraphs;

                // Create a new graph for each micropulse period to be plotted
                for (int m = 0; m<floor(2*procEvent.waveform->size()/MICRO_LENGTH); m++)
                {
                microGraphs.push_back(new TGraph());
                }

                // Fill each micropulse graph with waveform samples
                for (int l = 0; l<procEvent.waveform->size(); l++)
                {
                microGraphs[(int)floor(l/(double)MICRO_LENGTH)]->SetPoint(microGraphs[(int)floor(l/(double)MICRO_LENGTH)]->GetN(),fmod(2*l+WAVEFORM_OFFSET,MICRO_LENGTH),procEvent.waveform->at(l));
                //cout << "Adding value " << procEvent.waveform->at(l) << " to position " << fmod(l,MICRO_LENGTH) << " in microGraph " << floor(l/(double)MICRO_LENGTH) << endl;
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

                //triggerH->Write();
                //cout << "Finished processing waveform " << j << endl << endl;

                /*if(j==10)
                  {
                  break;
                  }*/

                monitorHisto->Fill(prevTargetPos);

                fill(fullWaveform.begin(),fullWaveform.end(),BASELINE);

                firstTimetagInSeries = procEvent.completeTime;
            }

            if (procEvent.targetPos == 0)
            {
                continue;
            }

            double timeOffset = procEvent.completeTime-firstTimetagInSeries;

            if(timeOffset<MACRO_LENGTH-2*procEvent.waveform->size())
            {
                for(int k=0; k<procEvent.waveform->size(); k++)
                {
                    fullWaveform[timeOffset/SAMPLE_PERIOD+k] = procEvent.waveform->at(k);
                }
            }
        }

        cout << "Triggers with gamma in wavelet: " << triggersWithGamma << endl;

        monitorHisto->Write();

        for(Plots p : targetPlots)
        {
            p.getTOFHisto()->Write();
            p.getEnergyHisto()->Write();
        }

        cout << "Number of good fits: " << numberGoodFits << endl;
        cout << "onePeak = " << numberOnePeakFits << endl; 
        cout << "onePeakExpBack = " << numberOnePeakExpBackFits << endl; 
        cout << "twoPeaks = " << numberTwoPeakFits << endl << endl; 
        cout << "Number of bad fits: " << numberBadFits << endl;
    }

    else
    {
        cerr << "Error: digitizer mode was " << mode << "; only possible options are 'DPP' or 'waveform'. Exiting..." << endl;
        exit(1);
    }
}

void waveform(string inFileName, string outFileName, vector<string>channelMap, string mode)
{
    TFile* inFile = new TFile(inFileName.c_str(),"READ");
    if(!inFile->IsOpen())
    {
        cerr << "Error: failed to open sorted.root" << endl;
        exit(1);
    }

    TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

    for(int i=0; i<channelMap.size(); i++)
    {
        if(channelMap[i]!="lowThresholdDet"/* && channelMap[i]!="monitor"*/)
        {
            continue;
        }

        string treeName;

        if(mode=="DPP")
        {
            treeName = channelMap[i];
        }

        else if(mode=="waveform")
        {
            treeName = channelMap[i]+"W";
        }

        else
        {
            cerr << "Error: digitizer mode was " << mode << "; only possible options are 'DPP' or 'waveform'. Exiting..." << endl;
            exit(1);
        }

        TTree* treeToSort = (TTree*)inFile->Get(treeName.c_str());
        
        if(!treeToSort)
        {
            cerr << "Error: couldn't find channel " << i << " tree when attempting to assign macropulses. Exiting... " << endl;
            exit(1);
        }

        outFile->cd("/");
        outFile->mkdir(channelMap[i].c_str(),channelMap[i].c_str());
        outFile->GetDirectory(channelMap[i].c_str())->cd();

        vector<Plots> targetPlots;

        // Extract triggers from waveforms
        for(unsigned int i=0; i<tarGates.size()-1; i++)
        {
            string name = POSITION_NAMES[i];
            if(name=="")
            {
                continue;
            }

            targetPlots.push_back(Plots(name.c_str()));
        }

        processWaveforms(treeToSort, targetPlots, outFile, mode);

        /*int totalEntries = ch4TreeWaveform->GetEntries();

        // total number of micropulses processed per target (for performing dead time
        // calculation)
        vector<long> microsPerTargetWaveform(6,0);

        for(int i=0; i<totalEntries; i++)
        {
        ch4TreeWaveform->GetEntry(i);

        if(procEvent.targetPos==0)
        {
        continue;
        }
        microsPerTargetWaveform[procEvent.targetPos-1] += SAMPLE_PERIOD*procEvent.waveform->size()/(double)MICRO_LENGTH;
        }

        // perform a manual dead-time correction
        //calculateDeadtime(microsPerTargetWaveform,plots);

        // Calculate cross-sections from waveforms' trigger time data
        //calculateCS(targets, waveformFile);

        int totalMicros = 0;
        for(int i=0; (size_t)i<microsPerTargetWaveform.size(); i++)
        {
        totalMicros += microsPerTargetWaveform[i];
        }

        for(Plots* p : plots)
        {
        numberTotalTriggers += p->getTOFHisto()->GetEntries();
        }

        cout << "Total micros on ch. : " << totalMicros << endl;
        cout << "Total number of triggers on ch. : " << numberTotalTriggers << endl;
        cout << "Triggers/micropulse: " << numberTotalTriggers/(double)totalMicros << endl;
        */
    }

    inFile->Close();
    outFile->Close();
    return;
}
