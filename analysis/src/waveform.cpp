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
#include "../include/analysisConstants.h"
#include "../include/plottingConstants.h"
#include "../include/targetConstants.h"
#include "../include/dataStructures.h"
#include "../include/target.h"
#include "../include/waveformFitting.h"
#include "../include/waveform.h"
#include "../include/plots.h"
#include "../include/branches.h"

using namespace std;

TH1I *relativeTriggerSampleHisto;
TH2I *triggerWalk;
TH1I *peakHisto;
TF1 *fittingFunc;

extern struct ProcessedEvent procEvent;

// Declare variables to be used for calculating neutron TOFs
int microNo;
double microTime, velocity, rKE;

// Keep track of the number of waveforms collected during each target's period
// in the beam
int targetCounts[6] = {0};

/*****************************************************************************/

// Keep track of sample number where triggers are found
vector<double> triggerList;
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
float calculateBaseline()
{
    vector<int> baselineWindow;

    // Load the first BASELINE_SAMPLES waveform samples into the baseline window
    for(int i=0; i<BASELINE_SAMPLES; i++)
    {
        baselineWindow.push_back(procEvent.waveform->at(i));
    }

    // Loop through the waveform to try to find a BASELINE_SAMPLES-long stretch
    // that can be used to calculate the baseline. Truncate the search at
    // BASELINE_LIMIT samples.

    for(int i=BASELINE_SAMPLES; i<BASELINE_LIMIT; i++)
    {
        BASELINE = baselineWithinThreshold(baselineWindow);

        // test to see if the baseline returned by baselineWithinThreshold is
        // non-zero.

        if(BASELINE==0)
        {
            // Baseline calculation failed, so move the baseline window forward
            // one step and try to calculate the baseline again.
            baselineWindow.erase(baselineWindow.begin());
            baselineWindow.push_back(procEvent.waveform->at(i));
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

    if((procEvent.waveform->at(i) <= BASELINE-THRESHOLD
       && (procEvent.waveform->at(i)-procEvent.waveform->at(i-1))/(double)SAMPLE_PERIOD <= DERIVATIVE_THRESHOLD)

       && (procEvent.waveform->at(i-1) > BASELINE-THRESHOLD
       || (procEvent.waveform->at(i-1)-procEvent.waveform->at(i-2))/(double)SAMPLE_PERIOD > DERIVATIVE_THRESHOLD))
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

fitData fitTrigger(int waveformNo, float triggerSample)
{
    // A trigger has been detected on the current waveform; fitTrigger attempts
    // to fit the region around this trigger using a series of progressively
    // more complicated fitting functions.
    // If any of these functions fits the waveform to within an error boundary,
    // fitTrigger accepts it as a good fit and exits.

    // reset fit data
    data.clear();

    // check to make sure we don't run off the end of the waveform
    if(triggerSample+PEAKFIT_WINDOW >= procEvent.waveform->size())
    {
        return data;
    }

    // extract the waveform chunk we'd like to fit
    stringstream temp;
    temp << "waveform" << waveformNo << "_peak" << triggerSample;

    peakHisto = new TH1I(temp.str().c_str(),temp.str().c_str(),PEAKFIT_WINDOW,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW));

    for (int i=0; i<PEAKFIT_WINDOW; i++)
    {
        peakHisto->SetBinContent(i,procEvent.waveform->at(triggerSample+PEAKFIT_OFFSET+i));
    }

    /*************************************************************************/
    // first, try fitting with one peak
    fittingFunc = new TF1("fittingFunc",onePeakForm,SAMPLE_PERIOD*PEAKFIT_OFFSET,SAMPLE_PERIOD*(PEAKFIT_OFFSET+PEAKFIT_WINDOW),nParamsOnePeak);
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
    /*************************************************************************/

    else
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
        /*************************************************************************/
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
            /*************************************************************************/
        }
    }

    //peakHisto->Write(); // uncomment to produce a fitted histo for each peak
    delete peakHisto;
    delete fittingFunc;

    /*if(data.chiSquare<5)
      {
      delete peakHisto;
      }*/

    // failed to fit
    return data;
}

double calculateGammaOffset()
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
    }

    gammaOffset-=pow(10,7)*FLIGHT_DISTANCE/C;
    //cout << "gammaOffset = " << gammaOffset << endl;

    return gammaOffset;
}

void fillTriggerHistos(double triggerTime, double gammaOffset, int waveformNo, vector<Plots*>& plots)
{
    // Calculate time of flight from trigger time
    microTime = fmod((triggerTime+WAVEFORM_OFFSET),MICRO_LENGTH)-gammaOffset;

    microNo = floor((triggerTime+WAVEFORM_OFFSET)/MICRO_LENGTH);

    // convert microTime into neutron velocity based on flight path distance
    velocity = pow(10.,7.)*FLIGHT_DISTANCE/microTime; // in meters/sec 

    // convert velocity to relativistic kinetic energy
    rKE = (pow((1.-pow((velocity/C),2.)),-0.5)-1.)*NEUTRON_MASS; // in MeV

    triggerWalk->Fill(microTime,waveformNo);

    if (procEvent.targetPos>0 && procEvent.targetPos<=NUMBER_OF_TARGETS)
    {
        TH1I* tof = plots[procEvent.targetPos-1]->getTOFHisto();
        TH1I* en = plots[procEvent.targetPos-1]->getEnergyHisto();

        tof->Fill(microTime);
        en->Fill(rKE);
    }
}

/*****************************************************************************/
void processTrigger(int waveformNo, float triggerSample)
{
    // Uncomment to use raw trigger sample as trigger time
    //float triggerTime = triggerSample;
    //fillTriggerHistos(triggerTime);
    //triggerList.push_back(triggerSample);

    // Uncomment to use fitted peak threshold-intercept as trigger time
    if(fitTrigger(waveformNo, triggerSample).goodFit)
    {
        triggerList.push_back(data.trigger1Time);
        triggerValues.push_back(procEvent.waveform->at(data.trigger1Time/2));

        if(data.peak2Amplitude && data.peak1Amplitude > 13000)
        {
            triggerList.push_back(data.trigger2Time);
            triggerValues.push_back(procEvent.waveform->at(data.trigger2Time/2));
        }

        numberGoodFits++;
    }

    else
    {
        numberBadFits++;
        triggerList.push_back(triggerSample*SAMPLE_PERIOD);
    }

}

void produceTriggerOverlay(int j)
{
    stringstream temp;
    temp << "waveform " << j;
    waveformH = new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,SAMPLE_PERIOD*procEvent.waveform->size());

    for(int k=0; (size_t)k<procEvent.waveform->size(); k++)
    {
        waveformH->SetBinContent(k,procEvent.waveform->at(k));
    }

    waveformH->Write();

    temp << "triggers";
    triggerH = new TH1I(temp.str().c_str(),temp.str().c_str(),procEvent.waveform->size(),0,SAMPLE_PERIOD*(procEvent.waveform->size()));

    for(int k=0; (size_t)k<triggerList.size(); k++)
    {
        triggerH->SetBinContent(triggerList[k]/SAMPLE_PERIOD,triggerValues[k]);
    }

    TCanvas *c1 = new TCanvas;
    c1->DrawFrame(0,0,procEvent.waveform->size()+10,16383);

    triggerH->SetOption("P");

    triggerH->SetMarkerStyle(29);
    triggerH->SetMarkerSize(3);
    triggerH->SetMarkerColor(2);

    triggerH->Write();
}

void processWaveforms(TTree* ch4TreeWaveform, vector<Plots*>& plots)
{
    triggerWalk = new TH2I("triggerWalk","trigger time vs. waveform chunk #",200,0,200,1000,0,1000);

    relativeTriggerSampleHisto = new TH1I("relativeTriggerSampleHisto","relative trigger time, from start of fitted wavelet",100,PEAKFIT_OFFSET*SAMPLE_PERIOD,(PEAKFIT_OFFSET+PEAKFIT_WINDOW)*SAMPLE_PERIOD);

    // Loop through all channel-specific trees
    for(int i=0; (size_t)i<1; i++)
    {
        // point event variables at the correct tree in preparation for reading
        // data
        setBranchesHistosW(ch4TreeWaveform);

        int totalEntries = ch4TreeWaveform->GetEntries();
        cout << "Total waveforms = " << totalEntries << endl;

        TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan");

        if(!mycan)
        {
            mycan = new TCanvas("mycan","mycan");
        }

        // EVENT LOOP for sorting through channel-specific waveforms
        for(int j=0; j<totalEntries; j++)
        {

            cout << "Processing triggers on waveform " << j << "\r";
            fflush(stdout);

            triggerList.clear();
            triggerValues.clear();

            // pull individual waveform event
            ch4TreeWaveform->GetEntry(j);

            if (procEvent.targetPos == 0)
            {
                continue;
            }

            //cout << "waveform chunk time = " << procEvent.completeTime << endl;

            // calculate the baseline for this waveform
            BASELINE = calculateBaseline();

            // add event to target-position tracker
            if (procEvent.targetPos > 0)
            {
                targetCounts[procEvent.targetPos-1]++;
            }

            // Loop through all points in the waveform and fit peaks
            for(int k=PEAKFIT_WINDOW; (size_t)k<procEvent.waveform->size(); k++)
            {
                // Check to see if this point creates a new trigger
                if(isTrigger(k))
                {
                    // trigger found - plot/fit/extract time
                    processTrigger(j, k);

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
            double gammaOffset = calculateGammaOffset();

            if(gammaOffset!=0)
            {
                for(int m=0; (size_t)m<triggerList.size(); m++)
                {
                    fillTriggerHistos(triggerList[m], gammaOffset, j, plots);
                }
            }

            //produceTriggerOverlay(j);

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

            /*if(j==4)
            {
                break;
            }*/
        }

        for(Plots* p : plots)
        {
            p->getTOFHisto()->Write();
            p->getEnergyHisto()->Write();
        }
    }

    for(int k=0; (size_t)k<NUMBER_OF_TARGETS; k++)
    {
        cout << "target position " << k+1 << " counts = " << targetCounts[k] << endl;
    }
    cout << endl;
}

void calculateDeadtime(TTree* ch4TreeWaveform, vector<Plots*>& plots)
{
    setBranchesHistosW(ch4TreeWaveform);

    int totalEntries = ch4TreeWaveform->GetEntries();

    for(int i=0; i<totalEntries; i++)
    {
        ch4TreeWaveform->GetEntry(i);

        if(procEvent.targetPos==0)
        {
            continue;
        }
        microsPerTargetWaveform[procEvent.targetPos-1] += 2*procEvent.waveform->size()/(double)MICRO_LENGTH;
    }
    
    vector<vector<double>> eventsPerBinPerMicro(NUMBER_OF_TARGETS,vector<double>(0));

    // "deadtimeFraction" records the fraction of time that the detector is dead, for
    // neutrons of a certain energy. It is target-dependent.
    vector<vector<double>> deadtimeFraction(NUMBER_OF_TARGETS, vector<double>(0));

    // Calculate deadtime:
    // for each target,
    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        if(microsPerTargetWaveform[i] <= 0)
        {
            break;
        }

        // for each bin,
        TH1I* tof = plots[i]->getTOFHisto();
        TH1I* dtH = plots[i]->getDeadtimeHisto();

        for(int j=0; j<TOF_BINS; j++)
        {
            eventsPerBinPerMicro[i].push_back(tof->GetBinContent(j+1)/(double)microsPerTargetWaveform[i]);
            deadtimeFraction[i].push_back(0);
        }

        // find the fraction of the time that the detector is dead for each bin in the micropulse
        // set deadtime fraction base case

        int deadtimeBins = (TOF_BINS/TOF_RANGE)*DEADTIME_PERIOD;

        // use deadtime base case to calculate deadtime for remaining bins
        for(int j=0; (size_t)j<TOF_BINS; j++)
        {
            for(int k=j-deadtimeBins; k<j; k++)
            {
                if(k<0)
                {
                    deadtimeFraction[i][j] += eventsPerBinPerMicro[i][k+TOF_BINS]*(1-deadtimeFraction[i][j]);
                    continue;
                }

                deadtimeFraction[i][j] += eventsPerBinPerMicro[i][k]*(1-deadtimeFraction[i][j]);
            }
        }

        double averageDeadtime = 0;

        for(int j=0; (size_t)j<deadtimeFraction[i].size(); j++)
        {
            dtH->SetBinContent(j,1000000*deadtimeFraction[i][j]);
            averageDeadtime += deadtimeFraction[i][j];
        }
        dtH->Write();

        averageDeadtime /= deadtimeFraction[i].size();
        cout << "Average deadtime for target " << i << ": " << averageDeadtime << endl;
    }
}

void waveform(string inFileName, string outFileName)
{
    TFile* inFile = new TFile(inFileName.c_str(),"READ");
    if(!inFile->IsOpen())
    {
        cerr << "Error: failed to open resort.root" << endl;
        exit(1);
    }

    TTree* ch4TreeWaveform = (TTree*)inFile->Get("ch4ProcessedTreeW");

    TFile* outFile;
    outFile = new TFile(outFileName.c_str(),"UPDATE");

    vector<Plots*> plots;

    // if plots are empty
    if(!(TH1I*)outFile->Get("blankWEnergy"))
    {
        // Extract triggers from waveforms

        for(int i=0; i<NUMBER_OF_TARGETS; i++)
        {
            string name = positionNames[i] + "W";
            plots.push_back(new Plots(name.c_str()));
        }

        processWaveforms(ch4TreeWaveform, plots);

        cout << "Number of good fits: " << numberGoodFits << endl;
        cout << "onePeak = " << numberOnePeakFits << endl; 
        cout << "onePeakExpBack = " << numberOnePeakExpBackFits << endl; 
        cout << "twoPeaks = " << numberTwoPeakFits << endl << endl; 
        cout << "Number of bad fits: " << numberBadFits << endl;

     //   outFile->Write();
    }

    else
    {
        for(int i=0; i<NUMBER_OF_TARGETS; i++)
        {
            string name = positionNames[i] + "W";
            plots.push_back(new Plots(name, outFile));
        }
    }

    // perform a manual dead-time correction
    calculateDeadtime(ch4TreeWaveform,plots);

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
 
    inFile->Close();
    outFile->Close();
}
