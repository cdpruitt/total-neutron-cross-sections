#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TDirectory.h"

#include "../include/physicalConstants.h"
#include "../include/correctForDeadtime.h"
#include "../include/config.h"

using namespace std;

extern Config config;

// logistic curve for deadtime response
// (calculated by fitting "time difference between events" histogram using a
// logistic curve)

// input x in ns
double logisticDeadtimeFunction(double x)
{
    if(x<config.deadtime.LOGISTIC_MU-15)
    {
        return 1;
    }

    return (1-1/(1+exp(-config.deadtime.LOGISTIC_K*(x-config.deadtime.LOGISTIC_MU))));
}

int generateDeadtimeCorrection(TH1D*& TOFtoCorrect, TH1D*& deadtimeHisto, const int& numberOfPeriods)
{
    int DEADTIME_BINS = (config.deadtime.LOGISTIC_MU+15)*config.plot.TOF_BINS_PER_NS;

    string targetName = TOFtoCorrect->GetName();
    int numberOfBins = TOFtoCorrect->GetNbinsX();

    vector<double> measuredRatePerBin(numberOfBins);

    for(int i=1; i<=numberOfBins; i++)
    {
        measuredRatePerBin[i-1] = TOFtoCorrect->GetBinContent(i)/(double)numberOfPeriods;
    }

    vector<long double> deadtimePerBin(numberOfBins);

    for(int i=0; i<measuredRatePerBin.size(); i++)
    {
        double cumulativeDeadtime = 0;
        for(int j=0; j<DEADTIME_BINS; j++)
        {
            if((i-j)<0)
            {
                deadtimePerBin[i] += measuredRatePerBin[(i-j)+numberOfBins]*logisticDeadtimeFunction(((double)j)/config.plot.TOF_BINS_PER_NS);
            }

            else
            {
                deadtimePerBin[i] += measuredRatePerBin[i-j]*logisticDeadtimeFunction(((double)j)/config.plot.TOF_BINS_PER_NS);
            }
        }

        deadtimeHisto->SetBinContent(i+1,deadtimePerBin[i]);
    }

    return 0;
}

int applyDeadtimeCorrection(string inputFileName, string deadtimeFileName, string macroFileName, string gammaCorrectionFileName, ofstream& logFile, string outputFileName)
{
    // test if output file already exists
    ifstream f(outputFileName);
    if(f.good())
    {
        cout << outputFileName << " already exists; skipping deadtime generation." << endl;
        logFile << outputFileName << " already exists; skipping deadtime generation." << endl;
        return 2;
    }
    f.close();

    cout << "Applying deadtime correction..." << endl;

    // open input file
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    // open deadtime file
    TFile* deadtimeFile = new TFile(deadtimeFileName.c_str(),"READ");
    if(!deadtimeFile->IsOpen())
    {
        cerr << "Error: failed to open " << deadtimeFileName << "  to fill histos." << endl;
        return 1;
    }

    // open macropulse file
    TFile* macroFile = new TFile(macroFileName.c_str(),"READ");
    if(!macroFile->IsOpen())
    {
        cerr << "Error: failed to open " << macroFileName << "  to fill histos." << endl;
        return 1;
    }

    // open gamma correction file
    TFile* gammaCorrectionFile = new TFile(gammaCorrectionFileName.c_str(),"READ");
    if(!gammaCorrectionFile->IsOpen())
    {
        cerr << "Error: failed to open " << gammaCorrectionFileName << "  to read gamma correction." << endl;
        return 1;
    }

    TDirectory* gammaDirectory = (TDirectory*)gammaCorrectionFile->Get(config.analysis.GAMMA_CORRECTION_TREE_NAME.c_str());
    if(!gammaDirectory)
    {
        cerr << "Error: failed to open " << config.analysis.GAMMA_CORRECTION_TREE_NAME << " directory in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;
        return 1;
    }

    gammaDirectory->cd();

    TH1D* gammaCorrectionHisto = (TH1D*)gammaDirectory->Get("gammaCorrection");
    if(!gammaCorrectionHisto)
    {
        cerr << "Error: failed to open gammaCorrections histo in " << gammaCorrectionFileName << " for reading gamma corrections." << endl;
        return 1;
    }

    double overallAverageGammaTime = 0;

    int gammaCorrectionBins = gammaCorrectionHisto->GetNbinsX();
    if(gammaCorrectionBins<=0)
    {
        cerr << "Cannot divide by 0 to calculate overall average gamma time (in applyDeadtimeCorrection)." << endl;
        return 1;
    }

    for(int i=1; i<=gammaCorrectionBins; i++)
    {
        overallAverageGammaTime += gammaCorrectionHisto->GetBinContent(i);
    }

    overallAverageGammaTime /= gammaCorrectionBins;

    int deadtimeBinOffset = floor(config.plot.TOF_BINS_PER_NS*(overallAverageGammaTime));

    // create outputFile
    TFile* outputFile = new TFile(outputFileName.c_str(),"CREATE");

    for(auto& channelName : config.cs.DETECTOR_NAMES)
    {
        TDirectory* detectorDirectory = inputFile->GetDirectory(channelName.c_str());
        if(!detectorDirectory)
        {
            cerr << "Error: failed to find " << channelName << " directory in " << inputFileName << "." << endl;
            return 1;
        }

        TDirectory* deadtimeDirectory = deadtimeFile->GetDirectory(channelName.c_str());
        if(!deadtimeDirectory)
        {
            cerr << "Error: failed to find " << channelName << " directory in " << deadtimeFileName << "." << endl;
            return 1;
        }

        TDirectory* macroDirectory = macroFile->GetDirectory(channelName.c_str());
        if(!macroDirectory)
        {
            cerr << "Error: failed to find " << channelName << " directory in " << macroFileName << "." << endl;
            return 1;
        }

        TDirectory* directory = outputFile->mkdir(channelName.c_str(),channelName.c_str());
        directory->cd();

        // find TOF histograms in input file
        for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
        {
            string targetName = config.target.TARGET_ORDER[i];

            string TOFHistoName = targetName + "TOF";
            TH1D* TOF = (TH1D*)detectorDirectory->Get(TOFHistoName.c_str());
            if(!TOF)
            {
                cerr << "Error: failed to find " << TOFHistoName << " in "
                    << channelName << " of " << inputFileName << "." << endl;

                outputFile->Close();
                macroFile->Close();
                inputFile->Close();
                return 1;
            }

            string deadtimeHistoName = targetName + "Deadtime";
            TH1D* deadtimeHisto = (TH1D*)deadtimeDirectory->Get(deadtimeHistoName.c_str());
            if(!deadtimeHisto)
            {
                cerr << "Error: failed to find " << deadtimeHistoName << " in "
                    << channelName << " of " << deadtimeFileName << "." << endl;

                outputFile->Close();
                macroFile->Close();
                inputFile->Close();
                return 1;
            }

            string macroHistoName = targetName + "MacroNumber";
            TH1I* macroHisto = (TH1I*)macroDirectory->Get(macroHistoName.c_str());
            if(!macroHisto)
            {
                cerr << "Error: failed to find " << macroHistoName << " in "
                    << channelName << " of " << macroFileName << "." << endl;

                outputFile->Close();
                macroFile->Close();
                inputFile->Close();
                return 1;
            }

            string correctedTOFName = TOFHistoName;
            TH1D* correctedTOF = (TH1D*)TOF->Clone(correctedTOFName.c_str());

            int numberOfMacroBins = macroHisto->GetNbinsX();
            int numberOfMacros = 0;

            for(int j=1; j<=numberOfMacroBins; j++)
            {
                if(macroHisto->GetBinContent(j)>=1)
                {
                    numberOfMacros++;
                }
            }

            double numberOfMicros = numberOfMacros
                *(config.facility.LAST_GOOD_MICRO-config.facility.FIRST_GOOD_MICRO);

            if(numberOfMicros <=0)
            {
                cerr << "Error: cannot apply deadtime for <= 0 periods." << endl;

                outputFile->Close();
                macroFile->Close();
                inputFile->Close();
                return 1;
            }

            int numberOfBins = TOF->GetNbinsX();

            for(int i=1; i<=numberOfBins; i++)
            {
                double originalBin = TOF->GetBinContent(i);

                double deadtime;

                if(i+deadtimeBinOffset > numberOfBins)
                {
                    deadtime = deadtimeHisto->GetBinContent(i+deadtimeBinOffset-numberOfBins);
                }

                else
                {
                    deadtime = deadtimeHisto->GetBinContent(i+deadtimeBinOffset);
                }

                correctedTOF->SetBinContent(i, -log(1-(originalBin/numberOfMicros)/(1-deadtime))*numberOfMicros);
            }

            correctedTOF->Write();
        }
    }

    outputFile->Close();
    macroFile->Close();
    inputFile->Close();

    return 0;
}

int generateDeadtimeCorrection(string inputFileName, ofstream& logFile, string outputFileName)
{
    // test if output file already exists
    ifstream f(outputFileName);
    if(f.good())
    {
        cout << outputFileName << " already exists; skipping deadtime generation." << endl;
        logFile << outputFileName << " already exists; skipping deadtime generation." << endl;
        return 2;
    }

    // open input file
    TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
    if(!inputFile->IsOpen())
    {
        cerr << "Error: failed to open " << inputFileName << "  to fill histos." << endl;
        return 1;
    }

    // create outputFile
    TFile* outputFile = new TFile(outputFileName.c_str(),"CREATE");

    for(auto& channelName : config.cs.DETECTOR_NAMES)
    {
        TDirectory* detectorDirectory = inputFile->GetDirectory(channelName.c_str());
        if(!detectorDirectory)
        {
            cerr << "Error: failed to find " << channelName << " directory in " << inputFileName << "." << endl;
            return 1;
        }

        TDirectory* directory = outputFile->mkdir(channelName.c_str(),channelName.c_str());

        directory->cd();

        // find TOF histograms in input file
        for(int i=0; i<config.target.TARGET_ORDER.size(); i++)
        {
            string targetName = config.target.TARGET_ORDER[i];

            string TOFHistoName = targetName + "TOF";
            TH1D* TOF = (TH1D*)detectorDirectory->Get(TOFHistoName.c_str());

            string macroHistoName = targetName + "MacroNumber";
            TH1I* macroHisto = (TH1I*)detectorDirectory->Get(macroHistoName.c_str());
            int numberOfBins = macroHisto->GetNbinsX();
            int numberOfMacros = 0;

            for(int j=1; j<=numberOfBins; j++)
            {
                if(macroHisto->GetBinContent(j)>=1)
                {
                    numberOfMacros++;
                }
            }

            double numberOfMicros = numberOfMacros
                *(config.facility.LAST_GOOD_MICRO-config.facility.FIRST_GOOD_MICRO);

            string deadtimeHistoName = targetName + "Deadtime";
            TH1D* deadtimeHisto = (TH1D*)TOF->Clone(deadtimeHistoName.c_str());
            generateDeadtimeCorrection(TOF, deadtimeHisto, numberOfMicros);

            deadtimeHisto->Write();
        }
    }

    outputFile->Close();
    inputFile->Close();

    return 0;
}
