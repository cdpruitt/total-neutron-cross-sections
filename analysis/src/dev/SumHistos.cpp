#include "TH1.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

// sub-run number of last run to sum over
int limit = 24;

int main(int argc, char *argv[])
{
    ifstream runList("runsToSort.txt");
    if(!runList.is_open())
    {
        cout << "Could not find runs list... exiting." << endl;
        exit(1);
    }

    string runDir = argv[1];
    string outpath = argv[2];

    stringstream fileInName;
    stringstream fileInWaveformName;
    stringstream fileOutName;

    stringstream treeName;
    treeName << "run" << runDir << "-0000"; 

    fileInName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_cross-sections.root";
    fileInWaveformName << outpath <<"/analysis/run" << runDir << "/" << treeName.str() << "_waveform.root";

    // Examine the first sub-run of runDir and find out the correct number of
    // bins to use for the summed histograms
    TFile* infile =  new TFile(fileInName.str().c_str(),"READ");
    TH1D* csBins = ((TH1D*)infile->Get("blankcs"));
    TH1D* relativeBins = ((TH1D*)infile->Get("relativeSnCS"));

    TFile* infileWaveform =  new TFile(fileInWaveformName.str().c_str(),"READ");
    TH1D* waveformBins = ((TH1D*)infileWaveform->Get("blankcs"));

    // each histo has N bins + 1 overflow + 1 underflow
    // thus, subtract two to get the number of 'normal' bins for summed histos
    int noCSBins = csBins->GetSize()-2;
    int noRelativeBins = relativeBins->GetSize()-2;

    int noWaveformBins = waveformBins->GetSize()-2;

    delete csBins;
    delete relativeBins;

    delete waveformBins;

    // Create output file to contain summed histos
    fileOutName << outpath << "/analysis/run" << runDir << "/" << "sum.root";

    TFile *outfile = new TFile(fileOutName.str().c_str(),"RECREATE");

    // Create summed histos
    // First, sum DPP-derived cross-sections
    TH1D *blankcsSum = new TH1D("blankcsSum","blankcsSum",noCSBins,0,700);
    TH1D *target1csSum = new TH1D("target1csSum","target1csSum",noCSBins,0,700);
    TH1D *target2csSum = new TH1D("target2csSum","target2csSum",noCSBins,0,700);
    TH1D *target3csSum = new TH1D("target3csSum","target3csSum",noCSBins,0,700);
    TH1D *target4csSum = new TH1D("target4csSum","target4csSum",noCSBins,0,700);
    TH1D *target5csSum = new TH1D("target5csSum","target5csSum",noCSBins,0,700);

    TH1D *blankcsSumLog = new TH1D("blankcsSumLog","blankcsSumLog",noCSBins,0,TMath::Log10(700));
    TH1D *target1csSumLog = new TH1D("target1csSumLog","target1csSumLog",noCSBins,0,TMath::Log10(700));
    TH1D *target2csSumLog = new TH1D("target2csSumLog","target2csSumLog",noCSBins,0,TMath::Log10(700));
    TH1D *target3csSumLog = new TH1D("target3csSumLog","target3csSumLog",noCSBins,0,TMath::Log10(700));
    TH1D *target4csSumLog = new TH1D("target4csSumLog","target4csSumLog",noCSBins,0,TMath::Log10(700));
    TH1D *target5csSumLog = new TH1D("target5csSumLog","target5csSumLog",noCSBins,0,TMath::Log10(700));

    // Next, sum waveform-derived cross-sections
    TH1D *blankcsSumWaveform = new TH1D("blankcsSumWaveform","blankcsSumWaveform",noWaveformBins,0,700);
    TH1D *target1csSumWaveform = new TH1D("target1csSumWaveform","target1csSumWaveform",noWaveformBins,0,700);
    TH1D *target2csSumWaveform = new TH1D("target2csSumWaveform","target2csSumWaveform",noWaveformBins,0,700);
    TH1D *target3csSumWaveform = new TH1D("target3csSumWaveform","target3csSumWaveform",noWaveformBins,0,700);
    TH1D *target4csSumWaveform = new TH1D("target4csSumWaveform","target4csSumWaveform",noWaveformBins,0,700);
    TH1D *target5csSumWaveform = new TH1D("target5csSumWaveform","target5csSumWaveform",noWaveformBins,0,700);

    TH1D *blankcsSumWaveformLog = new TH1D("blankcsSumWaveformLog","blankcsSumWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *target1csSumWaveformLog = new TH1D("target1csSumWaveformLog","target1csSumWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *target2csSumWaveformLog = new TH1D("target2csSumWaveformLog","target2csSumWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *target3csSumWaveformLog = new TH1D("target3csSumWaveformLog","target3csSumWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *target4csSumWaveformLog = new TH1D("target4csSumWaveformLog","target4csSumWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *target5csSumWaveformLog = new TH1D("target5csSumWaveformLog","target5csSumWaveformLog",noWaveformBins,0,TMath::Log10(700));

    // Next, sum DPP-derived relative cross-sections, where Sn112/Sn124
    TH1D *relativeSnCSSum = new TH1D("relativeSnCSSum","relativeSnCSSum",noRelativeBins,0,700);
    TH1D *relativeSnCSSumLog = new TH1D("relativeSnCSSumLog","relativeSnCSSumLog",noRelativeBins,0,TMath::Log10(700));

    // Loop through all sub-runs and add their histos together
    int runNo;

    while(!runList.eof())
    {
        runList >> runNo;

        for(int segment=0; segment<=limit; segment++)
        {
            if(segment < 10)
            {
                infile =  new TFile(Form("%s/analysis/run%i/run%i-000%i_cross-sections.root",outpath.c_str(),runNo,runNo,segment));
                infileWaveform =  new TFile(Form("%s/analysis/run%i/run%i-000%i_waveform.root",outpath.c_str(),runNo,runNo,segment));
            }

            else if(segment < 100)
            {
                infile =  new TFile(Form("%s/analysis/run%i/run%i-00%i_cross-sections.root",outpath.c_str(),runNo,runNo,segment));
                infileWaveform =  new TFile(Form("%s/analysis/run%i/run%i-00%i_waveform.root",outpath.c_str(),runNo,runNo,segment));
            }
            else
            {
                cout << "Segment number too large!" << endl;
                return 1;
            }

            if(!infile->IsOpen())
            {
                cout << "Can't open root file run " << runNo << " segment = " << segment << endl;
                outfile->Write();
                outfile->Close();
                exit(1);
            }

            cout << "Adding run " << runNo << " segment " << segment <<endl;

            // sum cross-section histos
            TH1D * blank = ((TH1D*)infile->Get("blankcs"));
            if (blank)
            {
                blankcsSum->Add(blank,1/(segment+1));
            }

            TH1D * target1 = (TH1D*)infile->Get("target1cs");
            if (target1)
            {
                target1csSum->Add(target1,1/(double)(limit+1));
            }

            TH1D * target2 = (TH1D*)infile->Get("target2cs");
            if (target2)
            {
                target2csSum->Add(target2,1/(double)(limit+1));
            }

            TH1D * target3 = (TH1D*)infile->Get("target3cs");
            if (target3)
            {
                target3csSum->Add(target3,1/(double)(limit+1));
            }

            TH1D * target4 = (TH1D*)infile->Get("target4cs");
            if (target4)
            {
                target4csSum->Add(target4,1/(double)(limit+1));
            }

            TH1D * target5 = (TH1D*)infile->Get("target5cs");
            if (target5)
            {
                target5csSum->Add(target5,1/(double)(limit+1));
            }

            // sum log-scale cross-section histos
            TH1D * blankLog = ((TH1D*)infile->Get("blankcsLog"));
            if (blankLog)
            {
                blankcsSumLog->Add(blankLog,1/(double)(limit+1));
            }

            TH1D * target1Log = (TH1D*)infile->Get("target1csLog");
            if (target1Log)
            {
                target1csSumLog->Add(target1Log,1/(double)(limit+1));
            }

            TH1D * target2Log = (TH1D*)infile->Get("target2csLog");
            if (target2Log)
            {
                target2csSumLog->Add(target2Log,1/(double)(limit+1));
            }

            TH1D * target3Log = (TH1D*)infile->Get("target3csLog");
            if (target3Log)
            {
                target3csSumLog->Add(target3Log,1/(double)(limit+1));
            }

            TH1D * target4Log = (TH1D*)infile->Get("target4csLog");
            if (target4Log)
            {
                target4csSumLog->Add(target4Log,1/(double)(limit+1));
            }

            TH1D * target5Log = (TH1D*)infile->Get("target5csLog");
            if (target5Log)
            {
                target5csSumLog->Add(target5Log,1/(double)(limit+1));
            }

            // Add relative cross-section plots
            TH1D * relativeSnCS = (TH1D*)infile->Get("relativeSnCS");
            if (relativeSnCS)
            {
                relativeSnCSSum->Add(relativeSnCS,1/(double)(limit+1));
            }

            TH1D * relativeSnCSLog = (TH1D*)infile->Get("relativeSnCSLog");
            if (relativeSnCSLog)
            {
                relativeSnCSSumLog->Add(relativeSnCSLog,1/(double)(limit+1));
            }

            // sum cross-section histos (from waveform fitting)
            TH1D * blankW = ((TH1D*)infileWaveform->Get("blankcs"));
            if (blankW)
            {
                blankcsSumWaveform->Add(blankW,1/(segment+1));
            }

            TH1D * target1W = (TH1D*)infileWaveform->Get("target1cs");
            if (target1W)
            {
                target1csSumWaveform->Add(target1W,1/(double)(limit+1));
            }

            TH1D * target2W = (TH1D*)infileWaveform->Get("target2cs");
            if (target2W)
            {
                target2csSumWaveform->Add(target2W,1/(double)(limit+1));
            }

            TH1D * target3W = (TH1D*)infileWaveform->Get("target3cs");
            if (target3W)
            {
                target3csSumWaveform->Add(target3W,1/(double)(limit+1));
            }

            TH1D * target4W = (TH1D*)infileWaveform->Get("target4cs");
            if (target4W)
            {
                target4csSumWaveform->Add(target4W,1/(double)(limit+1));
            }

            TH1D * target5W = (TH1D*)infileWaveform->Get("target5cs");
            if (target5W)
            {
                target5csSumWaveform->Add(target5W,1/(double)(limit+1));
            }

            // sum log-scale cross-section histos
            TH1D * blankLogW = ((TH1D*)infileWaveform->Get("blankcsLog"));
            if (blankLogW)
            {
                blankcsSumWaveformLog->Add(blankLogW,1/(double)(limit+1));
            }

            TH1D * target1LogW = (TH1D*)infileWaveform->Get("target1csLog");
            if (target1LogW)
            {
                target1csSumWaveformLog->Add(target1LogW,1/(double)(limit+1));
            }

            TH1D * target2LogW = (TH1D*)infileWaveform->Get("target2csLog");
            if (target2LogW)
            {
                target2csSumWaveformLog->Add(target2LogW,1/(double)(limit+1));
            }

            TH1D * target3LogW = (TH1D*)infileWaveform->Get("target3csLog");
            if (target3LogW)
            {
                target3csSumWaveformLog->Add(target3LogW,1/(double)(limit+1));
            }

            TH1D * target4LogW = (TH1D*)infileWaveform->Get("target4csLog");
            if (target4LogW)
            {
                target4csSumWaveformLog->Add(target4LogW,1/(double)(limit+1));
            }

            TH1D * target5LogW = (TH1D*)infileWaveform->Get("target5csLog");
            if (target5LogW)
            {
                target5csSumWaveformLog->Add(target5LogW,1/(double)(limit+1));
            }

            // Close cross-section.root and waveform.root of subrun
            infile->Close();
            infileWaveform->Close();

        }

        totalfile

        // Write run histograms to sum.root
        outfile->Write();

        outfile->Close();
    }
}
