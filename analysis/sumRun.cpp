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

const double CS_LOWER_BOUND = 1; // cross-section plots' lower bound, in MeV
const double CS_UPPER_BOUND = 700; // cross-section plots' upper bound, in MeV

int noCSBins = 50;

TH1D* logBins(TH1D *inputHisto)
{
    string newName;
    newName = inputHisto->GetName();
    newName += "Log";

    double newXMin = (((TAxis*)inputHisto->GetXaxis())->GetXmin());
    if (newXMin <= 0)
    {
        cout << "Error: can't take log of negative energy on cross-section plot" << endl;
        exit(1);
    }

    newXMin = TMath::Log10(newXMin);

    TH1D* outputHisto = new TH1D(newName.c_str(),newName.c_str(),noCSBins,
            TMath::Log10(((TAxis*)inputHisto->GetXaxis())->GetXmin()),
            TMath::Log10(((TAxis*)inputHisto->GetXaxis())->GetXmax()));

    // Pull bin data from input histo, and map to the log scale:
    TAxis* axis = outputHisto->GetXaxis();
    int nBins = axis->GetNbins();

    double xMin = axis->GetXmin();
    double xMax = axis->GetXmax();

    double binWidth = (xMax-xMin)/nBins;
    double *newBins = new double[nBins+1];

    for(int i=0; i<=nBins; i++)
    {
        newBins[i] = TMath::Power(10, xMin+i*binWidth);
    }

    // Assign the log-scale bins to the new histo
    ((TAxis*)outputHisto->GetXaxis())->Set(nBins,newBins);
    delete newBins;

    return outputHisto;
}

int main(int argc, char *argv[])
{
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
    TH1D* csBins = ((TH1D*)infile->Get("carbonSCS"));

    TFile* infileWaveform =  new TFile(fileInWaveformName.str().c_str(),"READ");
    TH1D* waveformBins = ((TH1D*)infileWaveform->Get("carbonSCS"));

    // each histo has N bins + 1 overflow + 1 underflow
    // thus, subtract two to get the number of 'normal' bins for summed histos
    noCSBins = csBins->GetSize()-2;

    int noWaveformBins = waveformBins->GetSize()-2;

    delete csBins;
    delete waveformBins;

    // Create output file to contain summed histos
    fileOutName << outpath << "/analysis/run" << runDir << "/" << "sum.root";

    TFile *outfile = new TFile(fileOutName.str().c_str(),"RECREATE");

    // Create summed histos
    // First, sum DPP-derived cross-sections
    TH1D *blankCSSum = new TH1D("blankCSSum","blankCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonSCSSum = new TH1D("carbonSCSSum","carbonSCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonLCSSum = new TH1D("carbonLCSSum","carbonLCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn112CSSum = new TH1D("Sn112CSSum","Sn112CSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *SnNatCSSum = new TH1D("SnNatCSSum","SnNatCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn124CSSum = new TH1D("Sn124CSSum","Sn124CSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    TH1D *blankCSSumLog = logBins(blankCSSum);
    TH1D *carbonSCSSumLog = logBins(carbonSCSSum);
    TH1D *carbonLCSSumLog = logBins(carbonLCSSum);
    TH1D *Sn112CSSumLog = logBins(Sn112CSSum);
    TH1D *SnNatCSSumLog = logBins(SnNatCSSum);
    TH1D *Sn124CSSumLog = logBins(Sn124CSSum);

    // Next, sum waveform-derived cross-sections
    TH1D *blankCSSumWaveform = new TH1D("blankCSSumWaveform","blankCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonSCSSumWaveform = new TH1D("carbonSCSSumWaveform","carbonSCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonLCSSumWaveform = new TH1D("carbonLCSSumWaveform","carbonLCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn112CSSumWaveform = new TH1D("Sn112CSSumWaveform","Sn112CSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *SnNatCSSumWaveform = new TH1D("SnNatCSSumWaveform","SnNatCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn124CSSumWaveform = new TH1D("Sn124CSSumWaveform","Sn124CSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    TH1D *blankCSSumWaveformLog = new TH1D("blankCSSumWaveformLog","blankCSSumWaveformLog",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonSCSSumWaveformLog = new TH1D("carbonSCSSumWaveformLog","carbonSCSSumWaveformLog",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *carbonLCSSumWaveformLog = new TH1D("carbonLCSSumWaveformLog","carbonLCSSumWaveformLog",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn112CSSumWaveformLog = new TH1D("Sn112CSSumWaveformLog","Sn112CSSumWaveformLog",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *SnNatCSSumWaveformLog = new TH1D("SnNatCSSumWaveformLog","SnNatCSSumWaveformLog",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TH1D *Sn124CSSumWaveformLog = new TH1D("Sn124CSSumWaveformLog","Sn124CSSumWaveformLog",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    vector<TH1D*> allHistos;

    allHistos.push_back(blankCSSum);
    allHistos.push_back(carbonSCSSum);
    allHistos.push_back(carbonLCSSum);
    allHistos.push_back(Sn112CSSum);
    allHistos.push_back(SnNatCSSum);
    allHistos.push_back(Sn124CSSum);

    allHistos.push_back(blankCSSumLog);
    allHistos.push_back(carbonSCSSumLog);
    allHistos.push_back(carbonLCSSumLog);
    allHistos.push_back(Sn112CSSumLog);
    allHistos.push_back(SnNatCSSumLog);
    allHistos.push_back(Sn124CSSumLog);

    allHistos.push_back(blankCSSumWaveform);
    allHistos.push_back(carbonSCSSumWaveform);
    allHistos.push_back(carbonLCSSumWaveform);
    allHistos.push_back(Sn112CSSumWaveform);
    allHistos.push_back(SnNatCSSumWaveform);
    allHistos.push_back(Sn124CSSumWaveform);

    allHistos.push_back(blankCSSumWaveformLog);
    allHistos.push_back(carbonSCSSumWaveformLog);
    allHistos.push_back(carbonLCSSumWaveformLog);
    allHistos.push_back(Sn112CSSumWaveformLog);
    allHistos.push_back(SnNatCSSumWaveformLog);
    allHistos.push_back(Sn124CSSumWaveformLog);

    // Loop through all sub-runs and add their histos together

    // sub-run number of last run to sum over
    int limit = 25;
    for(int segment = 0; segment<=limit; segment++)
    {
        // We need to form the proper name for the sub-run we want to open:
        if(segment < 10)
        {
            infile =  new TFile(Form("%s/analysis/run%s/run%s-000%i_cross-sections.root",outpath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            infileWaveform =  new TFile(Form("%s/analysis/run%s/run%s-000%i_waveform.root",outpath.c_str(),runDir.c_str(),runDir.c_str(),segment));
        }

        else if(segment < 100)
        {
            infile =  new TFile(Form("%s/analysis/run%s/run%s-00%i_cross-sections.root",outpath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            infileWaveform =  new TFile(Form("%s/analysis/run%s/run%s-00%i_waveform.root",outpath.c_str(),runDir.c_str(),runDir.c_str(),segment));
        }

        else
        {
            // There's some error in the sub-run file number; write outfile and exit
            cout << "Segment number too large!" << endl;
            outfile->Write();
            outfile->Close();
            exit(1);
        }

        // Attempt to open the sub-run to access its histograms
        if(!infile->IsOpen())
        {
            cout << "Can't open root file run " << runDir << " segment = " << segment << endl;

            if (segment>0)
            {

                /*TH1D* blankCSSumRaw = blankCSSum->Clone("blankCSSumRaw");
                TH1D* Sn112CSSumRaw = blankCSSum->Clone("Sn112CsSumRaw");
                TH1D* Sn124CSSumRaw = blankCSSum->Clone("Sn124CsSumRaw");

                TH1D* blankCSSumRaw = blankCSSum->Clone("blankCSSumRaw");
                TH1D* Sn112CSSumRaw = blankCSSum->Clone("Sn112CsSumRaw");
                TH1D* Sn124CSSumRaw = blankCSSum->Clone("Sn124CsSumRaw");
                */

                for (int i=0; (size_t)i<allHistos.size(); i++)
                {
                    allHistos[i]->Scale(1/(double)(segment));
                }
            }

            outfile->Write();
            outfile->Close();
            exit(0);
        }

        cout << "Adding run " << runDir << " segment " << segment <<endl;


        /****************** SUM HISTOGRAMS ******************/

        /************** normal cross-section histos ************/
        TH1D * blank = ((TH1D*)infile->Get("blankCS"));
        if (blank)
        {
            blankCSSum->Add(blank);
        }

        TH1D * carbonS = (TH1D*)infile->Get("carbonSCS");
        if (carbonS)
        {
            carbonSCSSum->Add(carbonS);
        }

        TH1D * carbonL = (TH1D*)infile->Get("carbonLCS");
        if (carbonL)
        {
            carbonLCSSum->Add(carbonL);
        }

        TH1D * Sn112 = (TH1D*)infile->Get("Sn112CS");
        if (Sn112)
        {
            Sn112CSSum->Add(Sn112);
        }

        TH1D * SnNat = (TH1D*)infile->Get("SnNatCS");
        if (SnNat)
        {
            SnNatCSSum->Add(SnNat);
        }

        TH1D * Sn124 = (TH1D*)infile->Get("Sn124CS");
        if (Sn124)
        {
            Sn124CSSum->Add(Sn124);
        }

        /************* LOG-SCALE cross-section histos **************/
        TH1D * blankLog = ((TH1D*)infile->Get("blankCSLog"));
        if (blankLog)
        {
            blankCSSumLog->Add(blankLog);
        }

        TH1D * carbonSLog = (TH1D*)infile->Get("carbonSCSLog");
        if (carbonSLog)
        {
            carbonSCSSumLog->Add(carbonSLog);
        }

        TH1D * carbonLLog = (TH1D*)infile->Get("carbonLCSLog");
        if (carbonLLog)
        {
            carbonLCSSumLog->Add(carbonLLog);
        }

        TH1D * Sn112Log = (TH1D*)infile->Get("Sn112CSLog");
        if (Sn112Log)
        {
            Sn112CSSumLog->Add(Sn112Log);
        }

        TH1D * SnNatLog = (TH1D*)infile->Get("SnNatCSLog");
        if (SnNatLog)
        {
            SnNatCSSumLog->Add(SnNatLog);
        }

        TH1D * Sn124Log = (TH1D*)infile->Get("Sn124CSLog");
        if (Sn124Log)
        {
            Sn124CSSumLog->Add(Sn124Log);
        }

        /************* WAVEFORM-DERIVED cross-section histos **************/
        TH1D * blankW = ((TH1D*)infileWaveform->Get("blankCS"));
        if (blankW)
        {
            blankCSSumWaveform->Add(blankW);
        }

        TH1D * carbonSW = (TH1D*)infileWaveform->Get("carbonSCS");
        if (carbonSW)
        {
            carbonSCSSumWaveform->Add(carbonSW);
        }

        TH1D * carbonLW = (TH1D*)infileWaveform->Get("carbonLCS");
        if (carbonLW)
        {
            carbonLCSSumWaveform->Add(carbonLW);
        }

        TH1D * Sn112W = (TH1D*)infileWaveform->Get("Sn112CS");
        if (Sn112W)
        {
            Sn112CSSumWaveform->Add(Sn112W);
        }

        TH1D * SnNatW = (TH1D*)infileWaveform->Get("SnNatCS");
        if (SnNatW)
        {
            SnNatCSSumWaveform->Add(SnNatW);
        }

        TH1D * Sn124W = (TH1D*)infileWaveform->Get("Sn124CS");
        if (Sn124W)
        {
            Sn124CSSumWaveform->Add(Sn124W);
        }

        /********** LOG-SCALE WAVEFORM-DERIVED cross-section histos ***********/
        TH1D * blankLogW = ((TH1D*)infileWaveform->Get("blankCSLog"));
        if (blankLogW)
        {
            blankCSSumWaveformLog->Add(blankLogW);
        }

        TH1D * carbonSLogW = (TH1D*)infileWaveform->Get("carbonSCSLog");
        if (carbonSLogW)
        {
            carbonSCSSumWaveformLog->Add(carbonSLogW);
        }

        TH1D * carbonLLogW = (TH1D*)infileWaveform->Get("carbonLCSLog");
        if (carbonLLogW)
        {
            carbonLCSSumWaveformLog->Add(carbonLLogW);
        }

        TH1D * Sn112LogW = (TH1D*)infileWaveform->Get("Sn112CSLog");
        if (Sn112LogW)
        {
            Sn112CSSumWaveformLog->Add(Sn112LogW);
        }

        TH1D * SnNatLogW = (TH1D*)infileWaveform->Get("SnNatCSLog");
        if (SnNatLogW)
        {
            SnNatCSSumWaveformLog->Add(SnNatLogW);
        }

        TH1D * Sn124LogW = (TH1D*)infileWaveform->Get("Sn124CSLog");
        if (Sn124LogW)
        {
            Sn124CSSumWaveformLog->Add(Sn124LogW);
        }

        // Close the sub-run input files
        infile->Close();
        infileWaveform->Close();

        // End of loop - move to next sub-run
    }

    // Write run histograms to sum.root
    outfile->Write();
    outfile->Close();
}
