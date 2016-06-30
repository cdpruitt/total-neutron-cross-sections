#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

/*TGraph* logBins(TGraph *inputHisto)
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

    int noBins = ((TAxis*)inputHisto->GetXaxis())->GetNbins();

    TGraph* outputHisto = new TGraph(newName.c_str(),newName.c_str(),noBins,
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
}*/

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

    TFile* infileWaveform =  new TFile(fileInWaveformName.str().c_str(),"READ");
    //int noWaveformBins = ((TGraph*)infileWaveform->Get("carbonSCS"))->GetSize()-2;

    vector<vector<double>*> SnNatCrossSections;
    vector<vector<double>*> SnNatEnergies;

    // Create summed histos
    // First, sum DPP-derived cross-sections
    /*TGraph *blankCSSum = new TGraph("blankCSSum","blankCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *carbonSCSSum = new TGraph("carbonSCSSum","carbonSCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *carbonLCSSum = new TGraph("carbonLCSSum","carbonLCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *Sn112CSSum = new TGraph("Sn112CSSum","Sn112CSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *SnNatCSSum = new TGraph("SnNatCSSum","SnNatCSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *Sn124CSSum = new TGraph("Sn124CSSum","Sn124CSSum",noCSBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    TGraph *blankCSSumLog = logBins(blankCSSum);
    TGraph *carbonSCSSumLog = logBins(carbonSCSSum);
    TGraph *carbonLCSSumLog = logBins(carbonLCSSum);
    TGraph *Sn112CSSumLog = logBins(Sn112CSSum);
    TGraph *SnNatCSSumLog = logBins(SnNatCSSum);
    TGraph *Sn124CSSumLog = logBins(Sn124CSSum);

    // Next, sum waveform-derived cross-sections
    TGraph *blankCSSumWaveform = new TGraph("blankCSSumWaveform","blankCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *carbonSCSSumWaveform = new TGraph("carbonSCSSumWaveform","carbonSCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *carbonLCSSumWaveform = new TGraph("carbonLCSSumWaveform","carbonLCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *Sn112CSSumWaveform = new TGraph("Sn112CSSumWaveform","Sn112CSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *SnNatCSSumWaveform = new TGraph("SnNatCSSumWaveform","SnNatCSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);
    TGraph *Sn124CSSumWaveform = new TGraph("Sn124CSSumWaveform","Sn124CSSumWaveform",noWaveformBins,CS_LOWER_BOUND,CS_UPPER_BOUND);

    TGraph *blankCSSumWaveformLog = logBins(blankCSSumWaveform);
    TGraph *carbonSCSSumWaveformLog = logBins(carbonSCSSumWaveform);
    TGraph *carbonLCSSumWaveformLog = logBins(carbonLCSSumWaveform);
    TGraph *Sn112CSSumWaveformLog = logBins(Sn112CSSumWaveform);
    TGraph *SnNatCSSumWaveformLog = logBins(SnNatCSSumWaveform);
    TGraph *Sn124CSSumWaveformLog = logBins(Sn124CSSumWaveform);

    vector<TGraph*> allHistos;

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
    */

    // Loop through all sub-runs and add their histos together

    // sub-run number of last run to sum over
    int limit = 24;
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
            break;
            exit(1);
        }

        // Attempt to open the sub-run to access its histograms
        if(!infile->IsOpen())
        {
            cout << "Can't open root file run " << runDir << " segment = " << segment << endl;
            continue;
        }

        // Successfully found the sub-run of interest
        cout << "Adding run " << runDir << " segment " << segment <<endl;

        TGraph * SnNat = (TGraph*)infile->Get("NatSn");
        if(SnNat)
        {
            int numPoints = SnNat->GetN();
            SnNatCrossSections.push_back(new vector<double>);
            SnNatEnergies.push_back(new vector<double>);

            for(int i=0; i<numPoints; i++)
            {
                SnNatCrossSections.back()->push_back(0);
                SnNatEnergies.back()->push_back(0);
                SnNat->GetPoint(i,SnNatEnergies[segment]->at(i),SnNatCrossSections[segment]->at(i));
            }
        }

        //SnNat->Write();

        /****************** SUM HISTOGRAMS ******************/

        /************** normal cross-section histos ************/
        /*TGraph * blank = ((TGraph*)infile->Get("blankCS"));
        if (blank)
        {
            blankCSSum->Add(blank);
        }

        TGraph * carbonS = (TGraph*)infile->Get("carbonSCS");
        if (carbonS)
        {
            carbonSCSSum->Add(carbonS);
        }

        TGraph * carbonL = (TGraph*)infile->Get("carbonLCS");
        if (carbonL)
        {
            carbonLCSSum->Add(carbonL);
        }

        TGraph * Sn112 = (TGraph*)infile->Get("Sn112CS");
        if (Sn112)
        {
            Sn112CSSum->Add(Sn112);
        }

        TGraph * SnNat = (TGraph*)infile->Get("SnNatCS");
        if (SnNat)
        {
            SnNatCSSum->Add(SnNat);
        }

        TGraph * Sn124 = (TGraph*)infile->Get("Sn124CS");
        if (Sn124)
        {
            Sn124CSSum->Add(Sn124);
        }
        */

        /************* LOG-SCALE cross-section histos **************/
        /*TGraph * blankLog = ((TGraph*)infile->Get("blankCSLog"));
        if (blankLog)
        {
            blankCSSumLog->Add(blankLog);
        }

        TGraph * carbonSLog = (TGraph*)infile->Get("carbonSCSLog");
        if (carbonSLog)
        {
            carbonSCSSumLog->Add(carbonSLog);
        }

        TGraph * carbonLLog = (TGraph*)infile->Get("carbonLCSLog");
        if (carbonLLog)
        {
            carbonLCSSumLog->Add(carbonLLog);
        }

        TGraph * Sn112Log = (TGraph*)infile->Get("Sn112CSLog");
        if (Sn112Log)
        {
            Sn112CSSumLog->Add(Sn112Log);
        }

        TGraph * SnNatLog = (TGraph*)infile->Get("SnNatCSLog");
        if (SnNatLog)
        {
            SnNatCSSumLog->Add(SnNatLog);
        }

        TGraph * Sn124Log = (TGraph*)infile->Get("Sn124CSLog");
        if (Sn124Log)
        {
            Sn124CSSumLog->Add(Sn124Log);
        }
        */

        /************* WAVEFORM-DERIVED cross-section histos **************/
        /*TGraph * blankW = ((TGraph*)infileWaveform->Get("blankCS"));
        if (blankW)
        {
            blankCSSumWaveform->Add(blankW);
        }

        TGraph * carbonSW = (TGraph*)infileWaveform->Get("carbonSCS");
        if (carbonSW)
        {
            carbonSCSSumWaveform->Add(carbonSW);
        }

        TGraph * carbonLW = (TGraph*)infileWaveform->Get("carbonLCS");
        if (carbonLW)
        {
            carbonLCSSumWaveform->Add(carbonLW);
        }

        TGraph * Sn112W = (TGraph*)infileWaveform->Get("Sn112CS");
        if (Sn112W)
        {
            Sn112CSSumWaveform->Add(Sn112W);
        }

        TGraph * SnNatW = (TGraph*)infileWaveform->Get("SnNatCS");
        if (SnNatW)
        {
            SnNatCSSumWaveform->Add(SnNatW);
        }

        TGraph * Sn124W = (TGraph*)infileWaveform->Get("Sn124CS");
        if (Sn124W)
        {
            Sn124CSSumWaveform->Add(Sn124W);
        }
        */

        /********** LOG-SCALE WAVEFORM-DERIVED cross-section histos ***********/
        /*TGraph * blankLogW = ((TGraph*)infileWaveform->Get("blankCSLog"));

        if (blankLogW)
        {
            blankCSSumWaveformLog->Add(blankLogW);
        }

        TGraph * carbonSLogW = (TGraph*)infileWaveform->Get("carbonSCSLog");
        if (carbonSLogW)
        {
            carbonSCSSumWaveformLog->Add(carbonSLogW);
        }

        TGraph * carbonLLogW = (TGraph*)infileWaveform->Get("carbonLCSLog");
        if (carbonLLogW)
        {
            carbonLCSSumWaveformLog->Add(carbonLLogW);
        }

        TGraph * Sn112LogW = (TGraph*)infileWaveform->Get("Sn112CSLog");
        if (Sn112LogW)
        {
            Sn112CSSumWaveformLog->Add(Sn112LogW);
        }

        TGraph * SnNatLogW = (TGraph*)infileWaveform->Get("SnNatCSLog");
        if (SnNatLogW)
        {
            SnNatCSSumWaveformLog->Add(SnNatLogW);
        }

        TGraph * Sn124LogW = (TGraph*)infileWaveform->Get("Sn124CSLog");
        if (Sn124LogW)
        {
            Sn124CSSumWaveformLog->Add(Sn124LogW);
        }
        */

        // Close the sub-run input files
        infile->Close();
        infileWaveform->Close();

        // End of loop - move to next sub-run
    }

    // Create output file to contain summed histos
    fileOutName << outpath << "/analysis/run" << runDir << "/" << "sum.root";

    TFile *outfile = new TFile(fileOutName.str().c_str(),"RECREATE");

    vector<double> SnNatSum(SnNatEnergies[0]->size(), 0);

    for(int i=0; i<SnNatCrossSections.size(); i++)
    {
        for(int j=0; j<SnNatEnergies[0]->size(); j++)
        {
            SnNatSum[j] += SnNatCrossSections[i]->at(j);
        }
    }

    // normalize by the total number of graphs
    for(int i=0; i<SnNatSum.size(); i++)
    {
        SnNatSum[i] /= SnNatCrossSections.size();
    }

    TGraph* SnNatGraph = new TGraph(SnNatSum.size(),&SnNatEnergies.back()->at(0),&SnNatSum[0]);
    SnNatGraph->Write();

    // Write run histograms to sum.root
    outfile->Write();
    outfile->Close();
}
