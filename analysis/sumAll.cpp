#include "TH1.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <vector>
#include <cmath>

using namespace std;

// Set number of bins for relative cross-section histograms
int noRelativeBins = 50;
int noCSBins = 0;
const double CLIFF_OFFSET = 0;

struct Error
{
} error;

long monCountsBlank;
long monCountsSn124;
long monCountsSn112;

double monCountsBlankError;
double monCountsSn124Error;
double monCountsSn112Error;

vector<long> detCountsBlank (200);
vector<long> detCountsSn124 (200);
vector<long> detCountsSn112 (200);

vector<double> detCountsBlankError;
vector<double> detCountsSn124Error;
vector<double> detCountsSn112Error;

vector<double> totalSn112StatError (200);
vector<double> totalSn124StatError (200);
vector<double> totalRelStatError (200);

double totalSn112SysError;
double totalSn124SysError;
double totalRelSysError;

vector<double> totalSn112Error (200);
vector<double> totalSn124Error (200);
vector<double> totalRelError (200);

TFile* infile;

void getStatisticalError(TH1D* Sn112CSTotal, TH1D* Sn124CSTotal)
{
    // Get monCounts, raw target counts from _histos.root files
    // Calculate statistical error from histogram counts

    int limit = 25;

    ifstream runList("runsToSort.txt");

    string runDir;
    string outPath;

    while (runList >> runDir)
    {

        if (stoi(runDir)<6)
        {
            outPath = "/data3/analysis/run";
        }

        else if (stoi(runDir)>128 && stoi(runDir)<160)
        {
            outPath = "/data2/analysis/run";
        }

        else if (stoi(runDir)>159 && stoi(runDir)<178)
        {
            outPath = "/data3/analysis/run";
        }

        else
        {
            cout << "Run directory outside bounds (runs 128-177) - check run number" << endl;
            exit(1);
        }

        // keep track of which order the targets are in, based on which run number we're
        // sorting
        vector<int> order;

        if(stoi(runDir)<=5)
        {
            cout << "Neon run - stopping sort." << endl;
            exit(0);
        }

        if(stoi(runDir)<=151)
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

        else if(stoi(runDir)>=153 && stoi(runDir)<=168)
        {
            // blank, Sn112, NatSn, Sn124
            order.push_back(0);
            order.push_back(3);
            order.push_back(4);
            order.push_back(5);
        }

        else if(stoi(runDir)>=169 && stoi(runDir)<=180)
        {
            // blank, Sn112, NatSn, Sn124, short carbon
            order.push_back(0);
            order.push_back(3);
            order.push_back(4);
            order.push_back(5);
            order.push_back(1);
        }

        int pos112 = find(order.begin(), order.end(), 3) - order.begin();
        int pos124 = find(order.begin(), order.end(), 5) - order.begin();

        for(int segment = 0; segment<=limit; segment++)
        {
            // We need to form the proper name for the sub-run we want to open:
            if(segment < 10)
            {
                infile =  new TFile(Form("%s%s/run%s-000%i_histos.root",outPath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            }

            else if(segment < 100)
            {
                infile =  new TFile(Form("%s%s/run%s-00%i_histos.root",outPath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            }

            else
            {
                // There's some error in the sub-run file number; write outfile and exit
                cout << "Segment number too large!" << endl;
                exit(1);
            }

            // Attempt to open the sub-run to access its histograms
            if(!infile->IsOpen())
            {
                cout << "Can't open root file run " << runDir << " segment = " << segment << endl;

                break;
            }

            gDirectory->cd("/");
            gDirectory->GetDirectory("monitor")->cd();

            monCountsBlank += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2);
            monCountsSn112 += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2+pos112);
            monCountsSn124 += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2+pos124);

            gDirectory->cd("/");
            gDirectory->GetDirectory("detS")->cd();

            TH1I* blank = ((TH1I*)gDirectory->Get("blankLog"));

            string name112 = "target" + to_string(pos112) + "Log";
            string name124 = "target" + to_string(pos124) + "Log";

            TH1I* Sn112 = ((TH1I*)gDirectory->Get(name112.c_str()));
            TH1I* Sn124 = ((TH1I*)gDirectory->Get(name124.c_str()));

            for(int i=0; i<blank->GetNbinsX(); i++)
            {
                detCountsBlank[i] += (blank->GetBinContent(i));
                detCountsSn112[i] += (Sn112->GetBinContent(i));
                detCountsSn124[i] += (Sn124->GetBinContent(i));
            }

            // Close the sub-run input files
            infile->Close();

            // End of loop - move to next sub-run
        }
    }

    /*monCountsBlankError = pow(monCountsBlank,0.5);
    monCountsSn112Error = pow(monCountsSn112,0.5);
    monCountsSn124Error = pow(monCountsSn124,0.5);
    */

    cout << "monCountsBlank = " <<  monCountsBlank << endl;
    cout << "monCountsSn112 = " <<  monCountsSn112 << endl;
    cout << "monCountsSn124 = " <<  monCountsSn124 << endl;

    cout << "detCountsBlank[10] = " << detCountsBlank[10] << endl;
    cout << "detCountsSn112[10] = " << detCountsSn112[10] << endl;
    cout << "detCountsSn124[10] = " << detCountsSn124[10] << endl;


    for(int i=0; i<totalSn112StatError.size(); i++)
    {
        totalSn112StatError[i] = pow(
                ((double)monCountsBlank/(pow(monCountsBlank,2))) +
                ((double)monCountsSn112/(pow(monCountsSn112,2))) +
                ((double)detCountsBlank[i]/(pow(detCountsBlank[i],2))) +
                ((double)detCountsSn112[i]/(pow(detCountsSn112[i],2)))
                ,0.5)/*/exp(abs(Sn112CSTotal->GetBinContent(i)))*/;

        totalSn124StatError[i] = pow(
                ((double)monCountsBlank/(pow(monCountsBlank,2))) +
                ((double)monCountsSn124/(pow(monCountsSn124,2))) +
                ((double)detCountsBlank[i]/(pow(detCountsBlank[i],2))) +
                ((double)detCountsSn124[i]/(pow(detCountsSn124[i],2)))
                ,0.5)/*/exp(abs(Sn124CSTotal->GetBinContent(i)))*/;

        /*totalRelStatError[i] = pow(
                pow(totalSn112Error[i],2) +
                pow(totalSn124Error[i],2)
                ,0.5);
                */
    }

    cout << totalSn112StatError[10] << endl;

        /*detCountsBlankpush_back(pow(detCountsBlank[i],0.5));
        detCountsSn112push_back(pow(detCountsSn112[i],0.5));
        detCountsSn124push_back(pow(detCountsSn124[i],0.5));
        */
}

void getSystemicError()
{
    // physical target data, listed in order:
    // {blank, Sn112, Natural Sn, Sn124, short carbon, long carbon} 

    // lengths of each target:
    double targetlength[6] = {0,1.37,2.74,1.365,1.370,1.370}; //cm
    double targetlengthUncertainty[6] = {0,0.01,0.01,0.005,0.005,0.005}; //cm

    // molar mass of each target:
    double targetMolMass[6] = {0,12.01,12.01,112,118.7,124}; //g/mol
    double targetMolMassUncertainty[6] = {0,0.01,0.01,0.01,0.01,0.01}; //g/mol

    // density of each target:
    double targetdensity[6] = {0,2.2,2.2,6.89,7.31,7.63}; //g/cm^3
    double targetdensityUncertainty[6] = {0,0.1,0.1,0.001,0.001,0.001}; //g/cm^3

    totalSn112SysError = pow(
            pow(targetlengthUncertainty[3]/targetlength[3],2) +
            pow(targetMolMassUncertainty[3]/targetMolMass[3],2) +
            pow(targetdensityUncertainty[3]/targetdensity[3],2)
            ,0.5);

    totalSn124SysError = pow(
            pow(targetlengthUncertainty[5]/targetlength[5],2) +
            pow(targetMolMassUncertainty[5]/targetMolMass[5],2) +
            pow(targetdensityUncertainty[5]/targetdensity[5],2)
            ,0.5);
}

void getTotalError()
{
    for(int i=0; i<totalSn112Error.size(); i++)
    {
        totalSn112Error[i] = pow(
                pow(totalSn112StatError[i],2) +
                pow(totalSn112SysError,2)
                ,0.5);
        totalSn124Error[i] = pow(
                pow(totalSn124StatError[i],2) +
                pow(totalSn124SysError,2)
                ,0.5);
        totalRelError[i] = pow(
                pow(totalSn112Error[i],2) +
                pow(totalSn124Error[i],2)
                ,0.5);
    }

    cout << totalSn112Error[10] << endl;
    cout << totalSn124Error[10] << endl;
    cout << totalRelError[10] << endl;
}

TH1D* logBins(TH1D *inputHisto)
{
    string newName;
    newName = inputHisto->GetName();
    newName += "Log";
    TH1D* outputHisto = new TH1D(newName.c_str(),newName.c_str(),noCSBins,
            TMath::Log10(((TAxis*)inputHisto->GetXaxis())->GetXmin()),
            TMath::Log10(((TAxis*)inputHisto->GetXaxis())->GetXmax()));

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

    ((TAxis*)outputHisto->GetXaxis())->Set(nBins,newBins);
    delete newBins;

    return outputHisto;
}

int main()
{
    // Declare stringstreams needed to open and access runs
    stringstream fileInName;
    stringstream fileInWaveformName;
    stringstream fileOutName;

    // Open the runList to see which runs should be added to the total

    ifstream runList("runsToSort.txt");
    if(!runList.is_open())
    {
        cout << "Could not find runs list... exiting." << endl;
        exit(1);
    }

    // Examine the first run of runlist and find out the correct number of
    // bins to use for the summed histograms
    if(runList.eof())
    {
        cout << "Error: empty runList file" << endl;
        exit(1);
    }

    string runDir;
    string outPath;

    runList >> runDir;

    if (stoi(runDir)<6)
    {
        outPath = "/data3/analysis/run";
    }

    else if (stoi(runDir)>128 && stoi(runDir)<160)
    {
        outPath = "/data2/analysis/run";
    }

    else if (stoi(runDir)>159 && stoi(runDir)<178)
    {
        outPath = "/data3/analysis/run";
    }

    else
    {
        cout << "Run directory outside bounds (runs 128-177) - check run number" << endl;
        exit(1);
    }

    fileInName << outPath << runDir << "/" << "sum.root";

    infile =  new TFile(fileInName.str().c_str(),"READ");
    TH1D* csBins = ((TH1D*)infile->Get("blankCSSum"));
    TH1D* waveformBins = ((TH1D*)infile->Get("blankCSSumWaveform"));

    // each histo has N bins + 1 overflow + 1 underflow
    // thus, subtract two to get the number of 'normal' bins for summed histos
    noCSBins = csBins->GetSize()-2;
    int noWaveformBins = waveformBins->GetSize()-2;

    delete csBins;
    delete waveformBins;

    // Find the total number of runs in runList
    int totalRuns = 1;

    while(runList >> runDir)
    {
        totalRuns++;
    }

    cout << "Total runs in runlist: " << totalRuns << endl << endl;

    runList.close();

    // Create output file to contain summed histos
    fileOutName << "/data3/analysis/total.root";

    TFile *outfile = new TFile(fileOutName.str().c_str(),"RECREATE");

    // Create summed histos
    // First, sum DPP-derived cross-sections
    TH1D *blankCSTotal = new TH1D("blankCSTotal","blankCSTotal",noCSBins,1,700);
    TH1D *carbonSCSTotal = new TH1D("carbonSCSTotal","carbonSCSTotal",noCSBins,1,700);
    TH1D *carbonLCSTotal = new TH1D("carbonLCSTotal","carbonLCSTotal",noCSBins,1,700);
    TH1D *Sn112CSTotal = new TH1D("Sn112CSTotal","Sn112CSTotal",noCSBins,1,700);
    TH1D *SnNatCSTotal = new TH1D("SnNatCSTotal","SnNatCSTotal",noCSBins,1,700);
    TH1D *Sn124CSTotal = new TH1D("Sn124CSTotal","Sn124CSTotal",noCSBins,1,700);

    TH1D *blankCSTotalLog = logBins(blankCSTotal);
    TH1D *carbonSCSTotalLog = logBins(carbonSCSTotal);
    TH1D *carbonLCSTotalLog = logBins(carbonLCSTotal);
    TH1D *Sn112CSTotalLog = logBins(Sn112CSTotal);
    TH1D *SnNatCSTotalLog = logBins(SnNatCSTotal);
    TH1D *Sn124CSTotalLog = logBins(Sn124CSTotal);

    // Next, sum waveform-derived cross-sections
    TH1D *blankCSTotalWaveform = new TH1D("blankCSTotalWaveform","blankCSTotalWaveform",noWaveformBins,0,700);
    TH1D *carbonSCSTotalWaveform = new TH1D("carbonSCSTotalWaveform","carbonSCSTotalWaveform",noWaveformBins,0,700);
    TH1D *carbonLCSTotalWaveform = new TH1D("carbonLCSTotalWaveform","carbonLCSTotalWaveform",noWaveformBins,0,700);
    TH1D *Sn112CSTotalWaveform = new TH1D("Sn112CSTotalWaveform","Sn112CSTotalWaveform",noWaveformBins,0,700);
    TH1D *SnNatCSTotalWaveform = new TH1D("SnNatCSTotalWaveform","SnNatCSTotalWaveform",noWaveformBins,0,700);
    TH1D *Sn124CSTotalWaveform = new TH1D("Sn124CSTotalWaveform","Sn124CSTotalWaveform",noWaveformBins,0,700);

    TH1D *blankCSTotalWaveformLog = new TH1D("blankCSTotalWaveformLog","blankCSTotalWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *carbonSCSTotalWaveformLog = new TH1D("carbonSCSTotalWaveformLog","carbonSCSTotalWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *carbonLCSTotalWaveformLog = new TH1D("carbonLCSTotalWaveformLog","carbonLCSTotalWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *Sn112CSTotalWaveformLog = new TH1D("Sn112CSTotalWaveformLog","Sn112CSTotalWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *SnNatCSTotalWaveformLog = new TH1D("SnNatCSTotalWaveformLog","SnNatCSTotalWaveformLog",noWaveformBins,0,TMath::Log10(700));
    TH1D *Sn124CSTotalWaveformLog = new TH1D("Sn124CSTotalWaveformLog","Sn124CSTotalWaveformLog",noWaveformBins,0,TMath::Log10(700));

    // Re-open runList from the start
    runList.open("runsToSort.txt");

    // Main loop over runList
    while (runList >> runDir)
    {
        // Open the run's summed histos file
        fileInName.str("");

        if (stoi(runDir)<6)
        {
            outPath = "/data3/analysis/run";
        }

        else if (stoi(runDir)>128 && stoi(runDir)<160)
        {
            outPath = "/data2/analysis/run";
        }

        else if (stoi(runDir)>159 && stoi(runDir)<178)
        {
            outPath = "/data3/analysis/run";
        }

        else
        {
            cout << "Run directory outside bounds (runs 128-177) - check run number" << endl;
            exit(1);
        }

        fileInName << outPath << runDir << "/" << "sum.root";
        infile =  new TFile(fileInName.str().c_str());

        if(!infile->IsOpen())
        {
            cout << "Can't open sum.root file of run" << runDir << endl;
            outfile->Write();
            outfile->Close();
            exit(1);
        }

        cout << "Adding run" << runDir << endl;


        /****************** SUM HISTOGRAMS ******************/

        /************** normal cross-section histos ************/
        TH1D * blank = ((TH1D*)infile->Get("blankCSSum"));
        if (blank)
        {
            blankCSTotal->Add(blank,1/(double)totalRuns);
        }

        TH1D * carbonS = (TH1D*)infile->Get("carbonSCSSum");
        if (carbonS)
        {
            carbonSCSTotal->Add(carbonS,1/(double)totalRuns);
        }

        TH1D * carbonL = (TH1D*)infile->Get("carbonLCSSum");
        if (carbonL)
        {
            carbonLCSTotal->Add(carbonL,1/(double)totalRuns);
        }

        TH1D * Sn112 = (TH1D*)infile->Get("Sn112CSSum");
        if (Sn112)
        {
            Sn112CSTotal->Add(Sn112,1/(double)totalRuns);
        }

        TH1D * SnNat = (TH1D*)infile->Get("SnNatCSSum");
        if (SnNat)
        {
            SnNatCSTotal->Add(SnNat,1/(double)totalRuns);
        }

        TH1D * Sn124 = (TH1D*)infile->Get("Sn124CSSum");
        if (Sn124)
        {
            Sn124CSTotal->Add(Sn124,1/(double)totalRuns);
        }

        /************* LOG-SCALE cross-section histos **************/
        TH1D * blankLog = ((TH1D*)infile->Get("blankCSSumLog"));
        if (blankLog)
        {
            blankCSTotalLog->Add(blankLog,1/(double)totalRuns);
        }

        TH1D * carbonSLog = (TH1D*)infile->Get("carbonSCSSumLog");
        if (carbonSLog)
        {
            carbonSCSTotalLog->Add(carbonSLog,1/(double)totalRuns);
        }

        TH1D * carbonLLog = (TH1D*)infile->Get("carbonLCSSumLog");
        if (carbonLLog)
        {
            carbonLCSTotalLog->Add(carbonLLog,1/(double)totalRuns);
        }

        TH1D * Sn112Log = (TH1D*)infile->Get("Sn112CSSumLog");
        if (Sn112Log)
        {
            Sn112CSTotalLog->Add(Sn112Log,1/(double)totalRuns);
        }

        TH1D * SnNatLog = (TH1D*)infile->Get("SnNatCSSumLog");
        if (SnNatLog)
        {
            SnNatCSTotalLog->Add(SnNatLog,1/(double)totalRuns);
        }

        TH1D * Sn124Log = (TH1D*)infile->Get("Sn124CSSumLog");
        if (Sn124Log)
        {
            Sn124CSTotalLog->Add(Sn124Log,1/(double)totalRuns);
        }

        /************* WAVEFORM-DERIVED cross-section histos **************/
        TH1D * blankW = ((TH1D*)infile->Get("blankCSSumWaveform"));
        if (blankW)
        {
            blankCSTotalWaveform->Add(blankW,1/(double)totalRuns);
        }

        TH1D * carbonSW = (TH1D*)infile->Get("carbonSCSSumWaveform");
        if (carbonSW)
        {
            carbonSCSTotalWaveform->Add(carbonSW,1/(double)totalRuns);
        }

        TH1D * carbonLW = (TH1D*)infile->Get("carbonLCSSumWaveform");
        if (carbonLW)
        {
            carbonLCSTotalWaveform->Add(carbonLW,1/(double)totalRuns);
        }

        TH1D * Sn112W = (TH1D*)infile->Get("Sn112CSSumWaveform");
        if (Sn112W)
        {
            Sn112CSTotalWaveform->Add(Sn112W,1/(double)totalRuns);
        }

        TH1D * SnNatW = (TH1D*)infile->Get("SnNatCSSumWaveform");
        if (SnNatW)
        {
            SnNatCSTotalWaveform->Add(SnNatW,1/(double)totalRuns);
        }

        TH1D * Sn124W = (TH1D*)infile->Get("Sn124CSSumWaveform");
        if (Sn124W)
        {
            Sn124CSTotalWaveform->Add(Sn124W,1/(double)totalRuns);
        }

        /********** LOG-SCALE WAVEFORM-DERIVED cross-section histos ***********/
        TH1D * blankLogW = ((TH1D*)infile->Get("blankCSSumWaveformLog"));
        if (blankLogW)
        {
            blankCSTotalWaveformLog->Add(blankLogW,1/(double)totalRuns);
        }

        TH1D * carbonSLogW = (TH1D*)infile->Get("carbonSCSSumWaveformLog");
        if (carbonSLogW)
        {
            carbonSCSTotalWaveformLog->Add(carbonSLogW,1/(double)totalRuns);
        }

        TH1D * carbonLLogW = (TH1D*)infile->Get("carbonLCSSumWaveformLog");
        if (carbonLLogW)
        {
            carbonLCSTotalWaveformLog->Add(carbonLLogW,1/(double)totalRuns);
        }

        TH1D * Sn112LogW = (TH1D*)infile->Get("Sn112CSSumWaveformLog");
        if (Sn112LogW)
        {
            Sn112CSTotalWaveformLog->Add(Sn112LogW,1/(double)totalRuns);
        }

        TH1D * SnNatLogW = (TH1D*)infile->Get("SnNatCSSumWaveformLog");
        if (SnNatLogW)
        {
            SnNatCSTotalWaveformLog->Add(SnNatLogW,1/(double)totalRuns);
        }

        TH1D * Sn124LogW = (TH1D*)infile->Get("Sn124CSSumWaveformLog");
        if (Sn124LogW)
        {
            Sn124CSTotalWaveformLog->Add(Sn124LogW,1/(double)totalRuns);
        }

        // Close the sub-run input files
        infile->Close();

        // End of loop - move to next sub-run
    }

    // create relative cross-section plot for 112Sn/124Sn
    TH1D *Sn124_plus_Sn112CS = (TH1D*)Sn124CSTotal->Clone("Sn124_plus_Sn112CS");
    Sn124_plus_Sn112CS->SetDirectory(outfile);
    Sn124_plus_Sn112CS->Add(Sn112CSTotal);

    TH1D *Sn124_minus_Sn112CS = (TH1D*)Sn124CSTotal->Clone("Sn124_minus_Sn112CS");
    Sn124_minus_Sn112CS->SetDirectory(outfile);
    Sn124_minus_Sn112CS->Add(Sn112CSTotal,-1);

    /*for(int i = 0; i<Sn124_plus_Sn112CS->GetXaxis()->GetNbins(); i++)
    {
        if(Sn124_plus_Sn112CS->GetBinContent(i) < 0)
        {
            Sn124_plus_Sn112CS->AddBinContent(i,CLIFF_OFFSET);
        }
    }*/

    //Sn124_plus_Sn112CS->Rebin(Sn124_plus_Sn112CS->GetSize()/(double)50);
    //Sn124_minus_Sn112CS->Rebin(Sn124_minus_Sn112CS->GetSize()/(double)50);

    //Sn124_plus_Sn112CS->Scale(1/(Sn124_plus_Sn112CS->GetSize()/(double)50));
    //Sn124_minus_Sn112CS->Scale(1/(Sn124_plus_Sn112CS->GetSize()/(double)50));

    TH1D *relativeSnCS = (TH1D*)Sn124_minus_Sn112CS->Clone("relativeSnCS");
    relativeSnCS->SetDirectory(outfile);
    relativeSnCS->SetTitle("#frac{#sigma_{^{124}Sn}-#sigma_{^{112}Sn}}{#sigma_{^{124}Sn}+#sigma_{^{112}Sn}}");

    relativeSnCS->Divide(Sn124_plus_Sn112CS);

    //getStatisticalError(Sn112CSTotalLog, Sn124CSTotalLog);
    //getSystemicError();
    //getTotalError();

    

    relativeSnCS->GetXaxis()->SetRangeUser(3,300);

    // create relative cross-section plot for 112Sn/124Sn
    TH1D *Sn124_plus_Sn112CSLog = (TH1D*)Sn124CSTotalLog->Clone("Sn124_plus_Sn112CSLog");
    Sn124_plus_Sn112CSLog->SetDirectory(outfile);
    Sn124_plus_Sn112CSLog->Add(Sn112CSTotalLog);

    TH1D *Sn124_minus_Sn112CSLog = (TH1D*)Sn124CSTotalLog->Clone("Sn124_minus_Sn112CSLog");
    Sn124_minus_Sn112CSLog->SetDirectory(outfile);
    Sn124_minus_Sn112CSLog->Add(Sn112CSTotalLog,-1);

    for(int i = 0; i<Sn124_plus_Sn112CS->GetXaxis()->GetNbins(); i++)
    {
        if(Sn124_plus_Sn112CS->GetBinContent(i) < 0)
        {
            Sn124_plus_Sn112CS->AddBinContent(i,CLIFF_OFFSET);
        }
    }

    for(int i = 0; i<Sn124_plus_Sn112CSLog->GetXaxis()->GetNbins(); i++)
    {
        if(Sn124_plus_Sn112CSLog->GetBinContent(i) < 0)
        {
            Sn124_plus_Sn112CSLog->AddBinContent(i,CLIFF_OFFSET);
        }
    }

    Sn124_plus_Sn112CS->Rebin(Sn124_plus_Sn112CS->GetSize()/noRelativeBins);
    Sn124_minus_Sn112CS->Rebin(Sn124_minus_Sn112CS->GetSize()/(double)noRelativeBins);

    Sn124_plus_Sn112CS->Scale(1/(Sn124_plus_Sn112CS->GetSize()/(double)noRelativeBins));
    Sn124_minus_Sn112CS->Scale(1/(Sn124_plus_Sn112CS->GetSize()/(double)noRelativeBins));

    Sn124_plus_Sn112CSLog->Rebin(Sn124_plus_Sn112CSLog->GetSize()/(double)noRelativeBins);
    Sn124_minus_Sn112CSLog->Rebin(Sn124_minus_Sn112CSLog->GetSize()/(double)noRelativeBins);

    Sn124_plus_Sn112CSLog->Scale(1/(Sn124_plus_Sn112CSLog->GetSize()/(double)noRelativeBins));
    Sn124_minus_Sn112CSLog->Scale(1/(Sn124_plus_Sn112CSLog->GetSize()/(double)noRelativeBins));

    TH1D *relativeSnCSLog = (TH1D*)Sn124_minus_Sn112CSLog->Clone("relativeSnCSLog");
    relativeSnCSLog->SetDirectory(outfile);
    relativeSnCSLog->SetTitle("#frac{#sigma_{^{124}Sn}-#sigma_{^{112}Sn}}{#sigma_{^{124}Sn}+#sigma_{^{112}Sn}}");

    relativeSnCSLog->Divide(Sn124_plus_Sn112CSLog);
    relativeSnCSLog->GetXaxis()->SetRangeUser(3,300);

    for (int nBins = 0; nBins<relativeSnCSLog->GetNbinsX(); nBins++)
    {
        relativeSnCSLog->SetBinError(nBins,totalRelError[nBins]);
        //cout << totalRelError[nBins] << endl;
    }

    vector<double> energy;
    vector<double> xsection;
    vector<double> error;

    for(int j=0; j<relativeSnCSLog->GetNbinsX(); j++)
    {
        if(j==30 || j==31 || j==32)
        {
            continue;
        }

        if(relativeSnCSLog->GetXaxis()->GetBinCenter(j)<3 || relativeSnCSLog->GetXaxis()->GetBinCenter(j)>300)
        {
            continue;
        }

        energy.push_back(relativeSnCSLog->GetXaxis()->GetBinCenter(j));
        xsection.push_back(relativeSnCSLog->GetBinContent(j));
        error.push_back((totalRelError[4*j]+totalRelError[4*j+1]+totalRelError[4*j+2]+totalRelError[4*j+3])/4);
    }

    TGraphErrors* relativeSnCSLogGraph = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
    //TGraph* relativeSnCSLogGraph = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);

    //cout << relativeSnCSLogGraph->
    outfile->cd();
    relativeSnCSLogGraph->Write("relativeGraph");

    //relativeSnCSLog->SetBinContent(23,0);
    //relativeSnCSLog->SetBinContent(24,0);

    // create relative cross-section plot for 112Sn/124Sn using waveform data
    TH1D *Sn124_plus_Sn112CS_W = (TH1D*)Sn124CSTotalWaveform->Clone("Sn124_plus_Sn112CS_W");
    Sn124_plus_Sn112CS_W->SetDirectory(outfile);
    Sn124_plus_Sn112CS_W->Add(Sn112CSTotalWaveform);

    TH1D *Sn124_minus_Sn112CS_W = (TH1D*)Sn124CSTotalWaveform->Clone("Sn124_minus_Sn112CS_W");
    Sn124_minus_Sn112CS_W->SetDirectory(outfile);
    Sn124_minus_Sn112CS_W->Add(Sn112CSTotalWaveform,-1);

    //Sn124_plus_Sn112CS_W->Rebin(Sn124_plus_Sn112CS_W->GetSize()/(double)50);
    //Sn124_minus_Sn112CS_W->Rebin(Sn124_minus_Sn112CS_W->GetSize()/(double)50);

    //Sn124_plus_Sn112CS_W->Scale(1/(Sn124_plus_Sn112CS_W->GetSize()/(double)50));
    //Sn124_minus_Sn112CS_W->Scale(1/(Sn124_plus_Sn112CS_W->GetSize()/(double)50));

    TH1D *relativeSnCS_W = (TH1D*)Sn124_minus_Sn112CS_W->Clone("relativeSnCS_W");
    relativeSnCS_W->SetDirectory(outfile);
    relativeSnCS_W->SetTitle("#frac{#sigma_{^{124}Sn}-#sigma_{^{112}Sn}}{#sigma_{^{124}Sn}+#sigma_{^{112}Sn}}");
    relativeSnCS_W->Divide(Sn124_plus_Sn112CS_W);

    relativeSnCS_W->GetXaxis()->SetRangeUser(3,100);


    // Write run histograms to sum.root
    outfile->Write();
    outfile->Close();
}
