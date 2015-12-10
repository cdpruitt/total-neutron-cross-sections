#include "TH1.h"
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

int main()
{

    const string direct = "/media/ExternalDrive1/analysis/run";

    ifstream runs("runsToSort.txt");
    if(!runs.is_open())
    {
        cout << "Could not find runs list...Exiting." << endl;
        return 0;
    }

    int runno = -1;

    const int noBins = 20000;
    const int noBinsscale = 200;

    TFile *outfile = new TFile("SummedHistos.root","RECREATE");

    TH1I *blankSum = new TH1I("blankSum","blankSum",noBins,0,700);
    TH1I *carbonSSum = new TH1I("carbonSSum","carbonSSum",noBins,0,700);
    TH1I *carbonLSum = new TH1I("carbonLSum","carbonLSum",noBins,0,700);
    TH1I *Sn112Sum = new TH1I("Sn112Sum","Sn112Sum",noBins,0,700);
    TH1I *NatSnSum = new TH1I("NatSnSum","NatSnSum",noBins,0,700);
    TH1I *Sn124Sum = new TH1I("Sn124Sum","Sn124Sum",noBins,0,700);
    TH1I *totalSum = new TH1I("totalSum","totalSum",noBins,0,700);

    TH1D *carbonScs = new TH1D("carbonScs","carbonScs",noBinsscale,0,700);
    TH1D *carbonLcs = new TH1D("carbonLcs","carbonLcs",noBinsscale,0,700);
    TH1D *Sn112cs = new TH1D("Sn112cs","Sn112cs",noBinsscale,0,700);
    TH1D *SnNatcs = new TH1D("SnNatcs","SnNatcs",noBinsscale,0,700);
    TH1D *Sn124cs = new TH1D("Sn124cs","Sn124cs",noBinsscale,0,700);

    TH1I *blankSumLog = new TH1I("blankSumLog","blankSum",noBins,0,TMath::Log10(700));
    TH1I *carbonSSumLog = new TH1I("carbonSSumLog","carbonSSum",noBins,0,TMath::Log10(700));
    TH1I *carbonLSumLog = new TH1I("carbonLSumLog","carbonLSum",noBins,0,TMath::Log10(700));
    TH1I *Sn112SumLog = new TH1I("Sn112SumLog","Sn112Sum",noBins,0,TMath::Log10(700));
    TH1I *NatSnSumLog = new TH1I("NatSnSumLog","NatSnSum",noBins,0,TMath::Log10(700));
    TH1I *Sn124SumLog = new TH1I("Sn124SumLog","Sn124Sum",noBins,0,TMath::Log10(700));
    TH1I *totalSumLog = new TH1I("totalSumLog","totalSum",noBins,0,TMath::Log10(700));

    TH1D *carbonScsLog = new TH1D("carbonScsLog","carbonScs",noBinsscale,0,TMath::Log10(700));
    TH1D *carbonLcsLog = new TH1D("carbonLcsLog","carbonLcs",noBinsscale,0,TMath::Log10(700));
    TH1D *Sn112csLog = new TH1D("Sn112csLog","Sn112cs",noBinsscale,0,TMath::Log10(700));
    TH1D *SnNatcsLog = new TH1D("SnNatcsLog","SnNatcs",noBinsscale,0,TMath::Log10(700));
    TH1D *Sn124csLog = new TH1D("Sn124csLog","Sn124cs",noBinsscale,0,TMath::Log10(700));

    TH1I *STOFSum = new TH1I("STOFSum","Summed-detector time of flight",10000,-1788.82*1.1,1788.82*1.1);

    TH1I *monSum = new TH1I("monSum","Monitor sum",noBins,0,700);

    double blankmonCounts = 0;
    double carbonSmonCounts = 0;
    double carbonLmonCounts = 0;
    double Sn112monCounts = 0;
    double SnNatmonCounts = 0;
    double Sn124monCounts = 0;

    double targetCounts[6][noBinsscale] = {{0.}};
    double targetCountsLog[6][noBinsscale] = {{0.}};

    for(;;) // loop over runnos
      {
        runs >> runno;
        if(runs.eof())break;

        int segment = 0;
        for(;;) //loop over run segments
        {
            TFile *myfile;
            if(segment <10)
                myfile =  new TFile(Form("/media/ExternalDrive1/analysis/run%i/run%i-000%i.root",runno,runno,segment));
            else if(segment <100)
                myfile =  new TFile(Form("/media/ExternalDrive1/analysis/run%i/run%i-00%i.root",runno,runno,segment));
            else if(segment < 1000)
                myfile =  new TFile(Form("/media/ExternalDrive1/analysis/run%i/run%i-0%i.root",runno,runno,segment));
            else
            {
                cout << "Segment number too large!! Edit here!" << endl;
                return 0;
            }

            if(!myfile->IsOpen())
            {
                cout << "Can't open root file run " << runno << " segment = " << segment << endl;
                break;
            }

            cout << "Adding run " << runno << " segment " << segment <<endl;

            TH1I * blank = ((TH1I*)myfile->Get("blank"));
            blankSum->Add(blank,1);

            TH1I * carbonS = (TH1I*)myfile->Get("carbonS");
            carbonSSum->Add(carbonS,1);

            TH1I * carbonL = (TH1I*)myfile->Get("carbonL");
            carbonLSum->Add(carbonL,1.);

            TH1I * Sn124 = (TH1I*)myfile->Get("Sn124");
            Sn124Sum->Add(Sn124,1.);

            TH1I * Sn112 = (TH1I*)myfile->Get("Sn112");
            Sn112Sum->Add(Sn112,1.);

            TH1I * NatSn = (TH1I*)myfile->Get("NatSn");
            NatSnSum->Add(NatSn,1.);

            TH1I * total = (TH1I*)myfile->Get("total");
            totalSum->Add(total,1.);

            TH1I * STOF = (TH1I*)myfile->Get("STOF");
            STOFSum->Add(STOF,1.);

            TH1I * monBlank = (TH1I*)myfile->Get("monBlank");
            monSum->Add(monBlank,1.);

            // log plots

            TH1I * blankLog = ((TH1I*)myfile->Get("blankLog"));
            blankSumLog->Add(blankLog,1);

            TH1I * carbonSLog = (TH1I*)myfile->Get("carbonSLog");
            carbonSSumLog->Add(carbonSLog,1);

            TH1I * carbonLLog = (TH1I*)myfile->Get("carbonLLog");
            carbonLSumLog->Add(carbonLLog,1.);

            TH1I * Sn124Log = (TH1I*)myfile->Get("Sn124Log");
            Sn124SumLog->Add(Sn124Log,1.);

            TH1I * Sn112Log = (TH1I*)myfile->Get("Sn112Log");
            Sn112SumLog->Add(Sn112Log,1.);

            TH1I * NatSnLog = (TH1I*)myfile->Get("NatSnLog");
            NatSnSumLog->Add(NatSnLog,1.);

            TH1I * totalLog = (TH1I*)myfile->Get("totalLog");
            totalSumLog->Add(totalLog,1.);

            // monitor paddle scalings

            blankmonCounts += ((TH1I*)myfile->Get("monBlank"))->GetEntries(); // scale x-sections by monitor counts
            carbonSmonCounts += ((TH1I*)myfile->Get("monCarbonS"))->GetEntries();
            carbonLmonCounts += ((TH1I*)myfile->Get("monCarbonL"))->GetEntries();
            Sn112monCounts += ((TH1I*)myfile->Get("monSn112"))->GetEntries();
            SnNatmonCounts += ((TH1I*)myfile->Get("monNatSn"))->GetEntries();
            Sn124monCounts += ((TH1I*)myfile->Get("monSn124"))->GetEntries();
            myfile->Close();

            if(segment == 0)
                break;

            segment++;

        }
    }
    blankSum->Rebin(100);
    carbonSSum->Rebin(100);
    carbonLSum->Rebin(100);
    Sn112Sum->Rebin(100);
    NatSnSum->Rebin(100);
    Sn124Sum->Rebin(100);

    blankSumLog->Rebin(100);
    carbonSSumLog->Rebin(100);
    carbonLSumLog->Rebin(100);
    Sn112SumLog->Rebin(100);
    NatSnSumLog->Rebin(100);
    Sn124SumLog->Rebin(100);

    cout << "blank Mon count = " << blankmonCounts << endl;
    cout << "carbonS Mon count = " << carbonSmonCounts << endl;
    cout << "carbonL Mon count = " << carbonLmonCounts << endl;
    cout << "Sn112 Mon count = " << Sn112monCounts << endl;
    cout << "SnNat Mon count = " << SnNatmonCounts << endl;
    cout << "Sn124 Mon count = " << Sn124monCounts << endl;

    double energy[noBinsscale] = {0.};
    double energyLog[noBinsscale] = {0.};

    for(int i=1;i<noBinsscale;i++)
      {
	energy[i] = blankSum->GetBinCenter(i);
	targetCounts[0][i] = blankSum->GetBinContent(i)/blankmonCounts;
	targetCounts[1][i] = carbonSSum->GetBinContent(i)/carbonSmonCounts;
	targetCounts[2][i] = carbonLSum->GetBinContent(i)/carbonLmonCounts;
	targetCounts[3][i] = Sn112Sum->GetBinContent(i)/Sn112monCounts;
	targetCounts[4][i] = NatSnSum->GetBinContent(i)/SnNatmonCounts;
	targetCounts[5][i] = Sn124Sum->GetBinContent(i)/Sn124monCounts;

        energyLog[i] = blankSumLog->GetBinCenter(i);
	targetCountsLog[0][i] = blankSumLog->GetBinContent(i)/blankmonCounts;
	targetCountsLog[1][i] = carbonSSumLog->GetBinContent(i)/carbonSmonCounts;
	targetCountsLog[2][i] = carbonLSumLog->GetBinContent(i)/carbonLmonCounts;
	targetCountsLog[3][i] = Sn112SumLog->GetBinContent(i)/Sn112monCounts;
	targetCountsLog[4][i] = NatSnSumLog->GetBinContent(i)/SnNatmonCounts;
	targetCountsLog[5][i] = Sn124SumLog->GetBinContent(i)/Sn124monCounts;

      }

    double sig[6][noBinsscale] ={{0.}};
    double sigLog[6][noBinsscale] ={{0.}};

    double targetlength[6] = {0,1.37,2.74,1.365,1.370,1.370}; //cm
    double targetMolMass[6] = {0,12.01,12.01,112,118.7,124}; //g/mol
    double targetdensity[6] = {0,2.3,2.3,6.89,7.31,7.63}; //g/cm^3
    double ava = 6.022*pow(10.,23.); //atoms/mol

    for(int j= 1;j<6;j++)
      {
	for(int i=1;i<noBinsscale;i++)
	  {
              if(targetCounts[0][i] <=0 || targetCounts[j][i] <=0)
              {
                  sig[j][i] =0.;
              }

              else
              {
                  sig[j][i] = -log(targetCounts[j][i]/targetCounts[0][i])/(targetlength[j]/targetMolMass[j]*targetdensity[j]*ava*pow(10.,-24));//barn
              }

              if(targetCountsLog[0][i] <=0 || targetCountsLog[j][i] <=0)
              {
                  sigLog[j][i] =0.;
              }

              else
              {
                  sigLog[j][i] = -log(targetCountsLog[j][i]/targetCountsLog[0][i])/(targetlength[j]/targetMolMass[j]*targetdensity[j]*ava*pow(10.,-24));//barn
              }
	  }
      }

    // TGraph *CScs = new TGraph(5000,energy,sig[0]);
    // TGraph *CLcs = new TGraph(5000,energy,sig[1]);
    // TGraph *Sn112cs = new TGraph(5000,energy,sig[2]);
    // TGraph *SnNatcs = new TGraph(5000,energy,sig[3]);
    // TGraph *Sn124cs = new TGraph(5000,energy,sig[4]);

    // carbonScs->Add(blankSum,carbonSSum,1.,-1);
    // carbonLcs->Add(blankSum,carbonLSum,1.,-1);
    // Sn112cs->Add(blankSum,Sn112Sum,1.,-1);
    // SnNatcs->Add(blankSum,NatSnSum,1.,-1);
    // Sn124cs->Add(blankSum,Sn124Sum,1.,-1);

    float scale = 0.01;

    for(int i=1; i<noBinsscale; i++)
     {
       carbonScs->SetBinContent(i,sig[1][i]);
       carbonLcs->SetBinContent(i,sig[2][i]);
       Sn112cs->SetBinContent(i,sig[3][i]);
       SnNatcs->SetBinContent(i,sig[4][i]);
       Sn124cs->SetBinContent(i,sig[5][i]);
       
     }

    for(int i=1; i<noBinsscale; i++)
    {
        carbonScsLog->SetBinContent(i,sigLog[1][i]);
        carbonLcsLog->SetBinContent(i,sigLog[2][i]);
        Sn112csLog->SetBinContent(i,sigLog[3][i]);
        SnNatcsLog->SetBinContent(i,sigLog[4][i]);
        Sn124csLog->SetBinContent(i,sigLog[5][i]);
    }

    carbonScsLog->GetXaxis()->SetTitle("Log E");
    carbonLcsLog->GetXaxis()->SetTitle("Log E");
    Sn112csLog->GetXaxis()->SetTitle("Log E");
    SnNatcsLog->GetXaxis()->SetTitle("Log E");
    Sn124csLog->GetXaxis()->SetTitle("Log E");

    carbonScsLog->GetXaxis()->CenterTitle();
    carbonLcsLog->GetXaxis()->CenterTitle();
    Sn112csLog->GetXaxis()->CenterTitle();
    SnNatcsLog->GetXaxis()->CenterTitle();
    Sn124csLog->GetXaxis()->CenterTitle();

    // outfile->cd();
    // CScs->Write();
    // CLcs->Write();
    // Sn112cs->Write();
    // SnNatcs->Write();
    // Sn124cs->Write();



    outfile->Write();


    outfile->Close();

}
