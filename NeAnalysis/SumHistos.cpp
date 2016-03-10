#include "TH1.h"
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

int limit = 2;

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

    const int noBins = 1000;
    const int noBinsscale = 500;

    TFile *outfile = new TFile("SummedHistos.root","RECREATE");

    TH1I *blankSum = new TH1I("blankSum","blankSum",noBins,0,700);
    TH1I *target1Sum = new TH1I("target1Sum","target1Sum",noBins,0,700);
    TH1I *target2Sum = new TH1I("target2Sum","target2Sum",noBins,0,700);
    TH1I *target3Sum = new TH1I("target3Sum","target3Sum",noBins,0,700);
    TH1I *target4Sum = new TH1I("target4Sum","target4Sum",noBins,0,700);
    TH1I *target5Sum = new TH1I("target5Sum","target5Sum",noBins,0,700);

    TH1D *target1csSum = new TH1D("target1csSum","target1csSum",noBins,0,700);
    TH1D *target2csSum = new TH1D("target2csSum","target2csSum",noBins,0,700);
    TH1D *target3csSum = new TH1D("target3csSum","target3csSum",noBins,0,700);
    TH1D *target4csSum = new TH1D("target4csSum","target4csSum",noBins,0,700);
    TH1D *target5csSum = new TH1D("target5csSum","target5csSum",noBins,0,700);

    TH1D *target1csSumLog = new TH1D("target1csSumLog","target1csSumLog",noBins,0,TMath::Log10(700));
    TH1D *target2csSumLog = new TH1D("target2csSumLog","target2csSumLog",noBins,0,TMath::Log10(700));
    TH1D *target3csSumLog = new TH1D("target3csSumLog","target3csSumLog",noBins,0,TMath::Log10(700));
    TH1D *target4csSumLog = new TH1D("target4csSumLog","target4csSumLog",noBins,0,TMath::Log10(700));
    TH1D *target5csSumLog = new TH1D("target5csSumLog","target5csSumLog",noBins,0,TMath::Log10(700));

    while(!runs.eof())
    {
        runs >> runno;

        for(int segment=0; segment<=limit; segment++)
        {
            TFile *myfile;
            if(segment <10)
                myfile =  new TFile(Form("/media/Drive3/analysis/run%i/run%i-000%i_cross-sections.root",runno,runno,segment));
            else if(segment <100)
                myfile =  new TFile(Form("/media/Drive3/analysis/run%i/run%i-00%i_cross-sections.root",runno,runno,segment));
            else if(segment < 1000)
                myfile =  new TFile(Form("/media/Drive3/analysis/run%i/run%i-0%i_cross-sections.root",runno,runno,segment));
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

            //TH1D * blank = ((TH1D*)myfile->Get("blankcs"));
            //blankcsSum->Add(blank,1/(segment+1));

            TH1D * target1 = (TH1D*)myfile->Get("target1cs");
            target1csSum->Add(target1,1/(double)(limit+1));

            TH1D * target2 = (TH1D*)myfile->Get("target2cs");
            target2csSum->Add(target2,1/(double)(limit+1));

            TH1D * target3 = (TH1D*)myfile->Get("target3cs");
            target3csSum->Add(target3,1/(double)(limit+1));

            TH1D * target4 = (TH1D*)myfile->Get("target4cs");
            target4csSum->Add(target4,1/(double)(limit+1));

            //TH1D * target5 = (TH1D*)myfile->Get("target5cs");
            //target5csSum->Add(target5,1/(double)(limit+1));


            //TH1D * blankLog = ((TH1D*)myfile->Get("blankcsLog"));
            //blankcsSumLog->Add(blank,1/(double)(limit+1));

            TH1D * target1Log = (TH1D*)myfile->Get("target1csLog");
            target1csSumLog->Add(target1Log,1/(double)(limit+1));

            TH1D * target2Log = (TH1D*)myfile->Get("target2csLog");
            target2csSumLog->Add(target2Log,1/(double)(limit+1));

            TH1D * target3Log = (TH1D*)myfile->Get("target3csLog");
            target3csSumLog->Add(target3Log,1/(double)(limit+1));

            TH1D * target4Log = (TH1D*)myfile->Get("target4csLog");
            target4csSumLog->Add(target4Log,1/(double)(limit+1));

            //TH1D * target5Log = (TH1D*)myfile->Get("target5csLog");
            //target5csSumLog->Add(target5,1/(double)(limit+1));

            myfile->Close();
        }
    }

/*    blankSum->Rebin(4);
    target1Sum->Rebin(4);
    target2Sum->Rebin(4);
    target3Sum->Rebin(4);
    target4Sum->Rebin(4);
    target5Sum->Rebin(4);
    */

    outfile->Write();
    outfile->Close();
}
