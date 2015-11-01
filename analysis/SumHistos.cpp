#include "TH1.h"
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int main()
{


  ifstream runs("runsToSort.txt");
  if(!runs.is_open())
    {
      cout << "Could not find runs list...Exiting." << endl;
      return 0;
    }

  int runno = -1;


  TFile *outfile = new TFile("SummedHistos.root","RECREATE");

  TH1I *blankSum = new TH1I("blankSum","blankSum",700000,0,700);
  TH1I *carbonSSum = new TH1I("carbonSSum","carbonSSum",700000,0,700);
  TH1I *carbonLSum = new TH1I("carbonLSum","carbonLSum",700000,0,700);
  TH1I *Sn112Sum = new TH1I("Sn112Sum","Sn112Sum",700000,0,700);
  TH1I *NatSnSum = new TH1I("NatSnSum","NatSnSum",700000,0,700);
  TH1I *Sn124Sum = new TH1I("Sn124Sum","Sn124Sum",700000,0,700);
  TH1I *totalSum = new TH1I("totalSum","totalSum",700000,0,700);

  TH1I *carbonSDiff = new TH1I("carbonSDiff","carbonSDiff",700000,0,700);
  TH1I *carbonLDiff = new TH1I("carbonLDiff","carbonLDiff",700000,0,700);
  
  TH1I *STOFSum = new TH1I("STOFSum","Summed-detector time of flight",100000,0,2000);
  
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
	      cout << "Can't open root file run " << runno << " segement = " << segment << endl;
	      break;
	    }

	  cout << "Adding run " << runno << " segment " << segment <<endl;
	  TH1I * blank = (TH1I*)myfile->Get("blank");
	  blankSum->Add(blank,1.);
	  TH1I * carbonS = (TH1I*)myfile->Get("carbonS");
	  carbonSSum->Add(carbonS,1.);
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




	  //	  cout << "Added " << endl;

	  myfile->Close();


	  //	  if(segment == )
	  // break;

	  segment++;

	}
    }


  carbonSDiff->Add(blankSum,carbonSSum,1.,-1.);
  carbonLDiff->Add(blankSum,carbonLSum,1.,-1.);



  outfile->Write();


  outfile->Close();

}
