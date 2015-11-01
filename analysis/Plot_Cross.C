void Plot_Cross()
{

  gROOT->SetStyle("Pub");

  TFile *myfile =new TFile("SummedHistos.root");
  if(!myfile->IsOpen())
    {
      cout << "Can't Open summed file..." << endl;
      return;
    }

  TCanvas *mycan = new TCanvas("mycan","mycan",1200,600);
  mycan->Divide(3,1);


  mycan->cd(1);

  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(0.03);

  TH1I *carbonS = (TH1I*)myfile->Get("carbonSDiff");

  carbonS->Draw();

  int maxi = carbonS->GetMaximumBin();
  float max = carbonS->GetBinContent(maxi);

  carbonS->GetXaxis()->SetRangeUser(0,10);
  carbonS->GetYaxis()->SetRangeUser(0,max+500);
  carbonS->GetXaxis()->SetTitle("Approx Energy [Mev]");
  carbonS->GetYaxis()->SetTitle("#sigma [arb]");
  carbonS->GetYaxis()->CenterTitle();
  carbonS->GetYaxis()->SetTitleOffset(1.6);
  carbonS->GetXaxis()->CenterTitle();

  TLatex *mytex = new TLatex(5.,max*0.7,"Carbon Short");
  mytex->Draw();

  mycan->cd(2);
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(0.03);

  TH1I *carbonL = (TH1I*)myfile->Get("carbonLDiff");

  carbonL->Draw();

  maxi = carbonL->GetMaximumBin();
  max = carbonL->GetBinContent(maxi);

  carbonL->GetXaxis()->SetRangeUser(0,10);
  carbonL->GetYaxis()->SetRangeUser(0,max+500);
  carbonL->GetXaxis()->SetTitle("Approx Energy [Mev]");
  carbonL->GetXaxis()->CenterTitle();
  carbonL->GetYaxis()->SetTitle("#sigma [arb]");
  carbonL->GetYaxis()->CenterTitle();
  carbonL->GetYaxis()->SetTitleOffset(1.6);


  mytex->DrawLatex(5.,max*0.7,"Carbon Long");



  ifstream Data("CarbonData.dat");
  if(!Data.is_open())
    {
      cout << "No Previous Data..." << endl;
      return;
    }

  char bogus[200];
  Data.getline(bogus,200);

  float energy[1977] = {0.};
  float xsection[1977] = {0.};
  float error[1977] ={0.};

  float dum,dum2,dum3;

  for(int i = 0;i<1977;i++)
    {
      Data >> dum >> dum2>> dum3;
      if(Data.eof())break;
      energy[i] =dum;
      xsection[i] = dum2;
      error[i] =dum3;
    }

  mycan->cd(3);
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(0.03);

  TGraphErrors *mygraph = new TGraphErrors(1977,energy,xsection,0,error);
  mygraph->Draw("AP");

  mygraph->GetXaxis()->SetTitle("Energy [MeV]");
  mygraph->GetXaxis()->CenterTitle();
  mygraph->GetXaxis()->SetRangeUser(0,10.);

  mygraph->GetYaxis()->SetTitle("#sigma [b]");
  mygraph->GetYaxis()->CenterTitle();
  mygraph->GetYaxis()->SetTitleOffSet(1.5);


  return;
}
