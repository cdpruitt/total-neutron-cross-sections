void Plot_Cross()
{

  gROOT->SetStyle("Pub");

  TFile *myfile =new TFile("/media/Drive3/analysis/run170/run170-0024_cross-sections.root","UPDATE");
  if(!myfile->IsOpen())
    {
      cout << "Can't Open summed file..." << endl;
      return;
    }
  /*

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

*/

/*mycan->cd(3);
  gPad->SetTopMargin(0.03);
  gPad->SetRightMargin(0.03);
  */

  ifstream SnData("SnNatData.dat");
  if(!SnData.is_open())
    {
      cout << "No Previous Data..." << endl;
      return;
    }

  char dummy[200];
  SnData.getline(dummy,200);

  vector<float> energy;
  vector<float> xsection;
  vector<float> error;

  float dum,dum2,dum3;

  while(!SnData.eof())
  {
      SnData >> dum >> dum2 >> dum3;

      energy.push_back(TMath::Log10(dum));
      xsection.push_back(dum2);
      error.push_back(dum3);
  }

  TGraphErrors *SnLitLog = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
  //SnLitLog->Draw("AP");

  SnLitLog->GetXaxis()->SetTitle("Energy [MeV] (log10)");
  SnLitLog->GetXaxis()->CenterTitle();
  SnLitLog->GetXaxis()->SetRangeUser(0,TMath::Log10(700.));

  SnLitLog->GetYaxis()->SetTitle("sigma [b]");
  SnLitLog->GetYaxis()->CenterTitle();
  //SnLitLog->GetYaxis()->SetTitleOffSet(1.5);
  SnLitLog->Write();

  // carbon literature data
  ifstream carbonData("CarbonData.dat");
  if(!carbonData.is_open())
  {
      cout << "No Previous Data..." << endl;
      return;
  }

  carbonData.getline(dummy,200);

  energy.clear();
  xsection.clear();
  error.clear();
  
  float dum,dum2,dum3;

  while(!carbonData.eof())
  {
      carbonData >> dum >> dum2 >> dum3;

      energy.push_back(TMath::Log10(dum));
      xsection.push_back(dum2);
      error.push_back(dum3);
  }

  TGraphErrors *carbonLitLog = new TGraphErrors(energy.size(),&energy[0],&xsection[0],0,&error[0]);
  carbonLitLog->Draw("AP");

  carbonLitLog->GetXaxis()->SetTitle("Energy [MeV] (log10)");
  carbonLitLog->GetXaxis()->CenterTitle();
  carbonLitLog->GetXaxis()->SetRangeUser(0,TMath::Log10(700.));

  carbonLitLog->GetYaxis()->SetTitle("sigma [b]");
  carbonLitLog->GetYaxis()->CenterTitle();
  //carbonLitLog->GetYaxis()->SetTitleOffSet(1.5);
  carbonLitLog->Write();

  myfile->Close();

  return;
}
