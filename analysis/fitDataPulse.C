

/*  Fitting function, six parameters
 *
 *  fitf = A * (t-t0)^n * exp[-(t-t0)/w] + [C + m*(t-t0)]
 *
 *  par[0] = A
 *  par[1] = t0
 *  par[2] = n
 *  par[3] = w
 *  par[4] = C
 *  par[5] = m
 *
 *  Hope to keep n, w, C? fixed
 *
 */

Double_t fitf(Double_t *x, Double_t *par)
{
  Double_t fitval;    // fitted value of function
  Double_t arg = 0;   // argument of exponential
  if (par[3]!=0) arg = pow((x[0]-par[1]),1)/par[3];
  fitval  = exp(-arg); 
  fitval *= par[0] * pow((x[0]-par[1]),par[2]);
  fitval += par[4] + par[5]*(x[0]-par[1]);
  if (x[0]<par[1]) fitval = par[4];
  return fitval;
}



// Main fitting function
void fitDataPulse()
{

  gROOT->Reset();
  gStyle->SetOptFit(1111);

  float ch2ns         = 5.;     // time spacing between data points
  int nptsV           = 256;    // number of voltage points 
  float neutFullScale = 256;    // voltage range +- fullscale/2
  float neutOffset    = 128;    // voltage offset

  int nparam = 6;   // number of parameters for fit
  float paramMin[6];
  float paramMax[6];
  float chisqMax = 1.5;
  float chisqThresh = 0.4;

  int npulses = 200; // number of pulses to analyze
  int peakTime[200][2];
  int peakValu[200][2];
  int nPeaks;
  int junk;

  //float parameters[nparam][npulses];
  float parameters[6][240];
  float chisq[240];

  int npts = 15;    // number of points per pulse
  //float realTime[npts];  // holds true tof
  //float time[npts];      // holds time values for plotting
  //float wave[npts];      // holds voltage values
  float realTime[15];  // holds true tof
  float time[15];      // holds time values for plotting
  float wave[15];      // holds voltage values
  float peak;

  
  float tmin = -20;
  float tmax =  tmin + ch2ns*npts; // max time on plot default set to one pulse
  float maxV =  neutFullScale/2.-neutOffset;
  float minV = -neutFullScale/2.-neutOffset;

  paramMin[0]=  0;   
  paramMax[0]= -200;
  paramMin[1]=  tmin;
  paramMax[1]=  tmax;
  paramMin[2]=  0.;
  paramMax[2]=  4.;
  paramMin[3]=  0;
  paramMax[3]=  (tmax-tmin);
  paramMin[4]=  119; // digitizer bits : 0 volts = 128
  paramMax[4]=  126;
  paramMin[5]=  0;
  paramMax[5]=  .1;

  //float pulseStart[npulses];
  //float pulseEnd[npulses];

  ifstream ifFile("pulseNeut.out");
  ostringstream titlestring;
  string title;

  TH1F *histChisq = new TH1F("chisq","",50,0,chisqMax);
  TH1F *histParam[6];
  for (int j=0; j<nparam; j++) {
    titlestring.str("");
    if (j==0)  titlestring << "A";
    if (j==1)  titlestring << "t0";
    if (j==2)  titlestring << "n";
    if (j==3)  titlestring << "w";
    if (j==4)  titlestring << "C";
    if (j==5)  titlestring << "m";
    title = titlestring.str();
    histParam[j]= new TH1F(title.c_str(),"",50,paramMin[j],paramMax[j]);
    titlestring << " value";
    title = titlestring.str();
    histParam[j]->GetXaxis()->SetTitle(title.c_str());
    histParam[j]->GetYaxis()->SetTitle("Counts");
  }


  TCanvas *pulsepad = new TCanvas("pulsepad");
  pulsepad->Divide(15,npulses/15);

  TH1F *plotData;

  for (int iP=0;iP<npulses;iP++) {

    // read number of peaks and positions
    ifFile >> nPeaks;
    for (int nP=0; nP<2; nP++) ifFile >> peakTime[iP][nP] >> peakValu[iP][nP];


    for (int i=0; i<npts; i++) {
      ifFile >> realTime[i] >> wave[i];
      if(ifFile.eof()) break;
      if(ifFile.bad()) break;

    // read pulse shape
      if (i==0) peak = wave[i];
      else if (wave[i]>peak) peak = wave[i];       
      time[i] = ch2ns*i;
    } // end loop over points in pulse

    pulsepad->cd(iP+1);  
    plotData = new TH1F("plotData","",npts-1,tmin,tmax);
    plotData->GetXaxis()->SetTitle("time (ns)");
    plotData->GetYaxis()->SetTitle("pulse (V)");
    plotData->SetStats(kFALSE);
    plotData->SetLineColor(2);
    plotData->SetMarkerStyle(7);

    for (int i=0; i<npts; i++) plotData->SetBinContent(i,wave[i]);



//------------------------ Fit pulse -------------------------//

    // create fitting object
    TF1 *func  = new TF1("func",fitf,tmin,tmax,6);
    //TF1 *func  = new TF1("func",gaus,tmin,tmax,3);

    // Set initial parameters:  a, t0, n, w, C
    float A_init  = 50;                // half the range
    float t0_init = tmin+2*tmax/npts;  // a couple points into pulse
    float n_init  = 2;                 // squared
    float w_init  = (tmax-tmin)/3;     // a third of whole pulse
    float C_init  = 122;               // background is zero
    float m_init  = 0.5;                 // background is zero
    func->SetParameters(A_init,t0_init,n_init,w_init,C_init,m_init);
    func->SetParNames("A","t0","n","w","C","m");

    // Set limits on parameter values (defined above)
    for (int j=0; j<1; j++) {
       func->SetParLimits(j,paramMin[j],paramMax[j]);
    }

    // Fix certain parameters
    func->FixParameter(2,0.5);         // for new runs
    //func->FixParameter(2,2);         // for new runs
    //func->FixParameter(2,1.457);
    //func->FixParameter(2,1.457);
    //func->FixParameter(2,1.12);
    func->FixParameter(3,10);          // for new runs, exp = 1
    //func->FixParameter(3,13.55); // for exp power = 1.
    //func->FixParameter(3,2.6);   // for exp power = 0.7
    //func->FixParameter(3,1);     // for exp power = 0.5
    //func->FixParameter(4,121);
    func->FixParameter(4,121);
    func->FixParameter(5,0);
    //func->SetParameters(A_init,t0_init,n_init,w_init,C_init);


    plotData->Fit("func","R Q ","p");


    chisq[iP] = func->GetChisquare();

    for (int j=0; j<nparam; j++) {
      parameters[j][iP] = func->GetParameter(j);
      //if (chisq[iP]<chisqThresh) histParam[j]->Fill(parameters[j][iP]);
      histParam[j]->Fill(parameters[j][iP]);
      //cout << parameters[j][iP]<<endl;
    }

    ostringstream drawParam;
    string drawPar;

    TLatex params;
    params.SetTextSize(0.02);
    for (int j=0; j<nparam; j++) {
      drawParam.str("");
      drawParam << func->GetParName(j) << " ";
      drawParam << parameters[j][iP];
      drawPar = drawParam.str();
      params.DrawLatex( 30,250-7*j,drawPar.c_str() );
    }


    histChisq->Fill(-chisq[iP]/parameters[0][iP]);
    if (-chisq[iP]/parameters[0][iP]>chisqThresh) plotData->Draw("same");
    
    //histChisq->Fill(chisq[iP]);
    //if (chisq[iP]>chisqThresh) plotData->Draw("same");

    TGraph *plotPeaks = new TGraph(nPeaks,peakTime[iP],peakValu[iP]);  
    plotPeaks->SetMarkerStyle(8);  // small dots
    plotPeaks->SetMarkerColor(4);
    plotPeaks->SetLineColor(1);
    plotPeaks->Draw("p");


  } // end loop over pulses

  TCanvas *paramPlot = new TCanvas("paramPlot");

  paramPlot->Divide(2,3);
  for (int j=0; j<nparam; j++) {
    paramPlot->cd(j+1);
    histParam[j]->Draw();
  }

  TLatex fitfunction;
  fitfunction.SetTextSize(0.08);
  fitfunction.SetTextColor(2);
  fitfunction.DrawLatex(0.1,30,
			"fitf = A*(t-t0)^n*exp[-(t-t0)/w] + [C+m*(t-t0)]");


  TCanvas *chisqPlot = new TCanvas("chisqPlot");
  histChisq->Draw();

}
