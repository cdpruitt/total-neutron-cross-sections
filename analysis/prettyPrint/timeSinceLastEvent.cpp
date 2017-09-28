void timeSinceLastEvent()
{
    TStyle * style = (TStyle*)gROOT->FindObject("histoStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* c = new TCanvas("c1","",1200,1200);

    style->SetOptStat(0);
    style->SetOptTitle(0);    
    style->SetPalette(1,0);
    style->SetCanvasColor(10);      
    style->SetCanvasBorderMode(0);    
    style->SetFrameLineWidth(3);
    style->SetFrameFillColor(10);
    style->SetPadColor(10);
    style->SetHistLineWidth(3);
    style->SetHistLineColor(kBlue);
    style->SetMarkerSize(0.9);
    style->SetMarkerStyle(8);
    style->SetFuncWidth(3);
    style->SetFuncColor(kRed);
    style->SetLabelColor(kBlack,"xyz");
    style->SetTitleSize(0.06,"xyz");
    style->SetTitleFillColor(10);
    style->SetTitleTextColor(kBlack);
    style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    string TOFFileName = "/data2/analysis/66/0001/histos.root";

    TFile* TOFFile = new TFile(TOFFileName.c_str(),"READ");

    if(!TOFFile)
    {
        cerr << "Error: failed to open " << TOFFileName << endl;
        exit(1);
    }

    string timeSinceHistoName = "time since last event";

    TOFFile->cd("lowThresholdDet");
     
    TH1D* timeSinceHisto = (TH1D*)gDirectory->Get(timeSinceHistoName.c_str());
    
    if(!timeSinceHisto)
    {
        cout << "Error: failed to find " << timeSinceHistoName << " in file " << TOFFileName << endl;
        exit(1);
    }

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // Set histo point and line characteristics
    timeSinceHisto->SetLineColor(kBlack);
    timeSinceHisto->SetLineWidth(5);
    timeSinceHisto->SetLineStyle(0);
    timeSinceHisto->SetMarkerColor(kBlack);

    // X-axis parameters
    timeSinceHisto->GetXaxis()->SetTitle("Consecutive event #Delta time (ns)");
    timeSinceHisto->GetXaxis()->SetTitleSize(0.05);
    timeSinceHisto->GetXaxis()->SetTitleFont(2);
    timeSinceHisto->GetXaxis()->SetTitleOffset(1.5);
    timeSinceHisto->GetXaxis()->CenterTitle();

    timeSinceHisto->GetXaxis()->SetLabelOffset(0.01);
    timeSinceHisto->GetXaxis()->SetLabelSize(0.05);
    timeSinceHisto->GetXaxis()->SetLabelFont(2);

    timeSinceHisto->GetXaxis()->SetNdivisions(10);
    timeSinceHisto->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    timeSinceHisto->GetYaxis()->SetTitle("Counts");
    timeSinceHisto->GetYaxis()->SetTitleSize(0.06);
    timeSinceHisto->GetYaxis()->SetTitleFont(2);
    timeSinceHisto->GetYaxis()->SetTitleOffset(1.3);
    timeSinceHisto->GetYaxis()->CenterTitle();

    timeSinceHisto->GetYaxis()->SetLabelOffset(0.01);
    timeSinceHisto->GetYaxis()->SetLabelSize(0.05);

    timeSinceHisto->GetYaxis()->SetLabelFont(2);
    timeSinceHisto->GetYaxis()->SetNdivisions(10);
    timeSinceHisto->GetYaxis()->SetTickLength(0.02);

    timeSinceHisto->Draw();
        
    timeSinceHisto->GetYaxis()->SetRangeUser(0,4000);
    //timeSinceHisto->GetXaxis()->SetRange(2000, 2420);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Sn");
    //latex.DrawLatex(0.32,0.4,"C");

    //gPad->SetLogy(1);

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.75,0.7,0.9,0.9);
    //legend->SetHeader("data","C");
    //legend->AddEntry(timeSinceHisto,"blank","l");
    //legend->Draw();

    //TOFFile->Close();
}
