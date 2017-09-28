void TOF_blankCSnPb_zoom()
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

    string blankTOFHistoName = "blankTOF";
    string PbTOFHistoName = "target1TOF";
    string CTOFHistoName = "target2TOF";
    string SnTOFHistoName = "target3TOF";

    TOFFile->cd("lowThresholdDet");
     
    TH1D* blankTOFHisto = (TH1D*)gDirectory->Get(blankTOFHistoName.c_str());
    TH1D* PbTOFHisto = (TH1D*)gDirectory->Get(PbTOFHistoName.c_str());
    TH1D* CTOFHisto = (TH1D*)gDirectory->Get(CTOFHistoName.c_str());
    TH1D* SnTOFHisto = (TH1D*)gDirectory->Get(SnTOFHistoName.c_str());

    if(!blankTOFHisto)
    {
        cout << "Error: failed to find " << blankTOFHistoName << " in file " << TOFFileName << endl;
        exit(1);
    }

    blankTOFHisto->Rebin(10);
    PbTOFHisto->Rebin(10);
    CTOFHisto->Rebin(10);
    SnTOFHisto->Rebin(10);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // Set histo point and line characteristics
    blankTOFHisto->SetLineColor(kBlack);
    blankTOFHisto->SetLineWidth(5);
    blankTOFHisto->SetLineStyle(0);
    blankTOFHisto->SetMarkerColor(kBlack);

    PbTOFHisto->SetLineColor(kRed+3);
    PbTOFHisto->SetLineWidth(5);
    PbTOFHisto->SetLineStyle(0);
    PbTOFHisto->SetMarkerColor(kRed+3);

    CTOFHisto->SetLineColor(kRed-9);
    CTOFHisto->SetLineWidth(5);
    CTOFHisto->SetLineStyle(0);
    CTOFHisto->SetMarkerColor(kRed-9);

    SnTOFHisto->SetLineColor(kRed);
    SnTOFHisto->SetLineWidth(5);
    SnTOFHisto->SetLineStyle(0);
    SnTOFHisto->SetMarkerColor(kRed);

    // X-axis parameters
    blankTOFHisto->GetXaxis()->SetTitle("Time of flight (ns)");
    blankTOFHisto->GetXaxis()->SetTitleSize(0.05);
    blankTOFHisto->GetXaxis()->SetTitleFont(2);
    blankTOFHisto->GetXaxis()->SetTitleOffset(1.5);
    blankTOFHisto->GetXaxis()->CenterTitle();

    blankTOFHisto->GetXaxis()->SetLabelOffset(0.01);
    blankTOFHisto->GetXaxis()->SetLabelSize(0.05);
    blankTOFHisto->GetXaxis()->SetLabelFont(2);

    blankTOFHisto->GetXaxis()->SetNdivisions(10);
    blankTOFHisto->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    blankTOFHisto->GetYaxis()->SetTitle("Counts");
    blankTOFHisto->GetYaxis()->SetTitleSize(0.06);
    blankTOFHisto->GetYaxis()->SetTitleFont(2);
    blankTOFHisto->GetYaxis()->SetTitleOffset(1.2);
    blankTOFHisto->GetYaxis()->CenterTitle();

    blankTOFHisto->GetYaxis()->SetLabelOffset(0.01);
    blankTOFHisto->GetYaxis()->SetLabelSize(0.05);

    blankTOFHisto->GetYaxis()->SetLabelFont(2);
    blankTOFHisto->GetYaxis()->SetNdivisions(10);
    blankTOFHisto->GetYaxis()->SetTickLength(0.02);

    blankTOFHisto->Draw();
    PbTOFHisto->Draw("same");
    CTOFHisto->Draw("same");
    SnTOFHisto->Draw("same");
    
    blankTOFHisto->GetYaxis()->SetRangeUser(0,11000);
    blankTOFHisto->GetXaxis()->SetRange(200, 242);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Sn");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.75,0.7,0.9,0.9);
    //legend->SetHeader("data","C");
    legend->AddEntry(blankTOFHisto,"blank","l");
    legend->AddEntry(CTOFHisto,"{}^{nat}C","l");
    legend->AddEntry(SnTOFHisto,"{}^{nat}Sn","l");
    legend->AddEntry(PbTOFHisto,"{}^{nat}Pb","l");
    legend->Draw();

    //TOFFile->Close();
}
