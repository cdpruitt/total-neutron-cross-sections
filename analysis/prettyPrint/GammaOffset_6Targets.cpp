void GammaOffset_6Targets()
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

    string blankTOFHistoName = "gammaOffsetHBlank";
    string target1TOFHistoName = "gammaOffsetHTarget1";
    string target2TOFHistoName = "gammaOffsetHTarget2";
    string target3TOFHistoName = "gammaOffsetHTarget3";
    string target4TOFHistoName = "gammaOffsetHTarget4";
    string target5TOFHistoName = "gammaOffsetHTarget5";

    TOFFile->cd("lowThresholdDet");
     
    TH1D* blankTOFHisto = (TH1D*)gDirectory->Get(blankTOFHistoName.c_str());
    TH1D* target1TOFHisto = (TH1D*)gDirectory->Get(target1TOFHistoName.c_str());
    TH1D* target2TOFHisto = (TH1D*)gDirectory->Get(target2TOFHistoName.c_str());
    TH1D* target3TOFHisto = (TH1D*)gDirectory->Get(target3TOFHistoName.c_str());
    TH1D* target4TOFHisto = (TH1D*)gDirectory->Get(target4TOFHistoName.c_str());
    TH1D* target5TOFHisto = (TH1D*)gDirectory->Get(target5TOFHistoName.c_str());

    if(!blankTOFHisto)
    {
        cout << "Error: failed to find " << blankTOFHistoName << " in file " << TOFFileName << endl;
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
    blankTOFHisto->SetLineColor(kBlack);
    blankTOFHisto->SetLineWidth(5);
    blankTOFHisto->SetLineStyle(0);
    blankTOFHisto->SetMarkerColor(kBlack);

    target1TOFHisto->SetLineColor(kRed+2);
    target1TOFHisto->SetLineWidth(5);
    target1TOFHisto->SetLineStyle(0);
    target1TOFHisto->SetMarkerColor(kRed+2);

    target2TOFHisto->SetLineColor(kRed-9);
    target2TOFHisto->SetLineWidth(5);
    target2TOFHisto->SetLineStyle(0);
    target2TOFHisto->SetMarkerColor(kRed-9);

    target3TOFHisto->SetLineColor(kMagenta-10);
    target3TOFHisto->SetLineWidth(5);
    target3TOFHisto->SetLineStyle(0);
    target3TOFHisto->SetMarkerColor(kMagenta-10);

    target4TOFHisto->SetLineColor(kMagenta-6);
    target4TOFHisto->SetLineWidth(5);
    target4TOFHisto->SetLineStyle(0);
    target4TOFHisto->SetMarkerColor(kMagenta-6);

    target5TOFHisto->SetLineColor(kMagenta-2);
    target5TOFHisto->SetLineWidth(5);
    target5TOFHisto->SetLineStyle(0);
    target5TOFHisto->SetMarkerColor(kMagenta-2);

    
    // X-axis parameters
    blankTOFHisto->GetXaxis()->SetTitle("Gamma time of flight deviation (ns)");
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
    blankTOFHisto->GetYaxis()->SetTitleOffset(1.1);
    blankTOFHisto->GetYaxis()->CenterTitle();

    blankTOFHisto->GetYaxis()->SetLabelOffset(0.01);
    blankTOFHisto->GetYaxis()->SetLabelSize(0.05);

    blankTOFHisto->GetYaxis()->SetLabelFont(2);
    blankTOFHisto->GetYaxis()->SetNdivisions(10);
    blankTOFHisto->GetYaxis()->SetTickLength(0.02);

    blankTOFHisto->Draw();
    target1TOFHisto->Draw("same");
    target2TOFHisto->Draw("same");
    target3TOFHisto->Draw("same");
    target4TOFHisto->Draw("same");
    target5TOFHisto->Draw("same");

    //blankTOFHisto->GetYaxis()->SetRangeUser(0,12000);
    blankTOFHisto->GetXaxis()->SetRange(300, 700);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Sn");
    //latex.DrawLatex(0.32,0.4,"C");

    //gPad->SetLogy(1);

    // Define legend format and contents
    TLegend *legend = new TLegend(0.21,0.63,0.35,0.93);
    //legend->SetHeader("data","C");
    legend->AddEntry(blankTOFHisto,"blank","l");
    legend->AddEntry(target1TOFHisto,"{}^{nat}Pb","l");
    legend->AddEntry(target2TOFHisto,"{}^{nat}C","l");
    legend->AddEntry(target3TOFHisto,"{}^{112}Sn","l");
    legend->AddEntry(target4TOFHisto,"{}^{nat}Sn","l");
    legend->AddEntry(target5TOFHisto,"{}^{124}Sn","l");

    legend->Draw();

    //TOFFile->Close();
}
