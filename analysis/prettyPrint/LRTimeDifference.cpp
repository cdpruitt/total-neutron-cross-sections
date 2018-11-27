{

    string fileName = "/data1/analysis/44/0000/detTimeCheck.root";

    TFile* file = new TFile(fileName.c_str(),"READ");

    file->ls();

    string LRTimeDifferenceName = "LR time difference";
    TH1D* LRTimeDifferenceHisto = (TH1D*)file->Get(LRTimeDifferenceName.c_str());

    if(!LRTimeDifferenceHisto)
    {
        cerr << "Error: couldn't find "
            << LRTimeDifferenceName << " in "
            << fileName << ". Exiting..." << endl;
        exit(1);
    }

    TStyle* style = (TStyle*)gROOT->FindObject("histoStyle");

    if(!style)      
    {
        style = new TStyle("histoStyle","histoStyle");
    }

    TCanvas* c = new TCanvas("c1","",1000,1000);

    style->SetOptStat(0);
    style->SetOptFit(0);
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

    gROOT->SetStyle("histoStyle");
    gROOT->ForceStyle();

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    //gPad->SetTicky(2);

    // Set histo point and line characteristics
    LRTimeDifferenceHisto->SetLineColor(kBlack);
    LRTimeDifferenceHisto->SetLineWidth(3);

    // X-axis parameters
    LRTimeDifferenceHisto->GetXaxis()->SetTitle("L-R Time Difference [ns]");
    LRTimeDifferenceHisto->GetXaxis()->SetTitleSize(0.05);
    LRTimeDifferenceHisto->GetXaxis()->SetTitleFont(2);
    LRTimeDifferenceHisto->GetXaxis()->SetTitleOffset(1.5);
    LRTimeDifferenceHisto->GetXaxis()->CenterTitle();

    LRTimeDifferenceHisto->GetXaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceHisto->GetXaxis()->SetLabelSize(0.05);
    LRTimeDifferenceHisto->GetXaxis()->SetLabelFont(2);

    LRTimeDifferenceHisto->GetXaxis()->SetNdivisions(6);
    //LRTimeDifferenceHisto->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    LRTimeDifferenceHisto->GetYaxis()->SetTitle("Counts");
    LRTimeDifferenceHisto->GetYaxis()->SetTitleSize(0.05);
    LRTimeDifferenceHisto->GetYaxis()->SetTitleFont(2);
    LRTimeDifferenceHisto->GetYaxis()->SetTitleOffset(1.7);
    LRTimeDifferenceHisto->GetYaxis()->CenterTitle();

    LRTimeDifferenceHisto->GetYaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceHisto->GetYaxis()->SetLabelSize(0.05);
    LRTimeDifferenceHisto->GetYaxis()->SetLabelFont(2);

    //LRTimeDifferenceHisto->GetYaxis()->SetNdivisions(10);
    //LRTimeDifferenceHisto->GetYaxis()->SetTickLength(0.02);

    LRTimeDifferenceHisto->GetXaxis()->SetRangeUser(-1.09,1.09);
    //LRTimeDifferenceHisto->GetYaxis()->SetRangeUser(35,55);

    LRTimeDifferenceHisto->Draw();

    TF1* fit = new TF1("fit","gaus",-2,2);
    fit->SetLineWidth(4);
    LRTimeDifferenceHisto->Fit("fit");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.SetTextAlign(33); // align at top

    latex.SetTextColor(kBlack);
    char FWHM [50];
    sprintf(FWHM, "%.3f", 2.355*fit->GetParameter(2));
    string text = "FWHM = " + string(FWHM) + " ns";
    latex.DrawLatex(0.92,0.49,text.c_str());

    char mean [50];
    sprintf(mean, "%.3f", fit->GetParameter(1));
    string muText = "Mean = " + string(mean) + " ns";
    latex.DrawLatex(0.92,0.42,muText.c_str());

    TLine *meanLine = new TLine(fit->GetParameter(1), 0, fit->GetParameter(1), fit->GetParameter(0));
    meanLine->SetLineStyle(7);
    meanLine->SetLineWidth(3);
    meanLine->SetLineColor(kGray+1);
    meanLine->Draw();

    //gPad->SetLogy(1);

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.21,0.63,0.35,0.93);
    //legend->SetHeader("data","C");
    //legend->AddEntry(LRTimeDifferenceHisto,"blank","l");
    //legend->AddEntry(target1TOFHisto,"{}^{nat}Pb","l");
    //legend->AddEntry(target2TOFHisto,"{}^{nat}C","l");
    //legend->AddEntry(target3TOFHisto,"{}^{112}Sn","l");
    //legend->AddEntry(target4TOFHisto,"{}^{nat}Sn","l");
    //legend->AddEntry(target5TOFHisto,"{}^{124}Sn","l");

    //legend->Draw();

    //TOFFile->Close();
}
