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
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    //gPad->SetTicky(2);

    // Set histo point and line characteristics
    LRTimeDifferenceHisto->SetLineColor(kBlack);
    LRTimeDifferenceHisto->SetLineWidth(2);

    // X-axis parameters
    LRTimeDifferenceHisto->GetXaxis()->SetTitle("Time difference, left - right (ns)");
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
    //LRTimeDifferenceHisto->GetYaxis()->SetTitle("ns)");
    //LRTimeDifferenceHisto->GetYaxis()->SetTitleSize(0.05);
    //LRTimeDifferenceHisto->GetYaxis()->SetTitleFont(2);
    //LRTimeDifferenceHisto->GetYaxis()->SetTitleOffset(1.5);
    //LRTimeDifferenceHisto->GetYaxis()->CenterTitle();

    LRTimeDifferenceHisto->GetYaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceHisto->GetYaxis()->SetLabelSize(0.05);
    LRTimeDifferenceHisto->GetYaxis()->SetLabelFont(2);

    //LRTimeDifferenceHisto->GetYaxis()->SetNdivisions(10);
    //LRTimeDifferenceHisto->GetYaxis()->SetTickLength(0.02);

    LRTimeDifferenceHisto->GetXaxis()->SetRangeUser(-3,3);
    //LRTimeDifferenceHisto->GetYaxis()->SetRangeUser(35,55);

    LRTimeDifferenceHisto->Draw();

    gPad->SetLogy();

    LRTimeDifferenceHisto->Fit("gaus","","",-2,2);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(33); // align at top
    latex.SetTextColor(kRed);
    latex.DrawLatex(0.94,0.52,"#mu_{fit} = -0.0929 ns");
    latex.DrawLatex(0.94,0.44,"FWHM_{fit} = 0.520 ns");

    TLine *meanLine = new TLine(-0.0929, 0, -0.0929, 11000);
    meanLine->SetLineStyle(7);
    meanLine->SetLineWidth(2);
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
