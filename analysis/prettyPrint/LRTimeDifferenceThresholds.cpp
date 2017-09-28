{
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
    style->SetHistLineWidth(4);
    style->SetHistLineColor(kBlack);
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

    string file0Name = "/data1/analysis/44/0000/detTimeCheck_threshold0.root";
    string file500Name = "/data1/analysis/44/0000/detTimeCheck_threshold500.root";
    string file1500Name = "/data1/analysis/44/0000/detTimeCheck_threshold1500.root";
    string file4500Name = "/data1/analysis/44/0000/detTimeCheck_threshold4500.root";

    TFile* file0 = new TFile(file0Name.c_str(),"READ");
    TFile* file500 = new TFile(file500Name.c_str(),"READ");
    TFile* file1500 = new TFile(file1500Name.c_str(),"READ");
    TFile* file4500 = new TFile(file4500Name.c_str(),"READ");

    string LRTimeDifferenceName = "LR time difference";

    TH1D* LRTimeDifferenceHistoT0 = (TH1D*)file0->Get(LRTimeDifferenceName.c_str());
    TH1D* LRTimeDifferenceHistoT500 = (TH1D*)file500->Get(LRTimeDifferenceName.c_str());
    TH1D* LRTimeDifferenceHistoT1500 = (TH1D*)file1500->Get(LRTimeDifferenceName.c_str());
    TH1D* LRTimeDifferenceHistoT4500 = (TH1D*)file4500->Get(LRTimeDifferenceName.c_str());

    if(
            !LRTimeDifferenceHistoT0 ||
            !LRTimeDifferenceHistoT500 ||
            !LRTimeDifferenceHistoT1500 ||
            !LRTimeDifferenceHistoT4500
      )

    {
        cerr << "Error: couldn't find "
            << LRTimeDifferenceName << " in one of time check threshold ROOT files."
            << ". Exiting..." << endl;
        exit(1);
    }



    // Set histo point and line characteristics
    LRTimeDifferenceHistoT0->SetLineColor(kBlack);
    LRTimeDifferenceHistoT500->SetLineColor(kGray+3);
    LRTimeDifferenceHistoT1500->SetLineColor(kGray+2);
    LRTimeDifferenceHistoT4500->SetLineColor(kGray+1);

    // X-axis parameters
    LRTimeDifferenceHistoT0->GetXaxis()->SetTitle("Time difference, left - right (ns)");
    LRTimeDifferenceHistoT0->GetXaxis()->SetTitleSize(0.05);
    LRTimeDifferenceHistoT0->GetXaxis()->SetTitleFont(2);
    LRTimeDifferenceHistoT0->GetXaxis()->SetTitleOffset(1.5);
    LRTimeDifferenceHistoT0->GetXaxis()->CenterTitle();

    LRTimeDifferenceHistoT0->GetXaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceHistoT0->GetXaxis()->SetLabelSize(0.05);
    LRTimeDifferenceHistoT0->GetXaxis()->SetLabelFont(2);

    LRTimeDifferenceHistoT0->GetXaxis()->SetNdivisions(6);
    //LRTimeDifferenceHistoT0->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    //LRTimeDifferenceHistoT0->GetYaxis()->SetTitle("ns)");
    //LRTimeDifferenceHistoT0->GetYaxis()->SetTitleSize(0.05);
    //LRTimeDifferenceHistoT0->GetYaxis()->SetTitleFont(2);
    //LRTimeDifferenceHistoT0->GetYaxis()->SetTitleOffset(1.5);
    //LRTimeDifferenceHistoT0->GetYaxis()->CenterTitle();

    LRTimeDifferenceHistoT0->GetYaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceHistoT0->GetYaxis()->SetLabelSize(0.05);
    LRTimeDifferenceHistoT0->GetYaxis()->SetLabelFont(2);

    //LRTimeDifferenceHistoT0->GetYaxis()->SetNdivisions(10);
    //LRTimeDifferenceHistoT0->GetYaxis()->SetTickLength(0.02);

    LRTimeDifferenceHistoT0->GetXaxis()->SetRangeUser(-1,1);
    //LRTimeDifferenceHistoT0->GetYaxis()->SetRangeUser(35,55);

    LRTimeDifferenceHistoT0->Draw();
    LRTimeDifferenceHistoT500->Draw("same");
    LRTimeDifferenceHistoT1500->Draw("same");
    LRTimeDifferenceHistoT4500->Draw("same");

    //gPad->SetLogy();

    TF1* fT0 = new TF1("fT0","gaus",-1,1);
    TF1* fT500 = new TF1("fT500","gaus",-1,1);
    TF1* fT1500 = new TF1("fT1500","gaus",-1,1);
    TF1* fT4500 = new TF1("fT4500","gaus",-1,1);

    fT0->SetLineWidth(3);
    fT500->SetLineWidth(3);
    fT1500->SetLineWidth(3);
    fT4500->SetLineWidth(3);

    /*fT0->SetLineStyle(7);
    fT500->SetLineStyle(7);
    fT1500->SetLineStyle(7);
    fT4500->SetLineStyle(7);
    */

    fT0->SetLineColor(kRed);
    fT500->SetLineColor(kRed);
    fT1500->SetLineColor(kRed);
    fT4500->SetLineColor(kRed);

    LRTimeDifferenceHistoT0->Fit("fT0","","",-1,1);
    LRTimeDifferenceHistoT500->Fit("fT500","","",-1,1);
    LRTimeDifferenceHistoT1500->Fit("fT1500","","",-1,1);
    LRTimeDifferenceHistoT4500->Fit("fT4500","","",-1,1);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top

    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.7,0.81,"Q_{sum}>300");

    TArrow *arrowT0 = new TArrow(0.35, 9500, -0.03, 10700, 0.015, "|>");
    arrowT0->SetAngle(30);
    arrowT0->SetLineWidth(1);
    arrowT0->SetLineColor(kBlack);
    arrowT0->Draw();


    latex.SetTextColor(kGray+3);
    latex.DrawLatex(0.7,0.71,"Q_{sum}>500");

    TArrow *arrowT500 = new TArrow(0.35, 7900, 0, 7800, 0.015, "|>");
    arrowT500->SetAngle(30);
    arrowT500->SetLineWidth(1);
    arrowT500->SetLineColor(kGray+3);
    arrowT500->SetFillColor(kGray+3);
    arrowT500->Draw();


    latex.SetTextColor(kGray+2);
    latex.DrawLatex(0.7,0.61,"Q_{sum}>1500");

    TArrow *arrowT1500 = new TArrow(0.35, 6200, 0, 4700, 0.015, "|>");
    arrowT1500->SetAngle(30);
    arrowT1500->SetLineWidth(1);
    arrowT1500->SetLineColor(kGray+2);
    arrowT1500->SetFillColor(kGray+2);
    arrowT1500->Draw();


    latex.SetTextColor(kGray+1);
    latex.DrawLatex(0.7,0.51,"Q_{sum}>4500");

    TArrow *arrowT4500 = new TArrow(0.35, 4700, 0.03, 1600, 0.015, "|>");
    arrowT4500->SetAngle(30);
    arrowT4500->SetLineWidth(1);
    arrowT4500->SetLineColor(kGray+1);
    arrowT4500->SetFillColor(kGray+1);
    arrowT4500->Draw();

    //latex.DrawLatex(0.94,0.52,"#mu_{fit} = -0.0929 ns");
    //latex.DrawLatex(0.94,0.44,"FWHM_{fit} = 0.520 ns");

    /*TLine *meanLine = new TLine(-0.0929, 0, -0.0929, 11000);
    meanLine->SetLineStyle(7);
    meanLine->SetLineWidth(2);
    meanLine->SetLineColor(kGray+2);
    meanLine->Draw();
    */

    //gPad->SetLogy(1);

    // Define legend format and contents
    /*TLegend *legend = new TLegend(0.60,0.60,0.9,0.75);
    legend->SetHeader("Integrated charge threshold");
    legend->AddEntry(LRTimeDifferenceHistoT0,">300 ","l");
    legend->AddEntry(LRTimeDifferenceHistoT500,">500 ","l");
    legend->AddEntry(LRTimeDifferenceHistoT1500,">1500 ","l");
    legend->AddEntry(LRTimeDifferenceHistoT4500,">4500 ","l");
    */

    //legend->Draw();
}
