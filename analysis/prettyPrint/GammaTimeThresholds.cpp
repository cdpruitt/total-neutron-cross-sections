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
    gPad->SetLeftMargin(0.20);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.20);
    //gPad->SetTicky(2);

    string file0Name = "/data1/analysis/44/0000/detTimeCheck_0to655.root";
    string file1Name = "/data1/analysis/44/0000/detTimeCheck_655to1250.root";
    string file2Name = "/data1/analysis/44/0000/detTimeCheck_1250to2500.root";
    string file3Name = "/data1/analysis/44/0000/detTimeCheck_2500to5000.root";
    string file4Name = "/data1/analysis/44/0000/detTimeCheck_5000to10000.root";
    //string file5Name = "/data1/analysis/44/0000/detTimeCheck_10000to14000.root";

    TFile* file0 = new TFile(file0Name.c_str(),"READ");
    TFile* file1 = new TFile(file1Name.c_str(),"READ");
    TFile* file2 = new TFile(file2Name.c_str(),"READ");
    TFile* file3 = new TFile(file3Name.c_str(),"READ");
    TFile* file4 = new TFile(file4Name.c_str(),"READ");

    string LRTimeDifferenceName = "LR time difference";

    TH1D* LRTimeDifferenceHisto0 = (TH1D*)file0->Get(LRTimeDifferenceName.c_str());
    TH1D* LRTimeDifferenceHisto1 = (TH1D*)file1->Get(LRTimeDifferenceName.c_str());
    TH1D* LRTimeDifferenceHisto2 = (TH1D*)file2->Get(LRTimeDifferenceName.c_str());
    TH1D* LRTimeDifferenceHisto3 = (TH1D*)file3->Get(LRTimeDifferenceName.c_str());
    TH1D* LRTimeDifferenceHisto4 = (TH1D*)file4->Get(LRTimeDifferenceName.c_str());

    if(
            !LRTimeDifferenceHisto0 ||
            !LRTimeDifferenceHisto1 ||
            !LRTimeDifferenceHisto2 ||
            !LRTimeDifferenceHisto3 ||
            !LRTimeDifferenceHisto4
      )

    {
        cerr << "Error: couldn't find "
            << LRTimeDifferenceName << " in one of time check threshold ROOT files."
            << ". Exiting..." << endl;
        exit(1);
    }


    // Set histo point and line characteristics
    LRTimeDifferenceHisto0->SetLineColor(kRed+4);
    LRTimeDifferenceHisto1->SetLineColor(kRed+3);
    LRTimeDifferenceHisto2->SetLineColor(kRed+2);
    LRTimeDifferenceHisto3->SetLineColor(kRed+1);
    LRTimeDifferenceHisto4->SetLineColor(kRed);

    LRTimeDifferenceHisto0->SetLineWidth(4);
    LRTimeDifferenceHisto1->SetLineWidth(4);
    LRTimeDifferenceHisto2->SetLineWidth(4);
    LRTimeDifferenceHisto3->SetLineWidth(4);
    LRTimeDifferenceHisto4->SetLineWidth(4);

    // X-axis parameters
    LRTimeDifferenceHisto0->GetXaxis()->SetTitle("Left-Right Time Difference (ns)");
    LRTimeDifferenceHisto0->GetXaxis()->SetTitleSize(0.05);
    LRTimeDifferenceHisto0->GetXaxis()->SetTitleFont(2);
    LRTimeDifferenceHisto0->GetXaxis()->SetTitleOffset(1.5);
    LRTimeDifferenceHisto0->GetXaxis()->CenterTitle();

    LRTimeDifferenceHisto0->GetXaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceHisto0->GetXaxis()->SetLabelSize(0.05);
    LRTimeDifferenceHisto0->GetXaxis()->SetLabelFont(2);

    LRTimeDifferenceHisto0->GetXaxis()->SetNdivisions(6);
    //LRTimeDifferenceHisto0->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    //LRTimeDifferenceHisto0->GetYaxis()->SetTitle("ns)");
    //LRTimeDifferenceHisto0->GetYaxis()->SetTitleSize(0.05);
    //LRTimeDifferenceHisto0->GetYaxis()->SetTitleFont(2);
    //LRTimeDifferenceHisto0->GetYaxis()->SetTitleOffset(1.5);
    //LRTimeDifferenceHisto0->GetYaxis()->CenterTitle();

    LRTimeDifferenceHisto0->GetYaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceHisto0->GetYaxis()->SetLabelSize(0.05);
    LRTimeDifferenceHisto0->GetYaxis()->SetLabelFont(2);

    //LRTimeDifferenceHisto0->GetYaxis()->SetNdivisions(10);
    //LRTimeDifferenceHisto0->GetYaxis()->SetTickLength(0.02);

    LRTimeDifferenceHisto0->GetXaxis()->SetRangeUser(-1,1);
    LRTimeDifferenceHisto0->GetYaxis()->SetRangeUser(0,10500);

    LRTimeDifferenceHisto0->Draw();
    LRTimeDifferenceHisto1->Draw("same");
    LRTimeDifferenceHisto2->Draw("same");
    LRTimeDifferenceHisto3->Draw("same");
    LRTimeDifferenceHisto4->Draw("same");

    //gPad->SetLogy();

    TF1* f0 = new TF1("f0","gaus",-1,1);
    TF1* f1 = new TF1("f1","gaus",-1,1);
    TF1* f2 = new TF1("f2","gaus",-1,1);
    TF1* f3 = new TF1("f3","gaus",-1,1);
    TF1* f4 = new TF1("f4","gaus",-1,1);

    f0->SetLineWidth(3);
    f1->SetLineWidth(3);
    f2->SetLineWidth(3);
    f3->SetLineWidth(3);
    f4->SetLineWidth(3);

    /*fT0->SetLineStyle(7);
    fT500->SetLineStyle(7);
    fT1500->SetLineStyle(7);
    fT4500->SetLineStyle(7);
    */

    f0->SetLineColor(kRed);
    f1->SetLineColor(kRed);
    f2->SetLineColor(kRed);
    f3->SetLineColor(kRed);
    f4->SetLineColor(kRed);

    LRTimeDifferenceHisto0->Fit("f0","0","",-1,1);
    LRTimeDifferenceHisto1->Fit("f1","0","",-1,1);
    LRTimeDifferenceHisto2->Fit("f2","0","",-1,1);
    LRTimeDifferenceHisto3->Fit("f3","0","",-1,1);
    LRTimeDifferenceHisto4->Fit("f4","0","",-1,1);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top

    latex.SetTextColor(kRed+4);
    //latex.DrawLatex(0.65,0.60,"<2.5 MeV");

    TArrow *arrow0 = new TArrow(0.23, 5700, -0.02, 5700, 0.015, "|>");
    arrow0->SetAngle(30);
    arrow0->SetLineWidth(1);
    arrow0->SetLineColor(kRed+4);
    arrow0->SetFillColor(kRed+4);
    //arrow0->Draw();


    latex.SetTextColor(kRed+3);
    //latex.DrawLatex(0.65,0.67,"2.5-5 MeV");

    TArrow *arrow1 = new TArrow(0.23, 6700, -0.02, 7500, 0.015, "|>");
    arrow1->SetAngle(30);
    arrow1->SetLineWidth(1);
    arrow1->SetLineColor(kRed+3);
    arrow1->SetFillColor(kRed+3);
    //arrow1->Draw();


    latex.SetTextColor(kRed+2);
    //latex.DrawLatex(0.65,0.74,"5-10 MeV");

    TArrow *arrow2 = new TArrow(0.23, 7700, -0.09, 8800, 0.015, "|>");
    arrow2->SetAngle(30);
    arrow2->SetLineWidth(1);
    arrow2->SetLineColor(kRed+2);
    arrow2->SetFillColor(kRed+2);
    //arrow2->Draw();


    latex.SetTextColor(kRed+1);
    //latex.DrawLatex(0.65,0.81,"10-20 MeV");

    TArrow *arrow3 = new TArrow(0.35, 4700, 0.03, 1600, 0.015, "|>");
    arrow3->SetAngle(30);
    arrow3->SetLineWidth(1);
    arrow3->SetLineColor(kRed+1);
    arrow3->SetFillColor(kRed+1);
    //arrow3->Draw();

    latex.SetTextColor(kRed);
    //latex.DrawLatex(0.65,0.88,"20-40 MeV");

    TArrow *arrow4 = new TArrow(0.35, 4700, 0.03, 1600, 0.015, "|>");
    arrow4->SetAngle(30);
    arrow4->SetLineWidth(1);
    arrow4->SetLineColor(kRed);
    arrow4->SetFillColor(kRed);
    //arrow4->Draw();

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
    TLegend *legend = new TLegend(0.7,0.55,0.93,0.93);
    legend->SetTextSize(0.04);

    legend->SetHeader("Event Energy");

    LRTimeDifferenceHisto4->SetLineWidth(5);
    LRTimeDifferenceHisto3->SetLineWidth(5);
    LRTimeDifferenceHisto2->SetLineWidth(5);
    LRTimeDifferenceHisto1->SetLineWidth(5);
    LRTimeDifferenceHisto0->SetLineWidth(5);

    TLegendEntry* leg4 = legend->AddEntry(LRTimeDifferenceHisto4,"20-40 MeV ","l");
    TLegendEntry* leg3 = legend->AddEntry(LRTimeDifferenceHisto3,"10-20 MeV ","l");
    TLegendEntry* leg2 = legend->AddEntry(LRTimeDifferenceHisto2,"5-10 MeV ","l");
    TLegendEntry* leg1 = legend->AddEntry(LRTimeDifferenceHisto1,"2.5-5 MeV ","l");
    TLegendEntry* leg0 = legend->AddEntry(LRTimeDifferenceHisto0,"0-2.5 MeV ","l");

    leg4->SetTextAlign(22);
    leg3->SetTextAlign(22);
    leg2->SetTextAlign(22);
    leg1->SetTextAlign(22);
    leg0->SetTextAlign(22);

    legend->Draw();

    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextSize(0.045);
    header->SetTextAlign(22);

    //legend->SetEntrySeparation(2);
}
