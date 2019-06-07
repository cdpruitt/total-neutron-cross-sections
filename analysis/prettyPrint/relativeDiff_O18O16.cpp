{
    string fileName = "/data2/analysis/relative.root";
    string ramsauerFileName = "../../theory/ramsauer.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(), "READ");
    
    string relGraphName = "O18O16, percent";
    string relGraphSEName = "O18O16SysErrors, percent";

    string SARelDiffGraphName = "RelDiff18_16";
    string RamsauerRelDiffGraphName = "RelDiffRamsauer18_16";
        
    TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)file->Get(relGraphName.c_str());
    TGraphAsymmErrors* relGraphSE = (TGraphAsymmErrors*)file->Get(relGraphSEName.c_str());
    TGraph* SARelDiffGraph = (TGraph*)ramsauerFile->Get(SARelDiffGraphName.c_str());
    TGraph* RamsauerRelDiffGraph = (TGraph*)ramsauerFile->Get(RamsauerRelDiffGraphName.c_str());

    TStyle* style = (TStyle*)gROOT->FindObject("graphStyle");

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

    // Set graph point and line characteristics
    relGraph->SetLineColor(kRed);
    relGraph->SetLineWidth(5);
    relGraph->SetLineStyle(0);
    relGraph->SetMarkerColor(kRed);
    relGraph->SetFillColor(kRed);
    relGraph->SetFillStyle(3002);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.2);
    gPad->SetTicky(2);

    SARelDiffGraph->SetLineStyle(9);
    SARelDiffGraph->SetLineWidth(3);
    SARelDiffGraph->SetLineColor(kBlack);

    RamsauerRelDiffGraph->SetLineStyle(7);
    RamsauerRelDiffGraph->SetLineWidth(5);
    RamsauerRelDiffGraph->SetLineColor(kGray+2);

    relGraphSE->SetFillColor(kBlue);
    relGraphSE->SetFillStyle(3002);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    TMultiGraph* mg = new TMultiGraph();

    mg->Add(relGraph,"3l");
    mg->Add(relGraphSE, "3");
    mg->Add(SARelDiffGraph, "l");
    //mg->Add(RamsauerRelDiffGraph, "l");

    mg->Draw("al");

    // X-axis parameters
    mg->GetXaxis()->SetTitle("Energy (MeV)");
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetTitleFont(2);
    mg->GetXaxis()->SetTitleOffset(1.4);
    mg->GetXaxis()->CenterTitle();

    mg->GetXaxis()->SetLabelOffset(0.01);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetXaxis()->SetLabelFont(2);

    mg->GetXaxis()->SetNdivisions(10);
    mg->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    mg->GetYaxis()->SetTitle("(#frac{#sigma_{18} - #sigma_{16}}{#sigma_{18} + #sigma_{16}}) [%]");
    mg->GetYaxis()->SetTitleSize(0.06);
    mg->GetYaxis()->SetTitleFont(2);
    mg->GetYaxis()->SetTitleOffset(1.0);
    mg->GetYaxis()->CenterTitle();

    mg->GetYaxis()->SetLabelOffset(0.01);
    mg->GetYaxis()->SetLabelSize(0.05);

    mg->GetYaxis()->SetLabelFont(2);
    mg->GetYaxis()->SetNdivisions(10);
    mg->GetYaxis()->SetTickLength(0.02);

    mg->GetYaxis()->SetRangeUser(-5,7.5);
    mg->GetXaxis()->SetLimits(8,600);

    gPad->SetLogx(1);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.47,0.52,"Ni");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.15, 0.83, 0.45, 0.95);
    legend->SetNColumns(2);
    legend->AddEntry(relGraph,"Exp data, sys + stat   ","f");
    legend->AddEntry(SARelDiffGraph,"SAS, r #alpha A^{1/3} ","l");
    legend->AddEntry(relGraphSE,"Exp data, sys only   ","f");
    //legend->AddEntry(RamsauerRelDiffGraph,"Ramsauer","l");
    legend->Draw();

    //TLine* zeroLine = new TLine(0, 0, 600, 0);
    //zeroLine->SetLineColor(kBlack);
    //zeroLine->SetLineWidth(3);
    //zeroLine->SetLineStyle(9);
    //zeroLine->Draw();

    /*TLine* SixthLine = new TLine(0, 0, 600, 0);
    zeroLine->SetLineColor(kBlack);
    zeroLine->SetLineWidth(3);
    zeroLine->SetLineStyle(9);
    zeroLine->Draw();

    TLine* ThirdLine = new TLine(0, 0, 600, 0);
    zeroLine->SetLineColor(kBlack);
    zeroLine->SetLineWidth(3);
    zeroLine->SetLineStyle(9);
    zeroLine->Draw();
    */

    file->Close();
}
