{
    string fileName = "/data2/analysis/relative.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    
    string relGraphName = "longCShortC, percent";
    string relGraphSEName = "longCShortCSysErrors, percent";
        
    TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)file->Get(relGraphName.c_str());
    TGraphAsymmErrors* relGraphSE = (TGraphAsymmErrors*)file->Get(relGraphSEName.c_str());

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

    relGraphSE->SetFillColor(kBlue);
    relGraphSE->SetFillStyle(3001);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.15);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    gPad->SetFrameLineWidth(3);

    TMultiGraph* mg = new TMultiGraph();

    mg->Add(relGraph,"3l");
    mg->Add(relGraphSE, "3");

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
    mg->GetYaxis()->SetTitle("(#frac{#sigma_{l} - #sigma_{s}}{#sigma_{l} + #sigma_{s}}) (%)");
    mg->GetYaxis()->SetTitleSize(0.06);
    mg->GetYaxis()->SetTitleFont(2);
    mg->GetYaxis()->SetTitleOffset(1.0);
    mg->GetYaxis()->CenterTitle();

    mg->GetYaxis()->SetLabelOffset(0.01);
    mg->GetYaxis()->SetLabelSize(0.05);

    mg->GetYaxis()->SetLabelFont(2);
    mg->GetYaxis()->SetNdivisions(5);
    mg->GetYaxis()->SetTickLength(0.02);

    gPad->SetLogx(1);
    
    mg->GetYaxis()->SetRangeUser(-2.49,2.49);
    mg->GetXaxis()->SetLimits(3,500);

    TLine* plusOneLine = new TLine(3, 1, 500, 1);
    plusOneLine->SetLineColor(kGray+4);
    plusOneLine->SetLineWidth(4);
    plusOneLine->SetLineStyle(9);
    plusOneLine->Draw();

    TLine* minusOneLine = new TLine(3, -1, 500, -1);
    minusOneLine->SetLineColor(kGray+4);
    minusOneLine->SetLineWidth(4);
    minusOneLine->SetLineStyle(9);
    minusOneLine->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13); // align at top
    latex.DrawLatex(0.47,0.52,"{}^{nat}C");

    file->Close();
}
