{

    TStyle* style = (TStyle*)gROOT->FindObject("histoStyle");

    if(!style)      
    {
        style = new TStyle("histoStyle","histoStyle");
    }

    TCanvas* c = new TCanvas("c1","",1000,1000);

    style->SetOptStat(0);
    style->SetOptTitle(0);    
    style->SetPalette(1,0);
    style->SetCanvasColor(10);      
    style->SetCanvasBorderMode(0);    
    style->SetFrameLineWidth(3);
    style->SetFrameFillColor(10);
    style->SetPadColor(10);
    style->SetHistLineWidth(3);
    style->SetHistLineColor(kBlack);
    style->SetMarkerSize(2);
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
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.15);
    //gPad->SetTicky(2);

    double x[4] = {300, 500, 1500, 4500};
    double y[4] = {0.2208, 0.1852, 0.1595, 0.1483};
    TGraph* LRTimeDifferenceFitGraph = new TGraph(4,x,y);

    // Set histo point and line characteristics
    LRTimeDifferenceFitGraph->SetMarkerColor(kBlack);

    // X-axis parameters
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitle("Q_{sum} threshold");
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitleSize(0.04);
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitleFont(2);
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitleOffset(1.5);
    LRTimeDifferenceFitGraph->GetXaxis()->CenterTitle();

    LRTimeDifferenceFitGraph->GetXaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceFitGraph->GetXaxis()->SetLabelSize(0.04);
    LRTimeDifferenceFitGraph->GetXaxis()->SetLabelFont(2);

    LRTimeDifferenceFitGraph->GetXaxis()->SetNdivisions(5);

    // Y-axis parameters
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitle("FWHM of fit to time difference (ns)");
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitleSize(0.04);
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitleFont(2);
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitleOffset(1.5);
    LRTimeDifferenceFitGraph->GetYaxis()->CenterTitle();

    LRTimeDifferenceFitGraph->GetYaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceFitGraph->GetYaxis()->SetLabelSize(0.05);
    LRTimeDifferenceFitGraph->GetYaxis()->SetLabelFont(2);

    LRTimeDifferenceFitGraph->GetXaxis()->SetRangeUser(0,5000);
    LRTimeDifferenceFitGraph->GetYaxis()->SetRangeUser(0.135,0.25);

    LRTimeDifferenceFitGraph->Draw("AP");

    TF1* fitToFWHM = new TF1("fitToFWHM","[0]+[1]/x",230,5000);
    LRTimeDifferenceFitGraph->Fit("fitToFWHM","","",230,5000);
    fitToFWHM->SetLineWidth(3);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top

    latex.SetTextColor(kGray+2);
    latex.DrawLatex(0.2,0.28,"lim = 0.143 ns");

    latex.SetTextSize(0.022);
    latex.DrawLatex(0.2,0.25,"x#rightarrow#infty");

    TArrow *arrow = new TArrow(1020, 0.1485, 1100, 0.1445, 0.015, "|>");
    arrow->SetAngle(30);
    arrow->SetLineWidth(1);
    arrow->SetLineColor(kGray+2);
    arrow->SetFillColor(kGray+2);
    arrow->Draw();

    //latex.SetTextColor(kBlack);
    //latex.DrawLatex(0.23,0.78,"300");

    //TLine *line300 = new TLine(300, 0.21, 300, 0.23);
    //line300->SetLineStyle(7);
    //line300->SetLineWidth(2);
    //line300->SetLineColor(kBlack);
    //line300->Draw();

    TLine *asymptote = new TLine(0, 0.1429, 4900, 0.1429);
    asymptote->SetLineStyle(7);
    asymptote->SetLineWidth(4);
    asymptote->SetLineColor(kGray+2);
    asymptote->Draw();

}
