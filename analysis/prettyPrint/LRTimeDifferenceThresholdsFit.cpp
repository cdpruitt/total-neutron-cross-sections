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
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.20);
    //gPad->SetTicky(2);

    double x[5] = {5, 10, 20, 40, 80};
    double y[5] = {0.67846, 0.47541, 0.40246, 0.36975, 0.35233};
    TGraph* LRTimeDifferenceFitGraph = new TGraph(5,x,y);

    // Set histo point and line characteristics
    LRTimeDifferenceFitGraph->SetMarkerColor(kBlack);

    // X-axis parameters
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitle("Energy of Event (MeV)");
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitleSize(0.04);
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitleFont(2);
    LRTimeDifferenceFitGraph->GetXaxis()->SetTitleOffset(1.5);
    LRTimeDifferenceFitGraph->GetXaxis()->CenterTitle();

    LRTimeDifferenceFitGraph->GetXaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceFitGraph->GetXaxis()->SetLabelSize(0.04);
    LRTimeDifferenceFitGraph->GetXaxis()->SetLabelFont(2);

    LRTimeDifferenceFitGraph->GetXaxis()->SetNdivisions(5);

    // Y-axis parameters
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitle("Time Difference FWHM (ns)");
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitleSize(0.04);
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitleFont(2);
    LRTimeDifferenceFitGraph->GetYaxis()->SetTitleOffset(1.7);
    LRTimeDifferenceFitGraph->GetYaxis()->CenterTitle();

    LRTimeDifferenceFitGraph->GetYaxis()->SetLabelOffset(0.01);
    LRTimeDifferenceFitGraph->GetYaxis()->SetLabelSize(0.04);
    LRTimeDifferenceFitGraph->GetYaxis()->SetLabelFont(2);

    LRTimeDifferenceFitGraph->GetXaxis()->SetRangeUser(0,85);
    LRTimeDifferenceFitGraph->GetYaxis()->SetRangeUser(0.33,0.70);

    LRTimeDifferenceFitGraph->Draw("AP");

    TF1* fitToFWHM = new TF1("fitToFWHM","[0]+[1]/x+[2]/(x*x)",4.4,84);
    LRTimeDifferenceFitGraph->Fit("fitToFWHM","","",4.4,84);
    fitToFWHM->SetLineWidth(3);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top

    latex.SetTextColor(kGray+2);
    latex.DrawLatex(0.20,0.28,"lim = 0.34 ns");

    latex.SetTextSize(0.022);
    latex.DrawLatex(0.20,0.25,"x#rightarrow#infty");

    TArrow *arrow = new TArrow(17, 0.3550, 19, 0.3460, 0.015, "|>");
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

    TLine *asymptote = new TLine(0, 0.3406, 85, 0.3406);
    asymptote->SetLineStyle(7);
    asymptote->SetLineWidth(4);
    asymptote->SetLineColor(kGray+2);
    asymptote->Draw();

}
