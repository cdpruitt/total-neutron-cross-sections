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
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.12);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    gPad->SetFrameLineWidth(3);

    double x[9] =
       {26.97, 27.00, 27.03, 27.06, 27.09,
        27.12, 27.15, 27.18, 27.21};

    double y[9] =
       {4.46156, 3.71723, 2.88004, 2.03289, 1.68967,
        2.03278, 2.68672, 3.39422, 4.30685};

    TGraph* DistanceStudyFitGraph = new TGraph(9,x,y);

    // Set histo point and line characteristics
    DistanceStudyFitGraph->SetMarkerColor(kBlack);

    // X-axis parameters
    DistanceStudyFitGraph->GetXaxis()->SetTitle("TOF Detector Distance (m)");
    DistanceStudyFitGraph->GetXaxis()->SetTitleSize(0.04);
    DistanceStudyFitGraph->GetXaxis()->SetTitleFont(2);
    DistanceStudyFitGraph->GetXaxis()->SetTitleOffset(1.5);
    DistanceStudyFitGraph->GetXaxis()->CenterTitle();

    DistanceStudyFitGraph->GetXaxis()->SetLabelOffset(0.01);
    DistanceStudyFitGraph->GetXaxis()->SetLabelSize(0.05);
    DistanceStudyFitGraph->GetXaxis()->SetLabelFont(2);

    DistanceStudyFitGraph->GetXaxis()->SetNdivisions(6);

    // Y-axis parameters
    DistanceStudyFitGraph->GetYaxis()->SetTitle("#sigma_{exp}-#sigma_{lit} RMS (%)");
    DistanceStudyFitGraph->GetYaxis()->SetTitleSize(0.05);
    DistanceStudyFitGraph->GetYaxis()->SetTitleFont(2);
    DistanceStudyFitGraph->GetYaxis()->SetTitleOffset(1.2);
    DistanceStudyFitGraph->GetYaxis()->CenterTitle();

    DistanceStudyFitGraph->GetYaxis()->SetLabelOffset(0.01);
    DistanceStudyFitGraph->GetYaxis()->SetLabelSize(0.05);
    DistanceStudyFitGraph->GetYaxis()->SetLabelFont(2);

    DistanceStudyFitGraph->GetYaxis()->SetNdivisions(4);

    DistanceStudyFitGraph->GetXaxis()->SetRangeUser(26.94, 27.24);
    DistanceStudyFitGraph->GetYaxis()->SetRangeUser(1.2,4.99);

    DistanceStudyFitGraph->Draw("AP");
    //gPad->SetGrid(1,0);

    //DistanceStudyFitGraph->GetXaxis()->ChangeLabel(0,-1,-1,-1,-1,-1,"27.00");

    TF1* fit = new TF1("fit","[0]*(x-[1])*(x-[1])*(x-[1])*(x-[1])+[2]*(x-[1])*(x-[1])+[3]",26.95,27.25);
    fit->SetParameter(1, 27.09);

    DistanceStudyFitGraph->Fit("fit", "", "", 26.97, 27.21);
    fit->SetLineWidth(3);

    /*TLatex latex;
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
    */

    /*TLine *line = new TLine(27.0939, 1.2, 27.0939, 5.0);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw();
    */
}
