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

    string fileName = "/data2/analysis/66/0000/gatedHistos.root";
    TFile* file = new TFile(fileName.c_str(),"READ");
    file->cd("lowThresholdDet");

    string histoName = "time since last event";
    TH1D* histo = (TH1D*)gDirectory->Get(histoName.c_str());

    if(!histo)
    {
        cerr << "Error: couldn't find "
            << histoName  << " in "
            << fileName << ". Exiting..." << endl;
        exit(1);
    }

    // Set histo point and line characteristics
    histo->SetMarkerColor(kBlack);

    // X-axis parameters
    histo->GetXaxis()->SetTitle("Time difference between consecutive events");
    histo->GetXaxis()->SetTitleSize(0.04);
    histo->GetXaxis()->SetTitleFont(2);
    histo->GetXaxis()->SetTitleOffset(1.5);
    histo->GetXaxis()->CenterTitle();

    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetXaxis()->SetLabelSize(0.04);
    histo->GetXaxis()->SetLabelFont(2);

    histo->GetXaxis()->SetNdivisions(6);

    // Y-axis parameters
    //histo->GetYaxis()->SetTitle("#sigma_{exp}-#sigma_{lit} (RMS, in %)");
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleFont(2);
    histo->GetYaxis()->SetTitleOffset(1.2);
    histo->GetYaxis()->CenterTitle();

    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetLabelFont(2);

    histo->GetXaxis()->SetRangeUser(200,250);
    histo->GetYaxis()->SetRangeUser(0,7000);

    histo->Draw("");

    TF1* timeDiffFit = new TF1("logisticFit","(1627.16-1.75548*x)/(1+exp(-[0]*(x-[1])))",
            200, 250);
    timeDiffFit->SetParameter(1,200);
    histo->Fit("logisticFit","", "", 200, 250);
    timeDiffFit->SetLineWidth(3);

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
