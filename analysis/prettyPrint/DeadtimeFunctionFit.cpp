{
    TStyle* style = (TStyle*)gROOT->FindObject("histoStyle");

    if(!style)      
    {
        style = new TStyle("histoStyle","histoStyle");
    }

    TCanvas* c = new TCanvas("c1","",1000,700);

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
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.12);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    gPad->SetFrameLineWidth(3);

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
    histo->GetXaxis()->SetTitle("Time difference between consecutive events (ns)");
    histo->GetXaxis()->SetTitleSize(0.04);
    histo->GetXaxis()->SetTitleFont(2);
    histo->GetXaxis()->SetTitleOffset(1.5);
    histo->GetXaxis()->CenterTitle();

    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetXaxis()->SetLabelSize(0.04);
    histo->GetXaxis()->SetLabelFont(2);

    histo->GetXaxis()->SetNdivisions(6);

    // Y-axis parameters
    histo->GetYaxis()->SetTitle("Adjacent events");
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleFont(2);
    histo->GetYaxis()->SetTitleOffset(1.2);
    histo->GetYaxis()->CenterTitle();

    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetLabelSize(0.04);
    histo->GetYaxis()->SetLabelFont(2);

    histo->GetXaxis()->SetRangeUser(200,280);
    histo->GetYaxis()->SetRangeUser(0.1,1325);

    histo->Draw("");

    TF1* timeDiffFit = new TF1("logisticFit","([0]-[1]*x)/(1+exp(-[2]*(x-[3])))",
            200, 300);
    timeDiffFit->SetParameter(0,1625);
    timeDiffFit->SetParameter(1,2);
    timeDiffFit->SetParameter(3,220);
    histo->Fit("logisticFit","", "", 200, 300);
    timeDiffFit->SetLineWidth(5);

    double deadtimeMean = timeDiffFit->GetParameter(3);
    double deadtimeLogisticWidth = timeDiffFit->GetParameter(2);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13); // align at top

    latex.SetTextColor(kBlack);
    char mean[50];
    sprintf(mean,"Mean deadtime = %3.1f ns", deadtimeMean);
    char logisticWidth[50];
    sprintf(logisticWidth,"Logistic width = %3.3f ns", deadtimeLogisticWidth);
    latex.SetTextColor(kBlue);
    latex.DrawLatex(0.50,0.40,mean);
    //latex.DrawLatex(0.50,0.38,logisticWidth);

    TArrow *arrow = new TArrow(236, 400, 230, 300, 0.025, "|>");
    arrow->SetAngle(30);
    arrow->SetLineWidth(1);
    arrow->SetLineColor(kBlue);
    arrow->SetFillColor(kBlue);
    arrow->Draw();

    latex.SetTextColor(kGray+2);
    latex.DrawLatex(0.20,0.87,"Linear trend");
    latex.DrawLatex(0.19,0.83,"(no deadtime)");

    //latex.SetTextSize(0.06);
    //latex.SetTextColor(kGray+2);
    //latex.DrawLatex(0.20,0.37,"Dead Zone");
    //latex.DrawLatex(0.70,0.57,"Live Zone");

    TLine *line = new TLine(deadtimeMean, 0, deadtimeMean, 1325);
    line->SetLineStyle(1);
    line->SetLineWidth(3);
    line->SetLineColor(kBlue);
    line->Draw();

    double lineStart = timeDiffFit->GetParameter(0) - timeDiffFit->GetParameter(1)*200;
    double lineStop = timeDiffFit->GetParameter(0) - timeDiffFit->GetParameter(1)*280;

    TLine *line2 = new TLine(200, lineStart, 280, lineStop);
    line2->SetLineStyle(7);
    line2->SetLineWidth(4);
    line2->SetLineColor(kGray+2);
    line2->Draw();

}
