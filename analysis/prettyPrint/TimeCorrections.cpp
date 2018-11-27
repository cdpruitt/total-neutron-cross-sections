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

    string uncorrectedFileName  = "/data1/analysis/diagnostics/timingCorrections/noSoftwareCFD.root";
    string CFDCorrectedFileName = "/data1/analysis/diagnostics/timingCorrections/withSoftwareCFDNoGammaC.root";
    string fullyCorrectedFileName = "/data1/analysis/diagnostics/timingCorrections/withGammaC.root";

    TFile* uncorrectedFile = new TFile(uncorrectedFileName.c_str(),"READ");
    TFile* CFDCorrectedFile = new TFile(CFDCorrectedFileName.c_str(),"READ");
    TFile* fullyCorrectedFile = new TFile(fullyCorrectedFileName.c_str(),"READ");

    string directoryName = "summedDet";
    string tofName = "blankTOF";

    const int SHIFT = 48;

    uncorrectedFile->cd(directoryName.c_str());
    TH1D* uncorrectedTOFUnshifted = (TH1D*)gDirectory->Get(tofName.c_str());
    TH1D* uncorrectedTOF = (TH1D*)uncorrectedTOFUnshifted->Clone();
    for(int i=SHIFT+1; i<uncorrectedTOFUnshifted->GetNbinsX(); i++)
    {
        uncorrectedTOF->SetBinContent(i-SHIFT, uncorrectedTOFUnshifted->GetBinContent(i));
    }

    CFDCorrectedFile->cd(directoryName.c_str());
    TH1D* CFDCorrectedTOFUnshifted = (TH1D*)gDirectory->Get(tofName.c_str());

    TH1D* CFDCorrectedTOF = (TH1D*)CFDCorrectedTOFUnshifted->Clone();
    for(int i=SHIFT+1; i<CFDCorrectedTOF->GetNbinsX(); i++)
    {
        CFDCorrectedTOF->SetBinContent(i-SHIFT, CFDCorrectedTOFUnshifted->GetBinContent(i));
    }

    fullyCorrectedFile->cd(directoryName.c_str());
    TH1D* fullyCorrectedTOF = (TH1D*)gDirectory->Get(tofName.c_str());

    if(
            !uncorrectedTOF ||
            !CFDCorrectedTOF ||
            !fullyCorrectedTOF
      )

    {
        cerr << "Error: couldn't open "
            << tofName << " in one of the ROOT files."
            << ". Exiting..." << endl;
        exit(1);
    }

    // Set histo point and line characteristics
    uncorrectedTOF->SetLineColor(kBlack);
    CFDCorrectedTOF->SetLineColor(kBlue);
    fullyCorrectedTOF->SetLineColor(kViolet);

    uncorrectedTOF->SetLineWidth(4);
    CFDCorrectedTOF->SetLineWidth(4);
    fullyCorrectedTOF->SetLineWidth(4);

    // X-axis parameters
    uncorrectedTOF->GetXaxis()->SetTitle("TOF (ns)");
    uncorrectedTOF->GetXaxis()->SetTitleSize(0.05);
    uncorrectedTOF->GetXaxis()->SetTitleFont(2);
    uncorrectedTOF->GetXaxis()->SetTitleOffset(1.5);
    uncorrectedTOF->GetXaxis()->CenterTitle();

    uncorrectedTOF->GetXaxis()->SetLabelOffset(0.01);
    uncorrectedTOF->GetXaxis()->SetLabelSize(0.05);
    uncorrectedTOF->GetXaxis()->SetLabelFont(2);

    uncorrectedTOF->GetXaxis()->SetNdivisions(10);

    // Y-axis parameters
    uncorrectedTOF->GetYaxis()->SetTitle("Counts");
    uncorrectedTOF->GetYaxis()->SetTitleSize(0.05);
    uncorrectedTOF->GetYaxis()->SetTitleFont(2);
    uncorrectedTOF->GetYaxis()->SetTitleOffset(1.9);
    uncorrectedTOF->GetYaxis()->CenterTitle();

    uncorrectedTOF->GetYaxis()->SetLabelOffset(0.01);
    uncorrectedTOF->GetYaxis()->SetLabelSize(0.05);
    uncorrectedTOF->GetYaxis()->SetLabelFont(2);

    uncorrectedTOF->GetYaxis()->SetNdivisions(5);
    //uncorrectedTOF->GetYaxis()->SetTickLength(0.02);

    uncorrectedTOF->GetXaxis()->SetRangeUser(87,94);
    uncorrectedTOF->GetYaxis()->SetRangeUser(0,80000);

    uncorrectedTOF->Draw("");
    CFDCorrectedTOF->Draw("same");
    fullyCorrectedTOF->Draw("same");

    //gPad->SetLogy();

    TF1* f0 = new TF1("f0","gaus",87.5,93.5);
    f0->SetParameter(1, 90.5);
    f0->SetParameter(0, 80000);
    //f0->SetLineWidth(2);
    //f0->SetLineStyle(9);

    TF1* f1 = new TF1("f1","gaus",87.5,93.5);
    f1->SetParameter(1, 90.5);
    //f1->SetLineWidth(2);

    TF1* f2 = new TF1("f2","gaus",87.5,93.5);
    f2->SetParameter(1, 90.5);
    //f2->SetLineWidth(2);

    uncorrectedTOF->Fit("f0","0","",87.5,93.5);
    CFDCorrectedTOF->Fit("f1","0","",87.5,93.5);
    fullyCorrectedTOF->Fit("f2","0","",87.5,93.5);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top

    // uncorrected label
    latex.SetTextColor(kBlack);
    char FWHM [50];
    sprintf(FWHM, "%.3f", 2.355*f0->GetParameter(2));
    string uncorrectedText = "FWHM = " + string(FWHM) + " ns";
    latex.DrawLatex(0.76,0.43,"uncorrected");
    latex.DrawLatex(0.73,0.39,uncorrectedText.c_str());

    TArrow *arrow0 = new TArrow(92.8, 16000, 92.6, 10000, 0.015, "|>");
    arrow0->SetAngle(30);
    arrow0->SetLineWidth(3);
    arrow0->SetLineColor(kBlack);
    arrow0->SetFillColor(kBlack);
    arrow0->Draw();

    // CFD-corrected label
    latex.SetTextColor(kBlue);
    char FWHM2 [50];
    sprintf(FWHM2, "%.3f", 2.355*f1->GetParameter(2));
    string CFDCorrectedText = "FWHM = " + string(FWHM2) + " ns";
    latex.DrawLatex(0.69,0.58,"after CFD correction");
    latex.DrawLatex(0.71,0.54,CFDCorrectedText.c_str());

    TArrow *arrow1 = new TArrow(91.75, 32500, 91.3, 24000, 0.015, "|>");
    arrow1->SetAngle(30);
    arrow1->SetLineWidth(3);
    arrow1->SetLineColor(kBlue);
    arrow1->SetFillColor(kBlue);
    arrow1->Draw();

    // Gamma-corrected label
    latex.SetTextColor(kViolet);
    char FWHM3 [50];
    sprintf(FWHM3, "%.3f", 2.355*f2->GetParameter(2));
    string gammaCorrectedText = "FWHM = " + string(FWHM3) + " ns";
    latex.DrawLatex(0.69,0.72,"after #gamma-correction");
    latex.DrawLatex(0.69,0.68,gammaCorrectedText.c_str());

    TArrow *arrow2 = new TArrow(91.5, 48000, 90.95, 43500, 0.015, "|>");
    arrow2->SetAngle(30);
    arrow2->SetLineWidth(3);
    arrow2->SetLineColor(kViolet);
    arrow2->SetFillColor(kViolet);
    arrow2->Draw();

    /*TLine *meanLine = new TLine(-0.0929, 0, -0.0929, 11000);
    meanLine->SetLineStyle(7);
    meanLine->SetLineWidth(2);
    meanLine->SetLineColor(kGray+2);
    meanLine->Draw();
    */

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.73,0.55,0.93,0.93);
    //legend->SetTextSize(0.04);

    //legend->SetHeader("Event Energy");

    //TLegendEntry* leg2 = legend->AddEntry(fullyCorrectedTOF,"CFD + GC","l");
    //TLegendEntry* leg1 = legend->AddEntry(CFDCorrectedTOF,"CFD","l");
    //TLegendEntry* leg0 = legend->AddEntry(uncorrectedTOF,"raw","l");

    //leg2->SetTextAlign(22);
    //leg1->SetTextAlign(22);
    //leg0->SetTextAlign(22);

    //legend->Draw();

    //TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    //header->SetTextSize(0.045);
    //header->SetTextAlign(22);

    //legend->SetEntrySeparation(2);
}
