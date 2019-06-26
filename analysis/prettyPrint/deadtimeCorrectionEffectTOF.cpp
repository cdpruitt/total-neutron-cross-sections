{
    TStyle* style = (TStyle*)gROOT->FindObject("TOFHistoStyle");

    if(!style)      
    {
        style = new TStyle("TOFHistoStyle","TOFHistoStyle");
    }

    TCanvas* c = new TCanvas("c1","",1300,900);

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

    gROOT->SetStyle("TOFHistoStyle");
    gROOT->ForceStyle();

    string fileName = "/data1/analysis/total.root";

    TFile* file = new TFile(fileName.c_str(),"READ");

    string TOFHistoName = "blankTOF";
    string uncorrectedTOFHistoName = "blankTOFUncorrected";

    TH1D* TOFHisto = (TH1D*)(file->Get(TOFHistoName.c_str()));

    if(!TOFHisto)
    {
        cerr << "Error: couldn't find "
            << TOFHistoName  << " in "
            << fileName << ". Exiting..." << endl;
        exit(1);
    }

    TH1D* uncorrectedTOFHisto = (TH1D*)(file->Get(uncorrectedTOFHistoName.c_str()));

    if(!uncorrectedTOFHisto)
    {
        cerr << "Error: couldn't find "
            << uncorrectedTOFHistoName  << " in "
            << fileName << ". Exiting..." << endl;
        exit(1);
    }

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.14);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    gPad->SetFrameLineWidth(3);

    // Set TOFHisto point and line characteristics
    TOFHisto->SetMarkerColor(kBlack);
    TOFHisto->SetLineWidth(4);
    TOFHisto->SetLineColor(kRed);

    uncorrectedTOFHisto->SetLineWidth(6);

    // X-axis parameters
    uncorrectedTOFHisto->GetXaxis()->SetTitle("TOF (ns)");
    uncorrectedTOFHisto->GetXaxis()->SetTitleSize(0.05);
    uncorrectedTOFHisto->GetXaxis()->SetTitleFont(2);
    uncorrectedTOFHisto->GetXaxis()->SetTitleOffset(1.5);
    uncorrectedTOFHisto->GetXaxis()->CenterTitle();

    uncorrectedTOFHisto->GetXaxis()->SetLabelOffset(0.01);
    uncorrectedTOFHisto->GetXaxis()->SetLabelSize(0.05);
    uncorrectedTOFHisto->GetXaxis()->SetLabelFont(2);

    //uncorrectedTOFHisto->GetXaxis()->SetNdivisions(10);
    //uncorrectedTOFHisto->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    uncorrectedTOFHisto->GetYaxis()->SetTitle("Counts");
    uncorrectedTOFHisto->GetYaxis()->SetTitleSize(0.05);
    uncorrectedTOFHisto->GetYaxis()->SetTitleFont(2);
    uncorrectedTOFHisto->GetYaxis()->SetTitleOffset(1.);
    uncorrectedTOFHisto->GetYaxis()->CenterTitle();

    uncorrectedTOFHisto->GetYaxis()->SetLabelOffset(0.01);
    uncorrectedTOFHisto->GetYaxis()->SetLabelSize(0.05);
    uncorrectedTOFHisto->GetYaxis()->SetLabelFont(2);

    //uncorrectedTOFHisto->GetYaxis()->SetNdivisions(10);
    //uncorrectedTOFHisto->GetYaxis()->SetTickLength(0.02);

    uncorrectedTOFHisto->GetXaxis()->SetRangeUser(100, 499);
    uncorrectedTOFHisto->GetYaxis()->SetRangeUser(3000, 110000);

    uncorrectedTOFHisto->Draw();
    TOFHisto->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.15, 0.82, "300 MeV");
    latex.DrawLatex(0.24, 0.66, "150 MeV");
    latex.DrawLatex(0.40, 0.56, "71 MeV");

    gPad->SetLogy(1);

    // Define legend format and contents
    TLegend *legend = new TLegend(0.58,0.72,0.88,0.89);
    //legend->SetHeader("data","C");
    legend->SetTextSize(0.04);
    legend->SetTextAlign(12);
    legend->AddEntry(TOFHisto,"Deadtime-corrected","l");
    legend->AddEntry(uncorrectedTOFHisto,"Uncorrected","l");

    legend->Draw();

    TArrow *arrow = new TArrow(245, 15500, 245, 10500, 0.015, "|>");
    arrow->SetAngle(30);
    arrow->SetLineWidth(1);
    arrow->SetLineColor(kBlack);
    arrow->SetFillColor(kBlack);
    arrow->Draw();

    TArrow *arrow2 = new TArrow(139, 45000, 139, 29000, 0.015, "|>");
    arrow2->SetAngle(30);
    arrow2->SetLineWidth(1);
    arrow2->SetLineColor(kBlack);
    arrow2->SetFillColor(kBlack);
    arrow2->Draw();

    TArrow *arrow3 = new TArrow(178.5, 23000, 178.5, 16000, 0.015, "|>");
    arrow3->SetAngle(30);
    arrow3->SetLineWidth(1);
    arrow3->SetLineColor(kBlack);
    arrow3->SetFillColor(kBlack);
    arrow3->Draw();


    //TOFFile->Close();
}
