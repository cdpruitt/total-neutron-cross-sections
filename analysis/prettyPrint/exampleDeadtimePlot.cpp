{
    TStyle* style = (TStyle*)gROOT->FindObject("blankHistoStyle");

    if(!style)      
    {
        style = new TStyle("blankHistoStyle","blankHistoStyle");
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

    gROOT->SetStyle("blankHistoStyle");
    gROOT->ForceStyle();

    string fileName = "/data1/analysis/15/0000/deadtime.root";

    TFile* file = new TFile(fileName.c_str(),"READ");

    string directoryName = "summedDet";

    string blankHistoName = "blankDeadtime";
    string CNatHistoName = "CNatDeadtime";

    TDirectory* directory = (TDirectory*)file->GetDirectory(directoryName.c_str());

    TH1D* blankHisto = (TH1D*)(directory->Get(blankHistoName.c_str()));

    if(!blankHisto)
    {
        cerr << "Error: couldn't find "
            << blankHistoName  << " in "
            << fileName << ". Exiting..." << endl;
        exit(1);
    }

    TH1D* CNatHisto = (TH1D*)(directory->Get(CNatHistoName.c_str()));

    if(!CNatHisto)
    {
        cerr << "Error: couldn't find "
            << CNatHistoName  << " in "
            << fileName << ". Exiting..." << endl;
        exit(1);
    }

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.20);
    //gPad->SetTicky(2);

    // Set blankHisto point and line characteristics
    blankHisto->SetMarkerColor(kBlack);
    blankHisto->SetLineWidth(4);
    blankHisto->SetLineColor(kRed);

    CNatHisto->SetLineWidth(4);

    // X-axis parameters
    blankHisto->GetXaxis()->SetTitle("TOF (ns)");
    blankHisto->GetXaxis()->SetTitleSize(0.05);
    blankHisto->GetXaxis()->SetTitleFont(2);
    blankHisto->GetXaxis()->SetTitleOffset(1.5);
    blankHisto->GetXaxis()->CenterTitle();

    blankHisto->GetXaxis()->SetLabelOffset(0.01);
    blankHisto->GetXaxis()->SetLabelSize(0.05);
    blankHisto->GetXaxis()->SetLabelFont(2);

    //blankHisto->GetXaxis()->SetNdivisions(10);
    //blankHisto->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    blankHisto->GetYaxis()->SetTitle("Fraction Dead");
    blankHisto->GetYaxis()->SetTitleSize(0.05);
    blankHisto->GetYaxis()->SetTitleFont(2);
    blankHisto->GetYaxis()->SetTitleOffset(1.5);
    blankHisto->GetYaxis()->CenterTitle();

    blankHisto->GetYaxis()->SetLabelOffset(0.01);
    blankHisto->GetYaxis()->SetLabelSize(0.05);
    blankHisto->GetYaxis()->SetLabelFont(2);

    //blankHisto->GetYaxis()->SetNdivisions(10);
    //blankHisto->GetYaxis()->SetTickLength(0.02);

    //blankHisto->GetXaxis()->SetRangeUser(-1,1);

    blankHisto->Draw();
    CNatHisto->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top
    latex.SetTextColor(kBlack);

    //gPad->SetLogy(1);

    // Define legend format and contents
    TLegend *legend = new TLegend(0.66,0.72,0.78,0.89);
    //legend->SetHeader("data","C");
    legend->AddEntry(blankHisto,"blank","l");
    legend->AddEntry(CNatHisto,"{}^{nat}C","l");

    legend->Draw();

    //TOFFile->Close();
}
