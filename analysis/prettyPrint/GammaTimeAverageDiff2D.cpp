{
    string fileName = "/data1/analysis/15/0000/histos.root";

    TFile* file = new TFile(fileName.c_str(),"READ");

    string histoName = "gamma average diff, 2D";

    file->cd("summedDet");

    TH2D* histo = (TH2D*)gDirectory->Get(histoName.c_str());

    if(!histo)
    {
        cerr << "Error: couldn't find "
            << histoName  << " in "
            << fileName << ". Exiting..." << endl;
        exit(1);
    }

    TStyle* style = (TStyle*)gROOT->FindObject("histoStyle");

    if(!style)      
    {
        style = new TStyle("histoStyle","histoStyle");
    }

    TCanvas* c = new TCanvas("c1","",900,900);

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

    gROOT->SetStyle("histoStyle");
    gROOT->ForceStyle();

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.20);
    //gPad->SetTicky(2);

    // Set histo point and line characteristics
    histo->SetMarkerColor(kBlack);

    // X-axis parameters
    histo->GetXaxis()->SetTitle("Time Average Difference (ns)");
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetTitleFont(2);
    histo->GetXaxis()->SetTitleOffset(1.5);
    histo->GetXaxis()->CenterTitle();

    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetXaxis()->SetLabelSize(0.05);
    histo->GetXaxis()->SetLabelFont(2);

    //histo->GetXaxis()->SetNdivisions(10);
    //histo->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    histo->GetYaxis()->SetTitle("Gammas in Each Partition");
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleFont(2);
    histo->GetYaxis()->SetTitleOffset(1.5);
    histo->GetYaxis()->CenterTitle();

    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->GetYaxis()->SetLabelFont(2);

    //histo->GetYaxis()->SetNdivisions(10);
    //histo->GetYaxis()->SetTickLength(0.02);

    histo->GetXaxis()->SetRangeUser(-2,2);
    histo->GetYaxis()->SetRangeUser(0,30);

    histo->Draw("colz");

    //gPad->SetLogz();

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Sn");
    //latex.DrawLatex(0.32,0.4,"C");

    //gPad->SetLogy(1);

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.21,0.63,0.35,0.93);
    //legend->SetHeader("data","C");
    //legend->AddEntry(histo,"blank","l");
    //legend->AddEntry(target1TOFHisto,"{}^{nat}Pb","l");
    //legend->AddEntry(target2TOFHisto,"{}^{nat}C","l");
    //legend->AddEntry(target3TOFHisto,"{}^{112}Sn","l");
    //legend->AddEntry(target4TOFHisto,"{}^{nat}Sn","l");
    //legend->AddEntry(target5TOFHisto,"{}^{124}Sn","l");

    //legend->Draw();

    //TOFFile->Close();
}
