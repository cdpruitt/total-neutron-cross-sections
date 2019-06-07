{

    string fileName = "/data1/analysis/44/0000/detTimeCheck.root";

    TFile* file = new TFile(fileName.c_str(),"READ");

    string LRCorrelationHistoName = "LR time correlation";

    file->ls();

    TH2D* LRCorrelationHisto = (TH2D*)file->Get(LRCorrelationHistoName.c_str());

    if(!LRCorrelationHisto)
    {
        cerr << "Error: couldn't find "
            << LRCorrelationHistoName  << " in "
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
    LRCorrelationHisto->SetMarkerColor(kBlack);

    // X-axis parameters
    LRCorrelationHisto->GetXaxis()->SetTitle("Left PMT fine time (ns)");
    LRCorrelationHisto->GetXaxis()->SetTitleSize(0.05);
    LRCorrelationHisto->GetXaxis()->SetTitleFont(2);
    LRCorrelationHisto->GetXaxis()->SetTitleOffset(1.5);
    LRCorrelationHisto->GetXaxis()->CenterTitle();

    LRCorrelationHisto->GetXaxis()->SetLabelOffset(0.01);
    LRCorrelationHisto->GetXaxis()->SetLabelSize(0.05);
    LRCorrelationHisto->GetXaxis()->SetLabelFont(2);

    //LRCorrelationHisto->GetXaxis()->SetNdivisions(10);
    //LRCorrelationHisto->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    LRCorrelationHisto->GetYaxis()->SetTitle("Right PMT fine time (ns)");
    LRCorrelationHisto->GetYaxis()->SetTitleSize(0.05);
    LRCorrelationHisto->GetYaxis()->SetTitleFont(2);
    LRCorrelationHisto->GetYaxis()->SetTitleOffset(1.5);
    LRCorrelationHisto->GetYaxis()->CenterTitle();

    LRCorrelationHisto->GetYaxis()->SetLabelOffset(0.01);
    LRCorrelationHisto->GetYaxis()->SetLabelSize(0.05);
    LRCorrelationHisto->GetYaxis()->SetLabelFont(2);

    //LRCorrelationHisto->GetYaxis()->SetNdivisions(10);
    //LRCorrelationHisto->GetYaxis()->SetTickLength(0.02);

    LRCorrelationHisto->GetXaxis()->SetRangeUser(30,60);
    LRCorrelationHisto->GetYaxis()->SetRangeUser(30,60);

    LRCorrelationHisto->Draw("colz");

    gPad->SetLogz();

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
    //legend->AddEntry(LRCorrelationHisto,"blank","l");
    //legend->AddEntry(target1TOFHisto,"{}^{nat}Pb","l");
    //legend->AddEntry(target2TOFHisto,"{}^{nat}C","l");
    //legend->AddEntry(target3TOFHisto,"{}^{112}Sn","l");
    //legend->AddEntry(target4TOFHisto,"{}^{nat}Sn","l");
    //legend->AddEntry(target5TOFHisto,"{}^{124}Sn","l");

    //legend->Draw();

    //TOFFile->Close();
}
