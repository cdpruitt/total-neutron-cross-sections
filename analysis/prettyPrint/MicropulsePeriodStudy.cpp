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

    string file0Name = "/data1/analysis/diagnostics/MicroPeriodStudy/15/0000/histos_1788.8118ns.root";
    string file1Name = "/data1/analysis/diagnostics/MicroPeriodStudy/15/0000/histos_1788.8128ns.root";
    string file2Name = "/data1/analysis/diagnostics/MicroPeriodStudy/15/0000/histos_1788.8138ns.root";
    string file3Name = "/data1/analysis/diagnostics/MicroPeriodStudy/15/0000/histos_1788.8148ns.root";
    string file4Name = "/data1/analysis/diagnostics/MicroPeriodStudy/15/0000/histos_1788.8158ns.root";
    string file5Name = "/data1/analysis/diagnostics/MicroPeriodStudy/15/0000/histos_1788.8168ns.root";
    string file6Name = "/data1/analysis/diagnostics/MicroPeriodStudy/15/0000/histos_1788.8178ns.root";

    TFile* file0 = new TFile(file0Name.c_str(),"READ");
    TFile* file1 = new TFile(file1Name.c_str(),"READ");
    TFile* file2 = new TFile(file2Name.c_str(),"READ");
    TFile* file3 = new TFile(file3Name.c_str(),"READ");
    TFile* file4 = new TFile(file4Name.c_str(),"READ");
    TFile* file5 = new TFile(file5Name.c_str(),"READ");
    TFile* file6 = new TFile(file6Name.c_str(),"READ");

    TDirectory* file0Dir = file0->GetDirectory("summedDet");
    TDirectory* file1Dir = file1->GetDirectory("summedDet");
    TDirectory* file2Dir = file2->GetDirectory("summedDet");
    TDirectory* file3Dir = file3->GetDirectory("summedDet");
    TDirectory* file4Dir = file4->GetDirectory("summedDet");
    TDirectory* file5Dir = file5->GetDirectory("summedDet");
    TDirectory* file6Dir = file6->GetDirectory("summedDet");

    string MicropulsePeriodName = "blankTOF";

    TH1D* MicropulsePeriodHisto0 = (TH1D*)file0Dir->Get(MicropulsePeriodName.c_str());
    TH1D* MicropulsePeriodHisto1 = (TH1D*)file1Dir->Get(MicropulsePeriodName.c_str());
    TH1D* MicropulsePeriodHisto2 = (TH1D*)file2Dir->Get(MicropulsePeriodName.c_str());
    TH1D* MicropulsePeriodHisto3 = (TH1D*)file3Dir->Get(MicropulsePeriodName.c_str());
    TH1D* MicropulsePeriodHisto4 = (TH1D*)file4Dir->Get(MicropulsePeriodName.c_str());
    TH1D* MicropulsePeriodHisto5 = (TH1D*)file5Dir->Get(MicropulsePeriodName.c_str());
    TH1D* MicropulsePeriodHisto6 = (TH1D*)file6Dir->Get(MicropulsePeriodName.c_str());

    if(
            !MicropulsePeriodHisto0 ||
            !MicropulsePeriodHisto1 ||
            !MicropulsePeriodHisto2 ||
            !MicropulsePeriodHisto3 ||
            !MicropulsePeriodHisto4 ||
            !MicropulsePeriodHisto5 ||
            !MicropulsePeriodHisto6
      )

    {
        cerr << "Error: couldn't find "
            << MicropulsePeriodName << " in one of time check threshold ROOT files."
            << ". Exiting..." << endl;
        exit(1);
    }

    // Set histo point and line characteristics
    MicropulsePeriodHisto0->SetLineColor(kRed+4);
    MicropulsePeriodHisto1->SetLineColor(kRed+3);
    MicropulsePeriodHisto2->SetLineColor(kRed+2);
    MicropulsePeriodHisto3->SetLineColor(kRed+1);
    MicropulsePeriodHisto4->SetLineColor(kRed);
    MicropulsePeriodHisto5->SetLineColor(kRed-1);
    MicropulsePeriodHisto6->SetLineColor(kRed-2);

    MicropulsePeriodHisto0->SetLineWidth(4);
    MicropulsePeriodHisto1->SetLineWidth(4);
    MicropulsePeriodHisto2->SetLineWidth(4);
    MicropulsePeriodHisto3->SetLineWidth(4);
    MicropulsePeriodHisto4->SetLineWidth(4);
    MicropulsePeriodHisto5->SetLineWidth(4);
    MicropulsePeriodHisto6->SetLineWidth(4);

    // X-axis parameters
    MicropulsePeriodHisto3->GetXaxis()->SetTitle("Time-of-Flight (ns)");
    MicropulsePeriodHisto3->GetXaxis()->SetTitleSize(0.05);
    MicropulsePeriodHisto3->GetXaxis()->SetTitleFont(2);
    MicropulsePeriodHisto3->GetXaxis()->SetTitleOffset(1.5);
    MicropulsePeriodHisto3->GetXaxis()->CenterTitle();

    MicropulsePeriodHisto3->GetXaxis()->SetLabelOffset(0.01);
    MicropulsePeriodHisto3->GetXaxis()->SetLabelSize(0.05);
    MicropulsePeriodHisto3->GetXaxis()->SetLabelFont(2);

    MicropulsePeriodHisto3->GetXaxis()->SetNdivisions(6);
    //MicropulsePeriodHisto3->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    //MicropulsePeriodHisto3->GetYaxis()->SetTitle("ns)");
    //MicropulsePeriodHisto3->GetYaxis()->SetTitleSize(0.05);
    //MicropulsePeriodHisto3->GetYaxis()->SetTitleFont(2);
    //MicropulsePeriodHisto3->GetYaxis()->SetTitleOffset(1.5);
    //MicropulsePeriodHisto3->GetYaxis()->CenterTitle();

    MicropulsePeriodHisto3->GetYaxis()->SetLabelOffset(0.01);
    MicropulsePeriodHisto3->GetYaxis()->SetLabelSize(0.05);
    MicropulsePeriodHisto3->GetYaxis()->SetLabelFont(2);

    //MicropulsePeriodHisto3->GetYaxis()->SetNdivisions(10);
    //MicropulsePeriodHisto3->GetYaxis()->SetTickLength(0.02);

    MicropulsePeriodHisto3->GetXaxis()->SetRangeUser(87,94);
    //MicropulsePeriodHisto3->GetYaxis()->SetRangeUser(0,10500);

    MicropulsePeriodHisto3->Draw();
    MicropulsePeriodHisto0->Draw("same");
    MicropulsePeriodHisto1->Draw("same");
    MicropulsePeriodHisto2->Draw("same");
    MicropulsePeriodHisto4->Draw("same");
    MicropulsePeriodHisto5->Draw("same");
    MicropulsePeriodHisto6->Draw("same");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.7,0.55,0.93,0.93);
    legend->SetTextSize(0.04);

    legend->SetHeader("Event Energy");

    MicropulsePeriodHisto6->SetLineWidth(5);
    MicropulsePeriodHisto5->SetLineWidth(5);
    MicropulsePeriodHisto4->SetLineWidth(5);
    MicropulsePeriodHisto3->SetLineWidth(5);
    MicropulsePeriodHisto2->SetLineWidth(5);
    MicropulsePeriodHisto1->SetLineWidth(5);
    MicropulsePeriodHisto0->SetLineWidth(5);

    TLegendEntry* leg6 = legend->AddEntry(MicropulsePeriodHisto6,"-3 ps","l");
    TLegendEntry* leg5 = legend->AddEntry(MicropulsePeriodHisto5,"-2 ps","l");
    TLegendEntry* leg4 = legend->AddEntry(MicropulsePeriodHisto4,"-1 ps","l");
    TLegendEntry* leg3 = legend->AddEntry(MicropulsePeriodHisto3,"0 ps","l");
    TLegendEntry* leg2 = legend->AddEntry(MicropulsePeriodHisto2,"1 ps","l");
    TLegendEntry* leg1 = legend->AddEntry(MicropulsePeriodHisto1,"2 ps","l");
    TLegendEntry* leg0 = legend->AddEntry(MicropulsePeriodHisto0,"3 ps","l");

    leg6->SetTextAlign(22);
    leg5->SetTextAlign(22);
    leg4->SetTextAlign(22);
    leg3->SetTextAlign(22);
    leg2->SetTextAlign(22);
    leg1->SetTextAlign(22);
    leg0->SetTextAlign(22);

    legend->Draw();

    TLegendEntry *header = (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextSize(0.045);
    header->SetTextAlign(22);
}
