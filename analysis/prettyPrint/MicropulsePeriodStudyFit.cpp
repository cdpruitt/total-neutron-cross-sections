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

    double x[5] = {-2, -1, 0, 1, 2};
    double y[5] = {0.9725, 0.8764, 0.8223, 0.8702, 0.9661};
    TGraph* MicropulsePeriodStudyFitGraph = new TGraph(5,x,y);

    double x2[5] = {-2, -1, 0, 1, 2};
    double y2[5] = {0.9491, 0.8562, 0.8254, 0.8626, 0.9627};
    TGraph* MicropulsePeriodStudyFitGraph2 = new TGraph(5,x2,y2);

    // Set histo point and line characteristics
    MicropulsePeriodStudyFitGraph->SetMarkerColor(kRed);
    MicropulsePeriodStudyFitGraph2->SetMarkerColor(kBlue);

    style->SetFuncColor(kRed);
    TF1* fitToFWHM = new TF1("fit","[0]*(x-[1])*(x-[1])+[2]",-2.5,2.5);
    MicropulsePeriodStudyFitGraph->Fit("fit","","",-2.2,2.2);
    fitToFWHM->SetLineColor(kRed);
    fitToFWHM->SetLineWidth(3);

    style->SetFuncColor(kBlue);
    TF1* fitToFWHM2 = new TF1("fit2","[0]*(x-[1])*(x-[1])+[2]",-2.5,2.5);
    MicropulsePeriodStudyFitGraph2->Fit("fit2","","",-2.2,2.2);
    fitToFWHM2->SetLineColor(kBlue);
    fitToFWHM2->SetLineWidth(3);

    TMultiGraph* allGraphs = new TMultiGraph();
    allGraphs->Add(MicropulsePeriodStudyFitGraph,"AP");
    allGraphs->Add(MicropulsePeriodStudyFitGraph2,"AP");

    allGraphs->Draw("AP");

    // X-axis parameters
    allGraphs->GetXaxis()->SetTitle("Micropulse period #Delta (ps)");
    allGraphs->GetXaxis()->SetTitleSize(0.04);
    allGraphs->GetXaxis()->SetTitleFont(2);
    allGraphs->GetXaxis()->SetTitleOffset(1.5);
    allGraphs->GetXaxis()->CenterTitle();

    allGraphs->GetXaxis()->SetLabelOffset(0.01);
    allGraphs->GetXaxis()->SetLabelSize(0.04);
    allGraphs->GetXaxis()->SetLabelFont(2);

    allGraphs->GetXaxis()->SetNdivisions(10);

    // Y-axis parameters
    allGraphs->GetYaxis()->SetTitle("Gamma peak FWHM (ns)");
    allGraphs->GetYaxis()->SetTitleSize(0.04);
    allGraphs->GetYaxis()->SetTitleFont(2);
    allGraphs->GetYaxis()->SetTitleOffset(1.7);
    allGraphs->GetYaxis()->CenterTitle();

    allGraphs->GetYaxis()->SetLabelOffset(0.01);
    allGraphs->GetYaxis()->SetLabelSize(0.04);
    allGraphs->GetYaxis()->SetLabelFont(2);

    allGraphs->GetXaxis()->SetRangeUser(-3,3);
    allGraphs->GetYaxis()->SetRangeUser(0.8,1.0);

    gPad->SetGrid(1,0);

    // Define legend format and contents
    TLegend *legend = new TLegend(0.57,0.71,0.75,0.83);
    legend->AddEntry(MicropulsePeriodStudyFitGraph,"13 Sept, 2017","p");
    legend->AddEntry(MicropulsePeriodStudyFitGraph2,"21 Sept, 2017","p");

    legend->Draw();
}
