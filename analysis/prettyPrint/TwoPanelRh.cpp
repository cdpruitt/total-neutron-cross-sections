{
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas", 850, 850);
    canvas->Divide(1,2, 0, 0);

    style->SetOptStat(0);
    style->SetOptTitle(0);    
    //style->SetPalette(1,0);
    style->SetCanvasColor(10);      
    style->SetCanvasBorderMode(0);    
    style->SetFrameLineWidth(3);
    style->SetFrameFillColor(10);
    style->SetPadColor(10);
    style->SetHistLineWidth(4);
    style->SetHistLineColor(kRed);
    style->SetMarkerSize(0.6);
    style->SetMarkerStyle(8);
    //style->SetFuncWidth(3);
    //style->SetFuncColor(kRed);
    style->SetLabelColor(kBlack,"xyz");
    style->SetTitleSize(0.08,"xyz");
    style->SetTitleFillColor(10);
    style->SetTitleTextColor(kBlack);
    style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    // read graphs
    string expFileName = "/data1/analysis/total.root";
    string relFileName = "/data1/analysis/relative.root";
    string litFileName = "/data1/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");

    if(!expFile || !relFile || !litFile)
    {
        cout << "Error: failed to open a required file (check filenames). Exiting..." << endl;
        exit(1);
    }

    string expRhNatGraphName = "RhNat";

    string relRhNatGraphName = "RhNat, expLit, percent";

    string litRhNatGraphName = "Natural Rh (n,tot)";
    
    TGraphAsymmErrors* expRhNatGraph = (TGraphAsymmErrors*)expFile->Get(expRhNatGraphName.c_str());

    TGraphAsymmErrors* relRhNatGraph = (TGraphAsymmErrors*)relFile->Get(relRhNatGraphName.c_str());

    TGraphAsymmErrors* litRhNatGraph = (TGraphAsymmErrors*)litFile->Get(litRhNatGraphName.c_str());

    if(!expRhNatGraph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    if(!relRhNatGraph)
    {
        cout << "Error: failed to open an relative diff to lit cross section graph." << endl;
        exit(1);
    }

    if(!litRhNatGraph)
    {
        cout << "Error: failed to open an lit cross section graph." << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expRhNatGraph->SetLineColor(kRed);
    expRhNatGraph->SetLineWidth(4);
    expRhNatGraph->SetLineStyle(1);
    expRhNatGraph->SetMarkerColor(kRed);

    relRhNatGraph->SetLineColor(kRed);
    relRhNatGraph->SetLineWidth(4);
    relRhNatGraph->SetLineStyle(0);
    relRhNatGraph->SetMarkerColor(kRed);

    litRhNatGraph->SetLineColor(kBlue);
    litRhNatGraph->SetLineWidth(4);
    litRhNatGraph->SetLineStyle(1);
    litRhNatGraph->SetMarkerColor(kBlue);

    // first panel
    {
        canvas->cd(1);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.25);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.005);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        mg->Add(expRhNatGraph, "l");
        mg->Add(litRhNatGraph, "l");

        mg->Draw("al");

        // X-axis parameters
        mg->GetXaxis()->SetTitle("Energy (MeV)");
        mg->GetXaxis()->SetTitleSize(0.08);
        mg->GetXaxis()->SetTitleFont(2);
        mg->GetXaxis()->SetTitleOffset(1.4);
        mg->GetXaxis()->CenterTitle();

        mg->GetXaxis()->SetLabelOffset(0.01);
        mg->GetXaxis()->SetLabelSize(0.08);
        mg->GetXaxis()->SetLabelFont(2);

        mg->GetXaxis()->SetNdivisions(10);
        mg->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        mg->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        mg->GetYaxis()->SetTitleSize(0.10);
        mg->GetYaxis()->SetTitleFont(2);
        mg->GetYaxis()->SetTitleOffset(0.6);
        mg->GetYaxis()->CenterTitle();

        mg->GetYaxis()->SetLabelOffset(0.01);
        mg->GetYaxis()->SetLabelSize(0.08);

        mg->GetYaxis()->SetLabelFont(2);
        mg->GetYaxis()->SetNdivisions(6);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetYaxis()->SetRangeUser(0.81,5.19);
        mg->GetXaxis()->SetLimits(2.9,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.15);
        latex.DrawLatex(0.02, 1, "a)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.70,0.69,0.93,0.91);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(1);
        legend->AddEntry(litRhNatGraph,"^{103}Rh (An)","l");
        legend->AddEntry(expRhNatGraph,"^{103}Rh (DSP)","l");

        legend->Draw();
    }

    // second panel
    {
        canvas->cd(2);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.25);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.);
        gPad->SetBottomMargin(0.25);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx();

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        mg->Add(relRhNatGraph, "l");

        mg->Draw("al");

        // X-axis parameters
        mg->GetXaxis()->SetTitle("Energy (MeV)");
        mg->GetXaxis()->SetTitleSize(0.08);
        mg->GetXaxis()->SetTitleFont(2);
        mg->GetXaxis()->SetTitleOffset(1.4);
        mg->GetXaxis()->CenterTitle();

        mg->GetXaxis()->SetLabelOffset(0.01);
        mg->GetXaxis()->SetLabelSize(0.08);
        mg->GetXaxis()->SetLabelFont(2);

        mg->GetXaxis()->SetNdivisions(10);
        mg->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        mg->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        mg->GetYaxis()->SetTitleSize(0.08);
        mg->GetYaxis()->SetTitleFont(2);
        mg->GetYaxis()->SetTitleOffset(0.9);
        mg->GetYaxis()->CenterTitle();

        mg->GetYaxis()->SetLabelOffset(0.01);
        mg->GetYaxis()->SetLabelSize(0.08);

        mg->GetYaxis()->SetLabelFont(2);
        mg->GetYaxis()->SetNdivisions(10);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetXaxis()->SetLimits(2.9,500);
        mg->GetYaxis()->SetRangeUser(-3.9,3.9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.13);
        latex.DrawLatex(0.02, 1, "b)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.84,0.74,0.96,0.96);
        legend->SetNColumns(1);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relRhNatGraph,"{}^{Nat}Rh","l");

        //legend->Draw();

        TLine* zeroLine = new TLine(3, 0, 500, 0);
        zeroLine->SetLineColor(kBlack);
        zeroLine->SetLineWidth(3);
        zeroLine->SetLineStyle(9);
        zeroLine->Draw();
    }

    // close it all up
    expFile->Close();
    relFile->Close();
    litFile->Close();
}
