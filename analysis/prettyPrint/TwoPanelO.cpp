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
    string expFileName = "/data2/analysis/total.root";
    string relFileName = "/data2/analysis/relative.root";
    string litFileName = "/data2/analysis/literatureData.root";
    string shiftedFileName = "/data2/analysis/shifted.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");
    TFile* shiftedFile = new TFile(shiftedFileName.c_str(),"READ");

    if(!expFile || !relFile || !litFile || !shiftedFile)
    {
        cout << "Error: failed to open a required file (check filenames). Exiting..." << endl;
        exit(1);
    }

    string expO16GraphName = "ONat_fromH2O";
    string expO18GraphName = "O18+1Barn";

    string relO16GraphName = "ONat_fromH2O, expLit, percent";
    string relO18GraphName = "O18, expLit, percent";

    string litO16GraphName = "NatO(n,tot)";
    string litO18GraphName = "18O(n,tot)+1Barn";
    
    TGraphAsymmErrors* expO16Graph = (TGraphAsymmErrors*)expFile->Get(expO16GraphName.c_str());
    TGraphAsymmErrors* expO18Graph = (TGraphAsymmErrors*)shiftedFile->Get(expO18GraphName.c_str());

    TGraphAsymmErrors* relO16Graph = (TGraphAsymmErrors*)relFile->Get(relO16GraphName.c_str());
    TGraphAsymmErrors* relO18Graph = (TGraphAsymmErrors*)relFile->Get(relO18GraphName.c_str());

    TGraphAsymmErrors* litO16Graph = (TGraphAsymmErrors*)litFile->Get(litO16GraphName.c_str());
    TGraphAsymmErrors* litO18Graph = (TGraphAsymmErrors*)litFile->Get(litO18GraphName.c_str());

    if(!expO16Graph || !expO18Graph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    if(!relO16Graph || !relO18Graph)
    {
        cout << "Error: failed to open an relative diff to lit cross section graph." << endl;
        exit(1);
    }

    if(!litO16Graph || !litO18Graph)
    {
        cout << "Error: failed to open a lit cross section graph." << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expO16Graph->SetLineColor(kRed);
    expO16Graph->SetLineWidth(4);
    expO16Graph->SetLineStyle(1);
    expO16Graph->SetMarkerColor(kRed);

    expO18Graph->SetLineColor(kRed+2);
    expO18Graph->SetLineWidth(4);
    expO18Graph->SetLineStyle(1);
    expO18Graph->SetMarkerColor(kRed+2);

    relO16Graph->SetLineColor(kRed);
    relO16Graph->SetLineWidth(4);
    relO16Graph->SetLineStyle(0);
    relO16Graph->SetMarkerColor(kRed);
    relO16Graph->SetMarkerSize(2);
    relO16Graph->SetMarkerStyle(22);

    relO18Graph->SetLineColor(kRed+2);
    relO18Graph->SetMarkerColor(kRed+2);
    relO18Graph->SetMarkerSize(2);
    relO18Graph->SetMarkerStyle(22);

    litO16Graph->SetLineColor(kBlue-7);
    litO16Graph->SetLineWidth(4);
    litO16Graph->SetLineStyle(1);
    litO16Graph->SetMarkerColor(kBlue-7);
    litO16Graph->SetMarkerSize(2);
    litO16Graph->SetMarkerStyle(22);

    litO18Graph->SetLineColor(kBlue);
    litO18Graph->SetLineWidth(4);
    litO18Graph->SetLineStyle(1);
    litO18Graph->SetMarkerColor(kBlue);
    litO18Graph->SetMarkerSize(2);
    litO18Graph->SetMarkerStyle(22);

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

        mg->Add(litO16Graph, "l");
        mg->Add(litO18Graph, "l");
        mg->Add(expO16Graph, "l");
        mg->Add(expO18Graph, "l");

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
        mg->GetYaxis()->SetNdivisions(8);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetYaxis()->SetRangeUser(0.1,4.59);
        mg->GetXaxis()->SetLimits(2.8,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top

        latex.DrawLatex(0.35, 0.68, "+1 barn");

        latex.SetTextSize(0.15);
        latex.DrawLatex(0.02, 1, "a)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.77,0.51,0.96,0.93);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(1);
        legend->AddEntry(litO16Graph,"^{16}O (An)","l");
        legend->AddEntry(expO16Graph,"^{16}O (DSP)","l");
        legend->AddEntry(litO18Graph,"^{18}O (An)","l");
        legend->AddEntry(expO18Graph,"^{18}O (DSP)","l");

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

        mg->Add(relO18Graph, "lp");
        mg->Add(relO16Graph, "lp");

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
        mg->GetYaxis()->SetNdivisions(8);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetXaxis()->SetLimits(2.8,500);
        mg->GetYaxis()->SetRangeUser(-14.9,14.9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.13);
        latex.DrawLatex(0.02, 1, "b)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.83,0.74,0.96,0.96);
        legend->SetNColumns(1);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relO16Graph,"{}^{16}O","lp");
        legend->AddEntry(relO18Graph,"{}^{18}O","lp");

        legend->Draw();

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
