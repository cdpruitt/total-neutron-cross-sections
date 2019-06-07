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
    //style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    // read graphs
    string expFileName = "/data2/analysis/total.root";
    string relFileName = "/data2/analysis/relative.root";
    string litFileName = "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");

    if(!expFile || !relFile || !litFile)
    {
        cout << "Error: failed to open a required file (check filenames). Exiting..." << endl;
        exit(1);
    }

    string expSn112GraphName = "Sn112";
    string expSn124GraphName = "Sn124";

    string relSn112GraphName = "Sn112, expLit, percent";
    string relSn124GraphName = "Sn124, expLit, percent";

    string litSn112GraphName = "Sn112 (n,tot)";
    string litSn124GraphName = "Sn124 (n,tot)";
    
    TGraphAsymmErrors* expSn112Graph = (TGraphAsymmErrors*)expFile->Get(expSn112GraphName.c_str());
    TGraphAsymmErrors* expSn124Graph = (TGraphAsymmErrors*)expFile->Get(expSn124GraphName.c_str());

    TGraphAsymmErrors* relSn112Graph = (TGraphAsymmErrors*)relFile->Get(relSn112GraphName.c_str());
    TGraphAsymmErrors* relSn124Graph = (TGraphAsymmErrors*)relFile->Get(relSn124GraphName.c_str());

    TGraphAsymmErrors* litSn112Graph = (TGraphAsymmErrors*)litFile->Get(litSn112GraphName.c_str());
    TGraphAsymmErrors* litSn124Graph = (TGraphAsymmErrors*)litFile->Get(litSn124GraphName.c_str());

    if(!expSn112Graph || !expSn124Graph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    if(!relSn112Graph || !relSn124Graph)
    {
        cout << "Error: failed to open an relative diff to lit cross section graph." << endl;
        exit(1);
    }

    if(!litSn112Graph || !litSn124Graph)
    {
        cout << "Error: failed to open an lit cross section graph." << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expSn112Graph->SetLineColor(kRed);
    expSn112Graph->SetLineWidth(4);
    expSn112Graph->SetLineStyle(1);
    expSn112Graph->SetMarkerColor(kRed);
    expSn112Graph->SetMarkerSize(2.5);

    expSn124Graph->SetLineColor(kRed+2);
    expSn124Graph->SetLineWidth(4);
    expSn124Graph->SetLineStyle(1);
    expSn124Graph->SetMarkerColor(kRed+2);
    expSn124Graph->SetMarkerSize(2.5);

    relSn112Graph->SetLineColor(kRed);
    relSn112Graph->SetLineWidth(4);
    relSn112Graph->SetLineStyle(0);
    relSn112Graph->SetMarkerColor(kRed);
    relSn112Graph->SetMarkerSize(2.5);
    relSn112Graph->SetMarkerStyle(22);

    relSn124Graph->SetLineColor(kRed+2);
    relSn124Graph->SetMarkerColor(kRed+2);
    relSn124Graph->SetMarkerSize(2.5);
    relSn124Graph->SetMarkerStyle(23);

    litSn112Graph->SetLineColor(kBlue-7);
    litSn112Graph->SetLineWidth(4);
    litSn112Graph->SetLineStyle(1);
    litSn112Graph->SetMarkerColor(kBlue-7);
    litSn112Graph->SetMarkerSize(2.5);
    litSn112Graph->SetMarkerStyle(22);

    litSn124Graph->SetLineColor(kBlue);
    litSn124Graph->SetLineWidth(4);
    litSn124Graph->SetLineStyle(1);
    litSn124Graph->SetMarkerColor(kBlue);
    litSn124Graph->SetMarkerSize(2.5);
    litSn124Graph->SetMarkerStyle(23);

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

        mg->Add(expSn112Graph, "l");
        mg->Add(expSn124Graph, "l");

        mg->Add(litSn112Graph, "p");
        mg->Add(litSn124Graph, "l");

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
        mg->GetYaxis()->SetNdivisions(4);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetYaxis()->SetRangeUser(1.4,5.49);
        mg->GetXaxis()->SetLimits(3,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.295,0.735,"Pb");
        //latex.DrawLatex(0.315,0.54,"Sn");
        //latex.DrawLatex(0.325,0.367,"Ni");
        //latex.DrawLatex(0.33,0.235,"C");

        latex.SetTextSize(0.15);
        latex.DrawLatex(0.02, 1, "a)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.52,0.73,0.97,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litSn112Graph,"^{112}Sn (An)","p");
        legend->AddEntry(litSn124Graph,"^{124}Sn (An)","l");
        legend->AddEntry(expSn112Graph,"^{112}Sn (DSP)","l");
        legend->AddEntry(expSn124Graph,"^{124}Sn (DSP)","l");

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

        mg->Add(relSn124Graph, "lp");
        mg->Add(relSn112Graph, "p");

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
        mg->GetYaxis()->SetNdivisions(5);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetXaxis()->SetLimits(3,500);
        mg->GetYaxis()->SetRangeUser(-4.9,4.9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.13);
        latex.DrawLatex(0.02, 1, "b)");

        TLine* zeroLine = new TLine(3, 0, 500, 0);
        zeroLine->SetLineColor(kBlack);
        zeroLine->SetLineWidth(3);
        zeroLine->SetLineStyle(9);
        zeroLine->Draw();

        // Define legend format and contents
        TLegend *legend = new TLegend(0.83,0.74,0.96,0.96);
        legend->SetNColumns(1);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relSn112Graph,"{}^{112}Sn","p");
        legend->AddEntry(relSn124Graph,"{}^{124}Sn","pl");

        legend->Draw();

    }

    // close it all up
    expFile->Close();
    relFile->Close();
    litFile->Close();
}
