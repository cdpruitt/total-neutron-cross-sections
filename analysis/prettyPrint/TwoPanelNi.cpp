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

    string expNi58GraphName = "Ni58";
    string expNi64GraphName = "Ni64PurityCorrected";

    string relNi58GraphName = "Ni58, expLit, percent";
    string relNi64GraphName = "Ni64, expLit, percent";

    string litNi58GraphName = "Ni58 (n,tot)";
    string litNi64GraphName = "Ni64 (n,tot)";
    
    TGraphAsymmErrors* expNi58Graph = (TGraphAsymmErrors*)expFile->Get(expNi58GraphName.c_str());
    TGraphAsymmErrors* expNi64Graph = (TGraphAsymmErrors*)expFile->Get(expNi64GraphName.c_str());

    TGraphAsymmErrors* relNi58Graph = (TGraphAsymmErrors*)relFile->Get(relNi58GraphName.c_str());
    TGraphAsymmErrors* relNi64Graph = (TGraphAsymmErrors*)relFile->Get(relNi64GraphName.c_str());

    TGraphAsymmErrors* litNi58Graph = (TGraphAsymmErrors*)litFile->Get(litNi58GraphName.c_str());
    TGraphAsymmErrors* litNi64Graph = (TGraphAsymmErrors*)litFile->Get(litNi64GraphName.c_str());

    if(!expNi58Graph || !expNi64Graph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    if(!relNi58Graph || !relNi64Graph)
    {
        cout << "Error: failed to open an relative diff to lit cross section graph." << endl;
        exit(1);
    }

    if(!litNi58Graph || !litNi64Graph)
    {
        cout << "Error: failed to open an lit cross section graph." << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expNi58Graph->SetLineColor(kRed);
    expNi58Graph->SetLineWidth(4);
    expNi58Graph->SetLineStyle(1);
    expNi58Graph->SetMarkerColor(kRed);

    expNi64Graph->SetLineColor(kRed+2);
    expNi64Graph->SetLineWidth(4);
    expNi64Graph->SetLineStyle(1);
    expNi64Graph->SetMarkerColor(kRed+2);

    relNi58Graph->SetLineColor(kRed);
    relNi58Graph->SetLineWidth(4);
    relNi58Graph->SetLineStyle(0);
    relNi58Graph->SetMarkerColor(kRed);

    relNi64Graph->SetLineColor(kRed+2);
    relNi64Graph->SetMarkerColor(kRed+2);
    relNi64Graph->SetMarkerSize(2);
    relNi64Graph->SetMarkerStyle(22);

    litNi58Graph->SetLineColor(kBlue-7);
    litNi58Graph->SetLineWidth(4);
    litNi58Graph->SetLineStyle(1);
    litNi58Graph->SetMarkerColor(kBlue-7);

    litNi64Graph->SetMarkerColor(kBlue);
    litNi64Graph->SetMarkerSize(3);
    litNi64Graph->SetMarkerStyle(22);

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

        mg->Add(litNi58Graph, "l");

        mg->Add(expNi58Graph, "l");
        mg->Add(expNi64Graph, "l");

        mg->Add(litNi64Graph, "p");

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

        mg->GetYaxis()->SetRangeUser(0.81,4.69);
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
        TLegend *legend = new TLegend(0.55,0.73,0.97,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litNi58Graph,"^{58}Ni (An)","l");
        legend->AddEntry(litNi64Graph,"^{64}Ni (An)","p");
        legend->AddEntry(expNi58Graph,"^{58}Ni (DSP)","l");
        legend->AddEntry(expNi64Graph,"^{64}Ni (DSP)","l");

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

        mg->Add(relNi58Graph, "l");
        mg->Add(relNi64Graph, "p");

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

        mg->GetXaxis()->SetLimits(3,500);
        mg->GetYaxis()->SetRangeUser(-14.9,14.9);

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
        legend->AddEntry(relNi58Graph,"{}^{58}Ni","l");
        legend->AddEntry(relNi64Graph,"{}^{64}Ni","p");

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
