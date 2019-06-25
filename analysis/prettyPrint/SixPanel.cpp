void CanvasPartition(TCanvas *C, const Int_t Nx=2, const Int_t Ny=2,
        Float_t lMargin=0.15, Float_t rMargin=0.05,
        Float_t bMargin=0.15, Float_t tMargin=0.05);

void SixPanel() {
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas", 2550, 850);
    //canvas->Divide(3, 2, 0.2, 0);

    int xPads = 3;
    int yPads = 2;

    double lMargin = 0.06;
    double rMargin = 0.0;
    double bMargin = 0.10;
    double tMargin = 0.0;

    CanvasPartition(canvas, xPads, yPads, lMargin, rMargin, bMargin, tMargin);

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
    string expFile2Name2 = "/data2/analysis/total.root";
    string relFile2Name2 = "/data2/analysis/relative.root";
    string litFile2Name2 = "/data2/analysis/literatureData.root";
    string shiftedFile2Name2 = "/data2/analysis/shifted.root";

    TFile* expFile2 = new TFile(expFile2Name2.c_str(),"READ");
    TFile* relFile2 = new TFile(relFile2Name2.c_str(),"READ");
    TFile* litFile2 = new TFile(litFile2Name2.c_str(),"READ");
    TFile* shiftedFile2 = new TFile(shiftedFile2Name2.c_str(),"READ");

    if(!expFile2 || !relFile2 || !litFile2 || !shiftedFile2)
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
    
    TGraphAsymmErrors* expO16Graph = (TGraphAsymmErrors*)expFile2->Get(expO16GraphName.c_str());
    TGraphAsymmErrors* expO18Graph = (TGraphAsymmErrors*)shiftedFile2->Get(expO18GraphName.c_str());

    TGraphAsymmErrors* relO16Graph = (TGraphAsymmErrors*)relFile2->Get(relO16GraphName.c_str());
    TGraphAsymmErrors* relO18Graph = (TGraphAsymmErrors*)relFile2->Get(relO18GraphName.c_str());

    TGraphAsymmErrors* litO16Graph = (TGraphAsymmErrors*)litFile2->Get(litO16GraphName.c_str());
    TGraphAsymmErrors* litO18Graph = (TGraphAsymmErrors*)litFile2->Get(litO18GraphName.c_str());

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
    relO18Graph->SetMarkerStyle(20);

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
    relNi58Graph->SetMarkerStyle(22);
    relNi58Graph->SetMarkerSize(2.5);

    relNi64Graph->SetLineColor(kRed+2);
    relNi64Graph->SetMarkerColor(kRed+2);
    relNi64Graph->SetMarkerSize(2.5);
    relNi64Graph->SetMarkerStyle(20);

    litNi58Graph->SetLineColor(kBlue-7);
    litNi58Graph->SetLineWidth(4);
    litNi58Graph->SetLineStyle(1);
    litNi58Graph->SetMarkerColor(kBlue-7);

    litNi64Graph->SetMarkerColor(kBlue);
    litNi64Graph->SetMarkerSize(3);
    litNi64Graph->SetMarkerStyle(22);

    string expSn112GraphName = "Sn112";
    string expSn124GraphName = "Sn124";

    string relSn112GraphName = "Sn112, expLit, percent";
    string relSn124GraphName = "Sn124, expLit, percent";

    string litSn112GraphName = "Sn112 (n,tot)";
    string litSn124GraphName = "Sn124 (n,tot)";
    
    TGraphAsymmErrors* expSn112Graph = (TGraphAsymmErrors*)expFile2->Get(expSn112GraphName.c_str());
    TGraphAsymmErrors* expSn124Graph = (TGraphAsymmErrors*)expFile2->Get(expSn124GraphName.c_str());

    TGraphAsymmErrors* relSn112Graph = (TGraphAsymmErrors*)relFile2->Get(relSn112GraphName.c_str());
    TGraphAsymmErrors* relSn124Graph = (TGraphAsymmErrors*)relFile2->Get(relSn124GraphName.c_str());

    TGraphAsymmErrors* litSn112Graph = (TGraphAsymmErrors*)litFile2->Get(litSn112GraphName.c_str());
    TGraphAsymmErrors* litSn124Graph = (TGraphAsymmErrors*)litFile2->Get(litSn124GraphName.c_str());

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
    relSn124Graph->SetMarkerStyle(20);

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
    litSn124Graph->SetMarkerStyle(20);

    // fourth panel
    {
        canvas->cd(0);

        TPad* pad = (TPad*)gROOT->FindObject("pad_0_1");
        pad->Draw();
        pad->cd();

        // Pad dimensions and margins
        pad->SetLeftMargin(0.18);
        pad->SetRightMargin(0.0);
        pad->SetTopMargin(0.03);
        pad->SetBottomMargin(0.0);
        pad->SetTicky(1);
        pad->SetTickx(1);
        pad->SetLogx(1);

        pad->SetFrameLineWidth(3);

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

        mg->GetYaxis()->SetRangeUser(0.1,5.59);
        mg->GetXaxis()->SetLimits(2.8,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top

        latex.DrawLatex(0.35, 0.68, "+1 barn");

        latex.SetTextSize(0.08);
        latex.DrawLatex(0.21, 0.93, "a)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.55,0.70,0.97,0.92);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litO16Graph,"^{16}O (An)","l");
        legend->AddEntry(expO16Graph,"^{16}O (DSP)","l");
        legend->AddEntry(litO18Graph,"^{18}O (An)","l");
        legend->AddEntry(expO18Graph,"^{18}O (DSP)","l");

        legend->Draw();
    }

    // first panel
    {
        canvas->cd(0);

        TPad* pad = (TPad*)gROOT->FindObject("pad_0_0");
        pad->Draw();
        pad->cd();

        // Pad dimensions and margins
        pad->SetLeftMargin(0.18);
        pad->SetRightMargin(0.0);
        pad->SetTopMargin(0.0);
        pad->SetBottomMargin(0.25);
        pad->SetTicky(1);
        pad->SetTickx(1);
        pad->SetLogx();

        pad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        double zeroX[2] = {3, 500};
        double zeroY[2] = {0, 0};
        TGraph* zeroGraph = new TGraph(2, zeroX, zeroY);
        zeroGraph->SetLineColor(kBlack);
        zeroGraph->SetLineWidth(3);
        zeroGraph->SetLineStyle(9);

        mg->Add(zeroGraph, "l");
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
        mg->GetYaxis()->SetTitleSize(0.07);
        mg->GetYaxis()->SetTitleFont(2);
        mg->GetYaxis()->SetTitleOffset(1.1);
        mg->GetYaxis()->CenterTitle();

        mg->GetYaxis()->SetLabelOffset(0.01);
        mg->GetYaxis()->SetLabelSize(0.07);

        mg->GetYaxis()->SetLabelFont(2);
        mg->GetYaxis()->SetNdivisions(8);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetXaxis()->SetLimits(2.8,500);
        mg->GetYaxis()->SetRangeUser(-14.9,14.9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.07);
        latex.DrawLatex(0.21, 0.93, "b)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.83,0.74,0.96,0.96);
        legend->SetNColumns(1);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relO16Graph,"{}^{16}O","lp");
        legend->AddEntry(relO18Graph,"{}^{18}O","lp");

        legend->Draw();
    }

    // fifth panel
    {
        canvas->cd(0);

        TPad* pad = (TPad*)gROOT->FindObject("pad_1_1");
        pad->Draw();
        pad->cd();

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.0);
        gPad->SetRightMargin(0.0);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.0);
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
        mg->GetYaxis()->SetNdivisions(8);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetYaxis()->SetRangeUser(0.1,5.59);
        mg->GetXaxis()->SetLimits(2.8,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.295,0.735,"Pb");
        //latex.DrawLatex(0.315,0.54,"Sn");
        //latex.DrawLatex(0.325,0.367,"Ni");
        //latex.DrawLatex(0.33,0.235,"C");

        latex.SetTextSize(0.08);
        latex.DrawLatex(0.03, 0.93, "c)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.45,0.70,0.97,0.92);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litNi58Graph,"^{58}Ni (An)","l");
        legend->AddEntry(expNi58Graph,"^{58}Ni (DSP)","l");
        legend->AddEntry(litNi64Graph,"^{64}Ni (An)","p");
        legend->AddEntry(expNi64Graph,"^{64}Ni (DSP)","l");

        legend->Draw();
    }

    // second panel
    {
        canvas->cd(0);

        TPad* pad = (TPad*)gROOT->FindObject("pad_1_0");
        pad->Draw();
        pad->cd();

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.0);
        gPad->SetRightMargin(0.0);
        gPad->SetTopMargin(0.0);
        gPad->SetBottomMargin(0.25);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx();

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        double zeroX[2] = {3, 500};
        double zeroY[2] = {0, 0};
        TGraph* zeroGraph = new TGraph(2, zeroX, zeroY);
        zeroGraph->SetLineColor(kBlack);
        zeroGraph->SetLineWidth(3);
        zeroGraph->SetLineStyle(9);

        mg->Add(zeroGraph, "l");
        mg->Add(relNi58Graph, "lp");
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
        mg->GetYaxis()->SetNdivisions(8);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetXaxis()->SetLimits(2.8,500);
        mg->GetYaxis()->SetRangeUser(-14.9,14.9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.07);
        latex.DrawLatex(0.03, 0.93, "d)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.78,0.74,0.96,0.96);
        legend->SetNColumns(1);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relNi58Graph,"{}^{58}Ni","lp");
        legend->AddEntry(relNi64Graph,"{}^{64}Ni","p");

        legend->Draw();
    }

    // sixth panel
    {
        canvas->cd(0);

        TPad* pad = (TPad*)gROOT->FindObject("pad_2_1");
        pad->Draw();
        pad->cd();


        // Pad dimensions and margins
        gPad->SetLeftMargin(0.0);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.0);
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
        mg->GetYaxis()->SetNdivisions(8);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetYaxis()->SetRangeUser(0.1,5.59);
        mg->GetXaxis()->SetLimits(2.8,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.295,0.735,"Pb");
        //latex.DrawLatex(0.315,0.54,"Sn");
        //latex.DrawLatex(0.325,0.367,"Ni");
        //latex.DrawLatex(0.33,0.235,"C");

        latex.SetTextSize(0.08);
        latex.DrawLatex(0.03, 0.93, "e)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.03,0.07,0.64,0.27);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litSn112Graph,"^{112}Sn (An)","p");
        legend->AddEntry(expSn112Graph,"^{112}Sn (DSP)","l");
        legend->AddEntry(litSn124Graph,"^{124}Sn (An)","l");
        legend->AddEntry(expSn124Graph,"^{124}Sn (DSP)","l");

        legend->Draw();
    }

    // third panel
    {
        canvas->cd(0);

        TPad* pad = (TPad*)gROOT->FindObject("pad_2_0");
        pad->Draw();
        pad->cd();

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.0);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.0);
        gPad->SetBottomMargin(0.25);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx();

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        double zeroX[2] = {3, 500};
        double zeroY[2] = {0, 0};
        TGraph* zeroGraph = new TGraph(2, zeroX, zeroY);
        zeroGraph->SetLineColor(kBlack);
        zeroGraph->SetLineWidth(3);
        zeroGraph->SetLineStyle(9);

        mg->Add(zeroGraph, "l");
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
        mg->GetYaxis()->SetNdivisions(8);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetXaxis()->SetLimits(2.8,500);
        mg->GetYaxis()->SetRangeUser(-14.9,14.9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.07);
        latex.DrawLatex(0.03, 0.93, "f)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.75,0.74,0.94,0.96);
        legend->SetNColumns(1);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relSn112Graph,"{}^{112}Sn","p");
        legend->AddEntry(relSn124Graph,"{}^{124}Sn","pl");

        legend->Draw();
    }

    // close it all up
    expFile2->Close();
    relFile2->Close();
    litFile2->Close();
    shiftedFile2->Close();

    expFile->Close();
    relFile->Close();
    litFile->Close();
}

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   for (Int_t i=0;i<Nx;i++) {

      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j=0;j<Ny;j++) {

         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);

         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);

         pad->Draw();
      }
   }
}
