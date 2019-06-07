{
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas",650, 1300);
    canvas->Divide(1, 3, 0, 0);

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
    string ramsauerFileName = "../../theory/ramsauer.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(), "READ");

    if(!expFile || !relFile || !ramsauerFile)
    {
        cout << "Error: failed to open a required file (check filenames). Exiting..." << endl;
        exit(1);
    }

    string expNi58GraphName = "Ni58";
    string SANi58GraphName = "SA_A=58";
    string RamsauerNi58GraphName = "Ramsauer_A=58";

    string expNi64GraphName = "Ni64PurityCorrected";
    string SANi64GraphName = "SA_A=64";
    string RamsauerNi64GraphName = "Ramsauer_A=64";

    string relGraphName = "Ni64Ni58, percent";
    string relGraphSEName = "Ni64Ni58SysErrors, percent";
    string SARelDiffGraphName = "RelDiff64_58";
    string RamsauerRelDiffGraphName = "RelDiffRamsauer64_58";

    TGraphAsymmErrors* expNi58Graph = (TGraphAsymmErrors*)expFile->Get(expNi58GraphName.c_str());
    TGraphAsymmErrors* expNi64Graph = (TGraphAsymmErrors*)expFile->Get(expNi64GraphName.c_str());

    TGraph* SANi58Graph = (TGraphAsymmErrors*)ramsauerFile->Get(SANi58GraphName.c_str());
    TGraph* SANi64Graph = (TGraphAsymmErrors*)ramsauerFile->Get(SANi64GraphName.c_str());

    TGraph* RamsauerNi58Graph = (TGraphAsymmErrors*)ramsauerFile->Get(RamsauerNi58GraphName.c_str());
    TGraph* RamsauerNi64Graph = (TGraphAsymmErrors*)ramsauerFile->Get(RamsauerNi64GraphName.c_str());

    TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)relFile->Get(relGraphName.c_str());
    TGraphAsymmErrors* relGraphSE = (TGraphAsymmErrors*)relFile->Get(relGraphSEName.c_str());

    TGraph* SARelDiffGraph = (TGraphAsymmErrors*)ramsauerFile->Get(SARelDiffGraphName.c_str());
    TGraph* RamsauerRelDiffGraph = (TGraphAsymmErrors*)ramsauerFile->Get(RamsauerRelDiffGraphName.c_str());

    if(!expNi58Graph || !expNi64Graph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    if(!SANi58Graph || !SANi64Graph || !RamsauerNi58Graph || !RamsauerNi64Graph)
    {
        cout << "Error: failed to open a SA or Ramsauer absolute cross section graph." << endl;
        exit(1);
    }

    if(!relGraph || !relGraphSE)
    {
        cout << "Error: failed to open an experimental relative difference cross section graph." << endl;
        exit(1);
    }

    if(!SARelDiffGraph || !RamsauerRelDiffGraph)
    {
        cout << "Error: failed to open SA or Ramsauer relative difference cross section graph." << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expNi58Graph->SetLineColor(kRed);
    expNi58Graph->SetLineWidth(4);
    expNi58Graph->SetLineStyle(1);
    expNi58Graph->SetMarkerColor(kRed);

    expNi64Graph->SetLineColor(kRed);
    expNi64Graph->SetLineWidth(4);
    expNi64Graph->SetLineStyle(1);
    expNi64Graph->SetMarkerColor(kRed);

    SANi58Graph->SetLineColor(kBlack);
    SANi58Graph->SetLineWidth(4);
    SANi58Graph->SetLineStyle(9);

    SANi64Graph->SetLineColor(kBlack);
    SANi64Graph->SetLineWidth(4);
    SANi64Graph->SetLineStyle(9);

    RamsauerNi58Graph->SetLineColor(kGray+2);
    RamsauerNi58Graph->SetLineWidth(4);
    RamsauerNi58Graph->SetLineStyle(4);

    RamsauerNi64Graph->SetLineColor(kGray+2);
    RamsauerNi64Graph->SetLineWidth(4);
    RamsauerNi64Graph->SetLineStyle(4);

    relGraph->SetLineColor(kRed);
    relGraph->SetLineWidth(5);
    relGraph->SetLineStyle(0);
    relGraph->SetMarkerColor(kRed);
    relGraph->SetFillColor(kRed);
    relGraph->SetFillStyle(3002);

    SARelDiffGraph->SetLineStyle(9);
    SARelDiffGraph->SetLineWidth(3);
    SARelDiffGraph->SetLineColor(kBlack);

    RamsauerRelDiffGraph->SetLineStyle(7);
    RamsauerRelDiffGraph->SetLineWidth(3);
    RamsauerRelDiffGraph->SetLineColor(kGray+2);

    relGraphSE->SetFillColor(kBlue);
    relGraphSE->SetFillStyle(3002);

    // first panel
    {
        canvas->cd(1);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.005);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        mg->Add(expNi58Graph, "l");
        mg->Add(SANi58Graph, "l");
        mg->Add(RamsauerNi58Graph, "l");

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
        latex.DrawLatex(0.27, 0.25, "Ni^{58}");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.55,0.73,0.97,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(expNi58Graph,"Exp data","l");
        legend->AddEntry(SANi58Graph,"SAS","l");
        legend->AddEntry(RamsauerNi58Graph,"Ramsauer","l");

        legend->Draw();
    }

    // second panel
    {
        canvas->cd(2);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.005);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        mg->Add(expNi64Graph, "l");
        mg->Add(SANi64Graph, "l");
        mg->Add(RamsauerNi64Graph, "l");

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
        latex.DrawLatex(0.27, 0.25, "Ni^{64}");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.55,0.73,0.97,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(expNi64Graph,"Exp data","l");
        legend->AddEntry(SANi64Graph,"SAS","l");
        legend->AddEntry(RamsauerNi64Graph,"Ramsauer","l");

        legend->Draw();
    }

    // third panel
    {
        canvas->cd(3);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.005);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        mg->Add(relGraph,"3l");
        mg->Add(relGraphSE, "3");
        mg->Add(SARelDiffGraph, "l");
        mg->Add(RamsauerRelDiffGraph, "l");

        mg->Draw("al");

        // X-axis parameters
        mg->GetXaxis()->SetTitle("Energy (MeV)");
        mg->GetXaxis()->SetTitleSize(0.05);
        mg->GetXaxis()->SetTitleFont(2);
        mg->GetXaxis()->SetTitleOffset(1.4);
        mg->GetXaxis()->CenterTitle();

        mg->GetXaxis()->SetLabelOffset(0.01);
        mg->GetXaxis()->SetLabelSize(0.05);
        mg->GetXaxis()->SetLabelFont(2);

        mg->GetXaxis()->SetNdivisions(10);
        mg->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        mg->GetYaxis()->SetTitle("(#frac{#sigma_{64} - #sigma_{58}}{#sigma_{64} + #sigma_{58}})");
        mg->GetYaxis()->SetTitleSize(0.06);
        mg->GetYaxis()->SetTitleFont(2);
        mg->GetYaxis()->SetTitleOffset(1.0);
        mg->GetYaxis()->CenterTitle();

        mg->GetYaxis()->SetLabelOffset(0.01);
        mg->GetYaxis()->SetLabelSize(0.05);

        mg->GetYaxis()->SetLabelFont(2);
        mg->GetYaxis()->SetNdivisions(5);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetYaxis()->SetRangeUser(0.9,5.1);
        mg->GetXaxis()->SetLimits(3,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.295,0.735,"Pb");
        //latex.DrawLatex(0.315,0.54,"Sn");
        //latex.DrawLatex(0.325,0.367,"Ni");
        //latex.DrawLatex(0.33,0.235,"C");

        latex.SetTextSize(0.10);
        latex.DrawLatex(0.23, 0.18, "Rel. Diff.");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.53, 0.83, 0.95, 0.95);
        legend->SetNColumns(2);
        legend->AddEntry(relGraph,"Exp data, sys + stat","f");
        legend->AddEntry(SARelDiffGraph,"SAS","l");
        legend->AddEntry(relGraphSE,"Exp data, sys only","f");
        legend->AddEntry(RamsauerRelDiffGraph,"Ramsauer","l");
        legend->Draw();
    }

    // close it all up
    expFile->Close();
    ramsauerFile->Close();
    relFile->Close();
}
