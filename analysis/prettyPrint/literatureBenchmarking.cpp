{
    TStyle* style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas",850,850);
    canvas->Divide(1, 2, 0, 0, 0);

    style->SetOptStat(0);
    style->SetOptTitle(0);    
    style->SetPalette(1,0);
    style->SetCanvasColor(10);      
    style->SetCanvasBorderMode(0);    
    style->SetFrameLineWidth(3);
    style->SetFrameFillColor(10);
    style->SetPadColor(10);
    style->SetHistLineWidth(3);
    style->SetHistLineColor(kRed-4);
    style->SetMarkerSize(0);
    style->SetMarkerStyle(8);
    style->SetFuncWidth(3);
    style->SetFuncColor(kRed);
    style->SetLabelColor(kBlack,"xyz");
    style->SetTitleSize(0.06,"xyz");
    style->SetTitleFillColor(10);
    style->SetTitleTextColor(kBlack);
    style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    string fileName = "/data1/analysis/literatureData.root";
    string expFileName = "/data1/analysis/total.root";
    string relFileName = "/data1/analysis/relative.root";
    string SnFileName = "/data2/analysis/total.root";
    string relSnFileName = "/data2/analysis/relative.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* SnFile = new TFile(SnFileName.c_str(),"READ");
    TFile* relSnFile = new TFile(relSnFileName.c_str(),"READ");

    if(!file || !expFile || !SnFile || !relFile || !relSnFile)
    {
        cout << "Error: couldn't open a root file needed for extracting graphs." << endl;
        exit(1);
    }
 
    string CGraphName = "Natural C (n,tot)";
    string NiGraphName = "Natural Ni (n,tot)";
    string SnGraphName = "Natural Sn (n,tot)";
    string PbGraphName = "Natural Pb (n,tot)";
        
    TGraphAsymmErrors* CGraph = (TGraphAsymmErrors*)file->Get(CGraphName.c_str());
    TGraphAsymmErrors* NiGraph = (TGraphAsymmErrors*)file->Get(NiGraphName.c_str());
    TGraphAsymmErrors* SnGraph = (TGraphAsymmErrors*)file->Get(SnGraphName.c_str());
    TGraphAsymmErrors* PbGraph = (TGraphAsymmErrors*)file->Get(PbGraphName.c_str());

    if(!CGraph || !NiGraph || !SnGraph || !PbGraph)
    {
        cout << "Error: couldn't open a literature graph." << endl;
        exit(1);
    }

    string expCGraphName = "CNat";
    string expNiGraphName = "NiNat";
    string expSnGraphName = "SnNat";
    string expPbGraphName = "PbNat";

    TGraphAsymmErrors* expCGraph = (TGraphAsymmErrors*)expFile->Get(expCGraphName.c_str());
    TGraphAsymmErrors* expNiGraph = (TGraphAsymmErrors*)expFile->Get(expNiGraphName.c_str());
    TGraphAsymmErrors* expSnGraph = (TGraphAsymmErrors*)SnFile->Get(expSnGraphName.c_str());
    TGraphAsymmErrors* expPbGraph = (TGraphAsymmErrors*)expFile->Get(expPbGraphName.c_str());

    if(!expCGraph || !expNiGraph || !expSnGraph || !expPbGraph)
    {
        cout << "Error: couldn't open an experimental graph." << endl;
        exit(1);
    }

    string relCGraphName = "CNat, expLit, percent";
    string relNiGraphName = "NiNat, expLit, percent";
    string relSnGraphName = "SnNat, expLit, percent";
    string relPbGraphName = "PbNat, expLit, percent";

    TGraphAsymmErrors* relCGraph = (TGraphAsymmErrors*)relFile->Get(relCGraphName.c_str());
    TGraphAsymmErrors* relNiGraph = (TGraphAsymmErrors*)relFile->Get(relNiGraphName.c_str());
    TGraphAsymmErrors* relSnGraph = (TGraphAsymmErrors*)relSnFile->Get(relSnGraphName.c_str());
    TGraphAsymmErrors* relPbGraph = (TGraphAsymmErrors*)relFile->Get(relPbGraphName.c_str());

    // Set graph point and line characteristics
    CGraph->SetLineWidth(4);
    CGraph->SetLineStyle(0);
    CGraph->SetLineColor(kBlue);
    expCGraph->SetLineWidth(4);
    expCGraph->SetLineStyle(0);
    expCGraph->SetLineColor(kRed);
    
    NiGraph->SetLineWidth(4);
    NiGraph->SetLineStyle(0);
    NiGraph->SetLineColor(kBlue);
    expNiGraph->SetLineWidth(4);
    expNiGraph->SetLineStyle(0);
    expNiGraph->SetLineColor(kRed);

    SnGraph->SetLineWidth(4);
    SnGraph->SetLineStyle(0);
    SnGraph->SetLineColor(kBlue);
    expSnGraph->SetLineWidth(4);
    expSnGraph->SetLineStyle(0);
    expSnGraph->SetLineColor(kRed);

    PbGraph->SetLineWidth(4);
    PbGraph->SetLineStyle(0);
    PbGraph->SetLineColor(kBlue);
    expPbGraph->SetLineWidth(4);
    expPbGraph->SetLineStyle(0);
    expPbGraph->SetLineColor(kRed);

    relCGraph->SetLineColor(kRed-9);
    relCGraph->SetLineWidth(4);
    relCGraph->SetLineStyle(0);
    relCGraph->SetMarkerColor(kRed-9);

    relNiGraph->SetLineColor(kRed-3);
    relNiGraph->SetLineWidth(4);
    relNiGraph->SetLineStyle(0);
    relNiGraph->SetMarkerColor(kRed-3);

    relSnGraph->SetLineColor(kRed);
    relSnGraph->SetLineWidth(4);
    relSnGraph->SetLineStyle(0);
    relSnGraph->SetMarkerColor(kRed);

    relPbGraph->SetLineColor(kRed+3);
    relPbGraph->SetLineWidth(4);
    relPbGraph->SetLineStyle(0);
    relPbGraph->SetMarkerColor(kRed+3);

    // first panel
    {
        canvas->cd(1);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.17);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.005);
        gPad->SetBottomMargin(0.0);
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        mg->Add(CGraph, "l");
        mg->Add(NiGraph, "l");
        mg->Add(SnGraph, "l");
        mg->Add(PbGraph, "l");

        mg->Add(expCGraph, "l");
        mg->Add(expNiGraph, "l");
        mg->Add(expSnGraph, "l");
        mg->Add(expPbGraph, "l");

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
        mg->GetYaxis()->SetNdivisions(10);
        mg->GetYaxis()->SetTickLength(0.02);

        mg->GetYaxis()->SetRangeUser(0.01,8.99);
        mg->GetXaxis()->SetLimits(3,500);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.295,0.735,"Pb");
        latex.DrawLatex(0.295,0.55,"Sn");
        latex.DrawLatex(0.315,0.347,"Ni");
        latex.DrawLatex(0.33,0.235,"C");

        latex.SetTextSize(0.10);
        latex.DrawLatex(0.20, 0.13, "(a)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.76,0.70,0.94,0.9);
        legend->SetTextSize(0.08);
        legend->SetTextAlign(12);
        legend->AddEntry(CGraph,"Analog","l");
        legend->AddEntry(expCGraph,"DSP","l");
        legend->Draw();
    }

    // second panel
    {
        canvas->cd(2);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.17);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.0);
        gPad->SetBottomMargin(0.25);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx();

        gPad->SetFrameLineWidth(3);

        TMultiGraph* mg = new TMultiGraph();

        mg->Add(relCGraph);
        mg->Add(relNiGraph);
        mg->Add(relSnGraph);
        mg->Add(relPbGraph);

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
        mg->GetYaxis()->SetRangeUser(-8,7.9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top

        latex.SetTextSize(0.09);
        latex.DrawLatex(0.20, 0.35, "(b)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.70,0.74,0.96,0.96);
        legend->SetNColumns(2);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relCGraph,"{}^{nat}C","l");
        legend->AddEntry(relNiGraph,"{}^{nat}Ni","l");
        legend->AddEntry(relSnGraph,"{}^{nat}Sn","l");
        legend->AddEntry(relPbGraph,"{}^{nat}Pb","l");

        legend->Draw();

        TLine* zeroLine = new TLine(3, 0, 500, 0);
        zeroLine->SetLineColor(kBlack);
        zeroLine->SetLineWidth(3);
        zeroLine->SetLineStyle(9);
        zeroLine->Draw();
    }

    file->Close();
    expFile->Close();
    relFile->Close();
    SnFile->Close();
    relSnFile->Close();
}
