{
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas",650,2000);
    canvas->SetCanvasSize(600,1500);
    canvas->Divide(1,4);

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
    style->SetTitleSize(0.06,"xyz");
    style->SetTitleFillColor(10);
    style->SetTitleTextColor(kBlack);
    style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    // read graphs
    string expFileName = "/data1/analysis/total.root";
    string relFileName = "/data1/analysis/relative.root";
    string corFileName = "/data1/analysis/corrected.root";
    string litFileName = "/data1/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* corFile = new TFile(corFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");
    
    string expCGraphName = "CNat";
    string expNiGraphName = "NiNat";
    string expPbGraphName = "PbNat";

    string relCGraphName = "NatC, expLit";
    string relNiGraphName = "NatNi, expLit";
    string relPbGraphName = "NatPb, expLit";

    string litCGraphName = "Natural carbon (n,tot)";
    string litNiGraphName = "Natural Ni (n,tot)";
    string litPbGraphName = "Natural Pb (n,tot)";
    
    TGraphErrors* expCGraph = (TGraphErrors*)expFile->Get(expCGraphName.c_str());
    TGraphErrors* expNiGraph = (TGraphErrors*)expFile->Get(expNiGraphName.c_str());
    TGraphErrors* expPbGraph = (TGraphErrors*)expFile->Get(expPbGraphName.c_str());

    TGraphErrors* relCGraph = (TGraphErrors*)relFile->Get(relCGraphName.c_str());
    TGraphErrors* relNiGraph = (TGraphErrors*)relFile->Get(relNiGraphName.c_str());
    TGraphErrors* relPbGraph = (TGraphErrors*)relFile->Get(relPbGraphName.c_str());

    TGraphErrors* corCGraph = (TGraphErrors*)corFile->Get(expCGraphName.c_str());
    TGraphErrors* corNiGraph = (TGraphErrors*)corFile->Get(expNiGraphName.c_str());
    TGraphErrors* corPbGraph = (TGraphErrors*)corFile->Get(expPbGraphName.c_str());

    TGraphErrors* litCGraph = (TGraphErrors*)litFile->Get(litCGraphName.c_str());
    TGraphErrors* litNiGraph = (TGraphErrors*)litFile->Get(litNiGraphName.c_str());
    TGraphErrors* litPbGraph = (TGraphErrors*)litFile->Get(litPbGraphName.c_str());

    // Set graph point and line characteristics
    expCGraph->SetLineColor(kRed);
    expCGraph->SetMarkerColor(kRed);

    expNiGraph->SetLineColor(kRed);
    expNiGraph->SetMarkerColor(kRed);

    expPbGraph->SetLineColor(kRed);
    expPbGraph->SetMarkerColor(kRed);

    relCGraph->SetLineColor(kRed-9);
    relCGraph->SetLineWidth(4);
    relCGraph->SetLineStyle(0);
    relCGraph->SetMarkerColor(kRed-9);

    relNiGraph->SetLineColor(kRed);
    relNiGraph->SetLineWidth(4);
    relNiGraph->SetLineStyle(0);
    relNiGraph->SetMarkerColor(kRed);

    relPbGraph->SetLineColor(kRed+3);
    relPbGraph->SetLineWidth(4);
    relPbGraph->SetLineStyle(0);
    relPbGraph->SetMarkerColor(kRed+3);

    corCGraph->SetLineColor(kRed);
    corCGraph->SetMarkerColor(kRed);

    corNiGraph->SetLineColor(kRed);
    corNiGraph->SetMarkerColor(kRed);

    corPbGraph->SetLineColor(kRed);
    corPbGraph->SetMarkerColor(kRed);

    litCGraph->SetLineColor(kBlack);
    litCGraph->SetLineWidth(3);
    litCGraph->SetLineStyle(2);
    litCGraph->SetMarkerColor(kBlack);

    litNiGraph->SetLineColor(kBlack);
    litNiGraph->SetLineWidth(3);
    litNiGraph->SetLineStyle(2);
    litNiGraph->SetMarkerColor(kBlack);

    litPbGraph->SetLineColor(kBlack);
    litPbGraph->SetLineWidth(3);
    litPbGraph->SetLineStyle(2);
    litPbGraph->SetMarkerColor(kBlack);

    {
        // X-axis parameters
        expCGraph->GetXaxis()->SetTitle("Energy (MeV)");
        expCGraph->GetXaxis()->SetTitleSize(0.05);
        expCGraph->GetXaxis()->SetTitleFont(2);
        expCGraph->GetXaxis()->SetTitleOffset(1.4);
        expCGraph->GetXaxis()->CenterTitle();

        expCGraph->GetXaxis()->SetLabelOffset(0.01);
        expCGraph->GetXaxis()->SetLabelSize(0.05);
        expCGraph->GetXaxis()->SetLabelFont(2);

        expCGraph->GetXaxis()->SetNdivisions(10);
        expCGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        expCGraph->GetYaxis()->SetTitle("#sigma_{tot} (barns)");
        expCGraph->GetYaxis()->SetTitleSize(0.06);
        expCGraph->GetYaxis()->SetTitleFont(2);
        expCGraph->GetYaxis()->SetTitleOffset(0.8);
        expCGraph->GetYaxis()->CenterTitle();

        expCGraph->GetYaxis()->SetLabelOffset(0.01);
        expCGraph->GetYaxis()->SetLabelSize(0.05);

        expCGraph->GetYaxis()->SetLabelFont(2);
        expCGraph->GetYaxis()->SetNdivisions(10);
        expCGraph->GetYaxis()->SetTickLength(0.02);

        // first panel
        canvas->cd(1);
        //expCGraph->Draw("");
        //expNiGraph->Draw("same");
        //expPbGraph->Draw("same");
        litCGraph->Draw("");
        litNiGraph->Draw("same");
        litPbGraph->Draw("same");
        expCGraph->Draw("same");
        expNiGraph->Draw("same");
        expPbGraph->Draw("same");

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.15);
        gPad->SetTicky(2);
        gPad->SetLogx(1);

        litCGraph->GetXaxis()->SetRangeUser(2,500);
        litCGraph->GetYaxis()->SetRangeUser(0,9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        latex.DrawLatex(0.35,0.52,"Ni (elem.)");
        latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
        legend->AddEntry(litCGraph,"lit data (analog)","l");
        legend->AddEntry(expCGraph,"new data (DSP)","l");
        legend->Draw();
    }

    {
        // X-axis parameters
        expCGraph->GetXaxis()->SetTitle("Energy (MeV)");
        expCGraph->GetXaxis()->SetTitleSize(0.05);
        expCGraph->GetXaxis()->SetTitleFont(2);
        expCGraph->GetXaxis()->SetTitleOffset(1.4);
        expCGraph->GetXaxis()->CenterTitle();

        expCGraph->GetXaxis()->SetLabelOffset(0.01);
        expCGraph->GetXaxis()->SetLabelSize(0.05);
        expCGraph->GetXaxis()->SetLabelFont(2);

        expCGraph->GetXaxis()->SetNdivisions(10);
        expCGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCGraph->GetYaxis()->SetTitle("(#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}})");
        relCGraph->GetYaxis()->SetTitleSize(0.06);
        relCGraph->GetYaxis()->SetTitleFont(2);
        relCGraph->GetYaxis()->SetTitleOffset(1.4);
        relCGraph->GetYaxis()->CenterTitle();

        relCGraph->GetYaxis()->SetLabelOffset(0.01);
        relCGraph->GetYaxis()->SetLabelSize(0.05);

        relCGraph->GetYaxis()->SetLabelFont(2);
        relCGraph->GetYaxis()->SetNdivisions(10);
        relCGraph->GetYaxis()->SetTickLength(0.02);

        // second panel
        canvas->cd(2);
        relCGraph->Draw("");
        relNiGraph->Draw("same");
        relPbGraph->Draw("same");

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.15);
        gPad->SetTicky(2);
        gPad->SetLogx(1);

        relCGraph->GetYaxis()->SetRangeUser(-0.12,0.05);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        latex.DrawLatex(0.35,0.52,"Ni (elem.)");
        latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.25,0.25,0.35,0.45);
        legend->SetTextSize(0.04);
        legend->AddEntry(relCGraph,"{}^{nat}C","l");
        legend->AddEntry(relNiGraph,"{}^{nat}Ni","l");
        legend->AddEntry(relPbGraph,"{}^{nat}Pb","l");

        legend->Draw();
    }

    {
        // third panel
        canvas->cd(3);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.15);
        gPad->SetTicky(2);
        gPad->SetLogx(1);

        litCGraph->Draw("");
        litNiGraph->Draw("same");
        litPbGraph->Draw("same");
        corCGraph->Draw("same");
        corNiGraph->Draw("same");
        corPbGraph->Draw("same");

        // Pad dimensions and margins
        corCGraph->GetYaxis()->SetRangeUser(0,9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        latex.DrawLatex(0.35,0.52,"Ni (elem.)");
        latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
        legend->AddEntry(litCGraph,"lit data (analog)","l");
        legend->AddEntry(corCGraph,"new data (DSP)","l");
        legend->Draw();
    }

    canvas->Update();

    TImage* image = TImage::Create();
    image->FromPad(canvas);
    image->WriteImage("../plots/crossSections/FourPanelNi.png");

    // close it all up
    expFile->Close();
    relFile->Close();
    corFile->Close();
    litFile->Close();
}
