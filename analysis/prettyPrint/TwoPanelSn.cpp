{
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas",750,850);
    //canvas->SetCanvasSize(750,750);
    canvas->Divide(1,2);

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
    
    string expCGraphName = "CNat";
    string expSnGraphName = "SnNat";
    string expPbGraphName = "PbNat";

    string relCGraphName = "CNat, expLit, percent";
    string relSnGraphName = "SnNat, expLit, percent";
    string relPbGraphName = "PbNat, expLit, percent";

    string litCGraphName = "Natural C (n,tot)";
    string litSnGraphName = "Natural Sn (n,tot)";
    string litPbGraphName = "Natural Pb (n,tot)";
    
    TGraphAsymmErrors* expCGraph = (TGraphAsymmErrors*)expFile->Get(expCGraphName.c_str());
    TGraphAsymmErrors* expSnGraph = (TGraphAsymmErrors*)expFile->Get(expSnGraphName.c_str());
    TGraphAsymmErrors* expPbGraph = (TGraphAsymmErrors*)expFile->Get(expPbGraphName.c_str());

    TGraphAsymmErrors* relCGraph = (TGraphAsymmErrors*)relFile->Get(relCGraphName.c_str());
    TGraphAsymmErrors* relSnGraph = (TGraphAsymmErrors*)relFile->Get(relSnGraphName.c_str());
    TGraphAsymmErrors* relPbGraph = (TGraphAsymmErrors*)relFile->Get(relPbGraphName.c_str());

    TGraphAsymmErrors* litCGraph = (TGraphAsymmErrors*)litFile->Get(litCGraphName.c_str());
    TGraphAsymmErrors* litSnGraph = (TGraphAsymmErrors*)litFile->Get(litSnGraphName.c_str());
    TGraphAsymmErrors* litPbGraph = (TGraphAsymmErrors*)litFile->Get(litPbGraphName.c_str());

    if(!expCGraph || !expSnGraph || !expPbGraph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    if(!relCGraph || !relSnGraph || !relPbGraph)
    {
        cout << "Error: failed to open an relative diff to lit cross section graph." << endl;
        exit(1);
    }

    if(!litCGraph || !litSnGraph || !litPbGraph)
    {
        cout << "Error: failed to open an lit cross section graph." << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expCGraph->SetLineColor(kRed);
    expCGraph->SetMarkerColor(kRed);
    expCGraph->SetLineWidth(4);
    expCGraph->SetMarkerSize(0.9);
    expCGraph->SetMarkerStyle(8);

    expSnGraph->SetLineColor(kRed);
    expSnGraph->SetMarkerColor(kRed);
    expSnGraph->SetLineWidth(4);

    expPbGraph->SetLineColor(kRed);
    expPbGraph->SetMarkerColor(kRed);
    expPbGraph->SetLineWidth(4);

    relCGraph->SetLineColor(kRed-9);
    relCGraph->SetLineWidth(4);
    relCGraph->SetLineStyle(0);
    relCGraph->SetMarkerColor(kRed-9);

    relSnGraph->SetLineColor(kRed);
    relSnGraph->SetLineWidth(4);
    relSnGraph->SetLineStyle(0);
    relSnGraph->SetMarkerColor(kRed);

    relPbGraph->SetLineColor(kRed+3);
    relPbGraph->SetLineWidth(4);
    relPbGraph->SetLineStyle(0);
    relPbGraph->SetMarkerColor(kRed+3);

    litCGraph->SetLineColor(kBlue);
    litCGraph->SetLineWidth(4);
    litCGraph->SetLineStyle(0);
    litCGraph->SetMarkerColor(kBlue);

    litSnGraph->SetLineColor(kBlue);
    litSnGraph->SetLineWidth(4);
    litSnGraph->SetLineStyle(0);
    litSnGraph->SetMarkerColor(kBlue);

    litPbGraph->SetLineColor(kBlue);
    litPbGraph->SetLineWidth(4);
    litPbGraph->SetLineStyle(0);
    litPbGraph->SetMarkerColor(kBlue);

    // first panel
    {
        canvas->cd(1);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.005);
        gPad->SetTicky(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        expCGraph->GetXaxis()->SetLabelOffset(0.01);
        expCGraph->GetXaxis()->SetLabelSize(0.0);
        expCGraph->GetXaxis()->SetLabelFont(2);

        expCGraph->GetXaxis()->SetNdivisions(10);
        expCGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        expCGraph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        expCGraph->GetYaxis()->SetTitleSize(0.10);
        expCGraph->GetYaxis()->SetTitleFont(2);
        expCGraph->GetYaxis()->SetTitleOffset(0.7);
        expCGraph->GetYaxis()->CenterTitle();

        expCGraph->GetYaxis()->SetLabelOffset(0.01);
        expCGraph->GetYaxis()->SetLabelSize(0.10);

        expCGraph->GetYaxis()->SetLabelFont(2);
        expCGraph->GetYaxis()->SetNdivisions(5);
        expCGraph->GetYaxis()->SetTickLength(0.02);

        expCGraph->GetXaxis()->SetRangeUser(2,600);
        expCGraph->GetYaxis()->SetRangeUser(0.01,9);

        expCGraph->Draw("AL");
        litCGraph->Draw("same");
        litSnGraph->Draw("same");
        litPbGraph->Draw("same");
        expCGraph->Draw("same");
        expSnGraph->Draw("same");
        expPbGraph->Draw("same");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.59,0.63,"^{nat}Pb");
        latex.DrawLatex(0.50,0.37,"^{nat}Sn");
        latex.DrawLatex(0.47,0.11,"^{nat}C");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.77,0.70,0.95,0.9);
        legend->SetTextSize(0.08);
        legend->SetTextAlign(12);
        legend->AddEntry(litCGraph,"Analog","l");
        legend->AddEntry(expCGraph,"DSP","l");
        legend->Draw();
    }

    // second panel
    {
        canvas->cd(2);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.01);
        gPad->SetTopMargin(0.);
        gPad->SetBottomMargin(0.25);
        gPad->SetTicky(1);
        gPad->SetTickx(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        // X-axis parameters
        relCGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relCGraph->GetXaxis()->SetTitleSize(0.08);
        relCGraph->GetXaxis()->SetTitleFont(2);
        relCGraph->GetXaxis()->SetTitleOffset(1.4);
        relCGraph->GetXaxis()->CenterTitle();

        relCGraph->GetXaxis()->SetLabelOffset(0.01);
        relCGraph->GetXaxis()->SetLabelSize(0.08);
        relCGraph->GetXaxis()->SetLabelFont(2);

        relCGraph->GetXaxis()->SetNdivisions(10);
        relCGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCGraph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relCGraph->GetYaxis()->SetTitleSize(0.08);
        relCGraph->GetYaxis()->SetTitleFont(2);
        relCGraph->GetYaxis()->SetTitleOffset(1.05);
        relCGraph->GetYaxis()->CenterTitle();

        relCGraph->GetYaxis()->SetLabelOffset(0.01);
        relCGraph->GetYaxis()->SetLabelSize(0.08);

        relCGraph->GetYaxis()->SetLabelFont(2);
        relCGraph->GetYaxis()->SetNdivisions(10);
        relCGraph->GetYaxis()->SetTickLength(0.02);

        relCGraph->Draw("");
        relSnGraph->Draw("same");
        relPbGraph->Draw("same");
        
        relCGraph->GetXaxis()->SetRangeUser(2,600);
        relCGraph->GetYaxis()->SetRangeUser(-8,7.9);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        //latex.DrawLatex(0.35,0.52,"Sn (elem.)");
        //latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.8,0.67,0.96,0.96);
        //legend->SetNColumns(3);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relCGraph,"{}^{nat}C","l");
        legend->AddEntry(relSnGraph,"{}^{nat}Sn","l");
        legend->AddEntry(relPbGraph,"{}^{nat}Pb","l");

        legend->Draw();

        TLine* zeroLine = new TLine(0, 0, 600, 0);
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
