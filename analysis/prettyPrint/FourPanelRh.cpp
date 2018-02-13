{
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* canvas = new TCanvas("canvas","canvas",1200,900);
    canvas->SetCanvasSize(1100,750);
    canvas->Divide(2,2,-0.01,-0.01,0);

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
    string corFileName = "/data1/analysis/corrected.root";
    string litFileName = "/data1/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* corFile = new TFile(corFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");

    if(!expFile || !relFile || !corFile || !litFile)
    {
        cout << "Error: failed to open a required file (check filenames). Exiting..." << endl;
        exit(1);
    }
    
    string expRhNatGraphName = "RhNat";

    string relRhNatGraphName = "RhNat, expLit, percent";

    string relCorRhNatGraphName = "RhNat, expLit, corrected, percent";

    string litRhNatGraphName = "Natural Rh (n,tot)";
    
    TGraphErrors* expRhNatGraph = (TGraphErrors*)expFile->Get(expRhNatGraphName.c_str());

    TGraphErrors* relRhNatGraph = (TGraphErrors*)relFile->Get(relRhNatGraphName.c_str());

    TGraphErrors* corRhNatGraph = (TGraphErrors*)corFile->Get(expRhNatGraphName.c_str());

    TGraphErrors* relCorRhNatGraph = (TGraphErrors*)relFile->Get(relCorRhNatGraphName.c_str());

    TGraphErrors* litRhNatGraph = (TGraphErrors*)litFile->Get(litRhNatGraphName.c_str());

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

    if(!corRhNatGraph)
    {
        cout << "Error: failed to open an corrected absolute cross section graph." << endl;
        exit(1);
    }

    if(!relCorRhNatGraph)
    {
        cout << "Error: failed to open an corrected relative diff to lit cross section graph." << endl;
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

    corRhNatGraph->SetLineColor(kRed);
    corRhNatGraph->SetMarkerColor(kRed);

    relCorRhNatGraph->SetLineColor(kRed);
    relCorRhNatGraph->SetLineWidth(4);
    relCorRhNatGraph->SetLineStyle(0);
    relCorRhNatGraph->SetMarkerColor(kRed);

    litRhNatGraph->SetLineColor(kBlue);
    litRhNatGraph->SetLineWidth(4);
    litRhNatGraph->SetLineStyle(1);
    litRhNatGraph->SetMarkerColor(kBlue);

    // first panel
    {
        canvas->cd(1);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.0);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.03);
        gPad->SetTicky(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        expRhNatGraph->GetXaxis()->SetLabelOffset(0.01);
        expRhNatGraph->GetXaxis()->SetLabelSize(0.06);
        expRhNatGraph->GetXaxis()->SetLabelFont(2);

        expRhNatGraph->GetXaxis()->SetNdivisions(10);
        expRhNatGraph->GetXaxis()->SetTickLength(0.03);

        //expRhNatGraph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        expRhNatGraph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        expRhNatGraph->GetYaxis()->SetTitleSize(0.10);
        expRhNatGraph->GetYaxis()->SetTitleFont(2);
        expRhNatGraph->GetYaxis()->SetTitleOffset(0.7);
        expRhNatGraph->GetYaxis()->CenterTitle();

        expRhNatGraph->GetYaxis()->SetLabelOffset(0.01);
        expRhNatGraph->GetYaxis()->SetLabelSize(0.08);

        expRhNatGraph->GetYaxis()->SetLabelFont(2);
        expRhNatGraph->GetYaxis()->SetNdivisions(5);
        expRhNatGraph->GetYaxis()->SetTickLength(0.02);

        expRhNatGraph->GetXaxis()->SetLimits(2,600);
        expRhNatGraph->GetYaxis()->SetRangeUser(1,5.5);

        gStyle->SetLineWidth(3);

        expRhNatGraph->Draw("");
        litRhNatGraph->Draw("same");
        
        /*TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.52,0.33,"^{nat}Ni");
        latex.DrawLatex(0.51,0.13,"^{nat}C");
        */

        // Define legend format and contents
        TLegend *legend = new TLegend(0.45,0.73,0.97,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litRhNatGraph,"^{nat}Rh (An)","l");
        legend->AddEntry(expRhNatGraph,"^{nat}Rh (DSP)","l");

        legend->Draw();
    }

    // second panel
    {
        canvas->cd(2);
        gPad->SetLeftMargin(0);
        gPad->SetRightMargin(0.005);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.03);
        gPad->SetTicky(2);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        corRhNatGraph->GetXaxis()->SetLabelOffset(0.01);
        corRhNatGraph->GetXaxis()->SetLabelSize(0.06);
        corRhNatGraph->GetXaxis()->SetLabelFont(2);

        corRhNatGraph->GetXaxis()->SetNdivisions(10);
        corRhNatGraph->GetXaxis()->SetTickLength(0.03);

        //corRhNatGraph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        corRhNatGraph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        corRhNatGraph->GetYaxis()->SetTitleSize(0.10);
        corRhNatGraph->GetYaxis()->SetTitleFont(2);
        corRhNatGraph->GetYaxis()->SetTitleOffset(0.7);
        corRhNatGraph->GetYaxis()->CenterTitle();

        corRhNatGraph->GetYaxis()->SetLabelOffset(0.01);
        corRhNatGraph->GetYaxis()->SetLabelSize(0.08);

        corRhNatGraph->GetYaxis()->SetLabelFont(2);
        corRhNatGraph->GetYaxis()->SetNdivisions(5);
        corRhNatGraph->GetYaxis()->SetTickLength(0.02);

        corRhNatGraph->GetXaxis()->SetLimits(2,600);
        corRhNatGraph->GetYaxis()->SetRangeUser(1,5.5);

        gStyle->SetLineWidth(3);

        corRhNatGraph->Draw("");
        litRhNatGraph->Draw("same");
        
        /*TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.52,0.33,"^{nat}Ni");
        latex.DrawLatex(0.51,0.13,"^{nat}C");
        */

        // Define legend format and contents
        TLegend *legend = new TLegend(0.39,0.73,0.96,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litRhNatGraph,"^{nat}Rh (An)","l");
        legend->AddEntry(corRhNatGraph,"^{nat}Rh (DSP)","l");

        legend->Draw();
    }

    // third panel
    {
        canvas->cd(3);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.20);
        gPad->SetRightMargin(0.0);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.25);
        gPad->SetTicky(2);
        gPad->SetTickx(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        // X-axis parameters
        relRhNatGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relRhNatGraph->GetXaxis()->SetTitleSize(0.08);
        relRhNatGraph->GetXaxis()->SetTitleFont(2);
        relRhNatGraph->GetXaxis()->SetTitleOffset(1.4);
        relRhNatGraph->GetXaxis()->CenterTitle();

        relRhNatGraph->GetXaxis()->SetLabelOffset(0.01);
        relRhNatGraph->GetXaxis()->SetLabelSize(0.08);
        relRhNatGraph->GetXaxis()->SetLabelFont(2);

        relRhNatGraph->GetXaxis()->SetNdivisions(10);
        relRhNatGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relRhNatGraph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relRhNatGraph->GetYaxis()->SetTitleSize(0.08);
        relRhNatGraph->GetYaxis()->SetTitleFont(2);
        relRhNatGraph->GetYaxis()->SetTitleOffset(1.05);
        relRhNatGraph->GetYaxis()->CenterTitle();

        relRhNatGraph->GetYaxis()->SetLabelOffset(0.01);
        relRhNatGraph->GetYaxis()->SetLabelSize(0.08);

        relRhNatGraph->GetYaxis()->SetLabelFont(2);
        relRhNatGraph->GetYaxis()->SetNdivisions(10);
        relRhNatGraph->GetYaxis()->SetTickLength(0.02);

        relRhNatGraph->Draw("");
        
        relRhNatGraph->GetXaxis()->SetLimits(2,600);
        relRhNatGraph->GetYaxis()->SetRangeUser(-20,20);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        //latex.DrawLatex(0.35,0.52,"Ni (elem.)");
        //latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.77,0.67,0.94,0.97);
        legend->SetTextSize(0.065);
        legend->SetTextAlign(12);
        legend->AddEntry(relRhNatGraph," ^{nat}Rh","l");

        legend->Draw();

        TLine* zeroLine = new TLine(0, 0, 600, 0);
        zeroLine->SetLineColor(kBlack);
        zeroLine->SetLineWidth(3);
        zeroLine->SetLineStyle(9);
        zeroLine->Draw();
    }

    // fourth panel
    {
        canvas->cd(4);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0.0);
        gPad->SetRightMargin(0.005);
        gPad->SetTopMargin(0.01);
        gPad->SetBottomMargin(0.25);
        gPad->SetTicky(2);
        gPad->SetTickx(1);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        // X-axis parameters
        relCorRhNatGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relCorRhNatGraph->GetXaxis()->SetTitleSize(0.08);
        relCorRhNatGraph->GetXaxis()->SetTitleFont(2);
        relCorRhNatGraph->GetXaxis()->SetTitleOffset(1.4);
        relCorRhNatGraph->GetXaxis()->CenterTitle();

        relCorRhNatGraph->GetXaxis()->SetLabelOffset(0.01);
        relCorRhNatGraph->GetXaxis()->SetLabelSize(0.08);
        relCorRhNatGraph->GetXaxis()->SetLabelFont(2);

        relCorRhNatGraph->GetXaxis()->SetNdivisions(10);
        relCorRhNatGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCorRhNatGraph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relCorRhNatGraph->GetYaxis()->SetTitleSize(0.08);
        relCorRhNatGraph->GetYaxis()->SetTitleFont(2);
        relCorRhNatGraph->GetYaxis()->SetTitleOffset(1.05);
        relCorRhNatGraph->GetYaxis()->CenterTitle();

        relCorRhNatGraph->GetYaxis()->SetLabelOffset(0.01);
        relCorRhNatGraph->GetYaxis()->SetLabelSize(0.08);

        relCorRhNatGraph->GetYaxis()->SetLabelFont(2);
        relCorRhNatGraph->GetYaxis()->SetNdivisions(10);
        relCorRhNatGraph->GetYaxis()->SetTickLength(0.02);

        relCorRhNatGraph->Draw("");
        
        relCorRhNatGraph->GetXaxis()->SetLimits(2,600);
        relCorRhNatGraph->GetYaxis()->SetRangeUser(-20,20);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        //latex.DrawLatex(0.35,0.52,"Ni (elem.)");
        //latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.72,0.66,0.92,0.97);
        legend->SetTextSize(0.065);
        legend->SetTextAlign(12);
        legend->AddEntry(relCorRhNatGraph," ^{nat}Rh","l");

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
    corFile->Close();
    litFile->Close();
}
