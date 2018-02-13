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
    string expFileName = "/data2/analysis/total.root";
    string relFileName = "/data2/analysis/relative.root";
    string corFileName = "/data2/analysis/corrected.root";
    string litFileName = "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* corFile = new TFile(corFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");

    if(!expFile || !relFile || !corFile || !litFile)
    {
        cout << "Error: failed to open a required file (check filenames). Exiting..." << endl;
        exit(1);
    }
    
    string expSn112GraphName = "Sn112";
    string expSn124GraphName = "Sn124";

    string relSn112GraphName = "Sn112, expLit, percent";
    string relSn124GraphName = "Sn124, expLit, percent";

    string relCorSn112GraphName = "Sn112, expLit, corrected, percent";
    string relCorSn124GraphName = "Sn124, expLit, corrected, percent";

    string litSn112GraphName = "Sn112 (n,tot)";
    string litSn124GraphName = "Sn124 (n,tot)";
    
    TGraphErrors* expSn112Graph = (TGraphErrors*)expFile->Get(expSn112GraphName.c_str());
    TGraphErrors* expSn124Graph = (TGraphErrors*)expFile->Get(expSn124GraphName.c_str());

    TGraphErrors* relSn112Graph = (TGraphErrors*)relFile->Get(relSn112GraphName.c_str());
    TGraphErrors* relSn124Graph = (TGraphErrors*)relFile->Get(relSn124GraphName.c_str());

    TGraphErrors* corSn112Graph = (TGraphErrors*)corFile->Get(expSn112GraphName.c_str());
    TGraphErrors* corSn124Graph = (TGraphErrors*)corFile->Get(expSn124GraphName.c_str());

    TGraphErrors* relCorSn112Graph = (TGraphErrors*)relFile->Get(relCorSn112GraphName.c_str());
    TGraphErrors* relCorSn124Graph = (TGraphErrors*)relFile->Get(relCorSn124GraphName.c_str());

    TGraphErrors* litSn112Graph = (TGraphErrors*)litFile->Get(litSn112GraphName.c_str());
    TGraphErrors* litSn124Graph = (TGraphErrors*)litFile->Get(litSn124GraphName.c_str());

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

    if(!corSn112Graph || !corSn124Graph)
    {
        cout << "Error: failed to open an corrected absolute cross section graph." << endl;
        exit(1);
    }

    if(!relCorSn112Graph || !relCorSn124Graph)
    {
        cout << "Error: failed to open an corrected relative diff to lit cross section graph." << endl;
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

    expSn124Graph->SetLineColor(kRed+3);
    expSn124Graph->SetLineWidth(4);
    expSn124Graph->SetLineStyle(1);
    expSn124Graph->SetMarkerColor(kRed+3);

    relSn112Graph->SetLineColor(kRed);
    relSn112Graph->SetLineWidth(4);
    relSn112Graph->SetLineStyle(0);
    relSn112Graph->SetMarkerColor(kRed);
    relSn112Graph->SetMarkerSize(3);
    relSn112Graph->SetMarkerStyle(34);

    relSn124Graph->SetLineColor(kRed+3);
    relSn124Graph->SetLineWidth(4);
    relSn124Graph->SetLineStyle(0);
    relSn124Graph->SetMarkerColor(kRed+3);
    relSn124Graph->SetMarkerSize(1.5);
    relSn124Graph->SetMarkerStyle(43);

    corSn112Graph->SetLineColor(kRed);
    corSn112Graph->SetMarkerColor(kRed);

    corSn124Graph->SetLineColor(kRed+3);
    corSn124Graph->SetMarkerColor(kRed+3);

    relCorSn112Graph->SetLineColor(kRed);
    relCorSn112Graph->SetLineWidth(4);
    relCorSn112Graph->SetLineStyle(0);
    relCorSn112Graph->SetMarkerColor(kRed);
    relCorSn112Graph->SetMarkerSize(3);
    relCorSn112Graph->SetMarkerStyle(34);

    relCorSn124Graph->SetLineColor(kRed+3);
    relCorSn124Graph->SetLineWidth(4);
    relCorSn124Graph->SetLineStyle(0);
    relCorSn124Graph->SetMarkerColor(kRed+3);
    relCorSn124Graph->SetMarkerSize(1.5);
    relCorSn124Graph->SetMarkerStyle(43);

    litSn124Graph->SetLineColor(kBlue);
    litSn124Graph->SetLineWidth(4);
    litSn124Graph->SetLineStyle(1);
    litSn124Graph->SetMarkerColor(kBlue);
    litSn124Graph->SetMarkerSize(1.5);
    litSn124Graph->SetMarkerStyle(43);

    litSn112Graph->SetLineColor(kViolet);
    litSn112Graph->SetLineWidth(4);
    litSn112Graph->SetLineStyle(1);
    litSn112Graph->SetMarkerColor(kViolet);
    litSn112Graph->SetMarkerSize(3);
    litSn112Graph->SetMarkerStyle(34);

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

        expSn112Graph->GetXaxis()->SetLabelOffset(0.01);
        expSn112Graph->GetXaxis()->SetLabelSize(0.06);
        expSn112Graph->GetXaxis()->SetLabelFont(2);

        expSn112Graph->GetXaxis()->SetNdivisions(10);
        expSn112Graph->GetXaxis()->SetTickLength(0.03);

        //expSn112Graph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        expSn112Graph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        expSn112Graph->GetYaxis()->SetTitleSize(0.10);
        expSn112Graph->GetYaxis()->SetTitleFont(2);
        expSn112Graph->GetYaxis()->SetTitleOffset(0.7);
        expSn112Graph->GetYaxis()->CenterTitle();

        expSn112Graph->GetYaxis()->SetLabelOffset(0.01);
        expSn112Graph->GetYaxis()->SetLabelSize(0.08);

        expSn112Graph->GetYaxis()->SetLabelFont(2);
        expSn112Graph->GetYaxis()->SetNdivisions(5);
        expSn112Graph->GetYaxis()->SetTickLength(0.02);

        expSn112Graph->GetXaxis()->SetLimits(3,600);
        expSn112Graph->GetYaxis()->SetRangeUser(1,6);

        gStyle->SetLineWidth(3);

        expSn112Graph->Draw("");
        expSn124Graph->Draw("same");
        litSn112Graph->Draw("p");
        litSn124Graph->Draw("p");

        
        /*TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.52,0.33,"^{nat}Sn");
        latex.DrawLatex(0.51,0.13,"^{nat}C");
        */

        // Define legend format and contents
        TLegend *legend = new TLegend(0.48,0.73,0.99,0.97);
        legend->SetTextSize(0.06);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litSn112Graph,"^{112}Sn (An)","p");
        legend->AddEntry(litSn124Graph,"^{124}Sn (An)","p");
        legend->AddEntry(expSn112Graph,"^{112}Sn (DSP)","l");
        legend->AddEntry(expSn124Graph,"^{124}Sn (DSP)","l");

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

        corSn112Graph->GetXaxis()->SetLabelOffset(0.01);
        corSn112Graph->GetXaxis()->SetLabelSize(0.06);
        corSn112Graph->GetXaxis()->SetLabelFont(2);

        corSn112Graph->GetXaxis()->SetNdivisions(10);
        corSn112Graph->GetXaxis()->SetTickLength(0.03);

        //corSn112Graph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        corSn112Graph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        corSn112Graph->GetYaxis()->SetTitleSize(0.10);
        corSn112Graph->GetYaxis()->SetTitleFont(2);
        corSn112Graph->GetYaxis()->SetTitleOffset(0.7);
        corSn112Graph->GetYaxis()->CenterTitle();

        corSn112Graph->GetYaxis()->SetLabelOffset(0.01);
        corSn112Graph->GetYaxis()->SetLabelSize(0.08);

        corSn112Graph->GetYaxis()->SetLabelFont(2);
        corSn112Graph->GetYaxis()->SetNdivisions(5);
        corSn112Graph->GetYaxis()->SetTickLength(0.02);

        corSn112Graph->GetXaxis()->SetLimits(3,600);
        corSn112Graph->GetYaxis()->SetRangeUser(1,6);

        gStyle->SetLineWidth(3);

        corSn112Graph->Draw("");
        corSn124Graph->Draw("same");
        litSn112Graph->Draw("p");
        litSn124Graph->Draw("p");

        
        /*TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.52,0.33,"^{nat}Sn");
        latex.DrawLatex(0.51,0.13,"^{nat}C");
        */

        // Define legend format and contents
        TLegend *legend = new TLegend(0.40,0.73,0.98,0.97);
        legend->SetTextSize(0.06);
        legend->SetTextAlign(12);
        legend->SetNColumns(2);
        legend->AddEntry(litSn112Graph,"^{112}Sn (An)","p");
        legend->AddEntry(litSn124Graph,"^{124}Sn (An)","p");
        legend->AddEntry(corSn112Graph,"^{112}Sn (DSP)","l");
        legend->AddEntry(corSn124Graph,"^{124}Sn (DSP)","l");

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
        relSn112Graph->GetXaxis()->SetTitle("Energy (MeV)");
        relSn112Graph->GetXaxis()->SetTitleSize(0.08);
        relSn112Graph->GetXaxis()->SetTitleFont(2);
        relSn112Graph->GetXaxis()->SetTitleOffset(1.4);
        relSn112Graph->GetXaxis()->CenterTitle();

        relSn112Graph->GetXaxis()->SetLabelOffset(0.01);
        relSn112Graph->GetXaxis()->SetLabelSize(0.08);
        relSn112Graph->GetXaxis()->SetLabelFont(2);

        relSn112Graph->GetXaxis()->SetNdivisions(10);
        relSn112Graph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relSn112Graph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relSn112Graph->GetYaxis()->SetTitleSize(0.08);
        relSn112Graph->GetYaxis()->SetTitleFont(2);
        relSn112Graph->GetYaxis()->SetTitleOffset(1.05);
        relSn112Graph->GetYaxis()->CenterTitle();

        relSn112Graph->GetYaxis()->SetLabelOffset(0.01);
        relSn112Graph->GetYaxis()->SetLabelSize(0.08);

        relSn112Graph->GetYaxis()->SetLabelFont(2);
        relSn112Graph->GetYaxis()->SetNdivisions(10);
        relSn112Graph->GetYaxis()->SetTickLength(0.02);

        relSn112Graph->Draw("");
        relSn124Graph->Draw("p");
        
        relSn112Graph->GetXaxis()->SetLimits(2,600);
        relSn112Graph->GetYaxis()->SetRangeUser(-5,5);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        //latex.DrawLatex(0.35,0.52,"Sn (elem.)");
        //latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.75,0.67,0.94,0.97);
        legend->SetTextSize(0.065);
        legend->SetTextAlign(12);
        legend->AddEntry(relSn112Graph," ^{112}Sn","p");
        legend->AddEntry(relSn124Graph," ^{124}Sn","p");

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
        relCorSn112Graph->GetXaxis()->SetTitle("Energy (MeV)");
        relCorSn112Graph->GetXaxis()->SetTitleSize(0.08);
        relCorSn112Graph->GetXaxis()->SetTitleFont(2);
        relCorSn112Graph->GetXaxis()->SetTitleOffset(1.4);
        relCorSn112Graph->GetXaxis()->CenterTitle();

        relCorSn112Graph->GetXaxis()->SetLabelOffset(0.01);
        relCorSn112Graph->GetXaxis()->SetLabelSize(0.08);
        relCorSn112Graph->GetXaxis()->SetLabelFont(2);

        relCorSn112Graph->GetXaxis()->SetNdivisions(10);
        relCorSn112Graph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCorSn112Graph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relCorSn112Graph->GetYaxis()->SetTitleSize(0.08);
        relCorSn112Graph->GetYaxis()->SetTitleFont(2);
        relCorSn112Graph->GetYaxis()->SetTitleOffset(1.05);
        relCorSn112Graph->GetYaxis()->CenterTitle();

        relCorSn112Graph->GetYaxis()->SetLabelOffset(0.01);
        relCorSn112Graph->GetYaxis()->SetLabelSize(0.08);

        relCorSn112Graph->GetYaxis()->SetLabelFont(2);
        relCorSn112Graph->GetYaxis()->SetNdivisions(10);
        relCorSn112Graph->GetYaxis()->SetTickLength(0.02);

        relCorSn112Graph->Draw("");
        relCorSn124Graph->Draw("p");
        
        relCorSn112Graph->GetXaxis()->SetLimits(2,600);
        relCorSn112Graph->GetYaxis()->SetRangeUser(-5,5);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        //latex.DrawLatex(0.35,0.52,"Sn (elem.)");
        //latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.75,0.67,0.94,0.97);
        legend->SetTextSize(0.065);
        legend->SetTextAlign(12);
        legend->AddEntry(relCorSn112Graph," ^{112}Sn","p");
        legend->AddEntry(relCorSn124Graph," ^{124}Sn","p");

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
