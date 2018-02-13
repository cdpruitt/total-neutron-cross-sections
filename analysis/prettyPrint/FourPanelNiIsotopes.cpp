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
    
    string expNi58GraphName = "Ni58";
    string expNi64GraphName = "Ni64";

    string relNi58GraphName = "Ni58, expLit, percent";
    string relNi64GraphName = "Ni64, expLit, percent";

    string relCorNi58GraphName = "Ni58, expLit, corrected, percent";
    string relCorNi64GraphName = "Ni64, expLit, corrected, percent";

    string litNi58GraphName = "Ni58 (n,tot)";
    string litNi64GraphName = "Ni64 (n,tot)";
    
    TGraphErrors* expNi58Graph = (TGraphErrors*)expFile->Get(expNi58GraphName.c_str());
    TGraphErrors* expNi64Graph = (TGraphErrors*)expFile->Get(expNi64GraphName.c_str());

    TGraphErrors* relNi58Graph = (TGraphErrors*)relFile->Get(relNi58GraphName.c_str());
    TGraphErrors* relNi64Graph = (TGraphErrors*)relFile->Get(relNi64GraphName.c_str());

    TGraphErrors* corNi58Graph = (TGraphErrors*)corFile->Get(expNi58GraphName.c_str());
    TGraphErrors* corNi64Graph = (TGraphErrors*)corFile->Get(expNi64GraphName.c_str());

    TGraphErrors* relCorNi58Graph = (TGraphErrors*)relFile->Get(relCorNi58GraphName.c_str());
    TGraphErrors* relCorNi64Graph = (TGraphErrors*)relFile->Get(relCorNi64GraphName.c_str());

    TGraphErrors* litNi58Graph = (TGraphErrors*)litFile->Get(litNi58GraphName.c_str());
    TGraphErrors* litNi64Graph = (TGraphErrors*)litFile->Get(litNi64GraphName.c_str());

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

    if(!corNi58Graph || !corNi64Graph)
    {
        cout << "Error: failed to open an corrected absolute cross section graph." << endl;
        exit(1);
    }

    if(!relCorNi58Graph || !relCorNi64Graph)
    {
        cout << "Error: failed to open an corrected relative diff to lit cross section graph." << endl;
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

    expNi64Graph->SetLineColor(kRed+3);
    expNi64Graph->SetLineWidth(4);
    expNi64Graph->SetLineStyle(1);
    expNi64Graph->SetMarkerColor(kRed+3);

    relNi58Graph->SetLineColor(kRed);
    relNi58Graph->SetLineWidth(4);
    relNi58Graph->SetLineStyle(0);
    relNi58Graph->SetMarkerColor(kRed);

    relNi64Graph->SetLineColor(kRed+3);
    relNi64Graph->SetLineWidth(4);
    relNi64Graph->SetLineStyle(0);
    relNi64Graph->SetMarkerColor(kRed+3);
    relNi64Graph->SetMarkerSize(3);
    relNi64Graph->SetMarkerStyle(34);

    corNi58Graph->SetLineColor(kRed);
    corNi58Graph->SetMarkerColor(kRed);

    corNi64Graph->SetLineColor(kRed+3);
    corNi64Graph->SetMarkerColor(kRed+3);

    relCorNi58Graph->SetLineColor(kRed);
    relCorNi58Graph->SetLineWidth(4);
    relCorNi58Graph->SetLineStyle(0);
    relCorNi58Graph->SetMarkerColor(kRed);

    relCorNi64Graph->SetLineColor(kRed+3);
    relCorNi64Graph->SetLineWidth(4);
    relCorNi64Graph->SetLineStyle(0);
    relCorNi64Graph->SetMarkerColor(kRed+3);
    relCorNi64Graph->SetMarkerSize(3);
    relCorNi64Graph->SetMarkerStyle(34);

    litNi58Graph->SetLineColor(kBlue);
    litNi58Graph->SetLineWidth(4);
    litNi58Graph->SetLineStyle(1);
    litNi58Graph->SetMarkerColor(kBlue);

    litNi64Graph->SetLineColor(kViolet);
    litNi64Graph->SetLineWidth(4);
    litNi64Graph->SetLineStyle(1);
    litNi64Graph->SetMarkerColor(kViolet);
    litNi64Graph->SetMarkerSize(3);
    litNi64Graph->SetMarkerStyle(34);

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

        expNi58Graph->GetXaxis()->SetLabelOffset(0.01);
        expNi58Graph->GetXaxis()->SetLabelSize(0.06);
        expNi58Graph->GetXaxis()->SetLabelFont(2);

        expNi58Graph->GetXaxis()->SetNdivisions(10);
        expNi58Graph->GetXaxis()->SetTickLength(0.03);

        //expNi58Graph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        expNi58Graph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        expNi58Graph->GetYaxis()->SetTitleSize(0.10);
        expNi58Graph->GetYaxis()->SetTitleFont(2);
        expNi58Graph->GetYaxis()->SetTitleOffset(0.7);
        expNi58Graph->GetYaxis()->CenterTitle();

        expNi58Graph->GetYaxis()->SetLabelOffset(0.01);
        expNi58Graph->GetYaxis()->SetLabelSize(0.08);

        expNi58Graph->GetYaxis()->SetLabelFont(2);
        expNi58Graph->GetYaxis()->SetNdivisions(5);
        expNi58Graph->GetYaxis()->SetTickLength(0.02);

        expNi58Graph->GetXaxis()->SetLimits(2,600);
        expNi58Graph->GetYaxis()->SetRangeUser(0.5,4.5);

        gStyle->SetLineWidth(3);

        expNi58Graph->Draw("");
        litNi58Graph->Draw("same");
        litNi64Graph->Draw("p");
        expNi58Graph->Draw("same");
        expNi64Graph->Draw("same");
        
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
        legend->AddEntry(litNi58Graph,"^{58}Ni (An)","l");
        legend->AddEntry(litNi64Graph,"^{64}Ni (An)","p");
        legend->AddEntry(expNi58Graph,"^{58}Ni (DSP)","l");
        legend->AddEntry(expNi64Graph,"^{64}Ni (DSP)","l");

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

        corNi58Graph->GetXaxis()->SetLabelOffset(0.01);
        corNi58Graph->GetXaxis()->SetLabelSize(0.06);
        corNi58Graph->GetXaxis()->SetLabelFont(2);

        corNi58Graph->GetXaxis()->SetNdivisions(10);
        corNi58Graph->GetXaxis()->SetTickLength(0.03);

        //corNi58Graph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        corNi58Graph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        corNi58Graph->GetYaxis()->SetTitleSize(0.10);
        corNi58Graph->GetYaxis()->SetTitleFont(2);
        corNi58Graph->GetYaxis()->SetTitleOffset(0.7);
        corNi58Graph->GetYaxis()->CenterTitle();

        corNi58Graph->GetYaxis()->SetLabelOffset(0.01);
        corNi58Graph->GetYaxis()->SetLabelSize(0.08);

        corNi58Graph->GetYaxis()->SetLabelFont(2);
        corNi58Graph->GetYaxis()->SetNdivisions(5);
        corNi58Graph->GetYaxis()->SetTickLength(0.02);

        corNi58Graph->GetXaxis()->SetLimits(2,600);
        corNi58Graph->GetYaxis()->SetRangeUser(0.5,4.5);

        gStyle->SetLineWidth(3);

        corNi58Graph->Draw("");
        litNi58Graph->Draw("same");
        litNi64Graph->Draw("p");
        corNi58Graph->Draw("same");
        corNi64Graph->Draw("same");
        
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
        legend->AddEntry(litNi58Graph,"^{58}Ni (An)","l");
        legend->AddEntry(litNi64Graph,"^{64}Ni (An)","p");
        legend->AddEntry(corNi58Graph,"^{58}Ni (DSP)","l");
        legend->AddEntry(corNi64Graph,"^{64}Ni (DSP)","l");

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
        relNi58Graph->GetXaxis()->SetTitle("Energy (MeV)");
        relNi58Graph->GetXaxis()->SetTitleSize(0.08);
        relNi58Graph->GetXaxis()->SetTitleFont(2);
        relNi58Graph->GetXaxis()->SetTitleOffset(1.4);
        relNi58Graph->GetXaxis()->CenterTitle();

        relNi58Graph->GetXaxis()->SetLabelOffset(0.01);
        relNi58Graph->GetXaxis()->SetLabelSize(0.08);
        relNi58Graph->GetXaxis()->SetLabelFont(2);

        relNi58Graph->GetXaxis()->SetNdivisions(10);
        relNi58Graph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relNi58Graph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relNi58Graph->GetYaxis()->SetTitleSize(0.08);
        relNi58Graph->GetYaxis()->SetTitleFont(2);
        relNi58Graph->GetYaxis()->SetTitleOffset(1.05);
        relNi58Graph->GetYaxis()->CenterTitle();

        relNi58Graph->GetYaxis()->SetLabelOffset(0.01);
        relNi58Graph->GetYaxis()->SetLabelSize(0.08);

        relNi58Graph->GetYaxis()->SetLabelFont(2);
        relNi58Graph->GetYaxis()->SetNdivisions(10);
        relNi58Graph->GetYaxis()->SetTickLength(0.02);

        relNi58Graph->Draw("");
        relNi64Graph->Draw("p");
        
        relNi58Graph->GetXaxis()->SetLimits(2,600);
        relNi58Graph->GetYaxis()->SetRangeUser(-20,20);

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
        legend->AddEntry(relNi58Graph," ^{58}Ni","l");
        legend->AddEntry(relNi64Graph," ^{64}Ni","p");

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
        relCorNi58Graph->GetXaxis()->SetTitle("Energy (MeV)");
        relCorNi58Graph->GetXaxis()->SetTitleSize(0.08);
        relCorNi58Graph->GetXaxis()->SetTitleFont(2);
        relCorNi58Graph->GetXaxis()->SetTitleOffset(1.4);
        relCorNi58Graph->GetXaxis()->CenterTitle();

        relCorNi58Graph->GetXaxis()->SetLabelOffset(0.01);
        relCorNi58Graph->GetXaxis()->SetLabelSize(0.08);
        relCorNi58Graph->GetXaxis()->SetLabelFont(2);

        relCorNi58Graph->GetXaxis()->SetNdivisions(10);
        relCorNi58Graph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCorNi58Graph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relCorNi58Graph->GetYaxis()->SetTitleSize(0.08);
        relCorNi58Graph->GetYaxis()->SetTitleFont(2);
        relCorNi58Graph->GetYaxis()->SetTitleOffset(1.05);
        relCorNi58Graph->GetYaxis()->CenterTitle();

        relCorNi58Graph->GetYaxis()->SetLabelOffset(0.01);
        relCorNi58Graph->GetYaxis()->SetLabelSize(0.08);

        relCorNi58Graph->GetYaxis()->SetLabelFont(2);
        relCorNi58Graph->GetYaxis()->SetNdivisions(10);
        relCorNi58Graph->GetYaxis()->SetTickLength(0.02);

        relCorNi58Graph->Draw("");
        relCorNi64Graph->Draw("p");
        
        relCorNi58Graph->GetXaxis()->SetLimits(2,600);
        relCorNi58Graph->GetYaxis()->SetRangeUser(-20,20);

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
        legend->AddEntry(relCorNi58Graph," ^{58}Ni","l");
        legend->AddEntry(relCorNi64Graph," ^{64}Ni","p");

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
