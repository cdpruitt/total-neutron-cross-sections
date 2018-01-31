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
    
    string expCGraphName = "CNat";
    string expNiGraphName = "NiNat";
    string expPbGraphName = "PbNat";

    string relCGraphName = "NatC, expLit, percent";
    string relNiGraphName = "NatNi, expLit, percent";
    string relPbGraphName = "NatPb, expLit, percent";

    string relCorCGraphName = "NatC, expLit, corrected, percent";
    string relCorNiGraphName = "NatNi, expLit, corrected, percent";
    string relCorPbGraphName = "NatPb, expLit, corrected, percent";

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

    TGraphErrors* relCorCGraph = (TGraphErrors*)relFile->Get(relCorCGraphName.c_str());
    TGraphErrors* relCorNiGraph = (TGraphErrors*)relFile->Get(relCorNiGraphName.c_str());
    TGraphErrors* relCorPbGraph = (TGraphErrors*)relFile->Get(relCorPbGraphName.c_str());

    TGraphErrors* litCGraph = (TGraphErrors*)litFile->Get(litCGraphName.c_str());
    TGraphErrors* litNiGraph = (TGraphErrors*)litFile->Get(litNiGraphName.c_str());
    TGraphErrors* litPbGraph = (TGraphErrors*)litFile->Get(litPbGraphName.c_str());

    if(!expCGraph || !expNiGraph || !expPbGraph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    if(!relCGraph || !relNiGraph || !relPbGraph)
    {
        cout << "Error: failed to open an relative diff to lit cross section graph." << endl;
        exit(1);
    }

    if(!corCGraph || !corNiGraph || !corPbGraph)
    {
        cout << "Error: failed to open an corrected absolute cross section graph." << endl;
        exit(1);
    }

    if(!relCorCGraph || !relCorNiGraph || !relCorPbGraph)
    {
        cout << "Error: failed to open an corrected relative diff to lit cross section graph." << endl;
        exit(1);
    }

    if(!litCGraph || !litNiGraph || !litPbGraph)
    {
        cout << "Error: failed to open an lit cross section graph." << endl;
        exit(1);
    }

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

    relCorCGraph->SetLineColor(kRed-9);
    relCorCGraph->SetLineWidth(4);
    relCorCGraph->SetLineStyle(0);
    relCorCGraph->SetMarkerColor(kRed-9);

    relCorNiGraph->SetLineColor(kRed);
    relCorNiGraph->SetLineWidth(4);
    relCorNiGraph->SetLineStyle(0);
    relCorNiGraph->SetMarkerColor(kRed);

    relCorPbGraph->SetLineColor(kRed+3);
    relCorPbGraph->SetLineWidth(4);
    relCorPbGraph->SetLineStyle(0);
    relCorPbGraph->SetMarkerColor(kRed+3);

    litCGraph->SetLineColor(kBlue);
    litCGraph->SetLineWidth(4);
    litCGraph->SetLineStyle(0);
    litCGraph->SetMarkerColor(kBlue);

    litNiGraph->SetLineColor(kBlue);
    litNiGraph->SetLineWidth(4);
    litNiGraph->SetLineStyle(0);
    litNiGraph->SetMarkerColor(kBlue);

    litPbGraph->SetLineColor(kBlue);
    litPbGraph->SetLineWidth(4);
    litPbGraph->SetLineStyle(0);
    litPbGraph->SetMarkerColor(kBlue);

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

        expCGraph->GetXaxis()->SetLabelOffset(0.01);
        expCGraph->GetXaxis()->SetLabelSize(0.06);
        expCGraph->GetXaxis()->SetLabelFont(2);

        expCGraph->GetXaxis()->SetNdivisions(10);
        expCGraph->GetXaxis()->SetTickLength(0.03);

        //expCGraph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        expCGraph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        expCGraph->GetYaxis()->SetTitleSize(0.10);
        expCGraph->GetYaxis()->SetTitleFont(2);
        expCGraph->GetYaxis()->SetTitleOffset(0.7);
        expCGraph->GetYaxis()->CenterTitle();

        expCGraph->GetYaxis()->SetLabelOffset(0.01);
        expCGraph->GetYaxis()->SetLabelSize(0.08);

        expCGraph->GetYaxis()->SetLabelFont(2);
        expCGraph->GetYaxis()->SetNdivisions(5);
        expCGraph->GetYaxis()->SetTickLength(0.02);

        expCGraph->GetXaxis()->SetRangeUser(2,600);
        expCGraph->GetYaxis()->SetRangeUser(0,9);

        gStyle->SetLineWidth(3);

        expCGraph->Draw("");
        litCGraph->Draw("same");
        litNiGraph->Draw("same");
        litPbGraph->Draw("same");
        expCGraph->Draw("same");
        expNiGraph->Draw("same");
        expPbGraph->Draw("same");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.62,0.57,"^{nat}Pb");
        latex.DrawLatex(0.52,0.33,"^{nat}Ni");
        latex.DrawLatex(0.51,0.13,"^{nat}C");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.7,0.70,0.9,0.9);
        legend->SetTextSize(0.08);
        legend->SetTextAlign(12);
        legend->AddEntry(litCGraph,"Analog","l");
        legend->AddEntry(corCGraph,"DSP","l");
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

        expCGraph->Draw("");
        litCGraph->Draw("same");
        litNiGraph->Draw("same");
        litPbGraph->Draw("same");
        corCGraph->Draw("same");
        corNiGraph->Draw("same");
        corPbGraph->Draw("same");

        // Pad dimensions and margins
        corCGraph->GetYaxis()->SetRangeUser(0,9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.08);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.515,0.57,"^{nat}Pb");
        latex.DrawLatex(0.40,0.33,"^{nat}Ni");
        latex.DrawLatex(0.38,0.13,"^{nat}C");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.7,0.70,0.9,0.9);
        legend->SetTextSize(0.08);
        legend->SetTextAlign(12);
        legend->AddEntry(litCGraph,"Analog","l");
        legend->AddEntry(corCGraph,"DSP","l");
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
        relNiGraph->Draw("same");
        relPbGraph->Draw("same");
        
        relCGraph->GetXaxis()->SetRangeUser(2,600);
        relCGraph->GetYaxis()->SetRangeUser(-8,8);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        //latex.DrawLatex(0.35,0.52,"Ni (elem.)");
        //latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.57,0.81,0.94,0.92);
        legend->SetNColumns(3);
        legend->SetTextSize(0.065);
        legend->SetTextAlign(12);
        legend->AddEntry(relCGraph,"{}^{nat}C","l");
        legend->AddEntry(relNiGraph,"{}^{nat}Ni","l");
        legend->AddEntry(relPbGraph,"{}^{nat}Pb","l");

        legend->Draw();
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
        relCorCGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relCorCGraph->GetXaxis()->SetTitleSize(0.08);
        relCorCGraph->GetXaxis()->SetTitleFont(2);
        relCorCGraph->GetXaxis()->SetTitleOffset(1.4);
        relCorCGraph->GetXaxis()->CenterTitle();

        relCorCGraph->GetXaxis()->SetLabelOffset(0.01);
        relCorCGraph->GetXaxis()->SetLabelSize(0.08);
        relCorCGraph->GetXaxis()->SetLabelFont(2);

        relCorCGraph->GetXaxis()->SetNdivisions(10);
        relCorCGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCorCGraph->GetYaxis()->SetTitle("#left(#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}}#right)");
        relCorCGraph->GetYaxis()->SetTitleSize(0.08);
        relCorCGraph->GetYaxis()->SetTitleFont(2);
        relCorCGraph->GetYaxis()->SetTitleOffset(1.4);
        relCorCGraph->GetYaxis()->CenterTitle();

        relCorCGraph->GetYaxis()->SetLabelOffset(0.01);
        relCorCGraph->GetYaxis()->SetLabelSize(0.08);

        relCorCGraph->GetYaxis()->SetLabelFont(2);
        relCorCGraph->GetYaxis()->SetNdivisions(10);
        relCorCGraph->GetYaxis()->SetTickLength(0.02);

        relCorCGraph->Draw("");
        relCorNiGraph->Draw("same");
        relCorPbGraph->Draw("same");

        relCorCGraph->GetXaxis()->SetRangeUser(2,600);
        relCorCGraph->GetYaxis()->SetRangeUser(-8,8);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
        //latex.DrawLatex(0.35,0.52,"Ni (elem.)");
        //latex.DrawLatex(0.32,0.4,"C (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.45,0.81,0.85,0.92);
        legend->SetNColumns(3);
        legend->SetTextSize(0.065);
        legend->SetTextAlign(12);
        legend->AddEntry(relCGraph,"{}^{nat}C","l");
        legend->AddEntry(relNiGraph,"{}^{nat}Ni","l");
        legend->AddEntry(relPbGraph,"{}^{nat}Pb","l");

        legend->Draw();
    }

    //TImage* image = TImage::Create();
    //image->FromPad(canvas);
    //image->WriteImage("../plots/crossSections/FourPanelNi.png");

    // close it all up
    expFile->Close();
    relFile->Close();
    corFile->Close();
    litFile->Close();
}
