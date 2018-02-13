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
    string expFileName = "/data2/analysis/total.root";
    string shiftedFileName = "/data2/analysis/shifted.root";
    string litFileName = "/data2/analysis/literatureData.root";

    string relFileName = "/data2/analysis/relative.root";
    string corFileName = "/data2/analysis/corrected.root";
    string shiftedCorFileName = "/data2/analysis/shifted_Corrected.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* shiftedFile = new TFile(shiftedFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");

    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* corFile = new TFile(corFileName.c_str(),"READ");
    TFile* shiftedCorFile = new TFile(shiftedCorFileName.c_str(),"READ");

    if(!expFile || !shiftedFile || !relFile || !corFile
            || !shiftedCorFile || !litFile)
    {
        cout << "Error: failed to open a required file (check filenames). Exiting..." << endl;
        exit(1);
    }

    string expOH2OGraphName = "ONat_fromH2O";
    string expO18GraphName = "O18+1Barn";

    string litOGraphName = "NatO(n,tot)";
    string litO18GraphName = "18O(n,tot)+1Barn";

    string corOH2OGraphName = "ONat_fromH2O";
    string corO18GraphName = "O18+1Barn";

    string relOH2OGraphName = "ONat_fromH2O, expLit, percent";
    string relO18GraphName = "O18, expLit, percent";

    string relCorOH2OGraphName = "ONat_fromH2O, expLit, corrected, percent";
    string relCorO18GraphName = "O18, expLit, corrected, percent";

    TGraphErrors* expOH2OGraph = (TGraphErrors*)expFile->Get(expOH2OGraphName.c_str());
    TGraphErrors* expO18Graph = (TGraphErrors*)shiftedFile->Get(expO18GraphName.c_str());

    if(!expOH2OGraph || !expO18Graph)
    {
        cout << "Error: failed to open an experimental absolute cross section graph." << endl;
        exit(1);
    }

    TGraphErrors* litOGraph = (TGraphErrors*)litFile->Get(litOGraphName.c_str());
    TGraphErrors* litO18Graph = (TGraphErrors*)litFile->Get(litO18GraphName.c_str());

    if(!litOGraph || !litO18Graph)
    {
        cout << "Error: failed to open a literature absolute cross section graph." << endl;
        exit(1);
    }

    TGraphErrors* relOH2OGraph = (TGraphErrors*)relFile->Get(relOH2OGraphName.c_str());
    TGraphErrors* relO18Graph = (TGraphErrors*)relFile->Get(relO18GraphName.c_str());

    if(!relOH2OGraph || !relO18Graph)
    {
        cout << "Error: failed to open a relative cross section graph." << endl;
        exit(1);
    }

    TGraphErrors* corOH2OGraph = (TGraphErrors*)corFile->Get(corOH2OGraphName.c_str());
    TGraphErrors* corO18Graph = (TGraphErrors*)shiftedFile->Get(corO18GraphName.c_str());

    if(!corOH2OGraph || !corO18Graph)
    {
        cout << "Error: failed to open a corrected absolute cross section graph." << endl;
        exit(1);
    }

    TGraphErrors* relCorOH2OGraph = (TGraphErrors*)relFile->Get(relCorOH2OGraphName.c_str());
    TGraphErrors* relCorO18Graph = (TGraphErrors*)relFile->Get(relCorO18GraphName.c_str());

    if(!relCorOH2OGraph || !relCorO18Graph)
    {
        cout << "Error: failed to open a relative corrected absolute cross section graph." << endl;
        exit(1);
    }

    // Set graph point and line characteristics
    expOH2OGraph->SetLineColor(kRed);
    expOH2OGraph->SetLineWidth(4);
    expOH2OGraph->SetLineStyle(1);
    expOH2OGraph->SetMarkerColor(kRed);

    expO18Graph->SetLineColor(kRed);
    expO18Graph->SetLineWidth(4);
    expO18Graph->SetLineStyle(2);
    expO18Graph->SetMarkerColor(kRed);

    /*********************************/

    relOH2OGraph->SetLineColor(kRed);
    relOH2OGraph->SetLineWidth(4);
    relOH2OGraph->SetLineStyle(1);
    relOH2OGraph->SetMarkerColor(kRed);

    relO18Graph->SetLineColor(kRed+2);
    relO18Graph->SetLineWidth(4);
    relO18Graph->SetLineStyle(1);
    relO18Graph->SetMarkerColor(kRed+2);

    /*********************************/

    corOH2OGraph->SetLineColor(kRed);
    corOH2OGraph->SetLineWidth(4);
    corOH2OGraph->SetLineStyle(1);
    corOH2OGraph->SetMarkerColor(kRed);

    corO18Graph->SetLineColor(kRed);
    corO18Graph->SetLineWidth(4);
    corO18Graph->SetLineStyle(2);
    corO18Graph->SetMarkerColor(kRed);

    /*********************************/

    relCorOH2OGraph->SetLineColor(kRed);
    relCorOH2OGraph->SetLineWidth(4);
    relCorOH2OGraph->SetLineStyle(1);
    relCorOH2OGraph->SetMarkerColor(kRed);

    relCorO18Graph->SetLineColor(kRed+2);
    relCorO18Graph->SetLineWidth(4);
    relCorO18Graph->SetLineStyle(1);
    relCorO18Graph->SetMarkerColor(kRed+2);

    /*********************************/

    litOGraph->SetLineColor(kBlue);
    litOGraph->SetLineWidth(4);
    litOGraph->SetLineStyle(0);
    litOGraph->SetMarkerColor(kBlue);

    litO18Graph->SetLineColor(kBlue);
    litO18Graph->SetLineWidth(4);
    litO18Graph->SetLineStyle(2);
    litO18Graph->SetMarkerColor(kBlue);

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

        expOH2OGraph->GetXaxis()->SetLabelOffset(0.01);
        expOH2OGraph->GetXaxis()->SetLabelSize(0.06);
        expOH2OGraph->GetXaxis()->SetLabelFont(2);

        expOH2OGraph->GetXaxis()->SetNdivisions(10);
        expOH2OGraph->GetXaxis()->SetTickLength(0.03);

        //expOH2OGraph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        expOH2OGraph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        expOH2OGraph->GetYaxis()->SetTitleSize(0.10);
        expOH2OGraph->GetYaxis()->SetTitleFont(2);
        expOH2OGraph->GetYaxis()->SetTitleOffset(0.7);
        expOH2OGraph->GetYaxis()->CenterTitle();

        expOH2OGraph->GetYaxis()->SetLabelOffset(0.01);
        expOH2OGraph->GetYaxis()->SetLabelSize(0.08);

        expOH2OGraph->GetYaxis()->SetLabelFont(2);
        expOH2OGraph->GetYaxis()->SetNdivisions(10);
        expOH2OGraph->GetYaxis()->SetTickLength(0.02);

        expOH2OGraph->GetXaxis()->SetRangeUser(3,600);
        expOH2OGraph->GetYaxis()->SetRangeUser(0,4);

        gStyle->SetLineWidth(3);

        expOH2OGraph->Draw("");
        litO18Graph->Draw("same");
        expO18Graph->Draw("same");

        litOGraph->Draw("same");
        expOH2OGraph->Draw("same");

 
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.07);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.71,0.26,"^{nat}O");
        latex.DrawLatex(0.71,0.50,"^{18}O (+1 barn)");

        // Define legend format and contents
        //TLegend *legend = new TLegend(0.43,0.78,1,0.968);
        //legend->SetTextSize(0.06);
        //legend->SetTextAlign(12);
        //legend->SetNColumns(2);
        //legend->AddEntry(litOGraph,"Analog,{}^{nat}O","l");
        //legend->AddEntry(litO18Graph,"Analog,{}^{18}O","l");
        //legend->AddEntry(expOH2OGraph,"DSP,{}^{nat}O","l");
        //legend->AddEntry(expO18Graph,"DSP,{}^{18}O","l");
        //legend->Draw();

        TLegend *legend = new TLegend(0.7,0.70,0.9,0.9);
        legend->SetTextSize(0.06);
        legend->SetTextAlign(12);
        legend->AddEntry(litOGraph,"Analog","l");
        legend->AddEntry(expOH2OGraph,"DSP","l");
        legend->Draw();

    }

    // second panel
    {
        canvas->cd(2);

        // Pad dimensions and margins
        gPad->SetLeftMargin(0);
        gPad->SetRightMargin(0.005);
        gPad->SetTopMargin(0.03);
        gPad->SetBottomMargin(0.03);
        gPad->SetTicky(2);
        gPad->SetLogx(1);

        gPad->SetFrameLineWidth(3);

        corOH2OGraph->GetXaxis()->SetLabelOffset(0.01);
        corOH2OGraph->GetXaxis()->SetLabelSize(0.06);
        corOH2OGraph->GetXaxis()->SetLabelFont(2);

        corOH2OGraph->GetXaxis()->SetNdivisions(10);
        corOH2OGraph->GetXaxis()->SetTickLength(0.03);

        //corOH2OGraph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        corOH2OGraph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        corOH2OGraph->GetYaxis()->SetTitleSize(0.10);
        corOH2OGraph->GetYaxis()->SetTitleFont(2);
        corOH2OGraph->GetYaxis()->SetTitleOffset(0.7);
        corOH2OGraph->GetYaxis()->CenterTitle();

        corOH2OGraph->GetYaxis()->SetLabelOffset(0.01);
        corOH2OGraph->GetYaxis()->SetLabelSize(0.08);

        corOH2OGraph->GetYaxis()->SetLabelFont(2);
        corOH2OGraph->GetYaxis()->SetNdivisions(10);
        corOH2OGraph->GetYaxis()->SetTickLength(0.02);

        corOH2OGraph->GetXaxis()->SetRangeUser(3,600);
        corOH2OGraph->GetYaxis()->SetRangeUser(0,4);

        gStyle->SetLineWidth(3);

        corOH2OGraph->Draw("");
        litO18Graph->Draw("same");
        corO18Graph->Draw("same");

        litOGraph->Draw("same");
        corOH2OGraph->Draw("same");

 
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.07);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.66,0.26,"^{nat}O");
        latex.DrawLatex(0.66,0.50,"^{18}O (+1 barn)");

        // Define legend format and contents
        //TLegend *legend = new TLegend(0.43,0.78,1,0.968);
        //legend->SetTextSize(0.06);
        //legend->SetTextAlign(12);
        //legend->SetNColumns(2);
        //legend->AddEntry(litOGraph,"Analog,{}^{nat}O","l");
        //legend->AddEntry(litO18Graph,"Analog,{}^{18}O","l");
        //legend->AddEntry(corOH2OGraph,"DSP,{}^{nat}O","l");
        //legend->AddEntry(corO18Graph,"DSP,{}^{18}O","l");
        //legend->Draw();

        TLegend *legend = new TLegend(0.7,0.70,0.9,0.9);
        legend->SetTextSize(0.06);
        legend->SetTextAlign(12);
        legend->AddEntry(litOGraph,"Analog","l");
        legend->AddEntry(corOH2OGraph,"DSP","l");
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
        relOH2OGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relOH2OGraph->GetXaxis()->SetTitleSize(0.08);
        relOH2OGraph->GetXaxis()->SetTitleFont(2);
        relOH2OGraph->GetXaxis()->SetTitleOffset(1.4);
        relOH2OGraph->GetXaxis()->CenterTitle();

        relOH2OGraph->GetXaxis()->SetLabelOffset(0.01);
        relOH2OGraph->GetXaxis()->SetLabelSize(0.08);
        relOH2OGraph->GetXaxis()->SetLabelFont(2);

        relOH2OGraph->GetXaxis()->SetNdivisions(20);
        relOH2OGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relOH2OGraph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relOH2OGraph->GetYaxis()->SetTitleSize(0.08);
        relOH2OGraph->GetYaxis()->SetTitleFont(2);
        relOH2OGraph->GetYaxis()->SetTitleOffset(1.05);
        relOH2OGraph->GetYaxis()->CenterTitle();

        relOH2OGraph->GetYaxis()->SetLabelOffset(0.01);
        relOH2OGraph->GetYaxis()->SetLabelSize(0.08);

        relOH2OGraph->GetYaxis()->SetLabelFont(2);
        relOH2OGraph->GetYaxis()->SetNdivisions(10);
        relOH2OGraph->GetYaxis()->SetTickLength(0.02);

        relOH2OGraph->Draw("");
        relO18Graph->Draw("same");
        relOH2OGraph->Draw("same");

        relOH2OGraph->GetXaxis()->SetRangeUser(3,600);
        relOH2OGraph->GetYaxis()->SetRangeUser(-15,15);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.32,0.4,"O (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.70,0.70,0.85,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relOH2OGraph,"{}^{nat}O","l");
        legend->AddEntry(relO18Graph,"{}^{18}O","l");

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
        relCorOH2OGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relCorOH2OGraph->GetXaxis()->SetTitleSize(0.08);
        relCorOH2OGraph->GetXaxis()->SetTitleFont(2);
        relCorOH2OGraph->GetXaxis()->SetTitleOffset(1.4);
        relCorOH2OGraph->GetXaxis()->CenterTitle();

        relCorOH2OGraph->GetXaxis()->SetLabelOffset(0.01);
        relCorOH2OGraph->GetXaxis()->SetLabelSize(0.08);
        relCorOH2OGraph->GetXaxis()->SetLabelFont(2);

        relCorOH2OGraph->GetXaxis()->SetNdivisions(20);
        relCorOH2OGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCorOH2OGraph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relCorOH2OGraph->GetYaxis()->SetTitleSize(0.08);
        relCorOH2OGraph->GetYaxis()->SetTitleFont(2);
        relCorOH2OGraph->GetYaxis()->SetTitleOffset(1.05);
        relCorOH2OGraph->GetYaxis()->CenterTitle();

        relCorOH2OGraph->GetYaxis()->SetLabelOffset(0.01);
        relCorOH2OGraph->GetYaxis()->SetLabelSize(0.08);

        relCorOH2OGraph->GetYaxis()->SetLabelFont(2);
        relCorOH2OGraph->GetYaxis()->SetNdivisions(10);
        relCorOH2OGraph->GetYaxis()->SetTickLength(0.02);

        relCorOH2OGraph->Draw("");
        relCorO18Graph->Draw("same");
        relCorOH2OGraph->Draw("same");

        relCorOH2OGraph->GetXaxis()->SetRangeUser(3,600);
        relCorOH2OGraph->GetYaxis()->SetRangeUser(-15,15);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.32,0.4,"O (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.63,0.70,0.80,0.95);
        legend->SetTextSize(0.07);
        legend->SetTextAlign(12);
        legend->AddEntry(relCorOH2OGraph,"{}^{nat}O","l");
        legend->AddEntry(relCorO18Graph,"{}^{18}O","l");

        legend->Draw();

        TLine* zeroLine = new TLine(0, 0, 600, 0);
        zeroLine->SetLineColor(kBlack);
        zeroLine->SetLineWidth(3);
        zeroLine->SetLineStyle(9);
        zeroLine->Draw();
    }

    // close it all up
    expFile->Close();
    shiftedFile->Close();
    relFile->Close();
    corFile->Close();
    shiftedCorFile->Close();
    litFile->Close();
}
