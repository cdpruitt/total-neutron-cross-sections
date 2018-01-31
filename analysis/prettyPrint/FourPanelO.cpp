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
    string relFileName = "/data2/analysis/relative.root";
    string corFileName = "/data2/analysis/corrected.root";
    string litFileName = "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* relFile = new TFile(relFileName.c_str(),"READ");
    TFile* corFile = new TFile(corFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");
    
    string expOGraphName = "ONat";
    string expD2OGraphName = "D2ONat";

    string relOGraphName = "NatO, expLit, percent";
    string relD2OGraphName = "NatD2O, expLit, percent";

    string relCorOGraphName = "NatO, expLit, corrected, percent";
    string relCorD2OGraphName = "NatD2O, expLit, corrected, percent";

    string litOGraphName = "Natural carbon (n,tot)";
    string litD2OGraphName = "Natural D2O (n,tot)";
    
    TGraphErrors* expOGraph = (TGraphErrors*)expFile->Get(expOGraphName.c_str());
    TGraphErrors* expD2OGraph = (TGraphErrors*)expFile->Get(expD2OGraphName.c_str());

    TGraphErrors* relOGraph = (TGraphErrors*)relFile->Get(relOGraphName.c_str());
    TGraphErrors* relD2OGraph = (TGraphErrors*)relFile->Get(relD2OGraphName.c_str());

    TGraphErrors* corOGraph = (TGraphErrors*)corFile->Get(expOGraphName.c_str());
    TGraphErrors* corD2OGraph = (TGraphErrors*)corFile->Get(expD2OGraphName.c_str());

    TGraphErrors* relCorOGraph = (TGraphErrors*)relFile->Get(relCorOGraphName.c_str());
    TGraphErrors* relCorD2OGraph = (TGraphErrors*)relFile->Get(relCorD2OGraphName.c_str());

    TGraphErrors* litOGraph = (TGraphErrors*)litFile->Get(litOGraphName.c_str());
    TGraphErrors* litD2OGraph = (TGraphErrors*)litFile->Get(litD2OGraphName.c_str());

    // Set graph point and line characteristics
    expOGraph->SetLineColor(kRed);
    expOGraph->SetMarkerColor(kRed);

    expD2OGraph->SetLineColor(kRed);
    expD2OGraph->SetMarkerColor(kRed);

    relOGraph->SetLineColor(kRed-9);
    relOGraph->SetLineWidth(4);
    relOGraph->SetLineStyle(0);
    relOGraph->SetMarkerColor(kRed-9);

    relD2OGraph->SetLineColor(kRed);
    relD2OGraph->SetLineWidth(4);
    relD2OGraph->SetLineStyle(0);
    relD2OGraph->SetMarkerColor(kRed);

    corOGraph->SetLineColor(kRed);
    corOGraph->SetMarkerColor(kRed);

    corD2OGraph->SetLineColor(kRed);
    corD2OGraph->SetMarkerColor(kRed);

    relCorOGraph->SetLineColor(kRed-9);
    relCorOGraph->SetLineWidth(4);
    relCorOGraph->SetLineStyle(0);
    relCorOGraph->SetMarkerColor(kRed-9);

    relCorD2OGraph->SetLineColor(kRed);
    relCorD2OGraph->SetLineWidth(4);
    relCorD2OGraph->SetLineStyle(0);
    relCorD2OGraph->SetMarkerColor(kRed);

    litOGraph->SetLineColor(kBlue);
    litOGraph->SetLineWidth(4);
    litOGraph->SetLineStyle(0);
    litOGraph->SetMarkerColor(kBlue);

    litD2OGraph->SetLineColor(kBlue);
    litD2OGraph->SetLineWidth(4);
    litD2OGraph->SetLineStyle(0);
    litD2OGraph->SetMarkerColor(kBlue);

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

        expOGraph->GetXaxis()->SetLabelOffset(0.01);
        expOGraph->GetXaxis()->SetLabelSize(0.06);
        expOGraph->GetXaxis()->SetLabelFont(2);

        expOGraph->GetXaxis()->SetNdivisions(10);
        expOGraph->GetXaxis()->SetTickLength(0.03);

        //expOGraph->GetXaxis()->SetLineWidth(3);

        // Y-axis parameters
        expOGraph->GetYaxis()->SetTitle("#sigma_{tot} [b]");
        expOGraph->GetYaxis()->SetTitleSize(0.10);
        expOGraph->GetYaxis()->SetTitleFont(2);
        expOGraph->GetYaxis()->SetTitleOffset(0.7);
        expOGraph->GetYaxis()->CenterTitle();

        expOGraph->GetYaxis()->SetLabelOffset(0.01);
        expOGraph->GetYaxis()->SetLabelSize(0.08);

        expOGraph->GetYaxis()->SetLabelFont(2);
        expOGraph->GetYaxis()->SetNdivisions(10);
        expOGraph->GetYaxis()->SetTickLength(0.02);

        expOGraph->GetXaxis()->SetRangeUser(2,500);
        expOGraph->GetYaxis()->SetRangeUser(0,9);

        gStyle->SetLineWidth(3);

        expOGraph->Draw("");
        litOGraph->Draw("same");
        litD2OGraph->Draw("same");
        expOGraph->Draw("same");
        expD2OGraph->Draw("same");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.06);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.52,0.33,"^{nat}D2O");
        latex.DrawLatex(0.51,0.13,"^{nat}O");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.7,0.70,0.9,0.9);
        legend->SetTextSize(0.06);
        legend->SetTextAlign(12);
        legend->AddEntry(litOGraph,"Analog","l");
        legend->AddEntry(corOGraph,"DSP","l");
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

        expOGraph->Draw("");
        litOGraph->Draw("same");
        litD2OGraph->Draw("same");
        corOGraph->Draw("same");
        corD2OGraph->Draw("same");

        // Pad dimensions and margins
        corOGraph->GetYaxis()->SetRangeUser(0,9);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.06);
        latex.SetTextAlign(13); // align at top
        latex.DrawLatex(0.40,0.33,"^{nat}D2O");
        latex.DrawLatex(0.38,0.13,"^{nat}O");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.7,0.70,0.9,0.9);
        legend->SetTextSize(0.06);
        legend->SetTextAlign(12);
        legend->AddEntry(litOGraph,"Analog","l");
        legend->AddEntry(corOGraph,"DSP","l");
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
        relOGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relOGraph->GetXaxis()->SetTitleSize(0.08);
        relOGraph->GetXaxis()->SetTitleFont(2);
        relOGraph->GetXaxis()->SetTitleOffset(1.4);
        relOGraph->GetXaxis()->CenterTitle();

        relOGraph->GetXaxis()->SetLabelOffset(0.01);
        relOGraph->GetXaxis()->SetLabelSize(0.08);
        relOGraph->GetXaxis()->SetLabelFont(2);

        relOGraph->GetXaxis()->SetNdivisions(10);
        relOGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relOGraph->GetYaxis()->SetTitle("#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}} [%]");
        relOGraph->GetYaxis()->SetTitleSize(0.08);
        relOGraph->GetYaxis()->SetTitleFont(2);
        relOGraph->GetYaxis()->SetTitleOffset(1.05);
        relOGraph->GetYaxis()->CenterTitle();

        relOGraph->GetYaxis()->SetLabelOffset(0.01);
        relOGraph->GetYaxis()->SetLabelSize(0.08);

        relOGraph->GetYaxis()->SetLabelFont(2);
        relOGraph->GetYaxis()->SetNdivisions(5);
        relOGraph->GetYaxis()->SetTickLength(0.02);

        relOGraph->Draw("");
        relD2OGraph->Draw("same");
        
        relOGraph->GetXaxis()->SetRangeUser(2,500);
        relOGraph->GetYaxis()->SetRangeUser(-10,10);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.35,0.52,"D2O (elem.)");
        //latex.DrawLatex(0.32,0.4,"O (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.57,0.81,0.94,0.92);
        legend->SetNColumns(3);
        legend->SetTextSize(0.055);
        legend->SetTextAlign(12);
        legend->AddEntry(relOGraph,"{}^{nat}O","l");
        legend->AddEntry(relD2OGraph,"{}^{nat}D2O","l");

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
        relCorOGraph->GetXaxis()->SetTitle("Energy (MeV)");
        relCorOGraph->GetXaxis()->SetTitleSize(0.08);
        relCorOGraph->GetXaxis()->SetTitleFont(2);
        relCorOGraph->GetXaxis()->SetTitleOffset(1.4);
        relCorOGraph->GetXaxis()->CenterTitle();

        relCorOGraph->GetXaxis()->SetLabelOffset(0.01);
        relCorOGraph->GetXaxis()->SetLabelSize(0.08);
        relCorOGraph->GetXaxis()->SetLabelFont(2);

        relCorOGraph->GetXaxis()->SetNdivisions(20);
        relCorOGraph->GetXaxis()->SetTickLength(0.03);

        // Y-axis parameters
        relCorOGraph->GetYaxis()->SetTitle("#left(#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}}#right)");
        relCorOGraph->GetYaxis()->SetTitleSize(0.08);
        relCorOGraph->GetYaxis()->SetTitleFont(2);
        relCorOGraph->GetYaxis()->SetTitleOffset(1.4);
        relCorOGraph->GetYaxis()->CenterTitle();

        relCorOGraph->GetYaxis()->SetLabelOffset(0.01);
        relCorOGraph->GetYaxis()->SetLabelSize(0.08);

        relCorOGraph->GetYaxis()->SetLabelFont(2);
        relCorOGraph->GetYaxis()->SetNdivisions(5);
        relCorOGraph->GetYaxis()->SetTickLength(0.02);

        relCorOGraph->Draw("");
        relCorD2OGraph->Draw("same");

        relCorOGraph->GetXaxis()->SetRangeUser(2,500);
        relCorOGraph->GetYaxis()->SetRangeUser(-10,10);

        //TLatex latex;
        //latex.SetNDC();
        //latex.SetTextSize(0.035);
        //latex.SetTextAlign(13); // align at top
        //latex.DrawLatex(0.35,0.52,"D2O (elem.)");
        //latex.DrawLatex(0.32,0.4,"O (elem.)");

        // Define legend format and contents
        TLegend *legend = new TLegend(0.45,0.81,0.85,0.92);
        legend->SetNColumns(3);
        legend->SetTextSize(0.055);
        legend->SetTextAlign(12);
        legend->AddEntry(relOGraph,"{}^{nat}O","l");
        legend->AddEntry(relD2OGraph,"{}^{nat}D2O","l");

        legend->Draw();
    }

    //TImage* image = TImage::Create();
    //image->FromPad(canvas);
    //image->WriteImage("../plots/crossSections/FourPanelD2O.png");

    // close it all up
    expFile->Close();
    relFile->Close();
    corFile->Close();
    litFile->Close();
}
