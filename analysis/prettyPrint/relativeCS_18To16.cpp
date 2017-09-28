{
    TStyle * style = (TStyle*)gROOT->FindObject("graphStyle");

    if(!style)      
    {
        style = new TStyle("graphStyle","graphStyle");
    }

    TCanvas* c = new TCanvas("c1","",1200,1200);

    style->SetOptStat(0);
    style->SetOptTitle(0);    
    style->SetPalette(1,0);
    style->SetCanvasColor(10);      
    style->SetCanvasBorderMode(0);    
    style->SetFrameLineWidth(3);
    style->SetFrameFillColor(10);
    style->SetPadColor(10);
    style->SetHistLineWidth(3);
    style->SetHistLineColor(kBlue);
    style->SetFuncWidth(3);
    style->SetFuncColor(kRed);
    style->SetLabelColor(kBlack,"xyz");
    style->SetTitleSize(0.06,"xyz");
    style->SetTitleFillColor(10);
    style->SetTitleTextColor(kBlack);
    style->SetEndErrorSize(0);

    gROOT->SetStyle("graphStyle");
    gROOT->ForceStyle();

    string graphFileName = "/data3/analysis/scaledown.root";
    TFile* graphFile = new TFile(graphFileName.c_str(),"READ");

    string graphName = "relDiffO_percent";
    TGraphErrors* graph = (TGraphErrors*)graphFile->Get(graphName.c_str());

    graph->SetLineColor(kRed);
    graph->SetLineWidth(5);
    graph->SetLineStyle(0);

    graph->SetMarkerColor(kRed);
    graph->SetMarkerSize(0.9);
    graph->SetMarkerStyle(20);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.2);
    gPad->SetTicky(2);

    // X-axis parameters
    graph->GetXaxis()->SetTitle("Energy (MeV)");
    graph->GetXaxis()->SetTitleSize(0.05);
    graph->GetXaxis()->SetTitleFont(2);
    graph->GetXaxis()->SetTitleOffset(1.4);
    graph->GetXaxis()->CenterTitle();

    graph->GetXaxis()->SetLabelOffset(0.01);
    graph->GetXaxis()->SetLabelSize(0.05);
    graph->GetXaxis()->SetLabelFont(2);

    graph->GetXaxis()->SetNdivisions(10);
    graph->GetXaxis()->SetTickLength(0.03);

    graph->GetXaxis()->SetRange(3,200);

    // Y-axis parameters
    graph->GetYaxis()->SetTitle("(#frac{#sigma_{18} - #sigma_{16}}{#sigma_{18} + #sigma_{16}}) x 100");
    graph->GetYaxis()->SetTitleSize(0.06);
    graph->GetYaxis()->SetTitleFont(2);
    graph->GetYaxis()->SetTitleOffset(1.35);
    graph->GetYaxis()->CenterTitle();

    graph->GetYaxis()->SetLabelOffset(0.01);
    graph->GetYaxis()->SetLabelSize(0.05);

    graph->GetYaxis()->SetLabelFont(2);
    graph->GetYaxis()->SetNdivisions(10);
    graph->GetYaxis()->SetTickLength(0.02);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.47,0.95,"#frac{#sigma^{^{124}Sn}-#sigma^{^{112}Sn}}{#sigma^{^{124}Sn}+#sigma^{^{112}Sn}}");;

    graph->Draw();

    gPad->SetLogx(1);

    TLine *line1 = new TLine(1.85,4.00,400,4.00);
    line1->SetLineStyle(7);
    line1->SetLineWidth(2);
    line1->SetLineColor(kBlack);
    line1->Draw();

    TLine *line2 = new TLine(1.85,1.98,400,1.98);
    line2->SetLineStyle(7);
    line2->SetLineWidth(2);
    line2->SetLineColor(kBlack);
    line2->Draw();

    graph->GetYaxis()->SetRangeUser(-35,35);

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.75,0.7,0.9,0.9);
    //legend->SetHeader("data","C");
    //legend->AddEntry(graph,"{}^{nat}Pb","l");
    //legend->Draw();

    graphFile->Close();
}
