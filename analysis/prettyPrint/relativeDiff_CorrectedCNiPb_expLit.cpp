{
    string fileName = "/data1/analysis/relative.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    
    string CGraphName = "NatC, expLit, corrected";
    string NiGraphName = "NatNi, expLit, corrected";
    string PbGraphName = "NatPb, expLit, corrected";
        
    TGraphErrors* CGraph = (TGraphErrors*)file->Get(CGraphName.c_str());
    TGraphErrors* NiGraph = (TGraphErrors*)file->Get(NiGraphName.c_str());
    TGraphErrors* PbGraph = (TGraphErrors*)file->Get(PbGraphName.c_str());

    TStyle* style = (TStyle*)gROOT->FindObject("graphStyle");

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
    style->SetMarkerSize(0.9);
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

    // Set graph point and line characteristics
    CGraph->SetLineColor(kRed-9);
    CGraph->SetLineWidth(4);
    CGraph->SetLineStyle(0);
    CGraph->SetMarkerColor(kRed-9);

    NiGraph->SetLineColor(kRed);
    NiGraph->SetLineWidth(4);
    NiGraph->SetLineStyle(0);
    NiGraph->SetMarkerColor(kRed);

    PbGraph->SetLineColor(kRed+3);
    PbGraph->SetLineWidth(4);
    PbGraph->SetLineStyle(0);
    PbGraph->SetMarkerColor(kRed+3);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.2);
    gPad->SetTicky(2);

    // X-axis parameters
    CGraph->GetXaxis()->SetTitle("Energy (MeV)");
    CGraph->GetXaxis()->SetTitleSize(0.05);
    CGraph->GetXaxis()->SetTitleFont(2);
    CGraph->GetXaxis()->SetTitleOffset(1.4);
    CGraph->GetXaxis()->CenterTitle();

    CGraph->GetXaxis()->SetLabelOffset(0.01);
    CGraph->GetXaxis()->SetLabelSize(0.05);
    CGraph->GetXaxis()->SetLabelFont(2);

    CGraph->GetXaxis()->SetNdivisions(10);
    CGraph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    CGraph->GetYaxis()->SetTitle("(#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}})");
    CGraph->GetYaxis()->SetTitleSize(0.06);
    CGraph->GetYaxis()->SetTitleFont(2);
    CGraph->GetYaxis()->SetTitleOffset(1.4);
    CGraph->GetYaxis()->CenterTitle();

    CGraph->GetYaxis()->SetLabelOffset(0.01);
    CGraph->GetYaxis()->SetLabelSize(0.05);

    CGraph->GetYaxis()->SetLabelFont(2);
    CGraph->GetYaxis()->SetNdivisions(10);
    CGraph->GetYaxis()->SetTickLength(0.02);

    CGraph->Draw("");
    NiGraph->Draw("same");
    PbGraph->Draw("same");
    
    gPad->SetLogx(1);
    
    CGraph->GetYaxis()->SetRangeUser(-0.05,0.05);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Ni");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.25,0.25,0.35,0.45);
    //legend->SetHeader("");
    legend->SetTextSize(0.04);
    legend->AddEntry(CGraph,"{}^{nat}C","l");
    legend->AddEntry(NiGraph,"{}^{nat}Ni","l");
    legend->AddEntry(PbGraph,"{}^{nat}Pb","l");

    legend->Draw();

    file->Close();
}
