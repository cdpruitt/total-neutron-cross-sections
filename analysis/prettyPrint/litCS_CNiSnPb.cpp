{
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
    style->SetHistLineColor(kRed-4);
    style->SetMarkerSize(0);
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

    string fileName = "/data1/analysis/literatureData.root";
    string ramsauerFileName = "~/neutronTCS/theory/ramsauer.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(),"READ");
 
    string CGraphName = "Natural C (n,tot)";
    string NiGraphName = "Natural Ni (n,tot)";
    string SnGraphName = "Natural Sn (n,tot)";
    string PbGraphName = "Natural Pb (n,tot)";
        
    TGraphAsymmErrors* CGraph = (TGraphAsymmErrors*)file->Get(CGraphName.c_str());
    TGraphAsymmErrors* NiGraph = (TGraphAsymmErrors*)file->Get(NiGraphName.c_str());
    TGraphAsymmErrors* SnGraph = (TGraphAsymmErrors*)file->Get(SnGraphName.c_str());
    TGraphAsymmErrors* PbGraph = (TGraphAsymmErrors*)file->Get(PbGraphName.c_str());

    string SACGraphName = "SA_A=12";
    string SANiGraphName = "SA_A=58.7";
    string SASnGraphName = "SA_A=118.7";
    string SAPbGraphName = "SA_A=207.2";

    TGraph* SACGraph = (TGraph*)ramsauerFile->Get(SACGraphName.c_str());
    TGraph* SANiGraph = (TGraph*)ramsauerFile->Get(SANiGraphName.c_str());
    TGraph* SASnGraph = (TGraph*)ramsauerFile->Get(SASnGraphName.c_str());
    TGraph* SAPbGraph = (TGraph*)ramsauerFile->Get(SAPbGraphName.c_str());

    string RamsauerCGraphName = "Ramsauer_A=12";
    string RamsauerNiGraphName = "Ramsauer_A=58.7";
    string RamsauerSnGraphName = "Ramsauer_A=118.7";
    string RamsauerPbGraphName = "Ramsauer_A=207.2";

    TGraph* RamsauerCGraph = (TGraph*)ramsauerFile->Get(RamsauerCGraphName.c_str());
    TGraph* RamsauerNiGraph = (TGraph*)ramsauerFile->Get(RamsauerNiGraphName.c_str());
    TGraph* RamsauerSnGraph = (TGraph*)ramsauerFile->Get(RamsauerSnGraphName.c_str());
    TGraph* RamsauerPbGraph = (TGraph*)ramsauerFile->Get(RamsauerPbGraphName.c_str());

    string OMSnGraphName = "OM_A=118.7";
    TGraph* OMSnGraph = (TGraph*)ramsauerFile->Get(OMSnGraphName.c_str());

    // Set graph point and line characteristics
    CGraph->SetLineWidth(4);
    CGraph->SetLineStyle(0);
    CGraph->SetLineColor(kRed);
    SACGraph->SetLineWidth(2);
    SACGraph->SetLineStyle(9);
    SACGraph->SetLineColor(kRed);

    NiGraph->SetLineWidth(4);
    NiGraph->SetLineStyle(0);
    NiGraph->SetLineColor(kPink-4);
    SANiGraph->SetLineWidth(2);
    SANiGraph->SetLineStyle(9);
    SANiGraph->SetLineColor(kPink-4);

    SnGraph->SetLineWidth(4);
    SnGraph->SetLineStyle(0);
    SnGraph->SetLineColor(kViolet-4);
    SASnGraph->SetLineWidth(2);
    SASnGraph->SetLineStyle(9);
    SASnGraph->SetLineColor(kViolet-4);

    OMSnGraph->SetLineWidth(3);
    OMSnGraph->SetLineStyle(2);
    OMSnGraph->SetLineColor(kViolet-4);

    PbGraph->SetLineWidth(4);
    PbGraph->SetLineStyle(0);
    PbGraph->SetLineColor(kBlue);
    SAPbGraph->SetLineWidth(2);
    SAPbGraph->SetLineStyle(9);
    SAPbGraph->SetLineColor(kBlue);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.15);
    //gPad->SetTicky(2);

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
    CGraph->GetYaxis()->SetTitle("#sigma_{tot} (Barns)");
    CGraph->GetYaxis()->SetTitleSize(0.06);
    CGraph->GetYaxis()->SetTitleFont(2);
    CGraph->GetYaxis()->SetTitleOffset(0.8);
    CGraph->GetYaxis()->CenterTitle();

    CGraph->GetYaxis()->SetLabelOffset(0.01);
    CGraph->GetYaxis()->SetLabelSize(0.05);

    CGraph->GetYaxis()->SetLabelFont(2);
    CGraph->GetYaxis()->SetNdivisions(10);
    CGraph->GetYaxis()->SetTickLength(0.02);

    CGraph->Draw("");
    NiGraph->Draw("same");
    SnGraph->Draw("same");
    PbGraph->Draw("same");

    SACGraph->Draw("same");
    SANiGraph->Draw("same");
    SASnGraph->Draw("same");
    OMSnGraph->Draw("same");
    SAPbGraph->Draw("same");

    //RamsauerCGraph->Draw("same");
    //RamsauerNiGraph->Draw("same");
    //RamsauerSnGraph->Draw("same");
    //RamsauerPbGraph->Draw("same");
 
    gPad->SetLogx(1);
    
    CGraph->GetXaxis()->SetRangeUser(2,500);
    CGraph->GetYaxis()->SetRangeUser(0,9);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.56,0.7,"Pb");
    //latex.DrawLatex(0.52,0.58,"Sn");
    //latex.DrawLatex(0.47,0.49,"Ni");
    //latex.DrawLatex(0.42,0.40,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.83,0.63,0.95,0.95);
    //legend->SetHeader("");
    legend->SetTextSize(0.05);
    //legend->SetNColumns(2);
    legend->AddEntry(PbGraph,"{}^{nat}Pb","l");
    legend->AddEntry(SnGraph,"{}^{nat}Sn","l");
    legend->AddEntry(NiGraph,"{}^{nat}Ni","l");
    legend->AddEntry(CGraph,"{}^{nat}C","l");

    legend->Draw();

    file->Close();
}
