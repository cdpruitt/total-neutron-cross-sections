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

    TFile* file = new TFile(fileName.c_str(),"READ");
    
    string CGraphName = "Natural C (n,tot)";
    string NiGraphName = "Natural Ni (n,tot)";
    string SnGraphName = "Natural Sn (n,tot)";
    string PbGraphName = "Natural Pb (n,tot)";
        
    TGraphErrors* CGraph = (TGraphErrors*)file->Get(CGraphName.c_str());
    TGraphErrors* NiGraph = (TGraphErrors*)file->Get(NiGraphName.c_str());
    TGraphErrors* SnGraph = (TGraphErrors*)file->Get(SnGraphName.c_str());
    TGraphErrors* PbGraph = (TGraphErrors*)file->Get(PbGraphName.c_str());

    // Set graph point and line characteristics
    CGraph->SetLineWidth(3);
    CGraph->SetLineStyle(0);

    NiGraph->SetLineWidth(3);
    NiGraph->SetLineStyle(0);

    SnGraph->SetLineWidth(3);
    SnGraph->SetLineStyle(0);

    PbGraph->SetLineWidth(3);
    PbGraph->SetLineStyle(0);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.2);
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
    
    gPad->SetLogx(1);
    
    CGraph->GetXaxis()->SetRangeUser(2,500);
    CGraph->GetYaxis()->SetRangeUser(0,7);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13); // align at top
    latex.DrawLatex(0.56,0.7,"Pb");
    latex.DrawLatex(0.52,0.58,"Sn");
    latex.DrawLatex(0.47,0.49,"Ni");
    latex.DrawLatex(0.42,0.40,"C");

    // Define legend format and contents
    //TLegend *legend = new TLegend(0.25,0.25,0.35,0.45);
    //legend->SetHeader("");
    //legend->SetTextSize(0.04);
    //legend->AddEntry(CGraph,"{}^{nat}C","l");
    //legend->AddEntry(NiGraph,"{}^{nat}Ni","l");
    //legend->AddEntry(PbGraph,"{}^{nat}Pb","l");

    //legend->Draw();

    file->Close();
}
