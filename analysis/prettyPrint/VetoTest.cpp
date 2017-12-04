{
    string fileName = "/data1/analysis/diagnostics/testVetoEffect/scaledown_noVeto.root";
    string fileName2 = "/data1/analysis/diagnostics/testVetoEffect/scaledown_WithVeto.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    TFile* file2 = new TFile(fileName2.c_str(),"READ");
 
    string PbGraphName = "relDiff_natPb_expLit_percent";
        
    TGraphErrors* PbGraph = (TGraphErrors*)file->Get(PbGraphName.c_str());
    TGraphErrors* PbGraph2 = (TGraphErrors*)file2->Get(PbGraphName.c_str());

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
    PbGraph->SetLineColor(kBlue);
    PbGraph->SetLineWidth(4);
    PbGraph->SetLineStyle(0);
    PbGraph->SetMarkerColor(kBlue);

    PbGraph2->SetLineColor(kRed);
    PbGraph2->SetLineWidth(4);
    PbGraph2->SetLineStyle(0);
    PbGraph2->SetMarkerColor(kRed);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.2);
    gPad->SetTicky(2);

    // X-axis parameters
    PbGraph->GetXaxis()->SetTitle("Energy (MeV)");
    PbGraph->GetXaxis()->SetTitleSize(0.05);
    PbGraph->GetXaxis()->SetTitleFont(2);
    PbGraph->GetXaxis()->SetTitleOffset(1.4);
    PbGraph->GetXaxis()->CenterTitle();

    PbGraph->GetXaxis()->SetLabelOffset(0.01);
    PbGraph->GetXaxis()->SetLabelSize(0.05);
    PbGraph->GetXaxis()->SetLabelFont(2);

    PbGraph->GetXaxis()->SetNdivisions(10);
    PbGraph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    PbGraph->GetYaxis()->SetTitle("(#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}}) x 100");
    PbGraph->GetYaxis()->SetTitleSize(0.06);
    PbGraph->GetYaxis()->SetTitleFont(2);
    PbGraph->GetYaxis()->SetTitleOffset(1.2);
    PbGraph->GetYaxis()->CenterTitle();

    PbGraph->GetYaxis()->SetLabelOffset(0.01);
    PbGraph->GetYaxis()->SetLabelSize(0.05);

    PbGraph->GetYaxis()->SetLabelFont(2);
    PbGraph->GetYaxis()->SetNdivisions(10);
    PbGraph->GetYaxis()->SetTickLength(0.02);

    PbGraph->Draw("");
    PbGraph2->Draw("same");
    
    gPad->SetLogx(1);
    
    PbGraph->GetYaxis()->SetRangeUser(-6,4);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb");
    //latex.DrawLatex(0.47,0.52,"Ni");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.25,0.25,0.55,0.45);
    //legend->SetHeader("");
    legend->SetTextSize(0.04);
    legend->AddEntry(PbGraph,"{}^{nat}Pb, no veto","l");
    legend->AddEntry(PbGraph2,"{}^{nat}Pb, with veto","l");

    legend->Draw();

    file->Close();
}
