void relativeDiff_O_expLit()
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

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.2);
    gPad->SetTicky(2);

    string fileName = "/data3/analysis/scaledown.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    
    string NatOGraphName = "relDiff_natO_expLit_percent";
    string O18GraphName = "relDiff_18O_expLit_percent";
        
    TGraphErrors* NatOGraph = (TGraphErrors*)file->Get(NatOGraphName.c_str());
    TGraphErrors* O18Graph = (TGraphErrors*)file->Get(O18GraphName.c_str());

    // Set graph point and line characteristics
    NatOGraph->SetLineColor(kRed-9);
    NatOGraph->SetLineWidth(4);
    NatOGraph->SetLineStyle(0);
    NatOGraph->SetMarkerColor(kRed-9);

    O18Graph->SetLineColor(kRed);
    O18Graph->SetLineWidth(4);
    O18Graph->SetLineStyle(0);
    O18Graph->SetMarkerColor(kRed);

    // X-axis parameters
    NatOGraph->GetXaxis()->SetTitle("Energy (MeV)");
    NatOGraph->GetXaxis()->SetTitleSize(0.05);
    NatOGraph->GetXaxis()->SetTitleFont(2);
    NatOGraph->GetXaxis()->SetTitleOffset(1.4);
    NatOGraph->GetXaxis()->CenterTitle();

    NatOGraph->GetXaxis()->SetLabelOffset(0.01);
    NatOGraph->GetXaxis()->SetLabelSize(0.05);
    NatOGraph->GetXaxis()->SetLabelFont(2);

    NatOGraph->GetXaxis()->SetNdivisions(10);
    NatOGraph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    NatOGraph->GetYaxis()->SetTitle("(#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}}) x 100");
    NatOGraph->GetYaxis()->SetTitleSize(0.06);
    NatOGraph->GetYaxis()->SetTitleFont(2);
    NatOGraph->GetYaxis()->SetTitleOffset(1.2);
    NatOGraph->GetYaxis()->CenterTitle();

    NatOGraph->GetYaxis()->SetLabelOffset(0.01);
    NatOGraph->GetYaxis()->SetLabelSize(0.05);

    NatOGraph->GetYaxis()->SetLabelFont(2);
    NatOGraph->GetYaxis()->SetNdivisions(10);
    NatOGraph->GetYaxis()->SetTickLength(0.02);

    NatOGraph->Draw("");
    O18Graph->Draw("same");
    
    gPad->SetLogx(1);
    
    NatOGraph->GetYaxis()->SetRangeUser(-15,15);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.8,0.30,0.93,0.5);
    //legend->SetHeader("");
    legend->SetTextSize(0.06);
    legend->AddEntry(NatOGraph,"{}^{nat}O","l");
    legend->AddEntry(O18Graph,"{}^{18}O","l");

    legend->Draw();

    file->Close();
}
