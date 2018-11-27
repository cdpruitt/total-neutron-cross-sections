{
    string fileName = "/data2/analysis/relative.root";
    string ramsauerFileName = "../../theory/ramsauer.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(), "READ");
    
    string relGraphName = "Sn124Sn112, percent";
    string SARelDiffGraphName = "RelDiff124_112";
        
    TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)file->Get(relGraphName.c_str());
    TGraph* SARelDiffGraph = (TGraph*)ramsauerFile->Get(SARelDiffGraphName.c_str());

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
    relGraph->SetLineColor(kRed);
    relGraph->SetLineWidth(4);
    relGraph->SetLineStyle(0);
    relGraph->SetMarkerColor(kRed);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // X-axis parameters
    relGraph->GetXaxis()->SetTitle("Energy (MeV)");
    relGraph->GetXaxis()->SetTitleSize(0.05);
    relGraph->GetXaxis()->SetTitleFont(2);
    relGraph->GetXaxis()->SetTitleOffset(1.4);
    relGraph->GetXaxis()->CenterTitle();

    relGraph->GetXaxis()->SetLabelOffset(0.01);
    relGraph->GetXaxis()->SetLabelSize(0.05);
    relGraph->GetXaxis()->SetLabelFont(2);

    relGraph->GetXaxis()->SetNdivisions(10);
    relGraph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    relGraph->GetYaxis()->SetTitle("(#frac{#sigma_{124} - #sigma_{112}}{#sigma_{124} + #sigma_{112}})");
    relGraph->GetYaxis()->SetTitleSize(0.06);
    relGraph->GetYaxis()->SetTitleFont(2);
    relGraph->GetYaxis()->SetTitleOffset(1.0);
    relGraph->GetYaxis()->CenterTitle();

    relGraph->GetYaxis()->SetLabelOffset(0.01);
    relGraph->GetYaxis()->SetLabelSize(0.05);

    relGraph->GetYaxis()->SetLabelFont(2);
    relGraph->GetYaxis()->SetNdivisions(5);
    relGraph->GetYaxis()->SetTickLength(0.02);

    SARelDiffGraph->SetLineStyle(9);
    SARelDiffGraph->SetLineWidth(3);
    SARelDiffGraph->SetLineColor(kGray);

    relGraph->Draw("AL");
    SARelDiffGraph->Draw("same");

    gPad->SetLogx(1);
    
    relGraph->GetYaxis()->SetRangeUser(0.2,4);
    relGraph->GetXaxis()->SetLimits(5,600);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.47,0.52,"Ni");
    //latex.DrawLatex(0.32,0.4,"C");

    file->Close();
}