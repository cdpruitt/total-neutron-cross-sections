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

    // read graphs
    string expFileName = "/data2/analysis/corrected.root";
    string litFileName = "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");
    
    string exp112GraphName = "Sn112_corrected";
    string expNatGraphName = "SnNat_corrected";
    string exp124GraphName = "Sn124_corrected";

    string litSnGraphName = "Natural Sn (n,tot)";
    
    TGraphAsymmErrors* exp112Graph = (TGraphAsymmErrors*)expFile->Get(exp112GraphName.c_str());
    TGraphAsymmErrors* expNatGraph = (TGraphAsymmErrors*)expFile->Get(expNatGraphName.c_str());
    TGraphAsymmErrors* exp124Graph = (TGraphAsymmErrors*)expFile->Get(exp124GraphName.c_str());

    TGraphAsymmErrors* litSnGraph = (TGraphAsymmErrors*)litFile->Get(litSnGraphName.c_str());

    // Set graph point and line characteristics
    exp112Graph->SetLineColor(kRed);
    exp112Graph->SetLineWidth(5);
    exp112Graph->SetLineStyle(0);
    exp112Graph->SetMarkerColor(kRed);

    expNatGraph->SetLineColor(kBlack);
    expNatGraph->SetLineWidth(5);
    expNatGraph->SetLineStyle(0);
    expNatGraph->SetMarkerColor(kBlack);

    exp124Graph->SetLineColor(kBlue);
    exp124Graph->SetLineWidth(5);
    exp124Graph->SetLineStyle(0);
    exp124Graph->SetMarkerColor(kBlue);

    litSnGraph->SetLineColor(kGray);
    litSnGraph->SetLineWidth(5);
    litSnGraph->SetLineStyle(2);
    litSnGraph->SetMarkerColor(kGray);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // X-axis parameters
    exp112Graph->GetXaxis()->SetTitle("Energy (MeV)");
    exp112Graph->GetXaxis()->SetTitleSize(0.05);
    exp112Graph->GetXaxis()->SetTitleFont(2);
    exp112Graph->GetXaxis()->SetTitleOffset(1.4);
    exp112Graph->GetXaxis()->CenterTitle();

    exp112Graph->GetXaxis()->SetLabelOffset(0.01);
    exp112Graph->GetXaxis()->SetLabelSize(0.05);
    exp112Graph->GetXaxis()->SetLabelFont(2);

    exp112Graph->GetXaxis()->SetNdivisions(10);
    exp112Graph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    exp112Graph->GetYaxis()->SetTitle("#sigma_{tot} (barns)");
    exp112Graph->GetYaxis()->SetTitleSize(0.06);
    exp112Graph->GetYaxis()->SetTitleFont(2);
    exp112Graph->GetYaxis()->SetTitleOffset(0.8);
    exp112Graph->GetYaxis()->CenterTitle();

    exp112Graph->GetYaxis()->SetLabelOffset(0.01);
    exp112Graph->GetYaxis()->SetLabelSize(0.05);

    exp112Graph->GetYaxis()->SetLabelFont(2);
    exp112Graph->GetYaxis()->SetNdivisions(10);
    exp112Graph->GetYaxis()->SetTickLength(0.02);

    exp112Graph->Draw("");
    expNatGraph->Draw("same");
    exp124Graph->Draw("same");
    litSnGraph->Draw("same");
    exp112Graph->Draw("same");
    expNatGraph->Draw("same");
    exp124Graph->Draw("same");

    gPad->SetLogx(1);
    
    exp112Graph->GetYaxis()->SetRangeUser(1.5,5);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.65,0.65,"Pb (elem.)");
    //latex.DrawLatex(0.35,0.52,"Sn (elem.)");
    //latex.DrawLatex(0.32,0.4,"C (elem.)");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.63, 0.68,0.9,0.9);
    legend->AddEntry(litSnGraph,"{}^{nat}Sn, literature","l");
    legend->AddEntry(exp112Graph,"{}^{112}Sn, lit corrected","l");
    legend->AddEntry(expNatGraph,"{}^{nat}Sn, lit corrected","l");
    legend->AddEntry(exp124Graph,"{}^{124}Sn, lit corrected","l");
    legend->Draw();

    expFile->Close();
    litFile->Close();
}
