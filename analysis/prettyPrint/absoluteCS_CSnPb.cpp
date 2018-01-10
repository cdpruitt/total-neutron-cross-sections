void absoluteCS_CSnPb()
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
    string expFileName = "/data2/analysis/total.root";
    string litFileName = "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileName.c_str(),"READ");
    TFile* litFile = new TFile(litFileName.c_str(),"READ");
    
    string expCGraphName = "CNat";
    string expSnGraphName = "SnNat";
    string expPbGraphName = "PbNat";

    string litCGraphName = "Natural C (n,tot)";
    string litSnGraphName = "Natural Sn (n,tot)";
    string litPbGraphName = "Natural Pb (n,tot)";
    
    TGraphErrors* expCGraph = (TGraphErrors*)expFile->Get(expCGraphName.c_str());
    TGraphErrors* expSnGraph = (TGraphErrors*)expFile->Get(expSnGraphName.c_str());
    TGraphErrors* expPbGraph = (TGraphErrors*)expFile->Get(expPbGraphName.c_str());

    TGraphErrors* litCGraph = (TGraphErrors*)litFile->Get(litCGraphName.c_str());
    TGraphErrors* litSnGraph = (TGraphErrors*)litFile->Get(litSnGraphName.c_str());
    TGraphErrors* litPbGraph = (TGraphErrors*)litFile->Get(litPbGraphName.c_str());

    // Set graph point and line characteristics
    expCGraph->SetLineColor(kRed);
    expCGraph->SetLineWidth(5);
    expCGraph->SetLineStyle(0);
    expCGraph->SetMarkerColor(kRed);

    expSnGraph->SetLineColor(kRed);
    expSnGraph->SetLineWidth(5);
    expSnGraph->SetLineStyle(0);
    expSnGraph->SetMarkerColor(kRed);

    expPbGraph->SetLineColor(kRed);
    expPbGraph->SetLineWidth(5);
    expPbGraph->SetLineStyle(0);
    expPbGraph->SetMarkerColor(kRed);

    litCGraph->SetLineColor(kBlack);
    litCGraph->SetLineWidth(3);
    litCGraph->SetLineStyle(2);
    litCGraph->SetMarkerColor(kBlack);

    litSnGraph->SetLineColor(kBlack);
    litSnGraph->SetLineWidth(3);
    litSnGraph->SetLineStyle(2);
    litSnGraph->SetMarkerColor(kBlack);

    litPbGraph->SetLineColor(kBlack);
    litPbGraph->SetLineWidth(3);
    litPbGraph->SetLineStyle(2);
    litPbGraph->SetMarkerColor(kBlack);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    // X-axis parameters
    expCGraph->GetXaxis()->SetTitle("Energy (MeV)");
    expCGraph->GetXaxis()->SetTitleSize(0.05);
    expCGraph->GetXaxis()->SetTitleFont(2);
    expCGraph->GetXaxis()->SetTitleOffset(1.4);
    expCGraph->GetXaxis()->CenterTitle();

    expCGraph->GetXaxis()->SetLabelOffset(0.01);
    expCGraph->GetXaxis()->SetLabelSize(0.05);
    expCGraph->GetXaxis()->SetLabelFont(2);

    expCGraph->GetXaxis()->SetNdivisions(10);
    expCGraph->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    expCGraph->GetYaxis()->SetTitle("#sigma_{tot} (barns)");
    expCGraph->GetYaxis()->SetTitleSize(0.06);
    expCGraph->GetYaxis()->SetTitleFont(2);
    expCGraph->GetYaxis()->SetTitleOffset(0.8);
    expCGraph->GetYaxis()->CenterTitle();

    expCGraph->GetYaxis()->SetLabelOffset(0.01);
    expCGraph->GetYaxis()->SetLabelSize(0.05);

    expCGraph->GetYaxis()->SetLabelFont(2);
    expCGraph->GetYaxis()->SetNdivisions(10);
    expCGraph->GetYaxis()->SetTickLength(0.02);

    expCGraph->Draw("");
    expSnGraph->Draw("same");
    expPbGraph->Draw("same");
    litCGraph->Draw("same");
    litSnGraph->Draw("same");
    litPbGraph->Draw("same");
    expCGraph->Draw("same");
    expSnGraph->Draw("same");
    expPbGraph->Draw("same");

    gPad->SetLogx(1);
    
    expCGraph->GetYaxis()->SetRangeUser(0,9);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextAlign(13); // align at top
    latex.DrawLatex(0.65,0.65,"Pb (elem.)");
    latex.DrawLatex(0.35,0.52,"Sn (elem.)");
    latex.DrawLatex(0.32,0.4,"C (elem.)");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
    legend->AddEntry(litCGraph,"lit data (analog)","l");
    legend->AddEntry(expCGraph,"new data (DSP)","l");
    legend->Draw();

    expFile->Close();
    litFile->Close();
}
