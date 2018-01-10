{
    string fileName = "/data2/analysis/relative.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    
    string relGraphName1 = "Natural O, expLit";
    string relGraphName2 = "O18, expLit";
    
    TGraphErrors* relGraph1 = (TGraphErrors*)file->Get(relGraphName1.c_str());
    TGraphErrors* relGraph2 = (TGraphErrors*)file->Get(relGraphName2.c_str());

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
    relGraph1->SetLineColor(kRed);
    relGraph1->SetLineWidth(4);
    relGraph1->SetLineStyle(0);
    relGraph1->SetMarkerColor(kRed);

    relGraph2->SetLineColor(kBlue);
    relGraph2->SetLineWidth(4);
    relGraph2->SetLineStyle(0);
    relGraph2->SetMarkerColor(kBlue);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.2);
    gPad->SetTicky(2);

    // X-axis parameters
    relGraph1->GetXaxis()->SetTitle("Energy (MeV)");
    relGraph1->GetXaxis()->SetTitleSize(0.05);
    relGraph1->GetXaxis()->SetTitleFont(2);
    relGraph1->GetXaxis()->SetTitleOffset(1.4);
    relGraph1->GetXaxis()->CenterTitle();

    relGraph1->GetXaxis()->SetLabelOffset(0.01);
    relGraph1->GetXaxis()->SetLabelSize(0.05);
    relGraph1->GetXaxis()->SetLabelFont(2);

    relGraph1->GetXaxis()->SetNdivisions(10);
    relGraph1->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    relGraph1->GetYaxis()->SetTitle("(#frac{#sigma_{exp} - #sigma_{lit}}{#sigma_{exp} + #sigma_{lit}})");
    relGraph1->GetYaxis()->SetTitleSize(0.06);
    relGraph1->GetYaxis()->SetTitleFont(2);
    relGraph1->GetYaxis()->SetTitleOffset(1.4);
    relGraph1->GetYaxis()->CenterTitle();

    relGraph1->GetYaxis()->SetLabelOffset(0.01);
    relGraph1->GetYaxis()->SetLabelSize(0.05);

    relGraph1->GetYaxis()->SetLabelFont(2);
    relGraph1->GetYaxis()->SetNdivisions(10);
    relGraph1->GetYaxis()->SetTickLength(0.02);

    relGraph1->Draw("");
    relGraph2->Draw("same");
    
    gPad->SetLogx(1);
    
    relGraph1->GetXaxis()->SetRangeUser(5,300);
    relGraph1->GetYaxis()->SetRangeUser(-0.25,0.15);

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.47,0.52,"Ni");
    //latex.DrawLatex(0.32,0.4,"C");

    // Define legend format and contents
    TLegend *legend = new TLegend(0.75,0.75,0.90,0.90);
    //legend->SetHeader("");
    legend->SetTextSize(0.06);
    legend->AddEntry(relGraph1,"{}^{nat}O","l");
    legend->AddEntry(relGraph2,"{}^{18}O","l");

    legend->Draw();

    file->Close();
}
