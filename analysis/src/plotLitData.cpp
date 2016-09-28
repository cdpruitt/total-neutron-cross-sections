
/*****************************************************************************/
// Style section: edit to change style of plotted datasets
const std::string colors[4] = {"kRed","kRed","kRed","kRed"};
const std::string markers[4] = {"kPlus","kPlus","kPlus","kPlus"};
/*****************************************************************************/

// ROOT macro for formatting and combining plots

void plotLitData()
{
    gStyle->SetOptStat(0);

    TStyle * Sty = (TStyle*)gROOT->FindObject("MyStyle");
    if(!Sty)      
    {
        Sty = new TStyle("MyStyle","MyStyle");
    }

    Sty->SetOptTitle(0);    
    Sty->SetOptStat(0);
    Sty->SetPalette(1,0);
    Sty->SetCanvasColor(10);      
    Sty->SetCanvasBorderMode(0);    
    Sty->SetFrameLineWidth(3);
    Sty->SetFrameFillColor(10);
    Sty->SetPadColor(10);
    Sty->SetPadTickX(1);
    Sty->SetPadTickY(1);
    Sty->SetPadBottomMargin(.15);
    Sty->SetPadTopMargin(.03);
    Sty->SetPadLeftMargin(.14);
    Sty->SetPadRightMargin(.06);
    Sty->SetHistLineWidth(3);
    Sty->SetHistLineColor(kBlue);
    Sty->SetFuncWidth(3);
    Sty->SetMarkerColor(kBlue);
    Sty->SetLineWidth(1);
    Sty->SetLabelSize(0.06,"xyz");
    Sty->SetLabelOffset(0.02,"y");
    Sty->SetLabelOffset(0.02,"x");
    Sty->SetLabelColor(kBlack,"xyz");
    Sty->SetMarkerSize(1);
    Sty->SetMarkerStyle(21);
    Sty->SetTitleSize(0.06,"xyz");
    Sty->SetTitleOffset(1.15,"y");
    Sty->SetTitleOffset(1.1,"x");
    Sty->SetTitleFillColor(10);
    Sty->SetTitleTextColor(kBlack);
    Sty->SetTickLength(.03,"xz");
    Sty->SetTickLength(.02,"y");
    Sty->SetNdivisions(5,"x");
    Sty->SetNdivisions(10,"yz");
    Sty->SetEndErrorSize(0);
    gROOT->SetStyle("MyStyle");
    gROOT->ForceStyle();

    // open calculated data file
    std::stringstream fileInName;
    fileInName << "/data2/analysis/literatureData.root";
    TFile* fileIn = new TFile(fileInName.str().c_str(),"READ");

    TGraphErrors* plot1 = (TGraphErrors*)gDirectory->Get("diff_sum");
    TGraphErrors* plot2 = (TGraphErrors*)gDirectory->Get("diff_sumLog");

    if(!plot1 || !plot2)
    {
        exit(1);
    }

    // open experimental data file
    std::stringstream expFileName;
    expFileName << "/data3/analysis/total.root";
    TFile* expFile = new TFile(expFileName.str().c_str(),"READ");

    TH1D* relativeCS = (TH1D*)gDirectory->Get("relativeSnCS");
    TH1D* relativeCSLog = (TH1D*)gDirectory->Get("relativeSnCSLog");

    TH1D* relativeCS_W = (TH1D*)gDirectory->Get("relativeSnCS_W");

    TGraphErrors* relativeCSLogGraph = (TGraphErrors*)gDirectory->Get("relativeGraph");
    relativeCSLogGraph->SetMarkerColor(kRed);
    //relativeCSLogGraph->SetMarkerSize(1);

    // create output file
    std::stringstream fileOutName;
    fileOutName << "plotLitData.root";
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"RECREATE");

    plot1->SetMarkerColor(kRed);
    plot2->SetMarkerColor(kRed);

    relativeCSLogGraph->GetXaxis()->SetTitle("Energy (MeV)");
    relativeCSLogGraph->GetXaxis()->CenterTitle();

    relativeCSLogGraph->GetYaxis()->SetTitle("#frac{#sigma^{124}_{tot}-#sigma^{112}_{tot}}{#sigma^{124}_{tot}+#sigma^{112}_{tot}}");
    relativeCSLogGraph->GetYaxis()->SetTitleSize(0.05);
    relativeCSLogGraph->GetYaxis()->SetTitleOffset(1.4);

    relativeCSLogGraph->GetYaxis()->CenterTitle();

    relativeCSLogGraph->SetName("relativeCS");

    //relativeCS->SetLineColor(kRed);
    relativeCSLog->SetLineColor(kRed);

    //relativeCS->SetTitle("Relative #sigma_{tot} difference: #frac{^{124}Sn-^{112}Sn}{^{124}Sn-^{112}Sn}");

    //relativeCS->Draw();

    relativeCSLogGraph->GetXaxis()->SetRangeUser(3,300);
    relativeCSLogGraph->GetYaxis()->SetRangeUser(0,0.06);
    relativeCSLogGraph->SetLineColor(kPink);

    relativeCSLogGraph->Draw("ap");
    //relativeCS_W->Draw("same");
    //plot2->Draw("same");

    plot1->GetXaxis()->SetRangeUser(3,300);

    plot1->Draw("same");

    //relativeCS->SetTitle("Relative");
    //relativeCS->GetTitle->CenterTitle();
    //relativeCSLog->SetTitle("Relative");
    //relativeCSLog->Clone("
    //relativeCSLog->GetTitle()->CenterTitle();


    // plot1->GetXaxis()->SetRangeUser(3,300);
    //relativeCS->GetXaxis()->SetRangeUser(3,300);

    TLine *line = new TLine(3,0.0345,300,0.0345);
    line->SetLineStyle(7);
    line->SetLineWidth(3);
    line->SetLineColor(kBlack);
    line->Draw();

    TLegend *legend = new TLegend(0.2,0.7,0.4,0.9);
    //legend->SetHeader("data");
    legend->AddEntry(relativeCSLogGraph,"experimental (preliminary)","p");
    legend->AddEntry(plot1,"DOM calculation","l");
    legend->AddEntry(line,"expected from size scaling","l");
    legend->Draw();


    /*    TMultiGraph *allPlots = new TMultiGraph();
    for(int i=0; i<allData.size(); i++)
    {
        allPlots->Add(((TGraphErrors*)allData[i].getPlot()),"p");
    }
    allPlots->Draw("a");

    allPlots->GetXaxis()->SetTitle("MeV");
    allPlots->GetXaxis()->CenterTitle();

    allPlots->GetYaxis()->SetTitle("Barns");
    allPlots->GetYaxis()->CenterTitle();

    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->SetHeader("Literature data");
    for(int i=0; i<allData.size(); i++)
    {
        legend->AddEntry(((TGraphErrors*)allData[i].getPlot()),allData[i].getReference(),"lep");
    }
    legend->Draw();

    allPlots->Write();
*/

    // clean up
    fileOut->Write();
    fileOut->Close();
}
