
/*****************************************************************************/
// Style section: edit to change style of plotted datasets
const std::string colors[4] = {"kRed","kRed","kRed","kRed"};
const std::string markers[4] = {"kPlus","kPlus","kPlus","kPlus"};
/*****************************************************************************/

// ROOT macro for formatting and combining plots

void plotRelativeCS()
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
    Sty->SetPadTopMargin(.18);
    Sty->SetPadLeftMargin(.18);
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

    // open data files
    std::stringstream expFileInName;
    expFileInName << "/data3/analysis/total.root";
    std::stringstream theoryFileInName;
    theoryFileInName << "/data2/analysis/literatureData.root";

    TFile* expFile = new TFile(expFileInName.str().c_str(),"READ");
    TGraphErrors* plot1 = (TGraphErrors*)gDirectory->Get("relativeCS");

    TFile* theoryFile = new TFile(theoryFileInName.str().c_str(),"READ");
    TGraphErrors* plot2 = (TGraphErrors*)gDirectory->Get("Sn124-Sn112DivSn124+Sn112");

    if(!plot1 || !plot2)
    {
        cout << "Failed to open one of the graphs for combining." << endl;
        exit(1);
    }

    // create output file
    std::stringstream fileOutName;
    fileOutName << "relativeGraph.root";
    TFile* fileOut = new TFile(fileOutName.str().c_str(),"RECREATE");

    plot1->SetMarkerColor(kBlue);
    plot2->SetMarkerColor(kRed);

    plot1->GetXaxis()->SetTitle("Energy (MeV)");
    plot1->GetXaxis()->CenterTitle();

    plot1->GetYaxis()->SetTitleSize(0.05);
    plot1->GetYaxis()->SetTitleOffset(1.6);

    plot1->GetYaxis()->CenterTitle();

    plot1->SetLineColor(kBlue);
    plot2->SetLineColor(kRed);

    plot1->GetXaxis()->SetRangeUser(6,300);
    plot2->GetXaxis()->SetRangeUser(6,300);
    //relativeCSLogGraph->GetYaxis()->SetRangeUser(0,0.06);

    plot1->Draw("alp");
    plot2->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(13); // align at top
    latex.DrawLatex(0.47,0.95,"#frac{#sigma^{^{124}Sn}-#sigma^{^{112}Sn}}{#sigma^{^{124}Sn}+#sigma^{^{112}Sn}}");;

    //plot1->SetTitle("#frac{#sigma^{124}_{tot}-#sigma^{112}_{tot}}{#sigma^{124}_{tot}+#sigma^{112}_{tot}}");

    gPad->SetLogx(1);

    //relativeCS->SetTitle("Relative");
    //relativeCS->GetTitle->CenterTitle();
    //relativeCSLog->SetTitle("Relative");
    //relativeCSLog->Clone("
    //relativeCSLog->GetTitle()->CenterTitle();

    // plot1->GetXaxis()->SetRangeUser(3,300);
    //relativeCS->GetXaxis()->SetRangeUser(3,300);

    //TLine *line = new TLine(3,0.0345,300,0.0345);
    //line->SetLineStyle(7);
    //line->SetLineWidth(3);
    //line->SetLineColor(kBlack);
    //line->Draw();

    TLegend *legend = new TLegend(0.6,0.2,0.9,0.35);
    //legend->SetHeader("data");
    legend->AddEntry(plot1,"experimental (preliminary)","lp");
    legend->AddEntry(plot2,"DOM calculation","l");
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
