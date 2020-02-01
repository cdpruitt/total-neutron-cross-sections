{
    string fileName = "/data2/analysis/relative.root";
    string ramsauerFileName = "../../theory/ramsauer.root";

    TFile* file = new TFile(fileName.c_str(),"READ");
    TFile* ramsauerFile = new TFile(ramsauerFileName.c_str(), "READ");
    
    string relGraphName = "Sn124Sn112, percent";
    string relGraphSEName = "Sn124Sn112SysErrors, percent";

    string SARelDiffGraphThirdName = "RelDiff124_112Third";
    string SARelDiffGraphSixthName = "RelDiff124_112Sixth";

    string RamsauerRelDiffGraphName = "RelDiffRamsauer124_112";
    string RamsauerRelDiffGraphSixthName = "RelDiffRamsauerSixth124_112";
       
    TGraphAsymmErrors* relGraph = (TGraphAsymmErrors*)file->Get(relGraphName.c_str());
    TGraphAsymmErrors* relGraphSE = (TGraphAsymmErrors*)file->Get(relGraphSEName.c_str());

    TGraph* SARelDiffGraphThird = (TGraph*)ramsauerFile->Get(SARelDiffGraphThirdName.c_str());
    TGraph* SARelDiffGraphSixth = (TGraph*)ramsauerFile->Get(SARelDiffGraphSixthName.c_str());
    TGraph* RamsauerRelDiffGraph = (TGraph*)ramsauerFile->Get(RamsauerRelDiffGraphName.c_str());
    TGraph* RamsauerRelDiffGraphSixth = (TGraph*)ramsauerFile->Get(RamsauerRelDiffGraphSixthName.c_str());

    if(!ramsauerFile)
    {
        cerr << "Error: couldn't open " << ramsauerFileName << endl;
        exit(1);
    }

    if(!relGraph)
    {
        cerr << "Error: failed to find " << relGraphName << endl;
        exit(1);
    }

    if(!SARelDiffGraphThird || !SARelDiffGraphSixth)
    {
        cerr << "Error: failed to find " << SARelDiffGraphThirdName << endl;
        exit(1);
    }

    if(!RamsauerRelDiffGraph)
    {
        cerr << "Error: failed to find " << RamsauerRelDiffGraphName << endl;
        exit(1);
    }

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

    // create DOM relative difference file
    vector<double> energies;
    vector<double> CSs;
    vector<double> energyErrors;
    vector<double> CSErrors;

    string DOMRelDiffFileName = "/home/cdpruitt/gDOM/output/sn124_sn112.txt";
    ifstream DOMRelDiffFile(DOMRelDiffFileName.c_str());
    string line;
    while(getline(DOMRelDiffFile,line))
    {
        // parse line into tokens
        vector<string> tokens;
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
                istream_iterator<string>(),
                back_inserter(tokens));

        // skip empty lines
        if(!tokens.size())
        {
            continue;
        }

        energies.push_back(stod(tokens[0]));
        energyErrors.push_back(0);
        CSs.push_back(100*stod(tokens[1]));
        CSErrors.push_back(0);
        //CSErrors.push_back(100*stod(tokens[2]));
    }
    
    TGraphAsymmErrors* DOMGraph = new TGraphAsymmErrors(energies.size(),
                                      &energies[0],
                                      &CSs[0],
                                      &energyErrors[0],
                                      &energyErrors[0],
                                      &CSErrors[0],
                                      &CSErrors[0]);

    // Set graph point and line characteristics
    relGraph->SetLineColor(kRed);
    relGraph->SetLineWidth(3);
    relGraph->SetLineStyle(0);
    relGraph->SetMarkerColor(kRed);
    relGraph->SetFillColor(kRed);
    relGraph->SetFillStyle(3002);

    SARelDiffGraphThird->SetLineStyle(9);
    SARelDiffGraphThird->SetLineWidth(3);
    SARelDiffGraphThird->SetLineColor(kGray+2);

    SARelDiffGraphSixth->SetLineStyle(7);
    SARelDiffGraphSixth->SetLineWidth(3);
    SARelDiffGraphSixth->SetLineColor(kGray+2);

    DOMGraph->SetLineStyle(3);
    DOMGraph->SetLineWidth(5);
    DOMGraph->SetLineColor(kBlack);
    //DOMGraph->SetMarkerColor(kBlue-8);
    //DOMGraph->SetFillColor(kBlue-8);
    //DOMGraph->SetFillStyle(3002);

    RamsauerRelDiffGraph->SetLineStyle(7);
    RamsauerRelDiffGraph->SetLineWidth(3);
    RamsauerRelDiffGraph->SetLineColor(kGray+2);

    RamsauerRelDiffGraphSixth->SetLineStyle(7);
    RamsauerRelDiffGraphSixth->SetLineWidth(3);
    RamsauerRelDiffGraphSixth->SetLineColor(kGray+2);

    relGraphSE->SetFillColor(kBlue);
    relGraphSE->SetFillStyle(3002);

    // Pad dimensions and margins
    gPad->SetPad(0.005, 0.995, 0.995, 0.005);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicky(2);

    TMultiGraph* mg = new TMultiGraph();

    mg->Add(relGraph,"3l");
    mg->Add(relGraphSE, "3");
    mg->Add(SARelDiffGraphThird, "l");
    mg->Add(SARelDiffGraphSixth, "l");
    mg->Add(DOMGraph, "l");
    //mg->Add(RamsauerRelDiffGraph, "l");
    //mg->Add(RamsauerRelDiffGraphSixth, "l");

    mg->Draw("al");

    // X-axis parameters
    mg->GetXaxis()->SetTitle("Energy (MeV)");
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetTitleFont(2);
    mg->GetXaxis()->SetTitleOffset(1.4);
    mg->GetXaxis()->CenterTitle();

    mg->GetXaxis()->SetLabelOffset(0.01);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetXaxis()->SetLabelFont(2);

    mg->GetXaxis()->SetNdivisions(10);
    mg->GetXaxis()->SetTickLength(0.03);

    // Y-axis parameters
    mg->GetYaxis()->SetTitle("(#frac{#sigma_{124} - #sigma_{112}}{#sigma_{124} + #sigma_{112}})");
    mg->GetYaxis()->SetTitleSize(0.06);
    mg->GetYaxis()->SetTitleFont(2);
    mg->GetYaxis()->SetTitleOffset(1.0);
    mg->GetYaxis()->CenterTitle();

    mg->GetYaxis()->SetLabelOffset(0.01);
    mg->GetYaxis()->SetLabelSize(0.05);

    mg->GetYaxis()->SetLabelFont(2);
    mg->GetYaxis()->SetNdivisions(10);
    mg->GetYaxis()->SetTickLength(0.02);

    gPad->SetLogx(1);
    
    mg->GetYaxis()->SetRangeUser(0.0,5.001);
    mg->GetXaxis()->SetLimits(5,600);

    // Define legend format and contents
    TLegend *legend = new TLegend(0.65, 0.20, 0.95, 0.34);
    legend->SetNColumns(2);
    legend->AddEntry(relGraph,"Exp data, sys + stat   ","f");
    legend->AddEntry(SARelDiffGraphThird,"SAS, r #alpha A^{1/3} ","l");
    legend->AddEntry(relGraphSE,"Exp data, sys only   ","f");
    legend->AddEntry(SARelDiffGraphSixth,"SAS, r #alpha A^{1/6} ","l");
    legend->AddEntry(DOMGraph,"DOM fit","l");
    //legend->AddEntry(RamsauerRelDiffGraph,"Ramsauer","l");
    legend->Draw();

    //TLatex latex;
    //latex.SetNDC();
    //latex.SetTextSize(0.05);
    //latex.SetTextAlign(13); // align at top
    //latex.DrawLatex(0.19,0.82,"A^{#frac{1}{3}}");
    //latex.DrawLatex(0.77,0.40,"A^{#frac{1}{6}}");

    file->Close();
}
