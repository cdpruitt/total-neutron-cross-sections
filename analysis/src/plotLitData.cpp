
/*****************************************************************************/
// Style section: edit to change style of plotted datasets
const std::string colors[4] = {"kRed","kRed","kRed","kRed"};
const std::string markers[4] = {"kPlus","kPlus","kPlus","kPlus"};
/*****************************************************************************/

// ROOT macro for formatting and combining plots

void plotLitData(string formatFileName, string outFileName)
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
    ifstream inFile(formatFileName);
    if(!inFile.is_open())
    {
        cout << "Failed to open " << formatFileName << endl;
        exit(1);
    }

    int numberOfFiles;
    vector<string> fileNames;

    int numberOfPlots;
    vector<string> plotTypes;
    vector<string> plotNames;
    vector<string> plotTitles;
    vector<string> plotLine;
    vector<string> plotPoints;
    vector<int> plotLineColor;
    vector<int> plotPointColor;
    vector<double> plotLineWidth;
    vector<double> plotPointSize;
    vector<string> plotXAxisTitles;
    vector<string> plotYAxisTitles;
    vector<string> XAxisTitles;
    vector<string> YAxisTitles;

    vector<TGraphErrors*> graphs;
    vector<TH1I*> histos;

    string dummy;
    string dummy2;
    inFile >> dummy >> dummy2;
    {
        if(dummy=="numberOfFiles")
        {
            numberOfFiles = stoi(dummy2);
        }
    }

    for(int i=0; i<numberOfFiles; i++)
    {
        getline(inFile,dummy); // discard blank line
        inFile >> dummy >> dummy2;
        if(dummy=="fileName")
        {
            fileNames.push_back(dummy2);
        }

        inFile >> dummy >> dummy2;

        if(dummy=="numberOfPlots")
        {
            numberOfPlots = stoi(dummy2);
        }

        while(inFile >> dummy >> dummy2 && dummy!="b")
        {
            //cout << dummy << " " << dummy2 << endl;
            if(dummy=="plotType")
            {
                plotTypes.push_back(dummy2);
            }

            else if(dummy=="plotName")
            {
                plotNames.push_back(dummy2);
            }

            else if(dummy=="plotTitle")
            {
                plotTitles.push_back(dummy2);
            }

            else if(dummy=="plotLine")
            {
                plotLine.push_back(dummy2);
            }

            else if(dummy=="plotPoints")
            {
                plotPoints.push_back(dummy2);
            }

            else if(dummy=="plotLineColor")
            {
                plotLineColor.push_back(stoi(dummy2));
            }

            else if(dummy=="plotPointColor")
            {
                plotPointColor.push_back(stoi(dummy2));
            }

            else if(dummy=="plotLineWidth")
            {
                plotLineWidth.push_back(stod(dummy2));
            }

            else if(dummy=="plotPointSize")
            {
                plotPointSize.push_back(stod(dummy2));
            }

            else if(dummy=="plotXAxisTitle")
            {
                plotXAxisTitles.push_back(dummy2);
            }
            else if(dummy=="plotYAxisTitle")
            {
                plotYAxisTitles.push_back(dummy2);
            }

            else if(dummy=="XAxisTitle")
            {
                XAxisTitles.push_back(dummy2);
            }
            else if(dummy=="YAxisTitle")
            {
                YAxisTitles.push_back(dummy2);
            }
        }
    }

    vector<TFile*> files;
    for(int i=0; i<numberOfFiles; i++)
    {
        files.push_back(new TFile(fileNames[i].c_str(),"READ"));
        for(int j=0; j<plotTypes.size(); j++)
        {
            if(plotTypes[j]=="graph")
            {
                graphs.push_back((TGraphErrors*)files[i]->Get(plotNames[j].c_str()));
                if(!graphs.back())
                {
                    cerr << "Error: failed to find " << plotNames[j] <<
                        " in " << fileNames[i] << ". Exiting..." << endl;
                    exit(1);
                }
            }

            if(plotTypes[j]=="histo")
            {
                histos.push_back((TH1I*)files[i]->Get(plotNames[j].c_str()));
                if(!histos.back())
                {
                    cerr << "Error: failed to find " << plotNames[j] <<
                        " in " << fileNames[i] << ". Exiting..." << endl;
                    exit(1);
                }
            }
        }
    }

    //relativeCSLogGraph->SetMarkerColor(kRed);
    //relativeCSLogGraph->SetMarkerSize(1);

    // create output file
    TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

    for(int i=0; i<histos.size(); i++)
    {
        if(plotTitles[i]=="true")
        {
            histos[i]->SetTitle(plotTitles[i].c_str());
        }

        if(plotXAxisTitles[i]=="true")
        {
            histos[i]->GetXaxis()->SetTitle(XAxisTitles[i].c_str());
            histos[i]->GetXaxis()->CenterTitle();
            histos[i]->GetXaxis()->SetTitleSize(0.05);
            histos[i]->GetXaxis()->SetTitleOffset(1.5);
        }

        if(plotYAxisTitles[i]=="true")
        {
            histos[i]->GetYaxis()->SetTitle(YAxisTitles[i].c_str());
            histos[i]->GetYaxis()->CenterTitle();
            histos[i]->GetYaxis()->SetTitleSize(0.04);
            histos[i]->GetYaxis()->SetTitleOffset(3);
        }

        if(plotLine[i]=="true")
        {
            histos[i]->SetLineColor(plotLineColor[i]);
            histos[i]->SetLineWidth(plotLineWidth[i]);
        }

        if(plotPoints[i]=="true")
        {
            histos[i]->SetMarkerColor(plotPointColor[i]);
            histos[i]->SetMarkerSize(plotPointSize[i]);
        }

        else
        {
            histos[i]->SetMarkerSize(0);
        }

        //histos[i]->Scale(1/1000.); // for deadtime plotting
        histos[i]->Write();
    }

    for(int i=0; i<graphs.size(); i++)
    {
        if(plotTitles[i]=="true")
        {
            graphs[i]->SetTitle(plotTitles[i].c_str());
        }

        if(plotXAxisTitles[i]=="true")
        {
            graphs[i]->GetXaxis()->SetTitle(XAxisTitles[i].c_str());
            graphs[i]->GetXaxis()->CenterTitle();
            graphs[i]->GetXaxis()->SetTitleSize(0.05);
            graphs[i]->GetXaxis()->SetTitleOffset(1.5);
        }

        if(plotYAxisTitles[i]=="true")
        {
            graphs[i]->GetYaxis()->SetTitle(YAxisTitles[i].c_str());
            graphs[i]->GetYaxis()->CenterTitle();
            graphs[i]->GetYaxis()->SetTitleSize(0.05);
            graphs[i]->GetYaxis()->SetTitleOffset(1.2);
        }

        if(plotLine[i]=="true")
        {
            graphs[i]->SetLineColor(plotLineColor[i]);
            graphs[i]->SetLineWidth(plotLineWidth[i]);
        }

        if(plotPoints[i]=="true")
        {
            graphs[i]->SetMarkerColor(plotPointColor[i]);
            graphs[i]->SetMarkerSize(plotPointSize[i]);
        }

        else
        {
            graphs[i]->SetMarkerSize(0);
        }

        graphs[i]->Write();
    }

    //relativeCSLogGraph->GetXaxis()->SetTitle("Energy (MeV)");
    //relativeCSLogGraph->GetXaxis()->CenterTitle();

    //relativeCSLogGraph->GetYaxis()->SetTitle("#frac{#sigma^{124}_{tot}-#sigma^{112}_{tot}}{#sigma^{124}_{tot}+#sigma^{112}_{tot}}");
    //relativeCSLogGraph->GetYaxis()->SetTitleSize(0.05);
    //relativeCSLogGraph->GetYaxis()->SetTitleOffset(1.4);

    //relativeCSLogGraph->GetYaxis()->CenterTitle();

    //relativeCSLogGraph->SetName("relativeCS");

    //relativeCS->SetLineColor(kRed);

    //relativeCS->SetTitle("Relative #sigma_{tot} difference: #frac{^{124}Sn-^{112}Sn}{^{124}Sn-^{112}Sn}");

    //relativeCS->Draw();

    //relativeCSLogGraph->GetXaxis()->SetRangeUser(3,300);
    //relativeCSLogGraph->GetYaxis()->SetRangeUser(0,0.06);

    //relativeCSLogGraph->Draw("ap");
    //relativeCS_W->Draw("same");
    //plot2->Draw("same");

    //TLine *line = new TLine(3,0.0345,300,0.0345);
    //line->SetLineStyle(7);
    //line->SetLineWidth(3);
    //line->SetLineColor(kBlack);
    //line->Draw();

    //TLegend *legend = new TLegend(0.2,0.7,0.4,0.9);
    //legend->SetHeader("data");
    //legend->AddEntry(relativeCSLogGraph,"experimental (preliminary)","p");
    //legend->AddEntry(plot1,"DOM calculation","l");
    //legend->AddEntry(line,"expected from size scaling","l");
    //legend->Draw();


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
    outFile->Write();
    outFile->Close();
}
