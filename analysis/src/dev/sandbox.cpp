for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        // test to see if plot vectors have already been linked to existing plots (i.e.,
        // fillHistos and fillAdvancedHistos were run).
        if(skippedHistoFilling)
        {
            string temp = positionNames[i] + "TOF";
            plots.TOFHistos.push_back((TH1I*)gDirectory->Get(temp.c_str()));

            temp = positionNames[i] + "Energy";
            plots.energyHistos.push_back((TH1I*)gDirectory->Get(temp.c_str()));
        }

        // prepare deadtime-corrected plots
        string temp = positionNames[i] + "TOFCorrected";
        plots.TOFHistosCorrected.push_back((TH1I*)plots.TOFHistos[i]->Clone(temp.c_str()));

        temp = positionNames[i] + "EnergyCorrected";
        plots.energyHistosCorrected.push_back((TH1I*)plots.energyHistos[i]->Clone(temp.c_str()));
        plots.energyHistosCorrected.back()->Reset();
    }

DEADTIME FROM HISTOS

/*        cout << "microsPerTarget[i] = " << microsPerTarget[i] << endl;

        //plots.energyHistosCorrected[i]->Sumw2();

        // loop through all bins
        for(int j=0; j<TOFCorrectedHistos[i]->GetNbinsX(); j++)
        {
            if(microsPerTarget[i] > 0)
            {
                eventsPerBinPerMicro[i].push_back(TOFCorrectedHistos[i]->GetBinContent(j)/(double)microsPerTarget[i]);
            }

            else
            {
                eventsPerBinPerMicro[i].push_back(0);
            }

            deadtimeFraction[i].push_back(0);
        }

        // find the fraction of the time that the detector is dead for each bin in the micropulse
        for(int j=0; (size_t)j<eventsPerBinPerMicro[i].size(); j++)
        {
            *//*int k = j-(FULL_DEADTIME+PARTIAL_DEADTIME)*(TOF_BINS/TOF_RANGE);
            while(k<j)
            {
                if(k>=0)
                {
                    // partially-dead region
                    if((j-k)>=(FULL_DEADTIME)*(TOF_BINS/TOF_RANGE))
                    {
                        deadtimeFraction[i][j] +=
                            eventsPerBinPerMicro[i][k]*
                            (k+(FULL_DEADTIME+PARTIAL_DEADTIME)*(TOF_BINS/TOF_RANGE)-j)/(PARTIAL_DEADTIME*(TOF_BINS/TOF_RANGE));
                    }

                    else
                    {
                        deadtimeFraction[i][j] += eventsPerBinPerMicro[i][k];
                    } 
                }
                k++;
            }*/
            
            /*deadtimeFraction[i][j] = deadtimeFraction[i][j-1]+eventsPerBinPerMicro[i][j];
            if(j>(TOF_BINS/TOF_RANGE)*DEADTIME_PERIOD)
            {
                deadtimeFraction[i][j] -= eventsPerBinPerMicro[i][j-(TOF_BINS/TOF_RANGE)*DEADTIME_PERIOD];
            }
        }
        */


PLOT DEADTIME FROM HISTOS
// Plot the calculated dead time fractions (for debugging)
        string temp;
        temp = "deadtime" + targetNames[i];
        for(int j=0; (size_t)j<deadtimeFraction[i].size(); j++)
        {
            target[i].deadtimeHistos.back()->SetBinContent(j,1000000*deadtimeFraction[i][j]);
        }
        plots.deadtimeHistos.back()->Write();




PER MACRO HISTOS FROM HISTOS

// Find number of macropulses for each target to use in error calculation
    gDirectory->cd("/");
    gDirectory->GetDirectory("targetChanger")->cd();

vector<long> tarCounts;

    for(int i=0; i<NUMBER_OF_TARGETS; i++)
    {
        tarCounts.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
    }
    
    /*
    // Scale perMacroHistos histogram by number of macropulses in that target
    for(int i=0; i<6; i++)
    {
    if(tarCounts[i]==0)
    {
    continue;
    }
    perMacroHistos[i]->Scale(pow(10,3)/tarCounts[i]);
    }
    */

SCALED CS HISTOS FROM HISTOS

    /*plots.CSGraphsScaledToLit.push_back(new TGraphErrors(output.energyBins.size(),&output.energyBins[0],&(output.crossSectionsScaledToLit[3]->at(0)),&xError[3]->at(0),&output.crossSectionsError[3]->at(0)));
    string temp = targetNames[3] + "Scaled";
    plots.CSGraphsScaledToLit.back()->SetNameTitle(temp.c_str(),temp.c_str());
    plots.CSGraphsScaledToLit.back()->Write();

    temp = "";
    temp = targetNames[5] + "Scaled";
    plots.CSGraphsScaledToLit.push_back(new TGraphErrors(output.energyBins.size(),&output.energyBins[0],&(output.crossSectionsScaledToLit[5]->at(0)),&xError[5]->at(0),&output.crossSectionsError[5]->at(0)));
    plots.CSGraphsScaledToLit.back()->SetNameTitle(temp.c_str(),temp.c_str());
    plots.CSGraphsScaledToLit.back()->Write();
    */

// read literature data for natural Sn
    TFile *litData = new TFile("/data2/analysis/literatureData.root","READ");
    TGraphErrors *SnNatLitData = (TGraphErrors*)litData->Get("Natural Sn (n,tot)");

    // scale Sn112 and Sn124 cross sections using the literature value for natural Sn
    for(int i=0; i<output.energyBins.size(); i++)
    {
        // avoid "divide by 0"
        if(output.crossSections[4]->at(i) != 0)
        {
            double litValue = SnNatLitData->Eval(output.energyBins[i]);
            output.crossSectionsScaledToLit[3]->push_back(output.crossSections[3]->at(i)*(litValue/output.crossSections[4]->at(i)));
            output.crossSectionsScaledToLit[5]->push_back(output.crossSections[5]->at(i)*(litValue/output.crossSections[4]->at(i)));
        }
    }

    litData->Close();


    FROM HISTOS

if(TOF_BINS%NUMBER_ENERGY_BINS!=0)
    {
        cout << "Error: number of TOF bins must be a multiple of the number of energy bins." << endl;
        exit(1);
    }

STATISTICS FROM SUMALL

/*void getStatisticalError(TH1D* Sn112CSTotal, TH1D* Sn124CSTotal)
{
    // Get monCounts, raw target counts from _histos.root files
    // Calculate statistical error from histogram counts

    int limit = 25;

    ifstream runList("runsToSort.txt");

    string runDir;
    string outPath;

    while (runList >> runDir)
    {
        if (stoi(runDir)<6)
        {
            outPath = "/data3/analysis/run";
        }

        else if (stoi(runDir)>128 && stoi(runDir)<160)
        {
            outPath = "/data2/analysis/run";
        }

        else if (stoi(runDir)>159 && stoi(runDir)<178)
        {
            outPath = "/data3/analysis/run";
        }

        else
        {
            cout << "Run directory outside bounds (runs 128-177) - check run number" << endl;
            exit(1);
        }

        // keep track of which order the targets are in, based on which run number we're
        // sorting
        vector<int> order;

        int pos112 = find(order.begin(), order.end(), 3) - order.begin();
        int pos124 = find(order.begin(), order.end(), 5) - order.begin();

        for(int segment = 0; segment<=limit; segment++)
        {
            // We need to form the proper name for the sub-run we want to open:
            if(segment < 10)
            {
                infile =  new TFile(Form("%s%s/run%s-000%i_histos.root",outPath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            }

            else if(segment < 100)
            {
                infile =  new TFile(Form("%s%s/run%s-00%i_histos.root",outPath.c_str(),runDir.c_str(),runDir.c_str(),segment));
            }

            else
            {
                // There's some error in the sub-run file number; write outfile and exit
                cout << "Segment number too large!" << endl;
                exit(1);
            }

            // Attempt to open the sub-run to access its histograms
            if(!infile->IsOpen())
            {
                cout << "Can't open root file run " << runDir << " segment = " << segment << endl;

                break;
            }

            gDirectory->cd("/");
            gDirectory->GetDirectory("monitor")->cd();

            monCountsBlank += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2);
            monCountsSn112 += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2+pos112);
            monCountsSn124 += ((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(2+pos124);

            gDirectory->cd("/");
            gDirectory->GetDirectory("detS")->cd();

            TH1I* blank = ((TH1I*)gDirectory->Get("blankLog"));

            string name112 = "target" + to_string(pos112) + "Log";
            string name124 = "target" + to_string(pos124) + "Log";

            TH1I* Sn112 = ((TH1I*)gDirectory->Get(name112.c_str()));
            TH1I* Sn124 = ((TH1I*)gDirectory->Get(name124.c_str()));

            for(int i=0; i<blank->GetNbinsX(); i++)
            {
                detCountsBlank[i] += (blank->GetBinContent(i));
                detCountsSn112[i] += (Sn112->GetBinContent(i));
                detCountsSn124[i] += (Sn124->GetBinContent(i));
            }

            // Close the sub-run input files
            infile->Close();

            // End of loop - move to next sub-run
        }
    }

    *//*monCountsBlankError = pow(monCountsBlank,0.5);
    monCountsSn112Error = pow(monCountsSn112,0.5);
    monCountsSn124Error = pow(monCountsSn124,0.5);
    *//*

    cout << "monCountsBlank = " <<  monCountsBlank << endl;
    cout << "monCountsSn112 = " <<  monCountsSn112 << endl;
    cout << "monCountsSn124 = " <<  monCountsSn124 << endl;

    cout << "detCountsBlank[10] = " << detCountsBlank[10] << endl;
    cout << "detCountsSn112[10] = " << detCountsSn112[10] << endl;
    cout << "detCountsSn124[10] = " << detCountsSn124[10] << endl;


    for(int i=0; (size_t)i<totalSn112StatError.size(); i++)
    {
        totalSn112StatError[i] = pow(
                ((double)monCountsBlank/(pow(monCountsBlank,2))) +
                ((double)monCountsSn112/(pow(monCountsSn112,2))) +
                ((double)detCountsBlank[i]/(pow(detCountsBlank[i],2))) +
                ((double)detCountsSn112[i]/(pow(detCountsSn112[i],2)))
                ,0.5)*//*/exp(abs(Sn112CSTotal->GetBinContent(i)))*//*;

        totalSn124StatError[i] = pow(
                ((double)monCountsBlank/(pow(monCountsBlank,2))) +
                ((double)monCountsSn124/(pow(monCountsSn124,2))) +
                ((double)detCountsBlank[i]/(pow(detCountsBlank[i],2))) +
                ((double)detCountsSn124[i]/(pow(detCountsSn124[i],2)))
                ,0.5)*//*/exp(abs(Sn124CSTotal->GetBinContent(i)))*//*;

        *//*totalRelStatError[i] = pow(
                pow(totalSn112Error[i],2) +
                pow(totalSn124Error[i],2)
                ,0.5);
                *//*
    }

    cout << totalSn112StatError[10] << endl;

        *//*detCountsBlankpush_back(pow(detCountsBlank[i],0.5));
        detCountsSn112push_back(pow(detCountsSn112[i],0.5));
        detCountsSn124push_back(pow(detCountsSn124[i],0.5));
        *//*
}

void getSystemicError()
{
    // physical target data, listed in order:
    // {blank, Sn112, Natural Sn, Sn124, short carbon, long carbon} 

    // lengths of each target:
    double targetlength[6] = {0,1.37,2.74,1.365,1.370,1.370}; //cm
    double targetlengthUncertainty[6] = {0,0.01,0.01,0.005,0.005,0.005}; //cm

    // molar mass of each target:
    double targetMolMass[6] = {0,12.01,12.01,112,118.7,124}; //g/mol
    double targetMolMassUncertainty[6] = {0,0.01,0.01,0.01,0.01,0.01}; //g/mol

    // density of each target:
    double targetdensity[6] = {0,2.2,2.2,6.89,7.31,7.63}; //g/cm^3
    double targetdensityUncertainty[6] = {0,0.1,0.1,0.001,0.001,0.001}; //g/cm^3

    totalSn112SysError = pow(
            pow(targetlengthUncertainty[3]/targetlength[3],2) +
            pow(targetMolMassUncertainty[3]/targetMolMass[3],2) +
            pow(targetdensityUncertainty[3]/targetdensity[3],2)
            ,0.5);

    totalSn124SysError = pow(
            pow(targetlengthUncertainty[5]/targetlength[5],2) +
            pow(targetMolMassUncertainty[5]/targetMolMass[5],2) +
            pow(targetdensityUncertainty[5]/targetdensity[5],2)
            ,0.5);
}

void getTotalError()
{
    for(int i=0; (size_t)i<totalSn112Error.size(); i++)
    {
        totalSn112Error[i] = pow(
                pow(totalSn112StatError[i],2) +
                pow(totalSn112SysError,2)
                ,0.5);
        totalSn124Error[i] = pow(
                pow(totalSn124StatError[i],2) +
                pow(totalSn124SysError,2)
                ,0.5);
        totalRelError[i] = pow(
                pow(totalSn112Error[i],2) +
                pow(totalSn124Error[i],2)
                ,0.5);
    }

    cout << totalSn112Error[10] << endl;
    cout << totalSn124Error[10] << endl;
    cout << totalRelError[10] << endl;
}
*/



