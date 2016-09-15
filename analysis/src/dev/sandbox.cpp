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



