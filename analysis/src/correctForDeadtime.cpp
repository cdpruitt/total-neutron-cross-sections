// experimentally-determined digitizer deadtime
const int DEADTIME_PERIOD = 150; // in ns
const int DEADTIME_TRANSITION_PERIOD = 15; // in ns

double applyDeadtimeCorrection(TH1D*& rawTOF, vector<double> deadtimesPerBin)
{
    // produce histo showing deadtime for each bin
    //string name = correctedTOF->GetName();
    //name = name + "deadtimeH";
    //TH1D* deadtimeH = (TH1D*)correctedTOF->Clone(name.c_str());

    double sumOfDeadtimes = 0; // for computing average deadtime per bin

    string name = rawTOF->GetName();
    name += "Corrected";
    TH1D* correctedTOFHisto = (TH1D*)rawTOF->Clone(name.c_str());

    for(int i=0; (size_t)i<deadtimesPerBin.size(); i++)
    {
        if(deadtimesPerBin[i]>=1)
        {
            cerr << "Error: attempted to correct for deadtime, but encountered deadtime >100% (deadtime was "
                << deadtimesPerBin[i] << "). Exiting." << endl;
            exit(1);
        }

        correctedTOFHisto->SetBinContent(i+1,rawTOF->GetBinContent(i+1)/(1-deadtimesPerBin[i]));
        sumOfDeadtimes += deadtimesPerBin[i];
    }

    correctedTOFHisto->Write();
    
    return sumOfDeadtimes/deadtimesPerBin.size();
}

vector<double> generateDeadtimeCorrection(TH1D* tof, long totalNumberOfMicros)
{
    vector<double> eventsPerMicroPerBin(config.plotConfig.TOF_BINS, 0);
    vector<double> deadtimesPerBin(config.plotConfig.TOF_BINS,0);

    string name = tof->GetName();
    string eventsPerMicroName = name + "EventsPerMicro";
    TH1D* eventsPerMicroPerBinH = new TH1D(eventsPerMicroName.c_str(), eventsPerMicroName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND);

    string deadtimeName = name + "Deadtime";
    TH1D* deadtimeHisto = new TH1D(deadtimeName.c_str(), deadtimeName.c_str(),config.plotConfig.TOF_BINS,config.plotConfig.TOF_LOWER_BOUND,config.plotConfig.TOF_UPPER_BOUND);

    if(totalNumberOfMicros==0)
    {
        cerr << "Error: tried to create deadtimes, but totalNumberOfMicros = 0" << endl;
        return deadtimesPerBin;
    }

    const int deadtimeBins = ((double)config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*DEADTIME_PERIOD;
    const int deadtimeTransitionBins = ((double)config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*DEADTIME_TRANSITION_PERIOD;

    for(int j=0; j<eventsPerMicroPerBin.size(); j++)
    {
        eventsPerMicroPerBin[j] = tof->GetBinContent(j+1)/((double)totalNumberOfMicros);
        eventsPerMicroPerBinH->SetBinContent(j+1,eventsPerMicroPerBin[j]);
    }

    eventsPerMicroPerBinH->Write();

    // find the fraction of the time that the detector is dead for each bin in the micropulse
    // set deadtime fraction base case

    // use deadtime base case to calculate deadtime for remaining bins

    for(int j=0; j<config.plotConfig.TOF_BINS; j++)
    {
        for(int k=j-(deadtimeBins+deadtimeTransitionBins); k<j; k++)
        {
            if(k<(j-deadtimeBins))
            {
                if(k<0)
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k+config.plotConfig.TOF_BINS]*(double)(1-deadtimesPerBin[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }

                else
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimesPerBin[j])*((deadtimeBins+deadtimeTransitionBins-(j-k))/(double)deadtimeTransitionBins);
                }
            }

            else
            {
                if(k<0)
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k+config.plotConfig.TOF_BINS]*(double)(1-deadtimesPerBin[j]);
                }

                else
                {
                    deadtimesPerBin[j] += eventsPerMicroPerBin[k]*(double)(1-deadtimesPerBin[j]);
                }
            }
        }

        deadtimesPerBin[j] += (eventsPerMicroPerBin[j]/2)*(1-deadtimesPerBin[j]); // last bin contributes 1/2 its value

        // scale up events per micro on-the-fly
        //eventsPerMicroPerBin[j] *= (1-deadtimesPerBin[j]);

        deadtimeHisto->SetBinContent(j+1,deadtimesPerBin[j]);
    }

    deadtimeHisto->Write();

    return deadtimesPerBin;
}

    vector<TH1D*> energyHistos;

    string energyName =  config.targetConfig.TARGET_ORDER[i] + "Energy";
    energyHistos.push_back(timeBinsToRKEBins(TOFHistos.back(),energyName));

    energyHistos[procEvent.targetPos-1]->Fill(rKE);

    for(auto histo : energyHistos)
    {
        histo->Write();
    }

    std::vector<long> macrosPerTarget;
    std::vector<long> microsPerTarget;

    outputFile->cd("/macroTime");

    for(unsigned int i=0; i<config.targetConfig.TARGET_ORDER.size(); i++)
    {
        macrosPerTarget.push_back(((TH1I*)gDirectory->Get("targetPosH"))->GetBinContent(i+2));
        //double badMacroRatio = ((TH1D*)gDirectory->Get("ratioBadMacrosH"))->GetBinContent(i+2);
        microsPerTarget.push_back(macrosPerTarget.back()/*(1-badMacroRatio)*/
                *(config.facilityConfig.MACRO_LENGTH/config.facilityConfig.MICRO_LENGTH));
    }

    outputFile->cd();
    outputFile->cd("/summedDet");

    // perform iterative deadtime correction, until average deadtime changes
    // <0.1%
    for(unsigned int i=0; (unsigned int)i<microsPerTarget.size(); i++)
    {
        vector<double> deadtimeBins = generateDeadtimeCorrection(TOFHistos[i], microsPerTarget[i]);

        double averageDeadtimeDiff = applyDeadtimeCorrection(TOFHistos[i], deadtimeBins) - 0;

        /*while(averageDeadtimeDiff>0.001)
          {
          rawTOF = correctedTOF;

          vector<double> prevDeadtimeBins = deadtimeBins;
          deadtimeBins = generateDeadtimeCorrection(rawTOF, microsPerTarget[i]);

          for(unsigned int j=0; j<deadtimeBins.size(); j++)
          {
          deadtimeBins[j] = deadtimeBins[j]-prevDeadtimeBins[j];
          }

          correctedTOF = ((TH1D*)plots[i]->getTOFHisto());
          deadtimeH = ((TH1D*)plots[i]->getDeadtimeHisto());
          averageDeadtimeDiff = applyDeadtimeCorrection(rawTOF, correctedTOF, deadtimeH, deadtimeBins) - averageDeadtimeDiff;
          }*/
    }

    // calculate likelihood of double peak in wavelet, for each TOF bin
    //const int waveletBins = (config.plotConfig.TOF_BINS/config.plotConfig.TOF_RANGE)*WAVELET_PERIOD;

    vector<double> eventsPerMicroPerBin(config.plotConfig.TOF_BINS);
    vector<double> doublePeakOddsPerBin(config.plotConfig.TOF_BINS);
