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

    //TH1D* tof = (TH1D*)plots[0]->getTOFHisto();
    /*for(int i=0; i<eventsPerMicroPerBin.size(); i++)
      {
      if(microsPerTarget[0]==0)
      {
      cerr << "Error: can't calculate double peak likelihood with 0 total micros" << endl;
      continue;
      }

      eventsPerMicroPerBin[i] = tof->GetBinContent(i)/(double)microsPerTarget[0];
      }

      for(int i=0; i<config.plotConfig.TOF_BINS; i++)
      {
      for(int j=i; j<i+waveletBins/2; j++)
      {
      if(j>config.plotConfig.TOF_BINS)
      {
      doublePeakOddsPerBin[i] += eventsPerMicroPerBin[j-config.plotConfig.TOF_BINS];
      continue;
      }

      doublePeakOddsPerBin[i] += eventsPerMicroPerBin[j];
      }
      }

      for(int i=0; i<doublePeakOddsPerBin.size(); i++)
      {
      doublePeakH->SetBinContent(i,pow(10,3)*doublePeakOddsPerBin[i]);
      }

      doublePeakH->Write();
      */


