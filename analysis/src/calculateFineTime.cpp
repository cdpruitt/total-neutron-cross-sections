// header
double calculateCFDTime(std::vector<int>* waveform,
        double baseline,
        double fraction,
        unsigned int delay,
        bool isPositiveSignal);


double calculateCFDTime(vector<int>* waveform, double baseline, double fraction, unsigned int delay, bool isPositiveSignal)
{
    if(delay<=0 || delay>=waveform->size())
    {
        cerr << "Error: cannot calculate CFD time with delay outside range of [0,waveform size] (" << delay << " was provided; waveform size() = " << waveform->size() << ")." << endl;
        return 0;
    }

    if(fraction<=0 || fraction>=1)
    {
        cerr << "Error: cannot calculate CFD time with CFD fraction outside range of [0,1] (" << fraction << " was provided)." << endl;
        return 0;
    }

    bool listenForZC = false;
    double CFDSample;
    double prevCFDSample = 0;

    if(isPositiveSignal)
    {
        for(unsigned int i=1; i<waveform->size()-(delay+1); i++)
        {
            // produce CFD sum: opposite-sign waveform*fraction + normal-sign waveform, centered at
            // baseline
            double CFDSample = waveform->at(i)-(fraction*waveform->at(i+delay)+baseline*(1-fraction));

            if(!listenForZC && CFDSample<-CFD_ZC_TRIGGER_THRESHOLD)
            {
                // approaching ZC - start looking for a ZC
                listenForZC = true;
            }

            if(listenForZC && CFDSample>0)
            {
                // found ZC: return time of crossing, i.e., (baseline-NZC)/(PZC-NZC)
                return i+(0-prevCFDSample)/(CFDSample-prevCFDSample)-TC_PRETRIGGER_SAMPLES;
            }

            prevCFDSample = CFDSample;
        }
    }

    if(!isPositiveSignal)
    {
        for(unsigned int i=1; i<waveform->size()-(delay+1); i++)
        {
            // produce CFD sum: opposite-sign waveform*fraction + normal-sign waveform, centered at
            // baseline
            double CFDSample = waveform->at(i)-(fraction*waveform->at(i+delay)+baseline*(1-fraction));

            if(!listenForZC && CFDSample>CFD_ZC_TRIGGER_THRESHOLD)
            {
                // approaching ZC - start looking for a ZC
                listenForZC = true;
            }

            if(listenForZC && CFDSample<0)
            {
                // found ZC: return time of crossing, i.e., (baseline-NZC)/(PZC-NZC)
                return i+(0-prevCFDSample)/(CFDSample-prevCFDSample)-DETECTOR_PRETRIGGER_SAMPLES;
            }

            prevCFDSample = CFDSample;

        }
    }
    
    //cerr << "Error: could not calculate fine time of waveform." << endl;

    return 0;
}

double calculateTCFineTime(vector<int>* waveform, double threshold)
{
    for(unsigned int i=1; i<waveform->size(); i++)
    {
        if(waveform->at(i) < threshold)
        {
            // found LED trigger; extract time via linear interpolation
            return i+(threshold-waveform->at(i-1))/(waveform->at(i)-waveform->at(i-1));
        }
    }

    cerr << "Error: failed to calculate macropulse start fine time. Exiting..." << endl;
    exit(1);
}

/*double calculateTCFineTime(vector<int>* waveform, double baseline, double fraction, unsigned int delay)
{
    if(delay<=0 || delay>=waveform->size())
    {
        cerr << "Error: cannot calculate CFD time with delay outside range of [0,waveform.size()] (" << delay << " was provided)." << endl;
        return 0;
    }

    if(fraction<=0 || fraction>=1)
    {
        cerr << "Error: cannot calculate CFD time with CFD fraction outside range of [0,1] (" << fraction << " was provided)." << endl;
        return 0;
    }

    bool listenForZC = false;
    double CFDSample;
    double prevCFDSample=0;

    for(unsigned int i=1; i<waveform->size()-(delay+1); i++)
    {
        // produce CFD sum: opposite-sign waveform*fraction + normal-sign waveform, centered at
        // baseline
        double CFDSample = waveform->at(i)-(fraction*waveform->at(i+delay)+baseline*(1-fraction));

        if(!listenForZC && CFDSample<-TC_ZC_TRIGGER_THRESHOLD)
        {
            // approaching ZC - start looking for a ZC
            listenForZC = true;
        }

        if(listenForZC && CFDSample>0)
        {
            // found ZC: return time of crossing, i.e., (baseline-NZC)/(PZC-NZC)
            return i+(0-prevCFDSample)/(CFDSample-prevCFDSample)-TC_PRETRIGGER_SAMPLES;
        }

        prevCFDSample = CFDSample;
    }

    return 0;
}*/


