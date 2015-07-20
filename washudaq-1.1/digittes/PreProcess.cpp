/******************************************************************************
*
* CAEN SpA - Front End Division
* Via Vetraia, 11 - 55049 - Viareggio ITALY
* +390594388398 - www.caen.it
*
***************************************************************************//**
* \note TERMS OF USE:
* This program is free software; you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free Software
* Foundation. This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The user relies on the
* software, documentation and results solely at his own risk.
******************************************************************************/

#include "digiTES.h"

// --------------------------------------------------------------------------------------------------------- 
// Global Variables
// --------------------------------------------------------------------------------------------------------- 
uint64_t PrevTime[MAX_NBRD][MAX_NCH];		// previous timestamps (used to check rollover and extend to 64 bit)
uint64_t ExtendedTime[MAX_NBRD][MAX_NCH];   // Extended timestamps (higher bits of the 64bit word)
uint32_t EventFormat[MAX_NBRD][MAX_NCH];
void *WaveformBuff;							// memory buffer for the waveform


// --------------------------------------------------------------------------------------------------------- 
// Description: Init Global Variables and allocate memory for the waveforms
// Inputs:		-
// Outputs:		AllocatedSize: total number of bytes allocated by the function
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int InitPreProcess(uint32_t *AllocatedSize)
{
	int b, ch, ret=0;
	for(b=0; b<MAX_NBRD; b++) {
		for (ch=0; ch<MAX_NCH; ch++) {
			PrevTime[b][ch] = 0;
			Qnev[b][ch] = 0;
			ExtendedTime[b][ch] = 0;
		}
	}
	if (!WDcfg.ReplayDataFile) {
		ret = CAEN_DGTZ_MallocDPPWaveforms(handle[0], &WaveformBuff, AllocatedSize); 
	}
	return ret;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: Allocate memory for the plotting waveform
// Inputs:		Waveform = pointer to the waveform to allocate
// Outputs:		AllocatedSize: total number of bytes allocated by the function
// Return: 0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int AllocateWaveform(Waveform_t *Waveform, uint32_t *AllocatedSize) 
{
	// HACK: free this memory buffers somewhere before closing
	Waveform->AnalogTrace[0]  = (uint16_t *)malloc(WDcfg.RecordLength * sizeof(uint16_t));
	Waveform->AnalogTrace[1]  = (uint16_t *)malloc(WDcfg.RecordLength * sizeof(uint16_t));
	Waveform->DigitalTrace[0] = (uint8_t *)malloc(WDcfg.RecordLength * sizeof(uint8_t));
	Waveform->DigitalTrace[1] = (uint8_t *)malloc(WDcfg.RecordLength * sizeof(uint8_t));
	Waveform->DigitalTrace[2] = (uint8_t *)malloc(WDcfg.RecordLength * sizeof(uint8_t));
	Waveform->DigitalTrace[3] = (uint8_t *)malloc(WDcfg.RecordLength * sizeof(uint8_t));
	*AllocatedSize = WDcfg.RecordLength * sizeof(uint16_t) * 2 + WDcfg.RecordLength * sizeof(uint8_t) * 4;
	return 0;
}



// --------------------------------------------------------------------------------------------------------- 
// Description:	Free memory buffer
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ClosePreProcess()
{
	CAEN_DGTZ_FreeDPPWaveforms(handle[0], WaveformBuff);
	return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: Convert data from the struct of the CAEN_DGTZ library to a generic struct valid for any DPP
// Inputs:		int b = Board Number
//				int ch = Channel Number
//				ev = Event number
//				Event = data source (event list given by the CAEN DGTZ library)
// Outputs:		EventData = data destination (generic event)
// Return: 0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int PreProcessEvent(int b, int ch, int ev, void *Event, GenericDPPEvent_t *EventData)
{
	uint64_t TSmask = 0x7FFFFFFF;

	EventData->FineTimeStamp = 0;
	if (WDcfg.DppType == DPP_CI) {
		CAEN_DGTZ_DPP_CI_Event_t *EventCI = &((CAEN_DGTZ_DPP_CI_Event_t *)Event)[ev];

		EventData->Waveforms = EventCI->Waveforms;
		EventFormat[b][ch] = EventCI->Format;
		if (EventCI->TimeTag < PrevTime[b][ch])
			ExtendedTime[b][ch] += (TSmask+1);
		EventData->TimeStamp = (ExtendedTime[b][ch] + (uint64_t)EventCI->TimeTag) * WDcfg.Tsampl + (uint64_t)WDcfg.DelayLine[b][ch];
		PrevTime[b][ch] = EventCI->TimeTag;
		EventData->TimeStamp = (uint64_t)EventCI->TimeTag;
		EventData->Energy = EventCI->Charge;
		EventData->Baseline = EventCI->Baseline;
		EventData->Flags = 0;
		EventData->psd = 0;

	} else if ((WDcfg.DppType == DPP_PSD_720) || (WDcfg.DppType == DPP_PSD_751) || (WDcfg.DppType == DPP_PSD_730)) {
		CAEN_DGTZ_DPP_PSD_Event_t *EventPSD = &((CAEN_DGTZ_DPP_PSD_Event_t *)Event)[ev];
		uint32_t ts;

		if (WDcfg.DppType == DPP_PSD_751) 
			TSmask = 0xFFFFFFFF;
		EventData->Waveforms = EventPSD->Waveforms;
		EventFormat[b][ch] = EventPSD->Format;
		ts = EventPSD->TimeTag;
		if (ts < PrevTime[b][ch])
			ExtendedTime[b][ch] += (TSmask+1);
		EventData->TimeStamp = (ExtendedTime[b][ch] + (uint64_t)ts) * WDcfg.Tsampl + (uint64_t)WDcfg.DelayLine[b][ch];
		PrevTime[b][ch] = ts;
		EventData->Energy = EventPSD->ChargeLong;
		EventData->Flags = EventPSD->Pur<<7;
		if (WDcfg.DppType == DPP_PSD_730) {
			uint16_t ZCneg, ZCpos;
			uint16_t dt = (1  +  WDcfg.CFDinterp[b][ch]*2) * 1000 * WDcfg.Tsampl;  // time interval between the interpolated points (in ps)
			uint16_t zc=0;

			//ZCneg=8192; ZCpos = 8193;
			ZCneg = EventPSD->Extras & 0xFFFF;
			ZCpos = (EventPSD->Extras>>16) & 0xFFFF;
			if ((ZCneg < 8192) && (ZCpos >= 8192))
				zc = dt * (8192 - ZCneg) / (ZCpos - ZCneg);  // fine tstamp is expressed in ps
			EventData->FineTimeStamp = zc;

			//fprintf(rolog, "N=%u, P=%u, ZC=%u, zce=%d\n", ZCneg, ZCpos, zc, zce);
			EventData->Baseline = 0;
		} else {
			EventData->Baseline = EventPSD->Baseline;
		}
		if ((EventPSD->ChargeLong != 0xFFFF) && (EventPSD->ChargeShort != 0x7FFF))
			EventData->psd = (float)(EventPSD->ChargeLong - EventPSD->ChargeShort)/(EventPSD->ChargeLong);
		else
			EventData->psd = 0;

	} else if ((WDcfg.DppType == DPP_PHA_724) || (WDcfg.DppType == DPP_PHA_730)) {
		CAEN_DGTZ_DPP_PHA_Event_t *EventPHA = &((CAEN_DGTZ_DPP_PHA_Event_t *)Event)[ev];

		EventData->Waveforms = EventPHA->Waveforms;
		EventData->Baseline = 0;
		EventFormat[b][ch] = EventPHA->Format;
		if (EventPHA->TimeTag < PrevTime[b][ch])
			ExtendedTime[b][ch] += (TSmask+1);
		EventData->TimeStamp = (ExtendedTime[b][ch] + (uint64_t)EventPHA->TimeTag) * WDcfg.Tsampl + (uint64_t)WDcfg.DelayLine[b][ch];
		PrevTime[b][ch] = EventPHA->TimeTag;
		EventData->Energy =  EventPHA->Energy;
		EventData->Flags = (char)EventPHA->Extras;


		if (WDcfg.DppType == DPP_PHA_730) {
			uint16_t ZCneg, ZCpos;
			uint16_t dt = (1  +  WDcfg.CFDinterp[b][ch]*2) * 1000 * WDcfg.Tsampl;  // time interval between the interpolated points (in ps)
			uint16_t zc=0;

			//ZCneg=8192; ZCpos = 8193;
			ZCneg = EventPHA->Extras & 0xFFFF;
			ZCpos = (EventPHA->Extras>>16) & 0xFFFF;
			if ((ZCneg < 8192) && (ZCpos >= 8192))
				zc = dt * (8192 - ZCneg) / (ZCpos - ZCneg);  // fine tstamp is expressed in ps
			EventData->FineTimeStamp = zc;

			//fprintf(rolog, "N=%u, P=%u, ZC=%u, zce=%d\n", ZCneg, ZCpos, zc, zce);
			EventData->Baseline = 0;
		}

	}

	return 0;
}



// --------------------------------------------------------------------------------------------------------- 
// Description: Decode Waveform data of an Event
// Inputs:		b = Board Number
//				Event = Event To decode
// Outputs		Wfm = Data destination 
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ConvertWaveform(int b, int ch, GenericDPPEvent_t Event, Waveform_t *Wfm)
{
	int i;
	if (Event.Waveforms == NULL)
		return -1;

	// Take trace settings from the global variable in the software and not from data reaad from the board (trace settings are present in the "format").
	// Trace settings are assigned here only in case of data from file.
	for(i=0; i<MAX_NTRACES; i++)
		Wfm->TraceSet[i] = TraceSet[i];


	if (WDcfg.ReplayDataFile) {
		uint16_t i, smask;
		uint32_t wfmt;

		Wfm->Ns = (*Event.Waveforms & 0xFFFF);
		wfmt = *(Event.Waveforms + 1);
		Wfm->TraceSet[0] = (wfmt >> 16) & 0xF;  
		Wfm->TraceSet[1] = (wfmt >> 20) & 0xF;
		Wfm->TraceSet[2] = (wfmt)       & 0xF;
		Wfm->TraceSet[3] = (wfmt >>  4) & 0xF;
		Wfm->TraceSet[4] = (wfmt >>  8) & 0xF;
		Wfm->TraceSet[5] = (wfmt >> 12) & 0xF;
		Wfm->DualTrace =   (wfmt >> 31) & 0x1;
		smask = (0xFFFF >> (16-WDcfg.Nbit));

		for(i=0; i<(Wfm->Ns/2); i++) {
			uint32_t w = *(Event.Waveforms + i + 2);
			Wfm->AnalogTrace[0][i*2] = w & smask;
			if (Wfm->DualTrace) {
				Wfm->AnalogTrace[0][i*2+1] = w & smask;
				Wfm->AnalogTrace[1][i*2] = (w>>16) & smask;
				Wfm->AnalogTrace[1][i*2+1] = (w>>16) & smask;
			} else {
				Wfm->AnalogTrace[0][i*2+1] = (w>>16) & smask;
				Wfm->AnalogTrace[1][i*2] = 0;
				Wfm->AnalogTrace[1][i*2+1] = 0;
			}
			Wfm->DigitalTrace[0][i*2]   = (w>>12) & 1;
			Wfm->DigitalTrace[0][i*2+1] = (w>>28) & 1;
			Wfm->DigitalTrace[1][i*2]   = (w>>13) & 1;
			Wfm->DigitalTrace[1][i*2+1] = (w>>29) & 1;
			Wfm->DigitalTrace[2][i*2]   = (w>>14) & 1;
			Wfm->DigitalTrace[2][i*2+1] = (w>>30) & 1;
			Wfm->DigitalTrace[3][i*2]   = (w>>15) & 1;
			Wfm->DigitalTrace[3][i*2+1] = (w>>31) & 1;
		}

	} else if (WDcfg.DppType == DPP_CI) {
		CAEN_DGTZ_DPP_CI_Event_t EventCI;
		CAEN_DGTZ_DPP_CI_Waveforms_t *WaveformCI = (CAEN_DGTZ_DPP_CI_Waveforms_t *)WaveformBuff;

		EventCI.Waveforms = Event.Waveforms;
		EventCI.Format = EventFormat[b][ch];
		CAEN_DGTZ_DecodeDPPWaveforms(handle[b], &EventCI, WaveformCI);
		Wfm->Ns = WaveformCI->Ns;
		Wfm->DualTrace = WaveformCI->dualTrace;

		if (TraceSet[0] == 7)
			WaveformProcessor(b, ch, Wfm->Ns, WaveformCI->Trace1, &Event, Wfm->AnalogTrace[0]);
		else
			memcpy(Wfm->AnalogTrace[0], WaveformCI->Trace1, Wfm->Ns * sizeof(uint16_t));
		if (TraceSet[1] == 7) 
			WaveformProcessor(b, ch, Wfm->Ns, WaveformCI->Trace1, &Event, Wfm->AnalogTrace[1]);
		else if (Wfm->DualTrace)
			memcpy(Wfm->AnalogTrace[1], WaveformCI->Trace2, Wfm->Ns * sizeof(uint16_t));

		memcpy(Wfm->DigitalTrace[0], WaveformCI->DTrace1, Wfm->Ns * sizeof(uint8_t));
		memcpy(Wfm->DigitalTrace[1], WaveformCI->DTrace2, Wfm->Ns * sizeof(uint8_t));
		memcpy(Wfm->DigitalTrace[2], WaveformCI->DTrace3, Wfm->Ns * sizeof(uint8_t));
		memcpy(Wfm->DigitalTrace[3], WaveformCI->DTrace4, Wfm->Ns * sizeof(uint8_t));

	} else if ((WDcfg.DppType == DPP_PSD_720) || (WDcfg.DppType == DPP_PSD_751) || (WDcfg.DppType == DPP_PSD_730)) {
		CAEN_DGTZ_DPP_PSD_Event_t EventPSD;
		CAEN_DGTZ_DPP_PSD_Waveforms_t *WaveformPSD = (CAEN_DGTZ_DPP_PSD_Waveforms_t *)WaveformBuff;

		EventPSD.Waveforms = Event.Waveforms;
		EventPSD.Format = EventFormat[b][ch];
		CAEN_DGTZ_DecodeDPPWaveforms(handle[b], &EventPSD, WaveformPSD);
		Wfm->Ns = WaveformPSD->Ns;
		if (Wfm->Ns > WDcfg.RecordLength)
			Wfm->Ns = WDcfg.RecordLength;  // in the 751 the nymber of samples can be higher than RecordLength
		Wfm->DualTrace = WaveformPSD->dualTrace;

		if (TraceSet[0] == 7) 
			WaveformProcessor(b, ch, Wfm->Ns, WaveformPSD->Trace1, &Event, Wfm->AnalogTrace[0]);
		else
			memcpy(Wfm->AnalogTrace[0], WaveformPSD->Trace1, Wfm->Ns * sizeof(uint16_t));

		if (TraceSet[1] == 7) 
			WaveformProcessor(b, ch, Wfm->Ns, WaveformPSD->Trace1, &Event, Wfm->AnalogTrace[1]);
		else if (Wfm->DualTrace)
			memcpy(Wfm->AnalogTrace[1], WaveformPSD->Trace2, Wfm->Ns * sizeof(uint16_t));

		if (WDcfg.DppType == DPP_PSD_720) {
			memcpy(Wfm->DigitalTrace[0], WaveformPSD->DTrace1, Wfm->Ns * sizeof(uint8_t));
			memcpy(Wfm->DigitalTrace[1], WaveformPSD->DTrace2, Wfm->Ns * sizeof(uint8_t));
			memcpy(Wfm->DigitalTrace[2], WaveformPSD->DTrace3, Wfm->Ns * sizeof(uint8_t));
			memcpy(Wfm->DigitalTrace[3], WaveformPSD->DTrace4, Wfm->Ns * sizeof(uint8_t));

		} else if (WDcfg.DppType == DPP_PSD_730) {
			memcpy(Wfm->DigitalTrace[0], WaveformPSD->DTrace1, Wfm->Ns * sizeof(uint8_t));
			memcpy(Wfm->DigitalTrace[1], WaveformPSD->DTrace2, Wfm->Ns * sizeof(uint8_t));
		} else if (WDcfg.DppType == DPP_PSD_751) {
			memcpy(Wfm->DigitalTrace[0], WaveformPSD->DTrace1, Wfm->Ns * sizeof(uint8_t));
			memcpy(Wfm->DigitalTrace[1], WaveformPSD->DTrace2, Wfm->Ns * sizeof(uint8_t));
			memcpy(Wfm->DigitalTrace[2], WaveformPSD->DTrace3, Wfm->Ns * sizeof(uint8_t));
		}


	} else if ((WDcfg.DppType == DPP_PHA_724) || (WDcfg.DppType == DPP_PHA_730)) {
		CAEN_DGTZ_DPP_PHA_Event_t EventPHA;
		CAEN_DGTZ_DPP_PHA_Waveforms_t *WaveformPHA = (CAEN_DGTZ_DPP_PHA_Waveforms_t *)WaveformBuff;

		EventPHA.Waveforms = Event.Waveforms;
		EventPHA.Format = EventFormat[b][ch];
		CAEN_DGTZ_DecodeDPPWaveforms(handle[b], &EventPHA, WaveformPHA);
		Wfm->Ns = WaveformPHA->Ns;
		Wfm->DualTrace = WaveformPHA->DualTrace;

		if (TraceSet[0] == 7)
			WaveformProcessor(b, ch, Wfm->Ns, (uint16_t *)WaveformPHA->Trace1, &Event, Wfm->AnalogTrace[0]);  // HACK: check sign here!
		else
			memcpy(Wfm->AnalogTrace[0], WaveformPHA->Trace1, Wfm->Ns * sizeof(uint16_t));
		if (TraceSet[1] == 7) 
			WaveformProcessor(b, ch, Wfm->Ns, (uint16_t *)WaveformPHA->Trace1, &Event, Wfm->AnalogTrace[1]);
		else if (Wfm->DualTrace)
			memcpy(Wfm->AnalogTrace[1], WaveformPHA->Trace2, Wfm->Ns * sizeof(uint16_t));

		if (TraceSet[0] == 3)   
			for (i=0; i<Wfm->Ns; i++) 
				Wfm->AnalogTrace[0][i] = (uint16_t)((int16_t)Wfm->AnalogTrace[0][i] / WDcfg.PHrescale[b][ch]);
		if ((TraceSet[1] == 2) || (TraceSet[1] == 3))  
			for (i=0; i<Wfm->Ns; i++) 
				Wfm->AnalogTrace[1][i] = (uint16_t)((int16_t)Wfm->AnalogTrace[1][i] / WDcfg.PHrescale[b][ch]);

		memcpy(Wfm->DigitalTrace[0], WaveformPHA->DTrace1, Wfm->Ns * sizeof(uint8_t));
		memcpy(Wfm->DigitalTrace[1], WaveformPHA->DTrace2, Wfm->Ns * sizeof(uint8_t));
	}
	return 0;
}



// --------------------------------------------------------------------------------------------------------- 
// Description: Off-line implementation of the discriminator (LED or CFD), time interpolation and PSD
//				The results are written into the Event (including discriminator waveform)
// Inputs:		b = Board Number
//				c = channel number
// InOuts:		Event = Event data
// Outputs		Waveout = waveform generated by the processor (e.g. CFD)
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int WaveformProcessor(int b, int ch, int ns, uint16_t *Wavein, GenericDPPEvent_t *Event, uint16_t *Waveout)
{
	int i, nsbl=0, ncross=0, armed=0, bslstop, sign;
	int LGwidth, SGwidth, PreGate;
	float bsl, Qs=0, Ql=0, atten;
	float ZCneg=0, ZCpos=0;
	uint32_t bslacc;
	uint16_t *wave;
	float *wavediscr;

	//ret = ConvertWaveform(b, ch, *Event, &PPwaveform);
	//wave = PPwaveform.AnalogTrace[0];
	wave = Wavein;
	LGwidth = WDcfg.GateWidth[b][ch]/WDcfg.Tsampl;
	SGwidth = WDcfg.GateWidth[b][ch]/WDcfg.Tsampl;
	PreGate = WDcfg.PreGate[b][ch]/WDcfg.Tsampl;

	if      (WDcfg.CFDatten[b][ch] == 0) atten = 0.25;
	else if (WDcfg.CFDatten[b][ch] == 1) atten = 0.50;
	else if (WDcfg.CFDatten[b][ch] == 2) atten = 0.75;
	else if (WDcfg.CFDatten[b][ch] == 3) atten = 1.00;

	// discriminator is intended for negative pulses; if positive, invert it (NOTE: in PHA, negaive signals are already inverted by the FPGA, so they are always positive)
	if ((WDcfg.PulsePolarity[b][ch] == 0) || (WDcfg.DppType == DPP_PHA_724))
		sign = -1;
	else
		sign = 1;

	wavediscr = (float *)malloc(ns * sizeof(*wavediscr));

	if (Event->Baseline == 0) {
		// Calculate Baseline (when not coming from on-line DPP)
		bslacc = 0;
		bsl = wave[0];
		for(i=0; i<ns; i++) {
			if (fabs(wave[i] - bsl) > (double)(WDcfg.TrgThreshold[b][ch]/2))
				break;
			bslacc += wave[i];
			nsbl++;
			bsl = (float)bslacc/nsbl;
		}
		bslstop = i-1;
		for(i=0; i<4; i++) {
			if ((bslstop - i) > 1) {
				bslacc -= wave[bslstop - i];
				nsbl--;
			}
		}
		bsl = (float)bslacc/nsbl;
		Event->Baseline = (uint16_t)bsl;
	} else {
		bsl = (float)Event->Baseline;
	}

	// calculate discriminator waveform (either LED or CFD)
	for(i=0; i<ns; i++) {
		if (WDcfg.DiscrMode[b][ch] == 1) {  // CFD
			if (i >= WDcfg.CFDdelay[b][ch])
				wavediscr[i] = sign * (atten  * (wave[i] - bsl) - (wave[i - WDcfg.CFDdelay[b][ch]] - bsl));
			else
				wavediscr[i] = sign * (wave[i] - bsl);
			if (!armed && (wavediscr[i] < -WDcfg.TrgThreshold[b][ch]))
				armed = 1;
			if (armed && (ncross==0) && (wavediscr[i] >= 0)) {
				ncross = i;
				if (WDcfg.CFDinterp[b][ch] == 0) {
					ZCneg = wavediscr[i-1];
					ZCpos = wavediscr[i];
				} else if (WDcfg.CFDinterp[b][ch] == 1) {
					ZCneg = (wavediscr[i-1] + wavediscr[i-2])/2;
					ZCpos = (wavediscr[i] + wavediscr[i+1])/2;
				}
			}
		} else {  // LED
			wavediscr[i] = sign * (wave[i] - bsl); 
			if ((ncross==0) && (wavediscr[i] < -WDcfg.TrgThreshold[b][ch])) {
				ncross = i;
				ZCneg = WDcfg.TrgThreshold[b][ch] - wavediscr[i-1];
				ZCpos = WDcfg.TrgThreshold[b][ch] - wavediscr[i];
			}
		}
		if (Waveout != NULL)
			Waveout[i] = (uint16_t)wavediscr[i] + (1<<WDcfg.Nbit)/2;


		if ((ncross>0) && (i <= (ncross + LGwidth - PreGate))) 
			Ql += (bsl - wave[i-PreGate]);
		if ((ncross>0) && (i <= (ncross + SGwidth - PreGate))) 
			Qs += (bsl - wave[i-PreGate]);
	}

	if (ncross > 0) {
		Event->TimeStamp += WDcfg.Tsampl * ncross;
		if ((ZCneg < 0) && (ZCpos >= 0)) {
			if (WDcfg.CFDinterp[b][ch] == 0)
				Event->FineTimeStamp = (uint16_t)(1000 * WDcfg.Tsampl * (-ZCneg) / (ZCpos - ZCneg));  // fine tstamp is expressed in ps
			else
				Event->FineTimeStamp = (uint16_t)(1000 * 2*WDcfg.Tsampl * ((-ZCneg) / (ZCneg - ZCpos) - 0.5));
		} else {
			Event->FineTimeStamp = 0;
		}
	} 

	//Event->Energy = (uint16_t)Ql;
	//Event->psd = (Ql-Qs)/Ql;

	free(wavediscr);

	return 0;
}

