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
// Description: Program the registers of the digitizer with the relevant parameters
// Inputs:		brd = board number
// Outputs:		-
// Return:		Error code (0=success) 
// --------------------------------------------------------------------------------------------------------- 
int ProgramDigitizer(int brd)
{
	int i, ret = 0, EnableMask=0;  
	char CfgStep[100];
	uint32_t reg;



	// ########################################################################################################
	sprintf(CfgStep, "General System Configuration and Acquisition Mode");
	// ########################################################################################################
	for(i=0; i<MAX_NCH; i++)
		EnableMask |= (WDcfg.EnableInput[brd][i]<<i);

	if ((WDcfg.DppType == DPP_PHA_730) || (WDcfg.DppType == DPP_PHA_724))
		WDcfg.Tsampl *= (1 << WDcfg.Decimation[0][0]);  // HACK: assuming the same for all channels

	// Set Fan speed (desktop versions)
	ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8168, (WDcfg.FanSpeed<<3) & 0x8);  

	// Acquisition mode
	ret |= CAEN_DGTZ_SetDPPAcquisitionMode(handle[brd], (CAEN_DGTZ_DPP_AcqMode_t)WDcfg.AcquisitionMode, CAEN_DGTZ_DPP_SAVE_PARAM_EnergyAndTime); 
	if (WDcfg.DigitizerModel == 730)
		CAEN_DGTZ_WriteRegister(handle[brd], 0x8004, 1<<17); // Enable Extra Word

	// Set the number of samples in the waveform
	ret |= CAEN_DGTZ_SetRecordLength(handle[brd], WDcfg.RecordLength, -1);
	// Set the enabled channels
	ret |= CAEN_DGTZ_SetChannelEnableMask(handle[brd], EnableMask);

	// some specific settings in the global CTRL register
	if (WDcfg.DppType == DPP_PHA_724)
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8000, 0x01000114);  // Channel Control Reg (indiv trg, seq readout) ??
	if (WDcfg.DppType == DPP_PSD_751)
		CAEN_DGTZ_WriteRegister(handle[brd], 0x8004, 1<<26); // Set new probe mode in x751 models
	if (WDcfg.DppType == DPP_PHA_730)
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8004, 1);  // Enable auto-flush
	/*if (WDcfg.DppType == DPP_PSD_720) {
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8004, 1<<17);  // Enable baseline or Extended Time stamp
		ret |= RegisterSetBits(handle[brd], 0x1080, 7, 7, 1);  // Enable ext time stamp
	}*/


	// Change direction for Trigger Validation signals (from mother board to daugther board)
	if (WDcfg.EnableCoinc) {
		printf(">>>>>> WARNING; coincidences are enabled\n\n");
		ret |= CAEN_DGTZ_ReadRegister (handle[brd], 0x8000, &reg);			 // Ch Config Register
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8000, reg | (1<<2));   // set trgout bidir
	}
	if (ret) goto abortcfg;


	// ########################################################################################################
	sprintf(CfgStep, "Front Panel I/O, Trigger Modes and Synchronization");
	// ########################################################################################################
	// Set the I/O level (CAEN_DGTZ_IOLevel_NIM or CAEN_DGTZ_IOLevel_TTL)
	ret |= CAEN_DGTZ_SetIOLevel(handle[brd], (CAEN_DGTZ_IOLevel_t)WDcfg.FPIOtype);

	// settings for the synchronization
	if (WDcfg.StartMode == START_MODE_INDEP_SW) {
		ret |= CAEN_DGTZ_SetAcquisitionMode(handle[brd], CAEN_DGTZ_SW_CONTROLLED);
		ret |= CAEN_DGTZ_SetRunSynchronizationMode(handle[brd], CAEN_DGTZ_RUN_SYNC_Disabled); 
	} else {
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], CAEN_DGTZ_ACQ_CONTROL_ADD, RUN_START_ON_TRGIN_RISING_EDGE);    // Arm acquisition (Run will start with 1st trigger)
		//ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8170, (WDcfg.NumBrd - brd - 1) * 4);   // Run Delay to deskew the start of acquisition
		// IGRIS setup
		if (brd==0) ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8170, 12);   
		if (brd==1) ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8170, 8);   
		if (brd==2) ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8170, 5);   
		if (brd==3) ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8170, 2);   
	}


	// Trigger modes
	if (WDcfg.TriggerMode == TRIGGER_MODE_SELF)
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], CAEN_DGTZ_TRIGGER_SRC_ENABLE_ADD, 0x80000000);   // accept SW trg only (ext trig is disabled)
	else
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], CAEN_DGTZ_TRIGGER_SRC_ENABLE_ADD, 0xC0000000);   // accept ext trg_in (from trg OR) or SW trg

	// TRGOUT/GPO mode
	if (WDcfg.TrgoutMode == TRGOUT_MODE_CH_TRG) {
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], CAEN_DGTZ_FP_TRIGGER_OUT_ENABLE_ADD, EnableMask);    // propagate self trg to TRGOUT
	} else {
		ret |= RegisterSetBits(handle[brd], CAEN_DGTZ_FRONT_PANEL_IO_CTRL_ADD, 15, 15, 1);
		ret |= RegisterSetBits(handle[brd], 0x8168, 0, 2, WDcfg.TrgoutMode);
	}
	if (ret) goto abortcfg;


	// ########################################################################################################
	sprintf(CfgStep, "Buffer Organization");
	// ########################################################################################################
	// Set Aggregation Mode
	ret |= CAEN_DGTZ_SetDPPEventAggregation(handle[brd], 0, 0);
	ret |= CAEN_DGTZ_SetMaxNumAggregatesBLT(handle[brd], 255);	// Number of buffers per BLT
	if (1) {  	// manual settings (for debug)
		printf(">>> WARNING: manual settings of memory segmentation\n");
		if (WDcfg.AcquisitionMode == ACQ_MODE_MIXED) {
			CAEN_DGTZ_WriteRegister(handle[brd], 0x800C, 0x3);		// Buffer organization (8 buffers)
			CAEN_DGTZ_WriteRegister(handle[brd], 0x8034, 1);		// Number of events per buffer
		} else {
			CAEN_DGTZ_WriteRegister(handle[brd], 0x800C, 0xA);		// Buffer organization (1024 buffers)
			CAEN_DGTZ_WriteRegister(handle[brd], 0x8034, 16);		// Number of events per buffer
		}
		CAEN_DGTZ_SetMaxNumAggregatesBLT(handle[brd], 1023);	// Number of buffers per BLT
	}
	if (0) {
		uint32_t d32;
		CAEN_DGTZ_ReadRegister(handle[brd], 0x8000, &d32);
		printf("0x8000: Ctrl Reg = %08X\n", d32);
		CAEN_DGTZ_ReadRegister(handle[brd], 0x1020, &d32);
		printf("0x1020: RecLength = %d\n", d32);
		CAEN_DGTZ_ReadRegister(handle[brd], 0x100C, &d32);
		printf("0x100C: Buffer Organization = %d\n", d32);
		CAEN_DGTZ_ReadRegister(handle[brd], 0x1034, &d32);
		printf("0x1034: Num Events x Buff= %d\n", d32);
		CAEN_DGTZ_ReadRegister(handle[brd], 0xEF1C, &d32);
		printf("0xEF1C: Num Aggr x BLT= %d\n", d32);
	}
	if (ret) goto abortcfg;




	// ########################################################################################################
	sprintf(CfgStep, "Channel individual settings");
	// ########################################################################################################
	for(i=0; i<MAX_NCH; i++) {
		if (WDcfg.EnableInput[brd][i]) { 
			// DC offset
			ret |= CAEN_DGTZ_SetChannelDCOffset(handle[brd], i, WDcfg.DCoffset[brd][i]);

			// Set the Pre-Trigger size (in nsec)
			ret |= CAEN_DGTZ_SetDPPPreTriggerSize(handle[brd], i, WDcfg.PreTrigger/WDcfg.Tsampl);

			if (WDcfg.DppType == DPP_PHA_730) {
				int MM, m, k, b, a, shf, ftd, thr, trgho, tu;
				uint32_t IPEperiod;

				tu = 8 * (1 << WDcfg.Decimation[brd][i]);
				MM =  (int)(WDcfg.TrapPoleZero[brd][i]/tu);
				k =   (int)(WDcfg.TrapRiseTime[brd][i]/tu);
				m =   (int)(WDcfg.TrapFlatTop[brd][i]/tu);
				ftd = (int)(WDcfg.PeakingTime[brd][i]/tu);
				b =   (int)(WDcfg.TTFdelay[brd][i]/tu);
				a =   (int)(WDcfg.TTFsmoothing[brd][i]);
				thr = (int)(WDcfg.TrgThreshold[brd][i]);
				trgho =  (int)(WDcfg.TrgHoldOff);
				for (shf=30; shf>0; shf--) {
					if ((1<<shf) < (k * MM / WDcfg.TrapDigitalGain[brd][i]))
						break;
				}
				WDcfg.PHrescale[brd][i] = (float)(k * MM / WDcfg.TrapDigitalGain[brd][i]) / (1 << shf);

				//ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1034 + (i<<8), 1);	// nevbuf
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1038 + (i<<8), (WDcfg.PreTrigger/WDcfg.Tsampl)/4);	// pretrg
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1054 + (i<<8), a);	// a
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1058 + (i<<8), b);	// b
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x105C + (i<<8), k);	// k
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1060 + (i<<8), m);	// m
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1064 + (i<<8), ftd);	// ftd
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1068 + (i<<8), MM);	// MM
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x106C + (i<<8), thr);  // thr
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1070 + (i<<8), 0);	// tww
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1074 + (i<<8), trgho);// trgho
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1078 + (i<<8), 8);	// pkho
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x107C + (i<<8), 4);	// blho
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1084 + (i<<8), 4);	// tvaw

				//ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1080 + (i<<8), 0);	
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 0, 5, shf);	// shf
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 12, 13, 0);	// nspk
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 16, 16, WDcfg.PulsePolarity[brd][i]);	// inv
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 18, 19, 0);	// trgmode
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 20, 22, WDcfg.NsBaseline[brd][i]);	// nsbl
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 8, 9, WDcfg.Decimation[brd][i]);	// decimation
				//ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 26, 26, 1);	// enable rollover tracing

		
				// Internal Pulse Emulator (for DPP_PHA_730 oly)
				if (WDcfg.IPErandom[brd][i])
					IPEperiod = 0xFFE0;  // HACK: make the calculation...
				else
					IPEperiod = (uint32_t)(10000000/WDcfg.IPEfrequency[brd][i]);

				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 14, 14, WDcfg.EnableIPE[brd][i] & 1);	// enable test pulse
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1044 + (i<<8), (WDcfg.IPErandom[brd][i] << 16) | (IPEperiod & 0xFFFF));  
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1048 + (i<<8), (uint32_t)(65536*exp(-1/(WDcfg.IPEdecay[brd][i]*125))) | (WDcfg.IPErisetime[brd][i]<<16) | (WDcfg.IPEaddnoise[brd][i]<<19));
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x104C + (i<<8), WDcfg.IPEamplitude[brd][i]);	  // amplitude

			} else if (WDcfg.DppType == DPP_PHA_724) {
				CAEN_DGTZ_DPP_PHA_Params_t DPPParams;
				memset(&DPPParams, 0, sizeof(CAEN_DGTZ_DPP_PHA_Params_t));
				DPPParams.k[i] = WDcfg.TrapRiseTime[brd][i];
				DPPParams.m[i] = WDcfg.TrapFlatTop[brd][i]; 
				DPPParams.M[i] = WDcfg.TrapPoleZero[brd][i];     
				DPPParams.ftd[i] = WDcfg.PeakingTime[brd][i];
				DPPParams.a[i] = WDcfg.TTFsmoothing[brd][i]; 
				DPPParams.b[i] = WDcfg.TTFdelay[brd][i];
				DPPParams.thr[i] = WDcfg.TrgThreshold[brd][i];
				DPPParams.trgho[i] = WDcfg.TrgHoldOff; 
				DPPParams.nsbl[i] = WDcfg.NsBaseline[brd][i];     
				DPPParams.nspk[i] = 0;
				DPPParams.pkho[i] = WDcfg.PeakHoldOff[brd][i];
				DPPParams.blho[i] = 1;
				DPPParams.enf[i] = WDcfg.TrapDigitalGain[brd][i];

				ret |= CAEN_DGTZ_SetDPPParameters(handle[brd], 1<<i, &DPPParams);

				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1038 + (i<<8), WDcfg.PreTrigger);		// pretrg
				ret |= CAEN_DGTZ_SetChannelPulsePolarity(handle[brd], i, (CAEN_DGTZ_PulsePolarity_t)WDcfg.PulsePolarity[brd][i]);

				// Set input range for DT5780
				if (0) {
					int drange;
					switch (WDcfg.InputDynamicRange[brd][i]) {
						case 0 : drange = 5; break;  // 0.6Vpp
						case 1 : drange = 6; break;  // 1.4Vpp
						case 2 : drange = 9; break;  // 3.7Vpp
						case 3 : drange = 10; break; // 9.5Vpp
					}
					ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x10B4 + (i<<8), drange & 0xF);  
				}

				// Set trapezoid+input
				ret |= RegisterSetBits(handle[brd], 0x8000, 11, 11, 1);	// dual trace
				ret |= RegisterSetBits(handle[brd], 0x8000, 12, 13, 3);	// trace1 = trapez
				ret |= RegisterSetBits(handle[brd], 0x8000, 14, 15, 0);	// trace2 = trapez

			} else {  // PSD & CI

				// Trigger Threshold
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1060 + (i<<8), WDcfg.TrgThreshold[brd][i]); 

				// Trigger Holdoff
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1074 + (i<<8), WDcfg.TrgHoldOff); 

				// Gate Settings
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1058 + (i<<8), WDcfg.GateWidth[brd][i]/WDcfg.Tsampl); 
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1054 + (i<<8), WDcfg.ShortGateWidth[brd][i]/WDcfg.Tsampl); 
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x105C + (i<<8), WDcfg.PreGate[brd][i]/WDcfg.Tsampl); 

				// Baseline
				ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1064 + (i<<8), WDcfg.FixedBaseline[brd][i]); 
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 20, 21, WDcfg.NsBaseline[brd][i]);

				// digital CFD and interpolated zero crossing (730 only)
				if (WDcfg.DppType == DPP_PSD_730) {
					ret |= RegisterSetBits(handle[brd], 0x103C + (i<<8), 0, 7, WDcfg.CFDdelay[brd][i]/WDcfg.Tsampl); 
					ret |= RegisterSetBits(handle[brd], 0x103C + (i<<8), 8, 9, WDcfg.CFDatten[brd][i]); 
					ret |= RegisterSetBits(handle[brd], 0x103C + (i<<8), 10, 11, WDcfg.CFDinterp[brd][i]); 
					// Set Discriminator Mode 0=LED, 1=CFD (x730 only)
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 6, 6, WDcfg.DiscrMode[brd][i]); 
					// Set Extra Word = CFD zero crossing samples 
					ret |= RegisterSetBits(handle[brd], 0x1084 + (i<<8), 8, 9, 1);   // Extra Word L = last sample before ZC
					ret |= RegisterSetBits(handle[brd], 0x1084 + (i<<8), 10, 11, 1); // Extra Word H = first sample after ZC
				}


				// Charge Pedestal
				if (WDcfg.EnablePedestal[brd][i])
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8),  4,  4, 1);  // enable pedestal

				// Charge Sensitivity
				if (WDcfg.DigitizerModel == 730) 
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8),  0,  2, WDcfg.ChargeSensitivity[brd][i]);
				else
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8),  0,  1, WDcfg.ChargeSensitivity[brd][i]);

				// Pulse Polarity
				ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 16, 16, WDcfg.PulsePolarity[brd][i]);

				// Internal Pulse Emulator
				if (WDcfg.EnableIPE[brd][i]) {
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 8, 8, 1);
					if (WDcfg.IPEfrequency[brd][i] <= 1)
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 9, 10, 0);
					else if (WDcfg.IPEfrequency[brd][i] <= 10)
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 9, 10, 1);
					else if (WDcfg.IPEfrequency[brd][i] <= 100)
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 9, 10, 2);
					else 
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 9, 10, 3);
				}

				// Set input range for x730 
				if (WDcfg.DigitizerModel == 730) {
					ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1028 + (i<<8), WDcfg.InputDynamicRange[brd][i] & 1);  // 0=2Vpp, 1=0.5Vpp
				}


				// Bit Downgrade
				if (WDcfg.DigitizerModel == 730) {
					int dgrade = 0; // 0=14bit, 1=13bit, 2=12bit, 3=10bit
					int dsamp = 0;  // 0=500MSps 1=250MSps
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 12, 13, dgrade);  
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 14, 14, dsamp);   
				}

				// Set PSD cut
				if (WDcfg.PSDcut[brd][i] > 0) {
					ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1078 + (i<<8),  WDcfg.PSDcut[brd][i]);   // Set threshold
					//ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 27, 28, 1);  // gammas only
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 27, 28, 2);  // neutrons only
				}

				// Set pile-up extension (x751 only) or Rejection (to be completed...)
				// When enabled, the gate are re-opened in case of pile up (two events are generated)
				if (WDcfg.PileUpMode[brd][i]) {
					// pile-up managent: 0=reject, 1=extend waveform (751 only)
					int pu_mode = 0;
					ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x107c + (i<<8), WDcfg.PurGap[brd][i]); 
					if (pu_mode==0)
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 26, 26, 1);   
					else if (WDcfg.DigitizerModel == 751)
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 25, 25, 1);   
				}

				// Coincidence 
				if (WDcfg.EnableCoinc) {
					// Trgout Width 
					ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x1070 + (i<<8), WDcfg.CoincWindow);

					// Set External TRGIN = Gate
					if (0)
						ret |= RegisterSetBits(handle[brd], 0x811C, 10, 11, 3);   

					// Set Extra Coinc Window Width
					//ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x106C + (i<<8), 10);   

					if (WDcfg.CoincMode[brd][i] == 1)  // coinc
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 18, 19, 1);
					else if (WDcfg.CoincMode[brd][i] == 2)  // anti-coinc
						ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 18, 19, 3);

					// Coinc between channels in the same couple (x730 only)
					//if ((WDcfg.DigitizerModel == 730) && ((i & 1) == 0)) {
					if ((WDcfg.DigitizerModel == 730)) {
						ret |= RegisterSetBits(handle[brd], 0x1084 + (i<<8), 0, 1, 3);  // couple trgout = OR
						ret |= RegisterSetBits(handle[brd], 0x1084 + (i<<8), 2, 2, 1);  // enable couple trgout
						//ret |= RegisterSetBits(handle[brd], 0x1084 + (i<<8), 4, 5, 2);  // validation from other ch in the couple
						ret |= RegisterSetBits(handle[brd], 0x1084 + (i<<8), 4, 5, 1);  // validation from other couples
						ret |= RegisterSetBits(handle[brd], 0x1084 + (i<<8), 6, 6, 1);  // enable trg validation
					}

					/*if ((WDcfg.DigitizerModel == 730) && ((i & 1) == 0)) {
						ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8180 + (i/2)*4, WDcfg.CoincMask[brd][i]);    // Coinc Mask
					}*/

				} else {  // coincidences disabled
					ret |= RegisterSetBits(handle[brd], 0x1080 + (i<<8), 18, 19, 0);
				}

			}
		}
	}
	if (ret) goto abortcfg;

	ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8180, 0x103);    // Coinc Mask
	ret |= CAEN_DGTZ_WriteRegister(handle[brd], 0x8184, 0x103);    // Coinc Mask


	// ########################################################################################################
	sprintf(CfgStep, "Generic Write accesses with mask");
	// ########################################################################################################
	for(i=0; i<WDcfg.GWn; i++) {
		ret |= CAEN_DGTZ_ReadRegister(handle[brd], WDcfg.GWaddr[brd][i], &reg);
		reg = (reg & ~WDcfg.GWmask[brd][i]) | (WDcfg.GWdata[brd][i] & WDcfg.GWmask[brd][i]);
		ret |= CAEN_DGTZ_WriteRegister(handle[brd], WDcfg.GWaddr[brd][i], reg);
	}
	if (ret) goto abortcfg;

	return 0;

	abortcfg:
	printf("Digitizer configuration failed on board %d:\n", brd);
	printf("Error at: %s. Exit Code = %d\n", CfgStep, ret);
	return ret;

}


// --------------------------------------------------------------------------------------------------------- 
// Description: Set virtual probe in the digitizer (analog and digital traces)
// Inputs:		brd = board number
// Outputs:		-
// Return:		Error code (0=success) 
// --------------------------------------------------------------------------------------------------------- 
int SetVirtualProbes(int brd)
{
	int ret=0;

	// Analog Probe 0
	if ((TraceSet[0] >= 0) && (TraceSet[0] < 4) && (TraceNames[0][TraceSet[0]][0] != '-')) {
		if (WDcfg.DppType == DPP_PSD_730)
			ret |= RegisterSetBits(handle[brd], 0x8000, 12, 12, TraceSet[0] & 0x1);
		else if (WDcfg.DppType == DPP_PHA_730)
			ret |= RegisterSetBits(handle[brd], 0x8000, 12, 13, TraceSet[0] & 0x3);
		else
			ret |= RegisterSetBits(handle[brd], 0x8000, 12, 13, TraceSet[0] & 0x3);
	}

	// Analog Probe 1
	if ((TraceSet[1] >= 0) && (TraceSet[0] < 4) && (TraceSet[1] < 4) && (TraceNames[1][TraceSet[1]][0] != '-')) {
		RegisterSetBits(handle[brd], 0x8000, 11, 11, 1);	// set dual trace
		if (WDcfg.DppType == DPP_PHA_724)
			ret |= RegisterSetBits(handle[brd], 0x8000, 14, 15, TraceSet[1] & 0x3);
		else if (WDcfg.DppType == DPP_PHA_730)
			ret |= RegisterSetBits(handle[brd], 0x8000, 14, 15, TraceSet[1] & 0x3);
		else
			ret |= RegisterSetBits(handle[brd], 0x8000, 13, 13, TraceSet[1] & 0x1);
	} else {
		RegisterSetBits(handle[brd], 0x8000, 11, 11, 0);	// unset dual trace
	}

	// Digital Probe 0
	if ((TraceSet[2] >= 0) && (TraceNames[2][TraceSet[2]][0] != '-')) {
		if (WDcfg.DigitizerModel == 730)
			ret |= RegisterSetBits(handle[brd], 0x8000, 23, 25, TraceSet[2] & 0x7);
		else
			ret |= RegisterSetBits(handle[brd], 0x8000, 20, 22, TraceSet[2] & 0x7);
	}

	// Digital Probe 1
	if ((TraceSet[3] >= 0) && (TraceNames[3][TraceSet[3]][0] != '-')) {
		if (WDcfg.DigitizerModel == 730)
			ret |= RegisterSetBits(handle[brd], 0x8000, 26, 28, TraceSet[3] & 0x7);
		else
			ret |= RegisterSetBits(handle[brd], 0x8000, 23, 25, TraceSet[3] & 0x7);
	}

	// Digital Probe 2
	if ((TraceSet[4] >= 0) && (TraceNames[4][TraceSet[4]][0] != '-')) {
		if (WDcfg.DigitizerModel == 720) 
			ret |= RegisterSetBits(handle[brd], 0x8000, 23, 25, TraceSet[4] & 0x7);
	}

	// Digital Probe 3
	if ((TraceSet[5] >= 0) && (TraceNames[4][TraceSet[4]][0] != '-')) {
		if (WDcfg.DigitizerModel == 720)
			ret |= RegisterSetBits(handle[brd], 0x8000, 26, 28, TraceSet[5] & 0x7);
	}

	return ret;
}



// --------------------------------------------------------------------------------------------------------- 
// Description: Set trace options and names for each type of digitizer and DPP
// Return:		Error code (0=success) 
// --------------------------------------------------------------------------------------------------------- 
int SetTraceNames()
{
	int i,j;
	for(i=0; i<MAX_NTRACES; i++) {
		for(j=0; j<8; j++)
			strcpy(TraceNames[i][j], "-");
	}

	switch(WDcfg.DppType) {
		// --------------------------------------------------------------------
		case DPP_PHA_724:
			TraceSet[0] = 3;
			TraceSet[1] = 0;
			TraceSet[2] = 0;
			TraceSet[3] = 0;
			TraceSet[4] = -1;
			TraceSet[5] = -1;
			strcpy(TraceNames[0][0], "Input");
			strcpy(TraceNames[0][1], "CR-RC");
			strcpy(TraceNames[0][2], "CR-RC2");
			strcpy(TraceNames[0][3], "Trapezoid");
			strcpy(TraceNames[1][0], "Input");
			strcpy(TraceNames[1][1], "Threshold");
			strcpy(TraceNames[1][2], "Trap-Baseline");
			strcpy(TraceNames[1][3], "Baseline");
			if (WDcfg.WaveformProcessor) {
				strcpy(TraceNames[0][7], "SW-WaveProcessor");
				strcpy(TraceNames[1][7], "SW-WaveProcessor");
			}


			strcpy(TraceNames[2][0], "Trigger");

			strcpy(TraceNames[3][0],  "ZCwindow");
			strcpy(TraceNames[3][1],  "Armed");
			strcpy(TraceNames[3][2],  "TrapRunning");
			strcpy(TraceNames[3][3],  "PileUp");
			strcpy(TraceNames[3][4],  "Peaking");
			strcpy(TraceNames[3][5],  "TrgValWindow");
			strcpy(TraceNames[3][6],  "BslHoldOff");
			strcpy(TraceNames[3][7],  "TrgHoldOff");
			strcpy(TraceNames[3][8],  "TrgValidation");
			strcpy(TraceNames[3][9],  "AcqVeto");
			strcpy(TraceNames[3][10], "BFMveto");
			strcpy(TraceNames[3][11], "ExtTrigger");
			break;

		// --------------------------------------------------------------------
		case DPP_CI:
			TraceSet[0] = 0;
			TraceSet[1] = -1;
			TraceSet[2] = 0;
			TraceSet[3] = 0;
			TraceSet[4] = 0;
			TraceSet[5] = 0;

			TraceSet[0] = 0;
			TraceSet[1] = -1;
			TraceSet[2] = 0;
			TraceSet[3] = 0;
			TraceSet[4] = 0;
			TraceSet[5] = 0;

			strcpy(TraceNames[0][0], "Input");
			strcpy(TraceNames[1][0], "Baseline");
			if (WDcfg.WaveformProcessor) {
				strcpy(TraceNames[0][7], "SW-WaveProcessor");
				strcpy(TraceNames[1][7], "SW-WaveProcessor");
			}


			strcpy(TraceNames[2][0], "Trigger");

			strcpy(TraceNames[3][0], "Gate");

			strcpy(TraceNames[4][0], "ExtTrg");
			strcpy(TraceNames[4][1], "OverThr");
			strcpy(TraceNames[4][2], "Trgout");
			strcpy(TraceNames[4][3], "TrgValWindow");
			strcpy(TraceNames[4][4], "PileUp");
			strcpy(TraceNames[4][5], "Coincidence");

			strcpy(TraceNames[5][1], "OverThr");
			strcpy(TraceNames[5][2], "TrgValidation");
			strcpy(TraceNames[5][3], "TrgHoldOff");
			strcpy(TraceNames[5][4], "PileUp");
			strcpy(TraceNames[5][5], "Coincidence");
			break;

		// --------------------------------------------------------------------
		case DPP_PSD_720:
			TraceSet[0] = 0;
			TraceSet[1] = -1;
			TraceSet[2] = 0;
			TraceSet[3] = 0;
			TraceSet[4] = 0;
			TraceSet[5] = 0;

			strcpy(TraceNames[0][0], "Input");
			strcpy(TraceNames[1][0], "Baseline");
			if (WDcfg.WaveformProcessor) {
				strcpy(TraceNames[0][7], "SW-WaveProcessor");
				strcpy(TraceNames[1][7], "SW-WaveProcessor");
			}


			strcpy(TraceNames[2][0], "Trigger");

			strcpy(TraceNames[3][0], "LongGate");

			strcpy(TraceNames[4][0], "ExtTrg");
			strcpy(TraceNames[4][1], "OverThr");
			strcpy(TraceNames[4][2], "Trgout");
			strcpy(TraceNames[4][3], "TrgValWindow");
			strcpy(TraceNames[4][4], "PileUp");
			strcpy(TraceNames[4][5], "Coincidence");

			strcpy(TraceNames[5][0], "ShortGate");
			strcpy(TraceNames[5][1], "OverThr");
			strcpy(TraceNames[5][2], "TrgValidation");
			strcpy(TraceNames[5][3], "TrgHoldOff");
			strcpy(TraceNames[5][4], "PileUp");
			strcpy(TraceNames[5][5], "Coincidence");
			break;

		// --------------------------------------------------------------------
		case DPP_PSD_751:
			TraceSet[0] = 0;
			TraceSet[1] = -1;
			TraceSet[2] = 0;
			TraceSet[3] = 0;
			TraceSet[4] = -1;
			TraceSet[5] = -1;

			strcpy(TraceNames[0][0], "Input");
			strcpy(TraceNames[1][0], "Baseline");
			if (WDcfg.WaveformProcessor) {
				strcpy(TraceNames[0][7], "SW-WaveProcessor");
				strcpy(TraceNames[1][7], "SW-WaveProcessor");
			}
			

			strcpy(TraceNames[2][0], "LongGate");
			strcpy(TraceNames[2][1], "OverThr");
			strcpy(TraceNames[2][2], "Trgout");
			strcpy(TraceNames[2][3], "TrgValWindow");
			strcpy(TraceNames[2][4], "PileUp");
			strcpy(TraceNames[2][5], "Coincidence");
			strcpy(TraceNames[2][6], "BslFreeze");

			strcpy(TraceNames[3][0], "ShortGate");
			strcpy(TraceNames[3][1], "OverThr");
			strcpy(TraceNames[3][2], "TrgValidation");
			strcpy(TraceNames[3][3], "TrgHoldOff");
			strcpy(TraceNames[3][4], "PileUp");
			strcpy(TraceNames[3][5], "Coincidence");
			break;

		// --------------------------------------------------------------------
		case DPP_PSD_730:
			TraceSet[0] = 0;
			TraceSet[1] = -1;
			TraceSet[2] = 0;
			TraceSet[3] = 0;
			TraceSet[4] = -1;
			TraceSet[5] = -1;

			strcpy(TraceNames[0][0], "Input");
			strcpy(TraceNames[0][1], "CFD");
			strcpy(TraceNames[1][0], "Baseline");
			strcpy(TraceNames[1][1], "CFD");
			if (WDcfg.WaveformProcessor) {
				strcpy(TraceNames[0][7], "SW-WaveProcessor");
				strcpy(TraceNames[1][7], "SW-WaveProcessor");
			}


			strcpy(TraceNames[2][0], "LongGate");
			strcpy(TraceNames[2][1], "OverThr");
			strcpy(TraceNames[2][2], "Trgout");
			strcpy(TraceNames[2][3], "TrgValWindow");
			strcpy(TraceNames[2][4], "PileUp");
			strcpy(TraceNames[2][5], "Coincidence");
			strcpy(TraceNames[2][6], "BslFreeze");
			strcpy(TraceNames[2][7], "Trigger");

			strcpy(TraceNames[3][0], "ShortGate");
			strcpy(TraceNames[3][1], "OverThr");
			strcpy(TraceNames[3][2], "TrgValidation");
			strcpy(TraceNames[3][3], "TrgHoldOff");
			strcpy(TraceNames[3][4], "PileUp");
			strcpy(TraceNames[3][5], "Coincidence");
			strcpy(TraceNames[3][7], "Trigger");
			break;

		// --------------------------------------------------------------------
		case DPP_PHA_730:
			TraceSet[0] = 0;
			TraceSet[1] = -1;
			TraceSet[2] = 0;
			TraceSet[3] = 0;
			TraceSet[4] = -1;
			TraceSet[5] = -1;

			strcpy(TraceNames[0][0], "Input");
			strcpy(TraceNames[0][1], "CR-RC2");
			strcpy(TraceNames[0][2], "CR-RC2");
			strcpy(TraceNames[0][3], "Trapezoid");
			strcpy(TraceNames[1][0], "Vin");
			strcpy(TraceNames[1][1], "Threshold");
			strcpy(TraceNames[1][2], "Trapez-Baseline");
			strcpy(TraceNames[1][3], "Baseline");
			if (WDcfg.WaveformProcessor) {
				strcpy(TraceNames[0][7], "SW-WaveProcessor");
				strcpy(TraceNames[1][7], "SW-WaveProcessor");
			}


			strcpy(TraceNames[2][0], "Trigger");

			strcpy(TraceNames[3][0], "DaFare");
			strcpy(TraceNames[3][1], "TTF_probe");
			strcpy(TraceNames[3][2], "PkRun");
			strcpy(TraceNames[3][3], "PileUp");
			strcpy(TraceNames[3][4], "Peaking");

			break;

		default:
			break;
	}
	return 0;
}
