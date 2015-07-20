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


#include <CAENDigitizer.h>
#include "digiTES.h"


int load_default=1;
int ch=-1, brd=-1;

int GetInt(FILE *f_ini)
{
	int ret;
	fscanf(f_ini, "%d", &ret);
	return ret;
}
int GetHex(FILE *f_ini)
{
	int ret;
	fscanf(f_ini, "%x", &ret);
	return ret;
}
float GetFloat(FILE *f_ini)
{
	float ret;
	fscanf(f_ini, "%f", &ret);
	return ret;
}

void SetChannelParam(int param[MAX_NBRD][MAX_NCH], int val)
{
	int i, b;
    if (ch == -1)
        for(b=0; b<MAX_NBRD; b++)
			for(i=0; i<MAX_NCH; i++)
				param[b][i] = val;
    else
        param[brd][ch] = val;
}

void SetChannelParamFloat(float param[MAX_NBRD][MAX_NCH], float val)
{
	int i, b;
    if (ch == -1)
        for(b=0; b<MAX_NBRD; b++)
			for(i=0; i<MAX_NCH; i++)
				param[b][i] = val;
    else
        param[brd][ch] = val;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: Read a config file, parse the parameters and set the relevant fields in the WDcfg structure
// Inputs:		f_ini: config file pinter
// Outputs:		WDcfg: struct with all parameters
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ParseConfigFile(FILE *f_ini, DPP_Config_t *WDcfg) 
{
	char str[1000], str1[1000];
	int i, b, val, Off=0, tr = -1;
	FILE *runnum, *fcfg;

	if (load_default) {
		/* Default settings */
		strcpy(WDcfg->InputDataFileName, "");
		strcpy(WDcfg->DataFilePath, ".\\DataFiles\\");
		WDcfg->ReplayDataFile = 0;
		WDcfg->NumBrd = 0;
		WDcfg->NumCh = 8;
		WDcfg->RunNumber = 0;
		WDcfg->OutFileFormat = OUTFILE_ASCII;
		WDcfg->HistoOutputFormat = HISTO_FILE_FORMAT_1COL;
		WDcfg->RunNumberInDataFiles = 0;
		WDcfg->AcquisitionMode = ACQ_MODE_MIXED;
		WDcfg->TriggerMode = TRIGGER_MODE_SELF;
		WDcfg->TrgoutMode = TRGOUT_MODE_CH_TRG;
		WDcfg->StartMode = START_MODE_INDEP_SW;
		WDcfg->StopOnTotalEvents = 0;
		WDcfg->StopOnEnergyEvents = 0;
		WDcfg->StopOnTimeEvents = 0;
		WDcfg->StopOnTime = 0;
		WDcfg->FPIOtype = CAEN_DGTZ_IOLevel_NIM;
		WDcfg->RecordLength = (1024*16);
		WDcfg->PreTrigger = 128;
		WDcfg->TrgHoldOff = 1;
		WDcfg->EnableCoinc = 0;
		WDcfg->CoincWindow = 10;
		WDcfg->FanSpeed = 0;
		WDcfg->EnableTimeCorrelFilter = 0;
		WDcfg->EnableEnergyFilter = 0;
		WDcfg->EnablePSDFilter = 0;
		WDcfg->WaveformProcessor = 0;
		WDcfg->TimeCorrelWindow = 100;
		WDcfg->GWn = 0;
		WDcfg->THnbin = TMAXNBITS;
		WDcfg->EHnbin = EMAXNBITS;
		WDcfg->MCSHnbin = 1024;
		WDcfg->DwellTime = 1000;
        for(b=0; b<MAX_NBRD; b++) {
			for(i=0; i<MAX_NCH; i++) {
				WDcfg->DCoffset[b][i]=8000;  
				WDcfg->EnableInput[b][i]=1;  
				WDcfg->TrgThreshold[b][i]=10;  
				WDcfg->PulsePolarity[b][i]=1;  
				WDcfg->NsBaseline[b][i]=2;  
				WDcfg->FixedBaseline[b][i]=0;  
				WDcfg->GateWidth[b][i]=20;  
				WDcfg->ShortGateWidth[b][i]=5;  
				WDcfg->PreGate[b][i]=2;  
				WDcfg->ChargeSensitivity[b][i]=0;  
				WDcfg->EnablePedestal[b][i]=0;
				WDcfg->InputDynamicRange[b][i]=0;  
				WDcfg->CFDatten[b][i]=0;  
				WDcfg->CFDdelay[b][i]=0;  
				WDcfg->CFDinterp[b][i]=0;  
				WDcfg->DiscrMode[b][i]=0;  
				WDcfg->EnableIPE[b][i]=0;  
				WDcfg->IPEaddnoise[b][i]=1;
				WDcfg->IPErandom[b][i]=1;
				WDcfg->IPEamplitude[b][i]=1000;
				WDcfg->IPEdecay[b][i]=50;
				WDcfg->IPEfrequency[b][i]=10;
				WDcfg->IPErisetime[b][i]=2;
				WDcfg->PSDcut[b][i]=0;
				WDcfg->PileUpMode[b][i]=0;
				WDcfg->PurGap[b][i]=10;
				WDcfg->TrapRiseTime[b][i]=5000;
				WDcfg->TrapFlatTop[b][i]=1000;
				WDcfg->TrapPoleZero[b][i]=50000;
				WDcfg->PeakingTime[b][i]=500;
				WDcfg->PeakHoldOff[b][i]=1000;
				WDcfg->TTFsmoothing[b][i]=8;
				WDcfg->TTFdelay[b][i]=100;
				WDcfg->Decimation[b][i]=0;
				WDcfg->TrapDigitalGain[b][i]=1.0;

				WDcfg->EnergyLLD[b][i]=0;  
				WDcfg->EnergyULD[b][i]=0;  
				WDcfg->PsdLLD[b][i]=0;  
				WDcfg->PsdULD[b][i]=0;  
				WDcfg->EnergyGain[b][i]=1.0;
				WDcfg->EnergyOffset[b][i]=0;
				WDcfg->ECalibration_m[b][i]=1.0;
				WDcfg->ECalibration_q[b][i]=0;
				WDcfg->DelayLine[b][i]=0;  
				WDcfg->CRboard[b][i]=0;  
				WDcfg->CRchannel[b][i]=0;

				WDcfg->THmin[b][i]=0;  
				WDcfg->THmax[b][i]=TMAXNBITS;  
				WDcfg->EHmin[b][i]=0;  
				WDcfg->EHmax[b][i]=EMAXNBITS;

				WDcfg->CoincMode[b][i]=0;  
				WDcfg->CoincMask[b][i]=0;  
			}
		}
	}

	// append cfg file to _cfg.txt
	if (load_default)
		fcfg = fopen("_cfg.txt", "w");
	else
		fcfg = fopen("_cfg.txt", "a");
	while ((fcfg != NULL) && (!feof(f_ini))) {
		char line[1000], str[100];
		fgets(line, 1000, f_ini);
		sscanf(line, "%s", str);
        if ((str[0] != '#') && (strcmp(str, "Load")!=0))
			fputs(line, fcfg);
	}
	fclose(fcfg);
	rewind(f_ini);


	// Read Run Number from file
	if (load_default) {
		runnum = fopen("RunNumber.txt", "r");
		if (runnum != NULL)   {
			fscanf(runnum, "%d", &WDcfg->RunNumber);
			fclose(runnum);
		}
		runnum = fopen("RunNumber.txt", "w");
		fprintf(runnum, "%d", WDcfg->RunNumber+1);
		fclose(runnum);
	}


	/* read config file and assign parameters */
	while(!feof(f_ini)) {
		int read;
        // read a word from the file
        read = fscanf(f_ini, "%s", str);
        if( !read || (read == EOF) || !strlen(str))
			continue;

        // skip comments
        if (str[0] == '#') {
			fgets(str, 1000, f_ini);
			continue;
		}

		if (strcmp(str, "@ON")==0) {
			Off = 0;
			continue;
		} else if (strcmp(str, "@OFF")==0) {
			Off = 1;
			continue;
		}
        if (Off)
            continue;

        // Section (COMMON or individual channel)
		if (str[0] == '[')								{	if (strstr(str, "COMMON")!=NULL) {
																ch = -1;
																brd = -1;
															} else if (strstr(str, "BOARD")!=NULL) {
																ch = -1;
																fscanf(f_ini, "%s", str1);
																sscanf(str1, "%d", &val);
																if (val < 0 || val >= MAX_NBRD) printf("%s: Invalid board number\n", str);
																else brd = val;
															} else if (strstr(str, "CHANNEL")!=NULL) {
																fscanf(f_ini, "%s", str1);
																sscanf(str1, "%d", &val);
																if (val < 0 || val >= MAX_NCH) printf("%s: Invalid channel number\n", str);
																else ch = val;
															}

        // LOAD: open another config file
		} else if (strcmp(str, "Load")==0)				{	FILE *cf;
															fscanf(f_ini, "%s", str1);
															cf = fopen(str1, "r");
															if (cf!=NULL) {
																load_default=0;
																ParseConfigFile(cf, WDcfg);
																fclose(cf);
																load_default=1;
															} else {
																printf("Can't open secondary config file %s\n", str1);
															}
 
		} else if (strcmp(str, "Open")==0)				{	if (brd==-1) {
																printf("%s: cannot be a common setting (must be in a [BOARD] section)\n", str); 
																fgets(str1, 100, f_ini);
															} else {
																fscanf(f_ini, "%s", str1);
																if	(strcmp(str1, "USB")==0)		WDcfg->LinkType[brd] = CAEN_DGTZ_USB;
																else if (strcmp(str1, "PCI")==0)	WDcfg->LinkType[brd] = CAEN_DGTZ_PCI_OpticalLink;
																else 	printf("%s: invalid setting for %s\n", str1, str);

																WDcfg->LinkNum[brd] = GetInt(f_ini);
																if (WDcfg->LinkType[brd] == CAEN_DGTZ_USB) WDcfg->ConetNode[brd] = 0;
																else WDcfg->ConetNode[brd] = GetInt(f_ini);
																WDcfg->BaseAddress[brd] = GetHex(f_ini);
																WDcfg->NumBrd++;
															}

		} else if (strcmp(str, "WriteRegister")==0)		{	for(i=0; i<MAX_NBRD-1; i++) {
																if ((brd == -1) || (i == brd)) {
																	if (WDcfg->GWn < MAX_GW) {
																		fscanf(f_ini, "%x", (int *)&WDcfg->GWaddr[i][WDcfg->GWn]);
																		fscanf(f_ini, "%x", (int *)&WDcfg->GWdata[i][WDcfg->GWn]);
																		fscanf(f_ini, "%x", (int *)&WDcfg->GWmask[i][WDcfg->GWn]);
																		WDcfg->GWn++;
																	} else {
																		printf("MAX_GW Generic Write exceeded (%d). Change MAX_GW and recompile\n", MAX_GW);
																	}
																}
															}

		} else if (strcmp(str, "InputDataFileName")==0)		{	fscanf(f_ini, "%s", WDcfg->InputDataFileName); WDcfg->ReplayDataFile = 1;
		} else if (strcmp(str, "DataFilePath")==0)			{	fscanf(f_ini, "%s", WDcfg->DataFilePath);
		} else if (strcmp(str, "RunNumberInDataFiles")==0)	{	WDcfg->RunNumberInDataFiles		= GetInt(f_ini);
		} else if (strcmp(str, "RecordLength")==0)			{	WDcfg->RecordLength				= GetInt(f_ini);
		} else if (strcmp(str, "PreTrigger")==0) 			{	WDcfg->PreTrigger				= GetInt(f_ini); 
		} else if (strcmp(str, "TrgHoldOff")==0) 			{	WDcfg->TrgHoldOff				= GetInt(f_ini);
		} else if (strcmp(str, "StopOnTime")==0) 			{	WDcfg->StopOnTime				= GetInt(f_ini);
		} else if (strcmp(str, "StopOnTotalEvents")==0)		{	WDcfg->StopOnTotalEvents		= GetInt(f_ini);
		} else if (strcmp(str, "StopOnTimeEvents")==0)		{	WDcfg->StopOnTimeEvents			= GetInt(f_ini);
		} else if (strcmp(str, "StopOnEnergyEvents")==0)	{	WDcfg->StopOnEnergyEvents		= GetInt(f_ini);
		} else if (strcmp(str, "THnbin")==0)				{	WDcfg->THnbin					= GetInt(f_ini);
		} else if (strcmp(str, "EHnbin")==0)				{	WDcfg->EHnbin					= GetInt(f_ini);
		} else if (strcmp(str, "MCSHnbin")==0)				{	WDcfg->MCSHnbin					= GetInt(f_ini);
		} else if (strcmp(str, "DwellTime")==0)				{	WDcfg->DwellTime				= GetInt(f_ini);
		} else if (strcmp(str, "EnableTimeCorrelFilter")==0){	WDcfg->EnableTimeCorrelFilter	= GetInt(f_ini);
		} else if (strcmp(str, "EnableEnergyFilter")==0)	{	WDcfg->EnableEnergyFilter		= GetInt(f_ini);
		} else if (strcmp(str, "EnablePSDFilter")==0)		{	WDcfg->EnablePSDFilter			= GetInt(f_ini);
		} else if (strcmp(str, "TimeCorrelWindow")==0) 		{	WDcfg->TimeCorrelWindow			= GetInt(f_ini);
		} else if (strcmp(str, "WaveformProcessor")==0) 	{	WDcfg->WaveformProcessor		= GetInt(f_ini);
		} else if (strcmp(str, "FanSpeed")==0)				{	WDcfg->FanSpeed					= GetInt(f_ini);
		} else if (strcmp(str, "TOFstartChannel")==0)		{	WDcfg->TOFstartChannel			= GetInt(f_ini);
		} else if (strcmp(str, "TOFstartBoard")==0)			{	WDcfg->TOFstartBoard			= GetInt(f_ini);
		} else if (strcmp(str, "SaveLists")==0)				{	WDcfg->SaveLists				= GetInt(f_ini);


		} else if (strcmp(str, "CoincWindow")==0)		{	WDcfg->CoincWindow		= GetInt(f_ini);
		} else if (strcmp(str, "CoincMask")==0)			{	int ch, c;
															int err=0;
															char str1[100];
			                                                fscanf(f_ini, "%d", &ch);
															if ((ch>=0) && (ch<9)) {
																WDcfg->CoincMask[brd][ch] = 0;
																for(i=0; i<9; i++) {
																	fscanf(f_ini, "%d", &c);
																	if (i == 8) WDcfg->CoincMask[brd][ch] |= (c<<30);
																	else        WDcfg->CoincMask[brd][ch] |= (c<<i);
																}
																fscanf(f_ini, "%s", &str1);
																if      (strcmp(str1, "AND")==0)	WDcfg->CoincMask[brd][ch] |= (1<<8);
																else if (strcmp(str1, "OR")==0)		WDcfg->CoincMask[brd][ch] |= (0<<8);
																else if (str1[0] == 'M')			WDcfg->CoincMask[brd][ch] |= (2<<8) | ((str1[1]-'0') << 10);
																else 								err = 1;
																if (ch < 8) {
																	fscanf(f_ini, "%s", &str1);
																	if (str1[0] == 'C')	WDcfg->CoincMode[brd][ch] = 1;
																	else				WDcfg->CoincMode[brd][ch] = 2;
																}
															} else {
																err = 1;
															}
															if (err) {
																WDcfg->CoincMask[brd][ch] = 0;
																WDcfg->CoincMode[brd][ch] = 0;
																printf("%s: invalid setting for %s\n", str1, str);
															}
		} else if (strcmp(str, "EnableCoinc")==0)		{	WDcfg->EnableCoinc		= GetInt(f_ini);
		} else if (strcmp(str, "EnableInput")==0)		{	SetChannelParam(WDcfg->EnableInput,			GetInt(f_ini));
		} else if (strcmp(str, "DCoffset")==0)			{	SetChannelParam(WDcfg->DCoffset,			(int)((GetFloat(f_ini)+50) * 65535 / 100));
		} else if (strcmp(str, "TriggerThreshold")==0)	{	SetChannelParam(WDcfg->TrgThreshold,		GetInt(f_ini));
		} else if (strcmp(str, "NSBaseline")==0) 		{	SetChannelParam(WDcfg->NsBaseline,			GetInt(f_ini));
		} else if (strcmp(str, "FixedBaseline")==0) 	{	SetChannelParam(WDcfg->FixedBaseline,		GetInt(f_ini));
		} else if (strcmp(str, "PreGate")==0) 			{	SetChannelParam(WDcfg->PreGate,				GetInt(f_ini));
		} else if (strcmp(str, "ChargeSensitivity")==0)	{	SetChannelParam(WDcfg->ChargeSensitivity,	GetInt(f_ini));
		} else if (strcmp(str, "InputDynamicRange")==0)	{	SetChannelParam(WDcfg->InputDynamicRange,	GetInt(f_ini));
		} else if (strcmp(str, "GateWidth")==0)			{	SetChannelParam(WDcfg->GateWidth,			GetInt(f_ini));
		} else if (strcmp(str, "ShortGateWidth")==0)	{	SetChannelParam(WDcfg->ShortGateWidth,		GetInt(f_ini));
		} else if (strcmp(str, "EnablePedestal")==0)	{	SetChannelParam(WDcfg->EnablePedestal,		GetInt(f_ini));
		} else if (strcmp(str, "TrapRiseTime")==0)		{	SetChannelParam(WDcfg->TrapRiseTime,		GetInt(f_ini));
		} else if (strcmp(str, "TrapFlatTop")==0)		{	SetChannelParam(WDcfg->TrapFlatTop,			GetInt(f_ini));
		} else if (strcmp(str, "TrapPoleZero")==0)		{	SetChannelParam(WDcfg->TrapPoleZero,		GetInt(f_ini));
		} else if (strcmp(str, "PeakingTime")==0)		{	SetChannelParam(WDcfg->PeakingTime,			GetInt(f_ini));
		} else if (strcmp(str, "PeakHoldOff")==0)		{	SetChannelParam(WDcfg->PeakHoldOff,			GetInt(f_ini));
		} else if (strcmp(str, "Decimation")==0)		{	SetChannelParam(WDcfg->Decimation,			GetInt(f_ini));
		} else if (strcmp(str, "TTFsmoothing")==0)		{	SetChannelParam(WDcfg->TTFsmoothing,		GetInt(f_ini));
		} else if (strcmp(str, "TTFdelay")==0)			{	SetChannelParam(WDcfg->TTFdelay,			GetInt(f_ini));
		} else if (strcmp(str, "TrapDigitalGain")==0)	{	SetChannelParamFloat(WDcfg->TrapDigitalGain,GetFloat(f_ini));
		} else if (strcmp(str, "PileUpMode")==0)		{	SetChannelParam(WDcfg->PileUpMode,			GetInt(f_ini));
		} else if (strcmp(str, "PurGap")==0)			{	SetChannelParam(WDcfg->PurGap,				GetInt(f_ini));
		} else if (strcmp(str, "CRboard")==0)			{	SetChannelParam(WDcfg->CRboard,				GetInt(f_ini));
		} else if (strcmp(str, "CRchannel")==0)			{	SetChannelParam(WDcfg->CRchannel,			GetInt(f_ini));
		} else if (strcmp(str, "EnergyLLD")==0)			{	SetChannelParamFloat(WDcfg->EnergyLLD,		GetFloat(f_ini));
		} else if (strcmp(str, "EnergyULD")==0)			{	SetChannelParamFloat(WDcfg->EnergyULD,		GetFloat(f_ini));
		} else if (strcmp(str, "PsdLLD")==0)			{	SetChannelParamFloat(WDcfg->PsdLLD,			GetFloat(f_ini));
		} else if (strcmp(str, "PsdULD")==0)			{	SetChannelParamFloat(WDcfg->PsdULD,			GetFloat(f_ini));
		} else if (strcmp(str, "DelayLine")==0)			{	SetChannelParamFloat(WDcfg->DelayLine,		GetFloat(f_ini));
		} else if (strcmp(str, "THmin")==0)				{	SetChannelParamFloat(WDcfg->THmin,			GetFloat(f_ini));
		} else if (strcmp(str, "THmax")==0)				{	SetChannelParamFloat(WDcfg->THmax,			GetFloat(f_ini));
		} else if (strcmp(str, "EHmin")==0)				{	SetChannelParamFloat(WDcfg->EHmin,			GetFloat(f_ini));
		} else if (strcmp(str, "EHmax")==0)				{	SetChannelParamFloat(WDcfg->EHmax,			GetFloat(f_ini));
		} else if (strcmp(str, "EnergyGain")==0)		{	SetChannelParamFloat(WDcfg->EnergyGain,		GetFloat(f_ini));
		} else if (strcmp(str, "EnergyOffset")==0)		{	SetChannelParamFloat(WDcfg->EnergyOffset,	GetFloat(f_ini));
		} else if (strcmp(str, "ECalibration_m")==0)	{	SetChannelParamFloat(WDcfg->ECalibration_m,	GetFloat(f_ini));
		} else if (strcmp(str, "ECalibration_q")==0)	{	SetChannelParamFloat(WDcfg->ECalibration_q,	GetFloat(f_ini));
		} else if (strcmp(str, "CFDdelay")==0)			{	SetChannelParam(WDcfg->CFDdelay,			GetInt(f_ini));
		} else if (strcmp(str, "CFDatten")==0)			{	SetChannelParam(WDcfg->CFDatten,			GetInt(f_ini));
		} else if (strcmp(str, "CFDinterp")==0)			{	SetChannelParam(WDcfg->CFDinterp,			GetInt(f_ini));
		} else if (strcmp(str, "DiscrMode")==0)			{	SetChannelParam(WDcfg->DiscrMode,			GetInt(f_ini));
		} else if (strcmp(str, "PSDcut")==0)			{	SetChannelParam(WDcfg->PSDcut,				GetInt(f_ini));
		} else if (strcmp(str, "EnableIPE")==0)			{	SetChannelParam(WDcfg->EnableIPE,			GetInt(f_ini));
		} else if (strcmp(str, "IPEamplitude")==0)		{	SetChannelParam(WDcfg->IPEamplitude,		GetInt(f_ini));
		} else if (strcmp(str, "IPErandom")==0)			{	SetChannelParam(WDcfg->IPErandom,			GetInt(f_ini));
		} else if (strcmp(str, "IPEaddnoise")==0)		{	SetChannelParam(WDcfg->IPEaddnoise,			GetInt(f_ini));
		} else if (strcmp(str, "IPErisetime")==0)		{	SetChannelParam(WDcfg->IPErisetime,			GetInt(f_ini));
		} else if (strcmp(str, "IPEfrequency")==0)		{	SetChannelParamFloat(WDcfg->IPEfrequency,	GetFloat(f_ini));
		} else if (strcmp(str, "IPEdecay")==0)			{	SetChannelParamFloat(WDcfg->IPEdecay,		GetFloat(f_ini));
		} else if (strcmp(str, "PulsePolarity")==0)		{	fscanf(f_ini, "%s", str1);
															if		(strcmp(str1, "NEGATIVE")==0)		SetChannelParam(WDcfg->PulsePolarity, CAEN_DGTZ_PulsePolarityNegative);
															else if (strcmp(str1, "POSITIVE")==0)		SetChannelParam(WDcfg->PulsePolarity, CAEN_DGTZ_PulsePolarityPositive);
															else 	printf("%s: invalid setting for %s\n", str1, str);
		} else if (strcmp(str, "AcquisitionMode")==0)	{	fscanf(f_ini, "%s", str1);
															if		(strcmp(str1, "LIST")==0)			WDcfg->AcquisitionMode = CAEN_DGTZ_DPP_ACQ_MODE_List;
															else if (strcmp(str1, "OSCILLOSCOPE")==0)	WDcfg->AcquisitionMode = CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope;
															else if (strcmp(str1, "MIXED")==0)			WDcfg->AcquisitionMode = CAEN_DGTZ_DPP_ACQ_MODE_Mixed;
															else 	printf("%s: invalid setting for %s\n", str1, str);
		} else if (strcmp(str, "StartMode")==0)			{	fscanf(f_ini, "%s", str1);
															if		(strcmp(str1, "INDEP_SW")==0)		WDcfg->StartMode = START_MODE_INDEP_SW;
															else if (strcmp(str1, "SYNC_SW")==0)		WDcfg->StartMode = START_MODE_SYNC_SW;
															else if (strcmp(str1, "SYNC_HW")==0) 		WDcfg->StartMode = START_MODE_SYNC_HW;
															else 	printf("%s: invalid setting for %s\n", str1, str);
		} else if (strcmp(str, "TriggerMode")==0)		{	fscanf(f_ini, "%s", str1);
															if		(strcmp(str1, "SELF_TRIGGER")==0)	WDcfg->TriggerMode = TRIGGER_MODE_SELF;
															else if (strcmp(str1, "EXT_TRIGGER")==0) 	WDcfg->TriggerMode = TRIGGER_MODE_EXTERNAL;
															else 	printf("%s: invalid setting for %s\n", str1, str);
		} else if (strcmp(str, "TrgoutMode")==0)		{	fscanf(f_ini, "%s", str1);
															if		(strcmp(str1, "CHANNEL_TRIGGERS")==0)	WDcfg->TriggerMode = TRGOUT_MODE_CH_TRG;
															else if (strcmp(str1, "SQR_WAVE_1KHZ")==0) 		WDcfg->TriggerMode = TRGOUT_MODE_SQR_1KHZ;
															else if (strcmp(str1, "PULSES_1KHZ")==0) 		WDcfg->TriggerMode = TRGOUT_MODE_PLS_1KHZ;
															else if (strcmp(str1, "SQR_WAVE_10KHZ")==0) 	WDcfg->TriggerMode = TRGOUT_MODE_SQR_10KHZ;
															else if (strcmp(str1, "PULSES_10KHZ")==0) 		WDcfg->TriggerMode = TRGOUT_MODE_PLS_10KHZ;
															else 	printf("%s: invalid setting for %s\n", str1, str);
		} else if (strcmp(str, "OutFileFormat")==0)		{	fscanf(f_ini, "%s", str1);
															if		(strcmp(str1, "BINARY")==0)			WDcfg->OutFileFormat = OUTFILE_BINARY;
															else if (strcmp(str1, "ASCII")==0)			WDcfg->OutFileFormat = OUTFILE_ASCII;
															else	printf("%s: invalid setting for %s\n", str1, str);
		} else if (strcmp(str, "HistoOutputFormat")==0)	{	fscanf(f_ini, "%s", str1);
															if		(strcmp(str1, "ASCII_1COL")==0)		WDcfg->HistoOutputFormat = HISTO_FILE_FORMAT_1COL;
															else if	(strcmp(str1, "ASCII_2COL")==0)		WDcfg->HistoOutputFormat = HISTO_FILE_FORMAT_2COL;
															else if (strcmp(str1, "ANSI_42")==0)		WDcfg->HistoOutputFormat = HISTO_FILE_FORMAT_ANSI42;
															else	printf("%s: invalid setting for %s\n", str1, str);
		} else if (strcmp(str, "FPIOtype")==0)			{	fscanf(f_ini, "%s", str1);
															if      (strcmp(str1, "TTL")==0)			WDcfg->FPIOtype = CAEN_DGTZ_IOLevel_TTL;
															else if (strcmp(str1, "NIM")==0)			WDcfg->FPIOtype = CAEN_DGTZ_IOLevel_NIM;
															else 	printf("%s: invalid setting for %s\n", str1, str);
        } else if (strcmp(str, "NamedPipe") == 0) {
			                                            fscanf(f_ini, "%s", WDcfg->pipeBaseName);
#ifndef WIN32
			fprintf(stderr, "*** Note the NamedPipe parameter is ignored for non windows operating systems\n");
#endif
		} else											{	printf("%s: unkwown parameter\n", str);
															fgets(str, (int)strlen(str)-1, f_ini);
		}

	}

#ifdef LINUX
	if (WDcfg->DataFilePath[strlen(WDcfg->DataFilePath)-1] != '/')	sprintf(WDcfg->DataFilePath, "%s/", WDcfg->DataFilePath);
#else
	if (WDcfg->DataFilePath[strlen(WDcfg->DataFilePath)-1] != '\\')	
		sprintf(WDcfg->DataFilePath, "%s\\", WDcfg->DataFilePath);
#endif


	return 0;
}


