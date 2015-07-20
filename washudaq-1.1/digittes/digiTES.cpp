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

/* To do list (15/10/14):
-------------------------------------------------------------------------------
- coincidence HW managed by the excel GUI (only between couples)
- option for external VETO/GATE (using TRG-IN)
. use extended time stamp given by the board (for those that have this option)
- read data until all baords are empty after a stop of acquisition 
- Compatibility with std firmware (non dpp)
- clear time stamps when resetting statistics (command "r" of the menu)
- read list files and raw data (also from Excel GUI)
- bindkay 'y' as autoscale on y-axix in gnuplot
- rise-time histogram (in DPP_PHA FW?)
- Options (channel by channel) for saving list and/or which histograms 
- automatic histo saving every N seconds (N programmable)
- calculate dead time and print it on the screen for each channel

-------------------------------------------------------------------------------
 Known bugs
-------------------------------------------------------------------------------
- Raw Data file saving not working

*/


#define digiTES_Revision    "2.94"  


#include <CAENDigitizer.h>

#include "digiTES.h"
#include "keyb.h"
#include "evb.h"
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#ifdef LINUX
	#define Sleep(x) usleep((x)*1000)
#endif

#ifdef WIN32 
#include "pipefunctions.h"
#endif

int DoStartSvc();

/* ###########################################################################
*  Global Variables
*  ########################################################################### */
DPP_Config_t	WDcfg;		// struct containing all acquisition parameters
Stats_t			Stats;		// struct containing variables for the statistics

int Quit=0;
int RestartAll=0;
int AcqRun=0;
int ChToPlot=0, BrdToPlot=0;
int ContinuousTrigger=0;
int SingleTrigger=0;
int EnableInternalPulseEmulator=0;
int DoRefresh=1;
int DoRefreshSingle=0;
int DoPlotWave=0;
int DoPlotHisto=0;
int DoSaveRawEvents=0;
int keVscale=0;
int handle[MAX_NBRD];            // board handles (for the CAEN_DGTZ library)
int TraceSet[MAX_NTRACES];
char TraceNames[MAX_NTRACES][MAX_NTRSETS][20];

FILE *hplot=NULL, *wplot=NULL;   // gnuplot pipes
FILE *cntlog=NULL;

#ifdef WIN32
HANDLE pipeHandle = INVALID_HANDLE_VALUE;
#endif

/* ###########################################################################
*  Functions
*  ########################################################################### */
int ClosePipeFile()
{
#ifdef WIN32
	if (pipeHandle != INVALID_HANDLE_VALUE) {
		ClosePipe(pipeHandle);
		pipeHandle = INVALID_HANDLE_VALUE;
	}
	return 1;
#else
	return 0;
#endif
}
static int WriteEventToPipe(int nFragments, pFragment_t pFrags)
{
#ifdef WIN32
	int i = sizeof(Fragment_t);
    int j = sizeof(GenericDPPEvent_t);

	if (pipeHandle == INVALID_HANDLE_VALUE) {
		pipeHandle = CreateServerPipe(WDcfg.pipeBaseName);
		if (pipeHandle == INVALID_HANDLE_VALUE) {
			return 0;
		}
		WaitForClient(pipeHandle);
	}
	/* Write the number of fragments and then all of the fragments  stop on any error */
	if (WriteToPipe(pipeHandle, &nFragments, sizeof(int)) == 0) return 0;
	for (i = 0; i < nFragments; i++) {
		if (WriteToPipe(pipeHandle, pFrags, sizeof(Fragment_t)) == 0) return 0;
		pFrags++;
	}
	return 1;                             /* Success */
#else
	return 0;                             
#endif
}

//--------------------------------------------------------------------------------------------------------
// Description: Fold the timestamp and fine time to get a picosecond trigger timestamp:
// Inputs:     builtEvents array of built events
//             f           fragment number we care about.

uint64_t ComputeFineTime(pFragment_t builtEvents, int f)
{
	uint64_t coarse = builtEvents[f].fragment.TimeStamp;    //ns
	uint64_t fine = builtEvents[f].fragment.FineTimeStamp;  //ps
	coarse = coarse * 1000;
	return (coarse  + fine);
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Read Board Information
// Inputs:		b=board index
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ReadBoardInfo(int b)
{
	int ret;
	uint32_t d32;
	int MajorNumber, MinorNumber;
	CAEN_DGTZ_BoardInfo_t BoardInfo;			// struct with the board information

	fprintf(cntlog, "Success, ");
	/* Once we have the handler to the digitizer, we use it to call the other functions */
	ret = CAEN_DGTZ_GetInfo(handle[b], &BoardInfo);
	if (ret) {
		printf("Can't read board info\n");
		return -1;
	}

	WDcfg.NumCh = BoardInfo.Channels;

	fprintf(cntlog, "%s, ", BoardInfo.ModelName);
	fprintf(cntlog, "%d, ", BoardInfo.Channels);
	fprintf(cntlog, "%d, ", BoardInfo.SerialNumber);
	if (BoardInfo.FamilyCode == 5) {
		WDcfg.DigitizerModel = 751;
		WDcfg.Tsampl = 1;
		WDcfg.Nbit = 10;
	} else if (BoardInfo.FamilyCode == 0) {
		WDcfg.DigitizerModel = 724;
		WDcfg.Tsampl = 10;
		WDcfg.Nbit = 14;
	} else if (BoardInfo.FamilyCode == 11) {
		WDcfg.DigitizerModel = 730;
		WDcfg.Tsampl = 2;
		WDcfg.Nbit = 14;
	} else {
		WDcfg.DigitizerModel = 720;
		WDcfg.Tsampl = 4;
		WDcfg.Nbit = 12;
	}



	/* Check firmware revision (only DPP firmware can be used with this Demo) */
	sscanf(BoardInfo.AMC_FirmwareRel, "%d", &MajorNumber);
	sscanf(&BoardInfo.AMC_FirmwareRel[4], "%d", &MinorNumber);
	if (MajorNumber == 128) {
		printf("This digitizer has a DPP_PHA firmware\n");
		WDcfg.DppType = DPP_PHA_724;
		fprintf(cntlog, "DPP_PHA_724, ");
	} else if (MajorNumber == 130) {
		printf("This digitizer has a DPP_CI firmware\n");
		WDcfg.DppType = DPP_CI;
		fprintf(cntlog, "DPP_CI_720, ");
	} else if (MajorNumber == 131) {
		printf("This digitizer has a DPP_PSD firmware\n");
		WDcfg.DppType = DPP_PSD_720;
		fprintf(cntlog, "DPP_PSD_720, ");
	} else if (MajorNumber == 132) {
		printf("This digitizer has a DPP_PSD firmware\n");
		WDcfg.DppType = DPP_PSD_751;
		fprintf(cntlog, "DPP_PSD_751, ");
	} else if (MajorNumber == 136) { 
		printf("This digitizer has a DPP_PSD firmware\n");
		WDcfg.DppType = DPP_PSD_730;
		fprintf(cntlog, "DPP_PSD_730, ");
	} else if (MajorNumber == 139) { 
		printf("This digitizer has a DPP_PHA firmware\n");
		WDcfg.DppType = DPP_PHA_730;
		fprintf(cntlog, "DPP_PHA_730, ");
	} else {
		printf("ERROR: This digitizer has not a DPP firmware\n");
		fprintf(cntlog, "Std FW, ");
		return -1;
	}
	printf("\nConnected to CAEN Digitizer Model %s; board num. %d\n", BoardInfo.ModelName, b);
	printf("ROC FPGA Release is %s\n", BoardInfo.ROC_FirmwareRel);
	printf("AMC FPGA Release is %s\n", BoardInfo.AMC_FirmwareRel);
	fprintf(cntlog, "%s, ", BoardInfo.ROC_FirmwareRel);
	fprintf(cntlog, "%s, ", BoardInfo.AMC_FirmwareRel);

	CAEN_DGTZ_ReadRegister(handle[b], 0x8158, &d32);
	if (d32 == 0x53D4) {
		printf("The DPP is licensed\n");
		fprintf(cntlog, "YES\n");
	} else {
		if (d32 > 0) {
			printf("\n>>> WARNING: DPP not licensed: %d minutes remaining\n\n", (int)((float)d32/0x53D4 * 30));
		} else {
			printf("\n>>> WARNING: DPP not licensed: time expired\n\n");
			return -1;
		}
		fprintf(cntlog, "NO\n");
	}

	SetTraceNames();
	SetVirtualProbes(b);

	return 0;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: force clock sync in one board
// Inputs:		b = board index
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ForceClockSync(int b)
{    
    int ret;
    Sleep(500);
    /* Force clock phase alignment */
    ret = CAEN_DGTZ_WriteRegister(handle[b], 0x813C, 1);
    /* Wait an appropriate time before proceeding */
    Sleep(2000);
    return ret;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: save run time settings
// --------------------------------------------------------------------------------------------------------- 
void SaveRTsettings()
{
	FILE *rt;
	int i;
	rt = fopen("rtset.txt", "w");
	fprintf(rt, "%d %d %d\n", ChToPlot, BrdToPlot, DoPlotHisto);
	for (i=0; i<MAX_NTRACES; i++)
		fprintf(rt, "%d ", TraceSet[i]);
	fclose(rt);
}

// --------------------------------------------------------------------------------------------------------- 
// Description: save run time settings
// --------------------------------------------------------------------------------------------------------- 
void ReadRTsettings()
{
	FILE *rt;
	int i;
	rt = fopen("rtset.txt", "r");
	if (rt==NULL)
		return;
	fscanf(rt, "%d %d %d\n", &ChToPlot, &BrdToPlot, &DoPlotHisto);
	for (i=0; i<MAX_NTRACES; i++)
		fscanf(rt, "%d ", &TraceSet[i]);
	SetVirtualProbes(BrdToPlot);
}



// --------------------------------------------------------------------------------------------------------- 
// Description: start the acquisition 
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int StartAcquisition() 
{
	int b;
	uint32_t d32a[MAX_NBRD], d32b[MAX_NBRD];

	if (WDcfg.StartMode == START_MODE_INDEP_SW) {
		for(b=0; b<WDcfg.NumBrd; b++) {
			CAEN_DGTZ_SWStartAcquisition(handle[b]);
		}
	} else if (WDcfg.StartMode == START_MODE_SYNC_SW) {
		for(b=0; b<WDcfg.NumBrd; b++) {
			uint32_t mask = (b==0) ? 0x80000000 : 0x40000000;
			CAEN_DGTZ_ReadRegister(handle[b], CAEN_DGTZ_TRIGGER_SRC_ENABLE_ADD, &d32a[b]);
			CAEN_DGTZ_ReadRegister(handle[b], CAEN_DGTZ_FP_TRIGGER_OUT_ENABLE_ADD, &d32b[b]);
			CAEN_DGTZ_WriteRegister(handle[b], CAEN_DGTZ_TRIGGER_SRC_ENABLE_ADD, mask);
			CAEN_DGTZ_WriteRegister(handle[b], CAEN_DGTZ_FP_TRIGGER_OUT_ENABLE_ADD, mask);
		}

		ForceClockSync(handle[0]);
		CAEN_DGTZ_SendSWtrigger(handle[0]); /* Send a software trigger to the 1st board to start the acquisition */

		/* set the registers back to the original settings */
		for(b=0; b<WDcfg.NumBrd; b++) {
			CAEN_DGTZ_WriteRegister(handle[b], CAEN_DGTZ_TRIGGER_SRC_ENABLE_ADD, d32a[b]);
			CAEN_DGTZ_WriteRegister(handle[b], CAEN_DGTZ_FP_TRIGGER_OUT_ENABLE_ADD, d32b[b]);
		}
	} else {
		// ???? HACK don't know what to write...
	}
	return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: stop the acquisition 
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int StopAcquisition()
{
	int b;
	for(b=0; b<WDcfg.NumBrd; b++) {
		CAEN_DGTZ_SWStopAcquisition(handle[b]); 
	}
	return 0;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Reset all the couters, histograms, etc...
// Return: 0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ResetStatistics()
{
	int b, ch;
	if (AcqRun)
		StopAcquisition();
	for(b=0; b<(WDcfg.NumBrd + NUM_VIRTUAL_BOARDS); b++) {
		for(ch=0; ch<WDcfg.NumCh; ch++) {
			if ((WDcfg.EnableInput[b][ch]) || (b >= WDcfg.NumBrd)) { 
				memset(Stats.EHisto[b][ch], 0, WDcfg.EHnbin * sizeof(uint32_t));
				memset(Stats.THisto[b][ch], 0, WDcfg.THnbin * sizeof(uint32_t));
				memset(Stats.MCSHisto[b][ch], 0, WDcfg.MCSHnbin * sizeof(uint32_t));
				memset(Stats.PSDHisto[b][ch], 0, 1024 * sizeof(uint32_t));
				memset(Stats.PSDvsEHisto[b][ch],  0, HISTO2D_NBINX * HISTO2D_NBINY * sizeof(uint32_t));
				Stats.EvCnt[b][ch] = 0;
				Stats.PrevEvCnt[b][ch] = 0;
				Stats.Matched[b][ch] = 0;
				Stats.ECnt[b][ch] = 0;
				Stats.TCnt[b][ch] = 0;
				Stats.MCSCnt[b][ch] = 0;
				Stats.PSDCnt[b][ch] = 0;
				Stats.Tmean[b][ch] = 0;
				Stats.Trms[b][ch] = 0;
			}
		}
	}
	if (AcqRun)
		StartAcquisition();
	return 0;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Allocate memory for the histograms
// Outputs:		AllocatedSize: total number of bytes allocated by the function
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int InitStatistics(uint32_t *AllocatedSize)
{
	int b, ch;

	*AllocatedSize = 0;
	for(b=0; b<(WDcfg.NumBrd + NUM_VIRTUAL_BOARDS); b++) {
		for(ch=0; ch<WDcfg.NumCh; ch++) {
			if ((WDcfg.EnableInput[b][ch]) || (b >= WDcfg.NumBrd)) { 
				Stats.EHisto[b][ch] = (uint32_t *)malloc( WDcfg.EHnbin*sizeof(uint32_t) );
				*AllocatedSize += WDcfg.EHnbin*sizeof(uint32_t);
				Stats.THisto[b][ch] = (uint32_t *)malloc( WDcfg.THnbin*sizeof(uint32_t) );
				*AllocatedSize += WDcfg.THnbin*sizeof(uint32_t);
				Stats.MCSHisto[b][ch] = (uint32_t *)malloc( WDcfg.MCSHnbin*sizeof(uint32_t) );
				*AllocatedSize += WDcfg.MCSHnbin*sizeof(uint32_t);
				Stats.PSDHisto[b][ch] = (uint32_t *)malloc( 1024 * sizeof(uint32_t) );
				*AllocatedSize +=  1024 * sizeof(uint32_t);
				Stats.PSDvsEHisto[b][ch] = (uint32_t *)malloc(HISTO2D_NBINX * HISTO2D_NBINY * sizeof(uint32_t));
				*AllocatedSize +=  HISTO2D_NBINX * HISTO2D_NBINY * sizeof(uint32_t);
			}
		}
	}
	ResetStatistics();
	return 0;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Open gnuplot for waveforms and histograms
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int OpenPlotter()
{
	/* open gnuplot in a pipe and the data file */
	if ((hplot == NULL) && ((WDcfg.AcquisitionMode != CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope)))
		hplot = popen(GNUPLOT_PATH, "w");  
	if (hplot==NULL) {
        printf("Can't open gnuplot\n");
        return -1;
    }
	fprintf(hplot, "set grid\n");
	fprintf(hplot, "set title 'Board %d - Channel %d'\n", BrdToPlot, ChToPlot);
	fflush(hplot);

	if ((wplot == NULL) && ((WDcfg.AcquisitionMode != CAEN_DGTZ_DPP_ACQ_MODE_List))) {
		wplot = popen(GNUPLOT_PATH, "w");  
		if (wplot==NULL) {
			printf("Can't open gnuplot\n");
			return -1;
		}
		fprintf(wplot, "set grid\n");
		fprintf(wplot, "set yrange [0:%d]\n", 1 << WDcfg.Nbit);
		fprintf(wplot, "set xlabel 'us'\n");
		fprintf(wplot, "set ylabel 'LSB'\n");
		fprintf(wplot, "set title 'Board %d - Channel %d'\n", BrdToPlot, ChToPlot);
		fflush(wplot);
		DoPlotWave=1;
	}
	return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: keyboard menu
// --------------------------------------------------------------------------------------------------------- 
void CheckKeyboardCommands()
{
	char c, c1;
	int ch, b, i, t;
	uint32_t rdata[MAX_NBRD];

	// Check keyboard
	if(kbhit()) {
		char enastring[2][10] = {"DISABLED", "ENABLED"};
		c = getch();

		menu:
		switch(c) {
		case 'q':  
			Quit = 1; 
			break;

		case 'i':
			for(b=0; b<WDcfg.NumBrd; b++)
				SaveRegImage(handle[b]);
			printf("Register Images saved to file");
			break;

		case 'c':  
			printf("Enter Channel number: ");
			scanf("%d", &ch);
			ChToPlot = ch;
			printf("Active Channel for plotting is now %d of Board %d\n", ChToPlot, BrdToPlot);
			break;

		case 'b':  
			printf("Enter board number: ");
			scanf("%d", &BrdToPlot);
			printf("Active Board for plotting is now %d\n", BrdToPlot);
			break;

		case 'r':
			printf("Use capital 'R' instead of 'r'\n");
			Sleep(200);
			break;
		case 'R':
			printf("Are you really sure (press 'y' for yes)?\n");
			c = getch();
			if (c == 'y')
				ResetStatistics();
			break;

		case 'w':
			/*
			DoPlotWave ^= 1;
			if (WDcfg.AcquisitionMode == CAEN_DGTZ_DPP_ACQ_MODE_List) {
				Quit = 1; 
				RestartAll = 1;
			}*/
			break;

		case 'g':
			EnableInternalPulseEmulator ^= 1;
			if (WDcfg.DppType == DPP_PHA_730) {
				for(b=0; b<WDcfg.NumBrd; b++) 
					for(ch=0; ch<WDcfg.NumCh; ch++)
						RegisterSetBits(handle[b], 0x1080 + (ch<<8), 14, 14, EnableInternalPulseEmulator);
				printf("Internal Test Pulse %s\n", enastring[EnableInternalPulseEmulator]);
			} else if ((WDcfg.DppType == DPP_PSD_730) || (WDcfg.DppType == DPP_PSD_720)) {
				for(b=0; b<WDcfg.NumBrd; b++) 
					for(ch=0; ch<WDcfg.NumCh; ch++)
						RegisterSetBits(handle[b], 0x1080 + (ch<<8), 8, 8, EnableInternalPulseEmulator);
				printf("Internal Test Pulse %s\n", enastring[EnableInternalPulseEmulator]);
			} else {
				printf("Option not available for this DPP type\n");
			}
			Sleep(300);
			break;


		case 'e':
			if (DoPlotHisto == HPLOT_ENERGY)
				DoPlotHisto = 0;
			else
				DoPlotHisto = HPLOT_ENERGY;
			break;

		case 'm':
			if (DoPlotHisto == HPLOT_MCS)
				DoPlotHisto = 0;
			else
				DoPlotHisto = HPLOT_MCS;
			break;

		case 'd':
			if (DoPlotHisto == HPLOT_TIME)
				DoPlotHisto = 0;
			else
				DoPlotHisto = HPLOT_TIME;
			break;

		case 'p':
			if (DoPlotHisto == HPLOT_PSD)
				DoPlotHisto = 0;
			else
				DoPlotHisto = HPLOT_PSD;
			break;

		case 'P':
			if (DoPlotHisto == HPLOT_ENERGYvsPSD)
				DoPlotHisto = 0;
			else
				DoPlotHisto = HPLOT_ENERGYvsPSD;
			break;

		case 'l':
			if (WDcfg.AcquisitionMode != CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope)
				WDcfg.SaveLists = WDcfg.SaveLists ? 0 : 1;
			printf("List Saving = %d\n", WDcfg.SaveLists);
			break;

		case 'D':
			DoSaveRawEvents = DoSaveRawEvents ? 0 : 1;
			printf("Raw Data Saving = %d\n", DoSaveRawEvents);
			break;

		case 'h':
			if (WDcfg.AcquisitionMode != CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope)
				SaveAllHistograms();
			break;

		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
			t = c - '1';
			if (c != '1')
				printf("x - OFF\n");
			for(i=0; i<MAX_NTRSETS; i++) {
				if (TraceNames[t][i][0] == '-')
					continue;
				printf("%d - %s\n", i, TraceNames[t][i]);
			}
			c1 = getch();
			if ((c1 == 'x') && (c != '1'))
				TraceSet[t] = -1;
			else
				TraceSet[t] = (int)(c1 - '0');
			SetVirtualProbes(BrdToPlot);
			Sleep(100);
			break;

		case 'T':
			ContinuousTrigger = ContinuousTrigger ? 0 : 1;
			printf("Continuous Trigger = %d\n", ContinuousTrigger);
			break;

		case 't':
			SingleTrigger = 1;
			break;

		case 'f':
			DoRefresh ^= 1;
			if (!DoRefresh) printf("Plots and Logs refresh is now disabled!\n");
			break;

		case 'x':
			keVscale ^= 1;
			break;

		case 'o':
			DoRefresh = 0;
			DoRefreshSingle = 1;
			break;

		case 's':
			// Start Acquisition
			Quit = 1; 
			RestartAll = 1;
			StartAcquisition();
			printf("Acquisition Started\n");
			AcqRun = 1;
			break;

		case 'S':
			StopAcquisition();
			SaveAllHistograms();
			printf("Acquisition Stopped\n");
			AcqRun = 0;
			break;

		case 'M':
			b = 0;
			while (1) {
				uint32_t ra, rd;
				printf("w : Write register\n");
				printf("r : Read register\n");
				printf("b : Change board (now is %d)\n", b);
				printf("q : Quit manual register access\n\n");
				c1 = getch();
				if (c1 == 'q') 
					break;
				printf("Enter Register Address (16 bit hex) : ");
				scanf("%x", &ra);
				if (c1 == 'b') {
					printf("Enter Board number : ");
					scanf("%d", &b);
				}
				if (c1 == 'w') {
					printf("Enter Register Data (32 bit hex) : ");
					scanf("%x", &rd);
					CAEN_DGTZ_WriteRegister(handle[b], ra, rd);
				}
				CAEN_DGTZ_ReadRegister(handle[b], ra, &rd);
				printf("\nRegister 0x%04X: 0x%08X\n\n\n", ra, rd);	
			}
			break;


		case 'k' :
			// propagate CLK to trgout on all boards
			for(b=0; b<WDcfg.NumBrd; b++) {
				CAEN_DGTZ_ReadRegister(handle[b], CAEN_DGTZ_FRONT_PANEL_IO_CTRL_ADD, &rdata[b]);
				CAEN_DGTZ_WriteRegister(handle[b], CAEN_DGTZ_FRONT_PANEL_IO_CTRL_ADD, 0x00050000);
			}
			printf("Trigger Clk is now output on TRGOUT.\n");
			printf("Press [r] to reload PLL config, any other key to quit clock monitor\n");
			while( (c=getch()) == 'r') {
				//CAEN_DGTZ_WriteRegister(handle[0], 0xEF34, 0);
				ForceClockSync(handle[0]);
				printf("PLL reloaded\n");
			}
			for(b=0; b<WDcfg.NumBrd; b++)
				CAEN_DGTZ_WriteRegister(handle[b], CAEN_DGTZ_FRONT_PANEL_IO_CTRL_ADD, rdata[b]);
			break;

		case 'u' : {
			char *buff;
			uint32_t nb;
			uint64_t tnb=0;
			long t, pt;

			ClearScreen();
			printf("Readout Bandwidth Test. Press a key to stop\n");
			buff = (char *)malloc(8*1024*1024);
			pt = get_time();
			while(!kbhit()) {
				for(b=0; b<WDcfg.NumBrd; b++) {
					CAEN_DGTZ_ReadData(handle[b], CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, buff, &nb);
					tnb += nb;
				}
				t=get_time();
				if ((t-pt) > 1000) {
					printf("Readout Rate = %.3f MB/s\n", ((float)tnb/(1024*1024)) / ((t-pt)/1000) );
					pt = t;
					tnb = 0;
				}
			}
			free(buff);
			ClearScreen();
			c = getch();
			}
			break;

		case ' ':
			printf("\n");
			printf("s )  Start acquisition (also valid for Reconfigure + Restart) \n");
			printf("S )  Stop acquisition\n");
			printf("M )  Manual Register Access\n");
			printf("1-6) Change Probe settings for oscilloscope traces\n");
			printf("k )  Propagate CLK to trgout on all boards\n");
			printf("t )  Force trigger (one shot)\n");
			printf("T )  Force trigger (enable/disable continuous) \n");
			printf("g )  Enable/Disable internal pulse emulator \n");
			printf("R )  Reset Histograms and Statistics\n");
			printf("h )  Save Histograms to files (one shot)\n");
			printf("l )  Enable/Disable List dump to files\n");
			printf("D )  Enable/Disable Raw Data dump to files\n");
			printf("b-c) Choose board/channel to plot\n");
			printf("f )  Enable/Disable automatic refersh (plots and logs)\n");
			printf("o )  One shot refresh (plots and logs)\n");
			printf("w )  Enable/Disable Waveform plot (switching between List-Mixed modes)\n");
			printf("e )  Enable/Disable Energy Histogram plot\n");
			printf("m )  Enable/Disable MCS Histogram plot\n");
			printf("d )  Enable/Disable DeltaT Histogram plot\n");
			printf("p )  Enable/Disable PSD Histogram plot\n");
			printf("P )  Enable/Disable Energy-PSD scatter 2D-plot\n");
			printf("x )  Toggle between keV and ADC channel in Energy Spectrum\n");
			printf("i )  Save Registers Image\n");
			printf("u )  Readout Bandwidth Test\n");
			printf("q )  Quit\n");
			printf("Space )  Print Menu\n\n");
			c = getch();
			goto menu;
			break;

		default: break;
		}
	}
	SaveRTsettings();
}



/* ########################################################################### */
/* MAIN                                                                        */
/* ########################################################################### */
int main(int argc, char *argv[])
{
	char ConfigFileName[100];					// Config file name
	CAEN_DGTZ_ErrorCode ret= CAEN_DGTZ_Success;					// return code of the CAEN_DGTZ library functions
	GenericDPPEvent_t evnt[MAX_NBRD][MAX_NCH];  // Events extracted from the queue
	Waveform_t Waveform;                        // Waveform struct for plotting
	Waveform_t WaveformRef;                     // Waveform struct for two channel plotting
	int WPch = 0, WPref = 0, WavePlotRdy = 0;   // Waveforms ready for plotting
	int SkipPlot = 0;							// don't do this plot (for event filtering)
	int EvRdy[MAX_NBRD][MAX_NCH];               // Event Ready for a specific channel 
	int TestConnection = 0;
	int ImmediateStart = 0;
	int LoadRTsettings = 0;
	int LineCmdRunNumber = -1;
	int nevFound;
	uint32_t StopCh[MAX_NBRD][MAX_NCH];
    int StartCh=0;  // HACK: read this two parameters from config file
	int StartBrd=0;
	// Other variables 
	int i, b, ch, rg;
	uint32_t AllocSize, TotAllocSize=0;
	int ExitSleep = 3000;
	uint64_t PrevNb=0;
	uint64_t CurrentTime, PrevRateTime, ElapsedTime, PrevPlotTime, PrevKeybTime, TotRunTime, StartTime;
	uint64_t PrevTimeStamp[MAX_NBRD][MAX_NCH];
	// input and output files 
	FILE *f_ini;                     // config file
	FILE *EventLog;
    pFragment_t builtEvents = (pFragment_t)malloc(MAX_NBRD*MAX_NCH*sizeof(Fragment_t));
	int nFrags;
    int servicestatus;

    printf("*************************************************\n");
    printf("CAEN SpA - digiTES Rev %s\n", digiTES_Revision);
    printf("*************************************************\n");
#ifdef WIN32
	//EmulateDataFile();
	//signal(SIGTERM, signalHandler);
	printf("Wait for the driver connection\n");
	servicestatus = DoStartSvc();
	if (servicestatus < 0) {
		printf("DigiTES is not able to communicate with A3818 Driver\n");
		printf("Only USB connections will be possible\n");
		printf("<Enter> to continue:");
		getch();
	}
	else {
		printf("Connection with the driver %dK\n", servicestatus);
	}
#endif
	Restart0:
	/* *************************************************************************************** */
	/* Init variables                                                                          */
	/* *************************************************************************************** */
	memset(&Stats, 0, sizeof(Stats));
	for(b=0; b<MAX_NBRD; b++) {
		for(ch=0; ch<MAX_NCH; ch++) {
			PrevTimeStamp[b][ch] = 0;
			EvRdy[b][ch] = 0;
			StopCh[b][ch] = 0;
		}
	}


	/* *************************************************************************************** */
	/* Open and parse configuration file                                                       */
	/* *************************************************************************************** */
	sprintf(ConfigFileName, "%s", "digiTES_Config_common.txt");
	for (i=1; i<argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 't') TestConnection = 1;
			if (argv[i][1] == 'i') ImmediateStart = 1;
			if (argv[i][1] == 'l') LoadRTsettings = 1;
			if (argv[i][1] == 'r') sscanf(&argv[i][2], "%d", &LineCmdRunNumber);
		} else {
		  sprintf(ConfigFileName, "%s", argv[i]);
		}
	}
		
	printf("Opening Configuration File %s\n", ConfigFileName);
	f_ini = fopen(ConfigFileName, "r");
	if (f_ini == NULL ) {
		printf("Can't open main config file %s\n", ConfigFileName);
		goto QuitProgram;
	}
	ParseConfigFile(f_ini, &WDcfg);
	printf("Found settings for %d board(s) in the config file\n", WDcfg.NumBrd);
	fclose(f_ini);
	if (DoPlotWave)
		WDcfg.AcquisitionMode = ACQ_MODE_MIXED;
	if (EnableInternalPulseEmulator) 
		for(b=0; b<WDcfg.NumBrd; b++) 
			for(ch=0; ch<WDcfg.NumCh; ch++)
				WDcfg.EnableIPE[b][ch] = 1;

	Restart:
	if (!WDcfg.ReplayDataFile) {
		/* *************************************************************************************** */
		/* Open the digitizers and read board information                                          */
		/* *************************************************************************************** */
		cntlog = fopen("ConnectionLog.txt", "w");
		for(b=0; b<WDcfg.NumBrd; b++) {
			char ct[2][10] = {"USB", "OPT_LINK"};

			if (!RestartAll) {
				fprintf(cntlog, "BOARD%d, ", b);
				printf("Opening digitizer n. %d through %s %d ", b, ct[WDcfg.LinkType[b]], WDcfg.LinkNum[b]);
				if (WDcfg.LinkType[b] != 0)    printf("%d ",  WDcfg.ConetNode[b]);
				if (WDcfg.BaseAddress[b] != 0) printf("BA=%08X ", WDcfg.BaseAddress[b]);
				ret = CAEN_DGTZ_OpenDigitizer((CAEN_DGTZ_ConnectionType)WDcfg.LinkType[b], WDcfg.LinkNum[b], WDcfg.ConetNode[b], WDcfg.BaseAddress[b], &handle[b]);
				if (ret) {
					int rn;
					printf(": Failed!!! Error Code = %d\n", ret);
					fprintf(cntlog, "Failed, -, -, -, -, -, -, -\n");
					if (TestConnection) 
						continue;
					printf("Press 'q' to quit, any other kay to read data from file (run off-line)\n");
					if (getch() == 'q') {
						ExitSleep = 0;
						goto QuitProgram;
					}
					printf("Enter RunNumber or -1 to select a specific file name: ");
					scanf("%d", &rn);
					if (rn>=0) {
						sprintf(WDcfg.InputDataFileName, "Run%d.dat", rn);
					} else {
						printf("Enter File Name: ");
						scanf("%s", WDcfg.InputDataFileName);
					}
					goto Restart;
				} else {
					printf(": Done\n");
				}
			}

			RestartAll = 0;
			// Reset the digitizer
			ret = CAEN_DGTZ_Reset(handle[b]);
			if (ret != 0) {
				printf("ERROR: can't reset the digitizer.\n");
				goto QuitProgram;    
			}
			CAEN_DGTZ_WriteRegister(handle[0], 0x8000, 0);  


			// Read Board Info
			if (ReadBoardInfo(b) < 0)
				goto QuitProgram;

		}
	} else {  // Replay from data file
		//sprintf(WDcfg.InputDataFileName, "List_0_2.txt");
		if (OpenInputDataFile() < 0)
			goto QuitProgram;
	}

	if (cntlog != NULL)
		fclose(cntlog);

	if (TestConnection) {
		ExitSleep = 0;
		goto QuitProgram;
	}

	// mask the channels not available for this model
	for(b=0; b<WDcfg.NumBrd; b++)
		for(i=WDcfg.NumCh; i<MAX_NCH; i++)
			WDcfg.EnableInput[b][i] = 0;

	// Set channel to plot as the first enabled channel in board 0
	BrdToPlot = 0;
	for(i=0; i<MAX_NCH; i++) {
		if (WDcfg.EnableInput[0][i]) {
			ChToPlot = i;
			break;
		}
	}

	// Recalculate max value of the histograms in order to have an integer bin
	for(b=0; b<MAX_NBRD; b++) {
		for(ch=0; ch<MAX_NCH; ch++) {
			double binsize;
			binsize = ceil(1000.0 * (WDcfg.THmax[b][ch] - WDcfg.THmin[b][ch]) / WDcfg.THnbin);  // in ps
			WDcfg.THmax[b][ch] = (float)(WDcfg.THmin[b][ch] + (binsize * WDcfg.THnbin) / 1000.0);
			binsize = ceil((WDcfg.EHmax[b][ch] - WDcfg.EHmin[b][ch]) / WDcfg.EHnbin); 
			WDcfg.EHmax[b][ch] = (float)(WDcfg.EHmin[b][ch] + (binsize * WDcfg.EHnbin));
		}
	}

	/* *************************************************************************************** */
	/* Program the digitizer (see function ProgramDigitizer)                                   */
	/* *************************************************************************************** */
	if (!WDcfg.ReplayDataFile) {
		for(b=0; b<WDcfg.NumBrd; b++) {
			ret = (CAEN_DGTZ_ErrorCode) ProgramDigitizer(b);
			if (ret) {
				printf("Failed to program the digitizer\n");
				goto QuitProgram;
			}
		}
	}

	/* *************************************************************************************** */
	/* Allocate memory buffers and init statistics                                             */
	/* *************************************************************************************** */
	// Allocate and initialize Histograms, Counters, memory buffers, waveforms, etc...
	if (InitStatistics(&AllocSize) < 0)	goto QuitProgram;
	TotAllocSize += AllocSize;
	if (InitReadout(&AllocSize) < 0) goto QuitProgram;
	TotAllocSize += AllocSize;
	if (InitPreProcess(&AllocSize) < 0)	goto QuitProgram;
	TotAllocSize += AllocSize;
	if (AllocateWaveform(&Waveform, &AllocSize) < 0) goto QuitProgram;
	TotAllocSize += AllocSize;
	if (AllocateWaveform(&WaveformRef, &AllocSize) < 0)	goto QuitProgram;
	TotAllocSize += AllocSize;
	printf("Allocated %.2f MB for histograms and counters\n", (float)TotAllocSize/(1024*1024));

	/* *************************************************************************************** */
	/* Open the plotters                                                                       */
	/* *************************************************************************************** */
	if (OpenPlotter() < 0)
		goto QuitProgram;


	/* *************************************************************************************** */
	/* Readout Loop                                                                            */
	/* *************************************************************************************** */
	if (LoadRTsettings)
		ReadRTsettings();
	if (LineCmdRunNumber != -1)
		WDcfg.RunNumber = LineCmdRunNumber;
	if (!AcqRun) {
		char c;
		if (ImmediateStart) {
			c = 's';
		} else {
			printf("press 's' to start the acquisition, any other key to continue with run stopped\n");
			c = getch();
		}
		if (c == 's') {
			StartAcquisition();
			AcqRun = 1;
		}
	}
	StartTime = get_time();
	PrevRateTime = StartTime;
	PrevPlotTime = StartTime;
	PrevKeybTime = StartTime;
	printf("Press [Space] for help\n\n");
	if (!DoRefresh) printf("Plots and Logs are OFF; press 'f' to enable or 'o' for a singles\n");
	if (ENEVLOG) EventLog = fopen("EventLog.txt", "w");

	while(!Quit) {
		CurrentTime = get_time();

		// Check from commands from the keyboard 
		if ((CurrentTime - PrevKeybTime) > 200) {
			CheckKeyboardCommands();
			PrevKeybTime = CurrentTime;
		}

		/* Send a software trigger to each board */
		if (!WDcfg.ReplayDataFile && (ContinuousTrigger || SingleTrigger)) {
			for(b=0; b<WDcfg.NumBrd; b++)
				CAEN_DGTZ_SendSWtrigger(handle[b]);
			if (SingleTrigger)				
				printf("Single Software Trigger issued\n");
			SingleTrigger = 0;
		}

		/* Log throughput and trigger rate to screen (every second) */
		ElapsedTime = CurrentTime - PrevRateTime; /* milliseconds */
		if (ElapsedTime > 1000) {
			int stoprun = 1;
			TotRunTime = CurrentTime - StartTime;
			if (DoRefresh || DoRefreshSingle) {
				ClearScreen();
				printf("Press [Space] for help\n\n");
				if (WDcfg.ReplayDataFile) {
					printf("\nData File readout rate=%.2f MB\n", (float)(Stats.TotNb-PrevNb)/((float)ElapsedTime*1048.576f));
				} else if (AcqRun) {
					printf("RUN n. %d: RealTime = %d sec\n", WDcfg.RunNumber, TotRunTime/1000);
					printf("Readout Rate=%.2f MB; AvgBltSize=%.2fKB; AvgNumEv=%.2f\n", (float)(Stats.TotNb-PrevNb)/((float)ElapsedTime*1048.576f), (float)((Stats.TotNb-PrevNb)/1024)/Stats.NumBlt, (float)Stats.TotNumEv/Stats.NumBlt);
				} else {
					printf("RUN n. %d: Stopped\n\n", WDcfg.RunNumber);
				}
				printf("Sum TOF: %7.2f cnts/s (%d Total cnts)\n", 1000.0 * (Stats.TCntSum - Stats.PrevTCntSum)/(float)ElapsedTime, Stats.TCntSum);
				Stats.PrevTCntSum = Stats.TCntSum;
				if (DoSaveRawEvents) {
					printf("Output Data Size = %.2f MB. ", (float)OutFileSize / (1024*1024));
					if (OutFileSize > (MAX_OUTPUT_DATA_SIZE*1024*1024)) {
						DoSaveRawEvents = 0;
						printf("Saving Stopped");
					}
					printf("\n");
				}
				for(b=0; b<WDcfg.NumBrd; b++) {
					printf("\nBrd  Ch |    Trg-Rate   Match-Rate   Match %%   QueueOccup   TotEvCnt\n");
					printf("-----------------------------------------------------------------------\n");
					for(ch=0; ch<WDcfg.NumCh; ch++) {
						double elapsed = (double)(Stats.CurrEventTime[b][ch] - Stats.PrevEventTime[b][ch]) / 1000000.0;  // elapsed time (in ms) for this channel
						uint32_t nev = Stats.EvCnt[b][ch]-Stats.PrevEvCnt[b][ch];
						if (!WDcfg.EnableInput[b][ch]) 
							printf("%3d %2d  |    Disabled\n", b, ch);
						else if (StopCh[b][ch])
							printf("%3d %2d  |    Stopped                                        %-10u\n", b, ch, Stats.EvCnt[b][ch]);
						else if (elapsed>0)
							printf("%3d %2d  | %7.2f KHz  %7.2f KHz  %6.2f %%  %6.2f%%       %-10u\n", b, ch, nev/elapsed, Stats.Matched[b][ch]/elapsed, 100.0*Stats.Matched[b][ch]/nev, 100.0*Qnev[b][ch]/EV_QUEUE_SIZE, Stats.EvCnt[b][ch]);
						else
							printf("%3d %2d  |    No Data                                        %-10u\n", b, ch, Stats.EvCnt[b][ch]);
						if ((WDcfg.StopOnTime > 0) && (TotRunTime > WDcfg.StopOnTime) && (Stats.PrevEvCnt[b][ch] == Stats.EvCnt[b][ch]))
							StopCh[b][ch] = 1;
						Stats.PrevEvCnt[b][ch] = Stats.EvCnt[b][ch];
						Stats.Matched[b][ch] = 0;
						Stats.PrevEventTime[b][ch] = Stats.CurrEventTime[b][ch];
					}
				}
				printf("\n\n");

				PrevNb = Stats.TotNb;
				Stats.NumBlt=0;
				Stats.TotNumEv=0;
				if (DoPlotHisto && (WDcfg.AcquisitionMode != CAEN_DGTZ_DPP_ACQ_MODE_Oscilloscope)) {
					if (!WDcfg.EnableInput[BrdToPlot][ChToPlot] && (BrdToPlot != (WDcfg.NumBrd + NUM_VIRTUAL_BOARDS - 1))) {
						printf("WARNING: the selected channel for plot is disabled (b=%d, ch=%d)\n", BrdToPlot, ChToPlot);
					} else {
						char title[500], xlabel[500];
						if (DoPlotHisto == HPLOT_ENERGY) {
							float a, b;
							if (keVscale) {
								a = WDcfg.EHmin[BrdToPlot][ChToPlot] * WDcfg.ECalibration_m[BrdToPlot][ChToPlot] + WDcfg.ECalibration_q[BrdToPlot][ChToPlot];
								b = WDcfg.EHmax[BrdToPlot][ChToPlot] * WDcfg.ECalibration_m[BrdToPlot][ChToPlot] + WDcfg.ECalibration_q[BrdToPlot][ChToPlot];
								sprintf(xlabel, "keV");
							} else {
								a = WDcfg.EHmin[BrdToPlot][ChToPlot];
								b = WDcfg.EHmax[BrdToPlot][ChToPlot];
								sprintf(xlabel, "ADC counts");
							}
							sprintf(title, "Brd%d Ch%d #Cnt=%d\n", BrdToPlot, ChToPlot, Stats.ECnt[BrdToPlot][ChToPlot]);
							PlotHisto(hplot, Stats.EHisto[BrdToPlot][ChToPlot], WDcfg.EHnbin, a, b, title, xlabel);
						} else if (DoPlotHisto == HPLOT_MCS) {
							sprintf(title, "Brd%d Ch%d #Cnt=%d\n", BrdToPlot, ChToPlot, Stats.MCSCnt[BrdToPlot][ChToPlot]);
							sprintf(xlabel, "Time (ms)");
							PlotHisto(hplot, Stats.MCSHisto[BrdToPlot][ChToPlot], WDcfg.MCSHnbin, 0, (float)WDcfg.MCSHnbin*WDcfg.DwellTime/1000, title, xlabel);
						} else if (DoPlotHisto == HPLOT_PSD) {
							sprintf(title, "Brd%d Ch%d #Cnt=%d\n", BrdToPlot, ChToPlot, Stats.PSDCnt[BrdToPlot][ChToPlot]);
							sprintf(xlabel, "PSD");
							PlotHisto(hplot, Stats.PSDHisto[BrdToPlot][ChToPlot], 1024, 0, 1, title, xlabel);
						} else if (DoPlotHisto == HPLOT_TIME) {
							double mean, rms;
							mean = Stats.Tmean[BrdToPlot][ChToPlot] / Stats.TCnt[BrdToPlot][ChToPlot];
							rms = sqrt(Stats.Trms[BrdToPlot][ChToPlot] / Stats.TCnt[BrdToPlot][ChToPlot] - mean*mean);
							sprintf(title, "Brd%d Ch%d #Cnt=%d M=%.3f ns, S=%.2f ps\n", BrdToPlot, ChToPlot, Stats.TCnt[BrdToPlot][ChToPlot], mean, rms*1000);
							sprintf(xlabel, "ns");
							PlotHisto(hplot, Stats.THisto[BrdToPlot][ChToPlot], WDcfg.THnbin, WDcfg.THmin[BrdToPlot][ChToPlot], WDcfg.THmax[BrdToPlot][ChToPlot], title, xlabel);
						} else if (DoPlotHisto == HPLOT_ENERGYvsPSD) {
							sprintf(title, "Brd%d Ch%d #Cnt=%d\n", BrdToPlot, ChToPlot, Stats.PSDCnt[BrdToPlot][ChToPlot]);
							PlotHisto2D(hplot, Stats.PSDvsEHisto[BrdToPlot][ChToPlot], HISTO2D_NBINX, HISTO2D_NBINY, title);
						}
					}
				}
			}
			for(b=0; b<WDcfg.NumBrd; b++)
				for(ch=0; ch<WDcfg.NumCh; ch++) 
					if (WDcfg.EnableInput[b][ch] && !StopCh[b][ch])
						stoprun = 0;
			if (stoprun) {
				StopAcquisition();
				SaveAllHistograms();
				AcqRun = 0;
			}
			PrevRateTime = CurrentTime;
		}

		if (!AcqRun) {
			Sleep(10);
			continue;
		}



		/* ----------------------------------------------------------------------------------- */
		/* Get events from the queues */   
		/* ----------------------------------------------------------------------------------- */
		//SelectEvents();  // apply filters (event selection) to the events in the queues
		if (WDcfg.EnableTimeCorrelFilter) {
			int StartCh=0;  // HACK: read this two parameters from config file
			int StartBrd=0;
			nevFound = GetEvent_Correlated(StartCh, StartBrd, evnt, EvRdy);
		} else {
			nevFound = 0;
			for(b=0; b<WDcfg.NumBrd; b++) {
				for(ch=0; ch<WDcfg.NumCh; ch++) {
					EvRdy[b][ch] = 0;
					if (!WDcfg.EnableInput[b][ch])
						continue;
					rg = GetEvent(b, ch, &evnt[b][ch]);
					if (rg < 0)
						goto QuitProgram;
					if (rg == 0) // no data from this channel
						continue;
					EvRdy[b][ch] = 1;
					nevFound++;
					Stats.Matched[b][ch]++;
					if (ENEVLOG) fprintf(EventLog, "Get Ch%d: TS=%u\n", ch, evnt[b][ch].TimeStamp);
				}
			}
		}

		if (DoPlotWave && ((CurrentTime-PrevPlotTime) > 200)) {
			if (EvRdy[BrdToPlot][ChToPlot]) {
				if (ConvertWaveform(BrdToPlot, ChToPlot, evnt[BrdToPlot][ChToPlot], &Waveform) == 0) {
					WavePlotRdy = 1;
				}
			}						
		}



		/* ----------------------------------------------------------------------------------- */
		/* Analyze Events data */
		/* ----------------------------------------------------------------------------------- */
		if (ENEVLOG) fprintf(EventLog, "Start Analysis\n");
		if (nevFound > 0) {
			GetBuiltEvent(StartBrd,StartCh,EvRdy,evnt,builtEvents, &nFrags);

            if ((strlen(WDcfg.pipeBaseName) > 0) && (nFrags > 0)) {
				WriteEventToPipe(nFrags, builtEvents);              /* Linux might make errors*/
			}

			for(b=0; b<WDcfg.NumBrd; b++) {
				for(ch=0; ch<WDcfg.NumCh; ch++) {
					double time;
					float  Energy;
					int BrdRef = 0; //WDcfg.CRboard[b][ch];
					int ChRef = 0; //WDcfg.CRchannel[b][ch];
					uint64_t Tbin;
					if (EvRdy[b][ch] != 1) {
						EvRdy[b][ch] = 0;
						continue;
					}
					EvRdy[b][ch] = 0;

					if (StopCh[b][ch]) {
						Stats.EvCnt[b][ch]--;  // single channel stop condition: data are still coming, but they are discarded
						continue;
					}
				
					/* Energy & PSD Spectra */
					Energy = evnt[b][ch].Energy * WDcfg.EnergyGain[b][ch] + WDcfg.EnergyOffset[b][ch];
					if ((Energy >= WDcfg.EHmin[b][ch]) && (Energy < WDcfg.EHmax[b][ch])) {

						uint16_t Ebin = (uint16_t)((Energy - WDcfg.EHmin[b][ch]) * WDcfg.EHnbin / (WDcfg.EHmax[b][ch] - WDcfg.EHmin[b][ch]));
						uint16_t psdbin = (int)(evnt[b][ch].psd*1024);
						uint16_t bx = Ebin * HISTO2D_NBINX / WDcfg.EHnbin;
						uint16_t by = (int)(evnt[b][ch].psd * HISTO2D_NBINY);

						Stats.EHisto[b][ch][Ebin]++;
						if ((NUM_VIRTUAL_BOARDS > 0) && (b > 0)) {  // HACK
							Stats.EHisto[WDcfg.NumBrd+NUM_VIRTUAL_BOARDS-1][0][Ebin]++;  // histogram sum
							Stats.ECntSum++;
						}
						Stats.ECnt[b][ch]++;

						if ((psdbin > 0) && (psdbin < 1024)) {
							Stats.PSDHisto[b][ch][psdbin]++;
							Stats.PSDCnt[b][ch]++;
						}
						if ((bx > 0) && (bx < HISTO2D_NBINX) && (by > 0) && (by < HISTO2D_NBINY))
							Stats.PSDvsEHisto[b][ch][bx * HISTO2D_NBINY + by]++;

					/* Time Spectrum */
					if ((b == BrdRef) && (ch == ChRef))  // self referred channel
						time = (double)(evnt[b][ch].TimeStamp - PrevTimeStamp[b][ch]);  // delta T between pulses on the same channel (in ps)
					else
						time = (double)evnt[b][ch].TimeStamp - (double)evnt[BrdRef][ChRef].TimeStamp + ((double)evnt[b][ch].FineTimeStamp - (double)evnt[BrdRef][ChRef].FineTimeStamp)/1000;  // delta T from Ref Channel (in ns)

					PrevTimeStamp[b][ch] = evnt[b][ch].TimeStamp;
					Tbin = (uint64_t)((time - WDcfg.THmin[b][ch]) * WDcfg.THnbin / (WDcfg.THmax[b][ch] - WDcfg.THmin[b][ch]));
					if ((Tbin>=0) && (Tbin<WDcfg.THnbin)) {
						Stats.THisto[b][ch][Tbin]++;
						if ((NUM_VIRTUAL_BOARDS > 0) && (b > 0)) {
							Stats.THisto[WDcfg.NumBrd+NUM_VIRTUAL_BOARDS-1][0][Tbin]++;  // histogram sum
							Stats.TCntSum++;
						}

						Stats.TCnt[b][ch]++;
						Stats.Tmean[b][ch] += time;
						Stats.Trms[b][ch] += (time*time);
					}

					/* Multi Channel Scaler */
					Tbin = evnt[b][ch].TimeStamp / ((uint64_t)WDcfg.DwellTime*1000);
					if ((Tbin >= 0) && (Tbin < (uint64_t)WDcfg.MCSHnbin)) {
						Stats.MCSHisto[b][ch][Tbin]++;
						Stats.MCSCnt[b][ch]++;
					}

						/* Write one event (T, E, S) into the relevant list file */
						if (WDcfg.SaveLists)	
							SaveList(b, ch, evnt[b][ch].TimeStamp * 1000 + evnt[b][ch].FineTimeStamp, (uint16_t)Energy, evnt[b][ch].psd);

						/* Check stopping criteria */
						if (((WDcfg.StopOnTotalEvents > 0)  && (Stats.EvCnt[b][ch] >= (uint32_t)WDcfg.StopOnTotalEvents))  ||
							((WDcfg.StopOnTimeEvents > 0)   && (Stats.TCnt[b][ch]  >= (uint32_t)WDcfg.StopOnTimeEvents))   ||
							((WDcfg.StopOnEnergyEvents > 0) && (Stats.ECnt[b][ch]  >= (uint32_t)WDcfg.StopOnEnergyEvents)) ||
							((WDcfg.StopOnTime > 0)         && (evnt[b][ch].TimeStamp >= ((uint64_t)WDcfg.StopOnTime * 1000000)))) {
							StopCh[b][ch] = 1;
						}
					}
				}

				// Waveform Plotting 
				if (WavePlotRdy && (DoRefresh || DoRefreshSingle)) {
					char title[500];
					if (!SkipPlot) {
						sprintf(title, "WAVEFORMS : Board %d - Channel %d", BrdToPlot, ChToPlot);
						PlotWaveforms(wplot, Waveform, title);
						PrevPlotTime = CurrentTime;
					}
					WavePlotRdy = 0;
				}
			}
			SkipPlot = 0;
			DoRefreshSingle = 0;

		}
	} // End of readout loop


	SaveAllHistograms();
	ExitSleep = 0;

QuitProgram:
	Quit = 0;
	/* stop the acquisition, close the device and free the buffers */
	for(b=0; b<WDcfg.NumBrd; b++) {
		CAEN_DGTZ_SWStopAcquisition(handle[b]);
		for (ch=0; ch<WDcfg.NumCh; ch++)
			if (Stats.EHisto[b][ch] != NULL)
				free(Stats.EHisto[b][ch]);
			if (Stats.PSDHisto[b][ch] != NULL)
				free(Stats.PSDHisto[b][ch]);
			if (Stats.PSDvsEHisto[b][ch] != NULL)
				free(Stats.PSDvsEHisto[b][ch]);
	}
	CloseReadout();
	ClosePreProcess();
	if (wplot != NULL)	
		fclose(wplot);
	if (hplot != NULL)
		fclose(hplot);
	hplot = NULL;
	wplot = NULL;
	if (RestartAll) 
		goto Restart0;
	for(b=0; b<WDcfg.NumBrd; b++)
		CAEN_DGTZ_CloseDigitizer(handle[b]);
	Sleep(ExitSleep);
	return ret;
}

