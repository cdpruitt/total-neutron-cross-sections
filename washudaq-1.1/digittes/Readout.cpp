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

#include "keyb.h"
#include "digiTES.h"

#include <stdio.h>
#ifdef WIN32
    #include <process.h>
#else
    #include <unistd.h>
    #include <stdint.h>   /* C99 compliant compilers: uint64_t */
#endif

// Event Queues and pointers 
GenericDPPEvent_t *EvQ[MAX_NBRD][MAX_NCH];    // Event queues
uint32_t EvWpnt[MAX_NBRD][MAX_NCH];           // Write pointer to EvQ
uint32_t EvRpnt[MAX_NBRD][MAX_NCH];           // Read pointer to EvQ
uint32_t Qnev[MAX_NBRD][MAX_NCH];             // Number of events in the queues
uint32_t AlmostFull[MAX_NBRD];				  // Queue Almost Full flags (one bit per channel)
int BoardEmpty[MAX_NBRD];					  // All the queues of a board are empty 
int AllEmpty=1;

// buffers for the readout
char *buffer = NULL;                // readout buffer
void *Events[MAX_NCH];              // memory buffers for the event decoding
uint32_t *wbuffer = NULL;           // buffer for the waveforms in off-line mode
uint32_t wbpnt = 0;                 // wbuffer pointer

// Data Files
uint32_t OutFileSize = 0;		// Size of the output data file (in bytes)
FILE *OutputDataFile = NULL;	// event data output file
FILE *InputDataFile = NULL;		// event data input file (for the data replay mode)
uint64_t emtstamp = 0;
uint64_t lastemtstamp = 0;

Waveform_t PPwaveform;			// waveform decoded for the Waveform processor (DPP off-line)


// Readout Log File
FILE *rolog;

#define IS_BAD_EVENT(b, ch)				((EvQ[b][ch][EvRpnt[b][ch]].Flags & 0xC000) == 0x8000)
#define IS_GOOD_EVENT(b, ch)			((EvQ[b][ch][EvRpnt[b][ch]].Flags & 0xC000) == 0xC000)
#define IS_SUSPENDED_EVENT(b, ch)		((EvQ[b][ch][EvRpnt[b][ch]].Flags & 0xC000) == 0x4000)
#define IS_UNPROCESSED_EVENT(b, ch)		((EvQ[b][ch][EvRpnt[b][ch]].Flags & 0xC000) == 0x0000)



// --------------------------------------------------------------------------------------------------------- 
// Description: Set a pointer to one event of the queue (without removing it from the queue)
// Inputs:	b = board number
//			ch = channel number
//			evnum = event number (0 is the 1st available)
// Outputs:	evnt = event data structure
// Return:	0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int PointEvent(int b, int ch, int evnum, GenericDPPEvent_t **evnt)
{
	*evnt = NULL;
	if (evnum >= (int)Qnev[b][ch])
		return 0;
	*evnt = &EvQ[b][ch][(EvRpnt[b][ch]+evnum) % EV_QUEUE_SIZE];
	return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// --------------------------------------------------------------------------------------------------------- 
static int PeekEvent(int b, int ch, GenericDPPEvent_t **evnt)
{
	*evnt = NULL;
	if (Qnev[b][ch] == 0)
		return 0;
	*evnt = &EvQ[b][ch][(EvRpnt[b][ch]) % EV_QUEUE_SIZE];
	return 1;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Set a pointer to one event of the queue (without removing it from the queue)
// Inputs:	b = board number
//			ch = channel number
//			evnum = event number (0 is the 1st available)
// Outputs:	evnt = event data structure
// Return:	0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
static int PopEvent(int b, int ch, GenericDPPEvent_t *evnt)
{
	if (Qnev[b][ch] == 0) return 0;
	memcpy(evnt, &EvQ[b][ch][EvRpnt[b][ch]], sizeof(GenericDPPEvent_t));
	Stats.CurrEventTime[b][ch] = EvQ[b][ch][EvRpnt[b][ch]].TimeStamp;
	EvRpnt[b][ch] = (EvRpnt[b][ch] + 1) % EV_QUEUE_SIZE;
	Stats.EvCnt[b][ch]++;
	if (!((b == WDcfg.TOFstartBoard) && (ch == WDcfg.TOFstartChannel)))
		Stats.EvCntSum++;
	Qnev[b][ch]--;
	if (ENROLOG) fprintf(rolog, "Read from ch %d: WP=%d RP=%d\n", ch, EvWpnt[b][ch], EvRpnt[b][ch]);
	return 1;
}





// --------------------------------------------------------------------------------------------------------- 
// Description: open output data file and write the header
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int OpenOutputDataFile() 
{
	char fname[200];
	char FileFormat = DATA_FILE_FORMAT;
	char *cfgimg;
	uint32_t nbcfg;
	uint32_t header[8];
	FILE *cfg;

	sprintf(fname, "%sRun%d.dat", WDcfg.DataFilePath, WDcfg.RunNumber);
	OutputDataFile = fopen(fname, "wb");
	if (OutputDataFile == NULL) {
		printf("Can't open Output Data File %s\n", fname);
		return -1;
	}

	// Write data file format
	fwrite(&FileFormat, 1, 1, OutputDataFile);

	// write data file header (1st word = header size)
	header[0] = 8;
	header[1] = WDcfg.RecordLength; 
	header[2] = WDcfg.DigitizerModel;
	header[3] = WDcfg.DppType;
	header[4] = WDcfg.NumBrd;
	header[5] = WDcfg.NumCh;
	header[6] = WDcfg.Tsampl;
	header[7] = WDcfg.Nbit;
	fwrite(header, sizeof(uint32_t), header[0], OutputDataFile);

	// write config files (appended) to the output data file
	cfgimg = (char *)malloc(1024*1024);
	cfg = fopen("_cfg.txt", "rb");
	if (cfg != NULL) {
		nbcfg = (int)fread(cfgimg, sizeof(char), 1024*1024, cfg);
		fwrite(&nbcfg, sizeof(uint32_t), 1, OutputDataFile);
		fwrite(cfgimg, sizeof(char), nbcfg, OutputDataFile);
	}
	fclose(cfg);
	free(cfgimg);
	return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: open input data file and read the header
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int OpenInputDataFile() 
{
	char fname[200];
	char DataFormat;
	char *cfgimg;
	uint32_t header[128];
	uint32_t nbcfg;
	FILE *cfg;

	sprintf(fname, "%s%s", WDcfg.DataFilePath, WDcfg.InputDataFileName);
	if (InputDataFile == NULL) {
		if (INPUT_DATA_FILE_TYPE == 0)
			InputDataFile = fopen(fname, "rb");
		else
			InputDataFile = fopen(fname, "r");
	} else {
		rewind(InputDataFile);
		emtstamp = lastemtstamp ;
	}

	if (InputDataFile == NULL) {
		printf("Can't open Input Data File %s\n", fname);
		return -1;
	}

	if (INPUT_DATA_FILE_TYPE == 0) {
		// read the header of the data file
		if (fread(&DataFormat, sizeof(char), 1, InputDataFile) < 1)  // Data format (version)
			return -1;
		if (DataFormat == 0) {
			fread(header, sizeof(uint32_t), 1, InputDataFile);
			fread(header+1, sizeof(uint32_t), header[0]-1, InputDataFile);
			WDcfg.RecordLength   = header[1]; 
			WDcfg.DigitizerModel = header[2];
			WDcfg.DppType        = header[3];
			WDcfg.NumBrd         = header[4];
			WDcfg.NumCh          = header[5];
			WDcfg.Tsampl         = header[6];
			WDcfg.Nbit           = header[7];
		}
		cfgimg = (char *)malloc(1024*1024);
		cfg = fopen("_cfg_datafile.txt", "wb");
		if (cfg != NULL) {
			fread(&nbcfg, sizeof(uint32_t), 1, InputDataFile);
			if (fread(cfgimg, sizeof(char), nbcfg, InputDataFile) < nbcfg) {
				free(cfgimg);
				return -1;
			}
			fwrite(cfgimg, sizeof(char), nbcfg, cfg);
		}
		fclose(cfg);
		free(cfgimg);
	} else {
		WDcfg.Tsampl = 2;
		WDcfg.Nbit = 14;
	}
	return 0;
}



// --------------------------------------------------------------------------------------------------------- 
// Description: Initialize variables and allocate memory buffers for the readout and the event queues
// Inputs:		-
// Outputs:		AllocatedMemSize = total num of bytes allocated
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int InitReadout(uint32_t *AllocatedMemSize)
{
	int ret=0, b, ch;
	int32_t AllocatedSize; 
	
	*AllocatedMemSize=0;

	if (ENROLOG) {
		rolog = fopen("ReadoutLog.txt", "w");
		printf("WARNING: Readout Log File enabled!!!\n");
	}

	// Init global variables
	for(b=0; b<MAX_NBRD; b++) {
		BoardEmpty[b] = 1;
		for (ch=0; ch<MAX_NCH; ch++) {
			EvWpnt[b][ch] = 0;
			EvRpnt[b][ch] = 0;
			EvRpnt[b][ch] = 0;
			EvWpnt[b][ch] = 0;
			Qnev[b][ch] = 0;
			AlmostFull[b] = 0;
		}
	}

	/* Allocate memory buffers and init statistics                                             */
	/* WARNING: The mallocs MUST be done after the digitizer programming, because the following functions 
	needs to know the digitizer configuration to allocate the right memory amount */
	if (!WDcfg.ReplayDataFile) {
		int tmpSize = 0;
		ret = CAEN_DGTZ_MallocReadoutBuffer(handle[0], &buffer, (uint32_t *)&AllocatedSize);
		*AllocatedMemSize += AllocatedSize;
		{
			int l,u;
			uint32_t EnableMask;
			void *tmpEvents[MAX_NCH];
			 Stats.AllocatedBoardHandleIndex = -1;
			for (l=0; l<WDcfg.NumBrd; l++) {
				ret |= CAEN_DGTZ_GetChannelEnableMask(handle[l], &EnableMask);
				ret |= CAEN_DGTZ_SetChannelEnableMask(handle[l], 0xFFFF);
				ret |= CAEN_DGTZ_MallocDPPEvents(handle[l], tmpEvents, (uint32_t *)&AllocatedSize);
				ret |= CAEN_DGTZ_SetChannelEnableMask(handle[l], EnableMask);
				if (AllocatedSize > tmpSize) {
					if (tmpSize != 0) CAEN_DGTZ_FreeDPPEvents(handle[Stats.AllocatedBoardHandleIndex], Events);
                     Stats.AllocatedBoardHandleIndex = l;
					for (u=0;u<MAX_NCH; u++) Events[u] =  tmpEvents[u];
					tmpSize = AllocatedSize;
				}

			}
		}
		*AllocatedMemSize += tmpSize;
	} else {
        wbuffer = (uint32_t *)malloc(WAVEFORM_BUFFER_SIZE);
        *AllocatedMemSize += WAVEFORM_BUFFER_SIZE;
    }
	if (WDcfg.WaveformProcessor) {
		AllocateWaveform(&PPwaveform, (uint32_t *)&AllocatedSize);
		*AllocatedMemSize += AllocatedSize;
	}

	// Allocate event queues
	for(b=0; b<MAX_NBRD; b++) {
		for (ch=0; ch<MAX_NCH; ch++) {
			EvQ[b][ch]=NULL; 
			EvQ[b][ch] = (GenericDPPEvent_t *)malloc(EV_QUEUE_SIZE * sizeof(GenericDPPEvent_t));
			*AllocatedMemSize += EV_QUEUE_SIZE * sizeof(GenericDPPEvent_t);
			if (EvQ[b][ch] == NULL)
				ret = -99;
		}
	}
	if (ret) {
		printf("Can't allocate memory buffers\n");
		return ret;
	}

	return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: free memory 
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int CloseReadout()
{
	int b, ch;
	if (!WDcfg.ReplayDataFile) {
		CAEN_DGTZ_FreeReadoutBuffer(&buffer);
		CAEN_DGTZ_FreeDPPEvents(handle[0], Events);
	} else {
		if (wbuffer != NULL)
			free(wbuffer);
    }

	for(b=0; b<MAX_NBRD; b++) {
		for (ch=0; ch<MAX_NCH; ch++) {
			if (EvQ[b][ch] != NULL)
				free(EvQ[b][ch]);
		}
	}
	if (InputDataFile != NULL)
		fclose(InputDataFile);
	if (OutputDataFile != NULL)
		fclose(OutputDataFile);
	return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: Read a data block from one board and put the events into the queue
// Inputs:		b = board to read
// Outputs		nb = num of bytes read
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ReadDataFromBoard(int b, uint32_t *nb)
{ 
	int ret = 0;
	int ch, ev;
	uint32_t NumEvents[MAX_NCH]; 

	memset(NumEvents, 0, sizeof(uint32_t)*MAX_NCH);
	ret = CAEN_DGTZ_ReadData(handle[b], CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, buffer, nb);
	if (ENROLOG) fprintf(rolog, "\n\nData Readout: %d bytes\n", *nb);
	if (ret) 
		return ret;
	if (*nb > 0) {
		// clear the waveform pointers to memory locations that have been overwritten by the ReadData
		for(ch=0; ch<MAX_NCH; ch++) {
			uint32_t rp = EvRpnt[b][ch];
			while(rp < EvWpnt[b][ch]) {
				if (EvQ[b][ch][rp].Waveforms < (uint32_t *)(buffer + *nb))
					EvQ[b][ch][rp].Waveforms = NULL;
				rp++;
			}
		}

		Stats.NumBlt++;
		ret = CAEN_DGTZ_GetDPPEvents(handle[b], buffer, *nb, Events, NumEvents);
		if (ret) 
			return ret;


		// Convert and copy events into the queue
		for(ch=0; ch<MAX_NCH; ch++) {
			if ((ENROLOG) && (NumEvents[ch]>0)) fprintf(rolog, "%d Events read for Brd %d - ch %d\n", NumEvents[ch], b, ch);
			Stats.TotNumEv += NumEvents[ch];
			for(ev=0; ev<(int)NumEvents[ch]; ev++) { 
				if (Qnev[b][ch] == EV_QUEUE_SIZE) {  // the queue is full
					if (ENROLOG) fprintf(rolog, "Lost event ");
					if (ev == (NumEvents[ch]-1)) {	
						GenericDPPEvent_t evnt;
						PreProcessEvent(b, ch, ev, Events[ch], &evnt);
						Stats.CurrEventTime[b][ch] = evnt.TimeStamp;
					}
				} else {
					PreProcessEvent(b, ch, ev, Events[ch], &EvQ[b][ch][EvWpnt[b][ch]]);
					if (WDcfg.WaveformProcessor) {
						ConvertWaveform(b, ch, EvQ[b][ch][EvWpnt[b][ch]], &PPwaveform);
						WaveformProcessor(b, ch, PPwaveform.Ns, PPwaveform.AnalogTrace[0], &EvQ[b][ch][EvWpnt[b][ch]], NULL);
					}
					if (ev == (NumEvents[ch]-1))
						Stats.CurrEventTime[b][ch] = EvQ[b][ch][EvWpnt[b][ch]].TimeStamp;

					if (DoSaveRawEvents) {
						char c[2];
						Waveform_t wfm;
						if (OutputDataFile == NULL) {
							OpenOutputDataFile();
						}
						if (OutputDataFile != NULL) {
							GenericDPPEvent_t oev;
							uint32_t wfmt;

							memcpy(&oev, &EvQ[b][ch][EvWpnt[b][ch]], sizeof(oev));
							if ((WDcfg.AcquisitionMode == ACQ_MODE_MIXED) || (WDcfg.AcquisitionMode == ACQ_MODE_OSCILLOSCOPE)) {  
								ConvertWaveform(b, ch, EvQ[b][ch][EvWpnt[b][ch]], &wfm);
								// HACK: change trace info...
								wfmt = (wfm.DualTrace << 31) || ((wfm.TraceSet[0] & 0xF) << 16) || ((wfm.TraceSet[1] & 0xF) << 20) ||
								       (wfm.TraceSet[2] & 0xF) || ((wfm.TraceSet[3] & 0xF) << 4) || ((wfm.TraceSet[4] & 0xF) << 8) || ((wfm.TraceSet[5] & 0xF) << 12);
								oev.Waveforms = (uint32_t *)(wfm.Ns);
							} else {
								wfmt = 0;
								oev.Waveforms = 0;
							}
							c[0] = (char)b; 
							c[1] = (char)ch;
							fwrite(c, sizeof(char), 2, OutputDataFile);
							fwrite(&oev, sizeof(GenericDPPEvent_t), 1, OutputDataFile);
							OutFileSize += 2 + sizeof(GenericDPPEvent_t);
							if ((WDcfg.AcquisitionMode == ACQ_MODE_MIXED) || (WDcfg.AcquisitionMode == ACQ_MODE_OSCILLOSCOPE)) { 
								fwrite(&wfmt, sizeof(wfmt), 1, OutputDataFile);
								fwrite(EvQ[b][ch][EvWpnt[b][ch]].Waveforms, sizeof(uint32_t), wfm.Ns/2, OutputDataFile);
								OutFileSize += sizeof(wfmt) + (wfm.Ns)*2;
							}
						}
					}
					if (ENROLOG) fprintf(rolog, "=> brd %d, ch %d: WP=%d RP=%d; TS=%u ", b, ch, EvWpnt[b][ch], EvRpnt[b][ch], EvQ[b][ch][EvWpnt[b][ch]].TimeStamp);
					EvWpnt[b][ch] = (EvWpnt[b][ch] + 1) % EV_QUEUE_SIZE;
					Qnev[b][ch]++;
					BoardEmpty[b] = 0;
				}
				if (ENROLOG) fprintf(rolog, "\n");
			}
		}
	}
	for(ch=0; ch<WDcfg.NumCh; ch++)
		if (Qnev[b][ch] > EV_QUEUE_ALMFULL_LEVEL)
			AlmostFull[b] |= (1<<ch);
	return 0;
}




// --------------------------------------------------------------------------------------------------------- 
// Description: Read a data block from data file and put the events into the queue
// Outputs:		nb = num of bytes read
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int ReadDataFromFile(uint32_t *nb)
{
	int ret = 0;
	int ch, b;
	uint32_t reclen;
	char c[2];
	int naggr = 0;
	GenericDPPEvent_t ev;

	*nb=0;
	if (InputDataFile == NULL) 	
		return 0;

	while(1) {
		int nbr;

		nbr = (int)fread(c, sizeof(char), 2, InputDataFile);  // read baord and channel number of the event (2 bytes)
		if (feof(InputDataFile) || (c[0]==-1)) {  // reached end of data file
			if (CONTINUOUS_REPLAY) {
				OpenInputDataFile();  // restart from begin
			} else {
				printf("End of file\n");
				break;  // quit loop
			}
			nbr = (int)fread(c, sizeof(char), 2, InputDataFile);  // read baord and channel number of the event (2 bytes)
		}
		if (INPUT_DATA_FILE_TYPE == 0) {  // raw data
			// Read board and channel number
			b = (int)c[0];
			ch = (int)c[1];
			if ((nbr < 2) || (b<0) || (b>MAX_NBRD) || (ch<0) || (ch>MAX_NCH)) 
				return -1;
			// Read one event from the file
			if ((int)fread(&ev, sizeof(GenericDPPEvent_t), 1, InputDataFile) < 1)
				return -1;
			reclen = (int)ev.Waveforms; // Num of 32 bit words for the waveform (0 if waveforms are disabled)
			ev.TimeStamp += emtstamp;
			if (lastemtstamp < ev.TimeStamp)
				lastemtstamp = ev.TimeStamp;
			if (reclen > 0) {  // read the waveform and format and write them to the buffer 
				wbuffer[wbpnt] = reclen;  // Num word of the waveform
				if ((int)fread(wbuffer + wbpnt + 1, sizeof(uint32_t), (reclen/2) + 1, InputDataFile) < ((reclen/2) + 1)) // waveform format (1 word) and data (nww words)
					return -1;
				ev.Waveforms = wbuffer + wbpnt;
			}
			*nb += sizeof(GenericDPPEvent_t);
		} else { // list file
			b = 0;
			ch = 0;
			wbpnt = 0;
			reclen = 0;
			*nb += fscanf(InputDataFile, "%ull", &ev.TimeStamp);
			*nb += fscanf(InputDataFile, "%u",	&ev.Energy);
			*nb += fscanf(InputDataFile, "%f", &ev.psd);
			ev.Flags = 0;
		}
		if (EvQ[b][ch] != NULL) {  // queue allocated (channel is not disabled)
			memcpy(&EvQ[b][ch][EvWpnt[b][ch]], &ev, sizeof(GenericDPPEvent_t));  // copy event to the queue
			Qnev[b][ch]++;
			Stats.TotNumEv++;
			naggr++;
			if (ENROLOG) fprintf(rolog, "=> brd %d, ch %d: WP=%d RP=%d; TS=%u\n", b, ch, EvWpnt[b][ch], EvRpnt[b][ch], EvQ[b][ch][EvWpnt[b][ch]].TimeStamp);
			EvWpnt[b][ch] = (EvWpnt[b][ch] + 1) % EV_QUEUE_SIZE;
			if ((Qnev[b][ch] == EV_QUEUE_SIZE) || (naggr == 64) || (wbpnt > (WAVEFORM_BUFFER_SIZE/4 - (reclen/2) - 2)))
				break;
		}
	}
	for(b=0; b<WDcfg.NumBrd; b++)
		for(ch=0; ch<WDcfg.NumCh; ch++)
			if (Qnev[b][ch] > EV_QUEUE_ALMFULL_LEVEL)
				AlmostFull[b] |= (1<<ch);
	wbpnt = 0;
	return 0;
}



#if 0
// --------------------------------------------------------------------------------------------------------- 
// Description: Get one event from the queue and advance read pointer
// Inputs:		b = board number
//				ch = channel number
// Outputs:		evnt = event data structure
// Return:		1=one event returned, 0=no event returned (empty), -1=error
// --------------------------------------------------------------------------------------------------------- 
int GetEvent(int b, int ch, GenericDPPEvent_t *evnt)
{
	int ret = 0, ret1, NeedNewData = 0;

	while(Qnev[b][ch] > 0) {
		if (!WDcfg.EnableTimeCorrelFilter && !WDcfg.EnableEnergyFilter && !WDcfg.EnablePSDFilter ) {
			ret = 1;
			break;
		} else if (IS_GOOD_EVENT(b, ch)) {  // Found good event; quit the loop and get it
			ret = 1;
			break;
		} else if (IS_BAD_EVENT(b, ch) || (IS_SUSPENDED_EVENT(b, ch) && AlmostFull[b])) {  // bad event (or suspended in an almost full condition) => discard it and try next one
			Stats.EvCnt[b][ch]++;
			Stats.CurrEventTime[b][ch] = EvQ[b][ch][EvRpnt[b][ch]].TimeStamp;
			EvRpnt[b][ch] = (EvRpnt[b][ch] + 1) % EV_QUEUE_SIZE;
			Qnev[b][ch]--;
			if (Qnev[b][ch] <= EV_QUEUE_ALMFULL_LEVEL)
				AlmostFull[b] &= ~(1<<ch);
			if (ENROLOG) fprintf(rolog, "Discarded event from queue in ch %d: WP=%d RP=%d\n", ch, EvWpnt[b][ch], EvRpnt[b][ch]);
			continue;
		} else if (IS_SUSPENDED_EVENT(b, ch))  {  // can't use this event; need to read more data to unblock it
			NeedNewData = 1;
			break;
		} else if (IS_UNPROCESSED_EVENT(b, ch)) {  // can't use this event; just wait for the next EventSelection call
			break;
		}
	}

	// If good event is found, get it from the queue
	if (ret == 1) {
		memcpy(evnt, &EvQ[b][ch][EvRpnt[b][ch]], sizeof(GenericDPPEvent_t));
		Stats.CurrEventTime[b][ch] = EvQ[b][ch][EvRpnt[b][ch]].TimeStamp;
		EvRpnt[b][ch] = (EvRpnt[b][ch] + 1) % EV_QUEUE_SIZE;
		Stats.EvCnt[b][ch]++;
		Qnev[b][ch]--;
		if (ENROLOG) fprintf(rolog, "Read from ch %d: WP=%d RP=%d\n", ch, EvWpnt[b][ch], EvRpnt[b][ch]);
	}

	// if the board is empty or need fresh data to unblock suspended events, then read more data and try to refill the queue;
	// Don't do that if the boards has at least one queue almost full (to avoid overflow)
	if ((BoardEmpty[b] || NeedNewData) && !AlmostFull[b]) { 
		uint32_t nb;
		if (!WDcfg.ReplayDataFile) 
			ret1 = ReadDataFromBoard(b, &nb);
		else
			ret1 = ReadDataFromFile(&nb);
		if (ret1 < 0) {
			printf("Data Readout Error\n");
			return -1;    
		}
		Stats.TotNb += nb;
	}

	// check for board empty (all queues must be empty)
	if (Qnev[b][ch] == 0) { 
		BoardEmpty[b] = 1;
		if (ENROLOG) fprintf(rolog, "Ch %d is empty\n", ch);
		for(ch=0; ch<WDcfg.NumCh; ch++)
			if (Qnev[b][ch] > 0)
				BoardEmpty[b] = 0;
	}
	if (Qnev[b][ch] <= EV_QUEUE_ALMFULL_LEVEL)
		AlmostFull[b] &= ~(1<<ch);
	return ret;
}
#endif



// --------------------------------------------------------------------------------------------------------- 
// Description: Get one event from the queue and advance read pointer
// Inputs:		b = board number
//				ch = channel number
// Outputs:		evnt = event data structure
// Return:		1=one event returned, 0=no event returned (empty), -1=error
// --------------------------------------------------------------------------------------------------------- 
int GetEvent(int b, int ch, GenericDPPEvent_t *evnt)
{
	int ret = 0, ret1, NeedNewData = 0;

	while(Qnev[b][ch] > 0) {
		if (!WDcfg.EnableTimeCorrelFilter && !WDcfg.EnableEnergyFilter && !WDcfg.EnablePSDFilter ) {
			ret = 1;
			break;
		} else if (IS_GOOD_EVENT(b, ch)) {  // Found good event; quit the loop and get it
			ret = 1;
			break;
		} else if (IS_BAD_EVENT(b, ch) || (IS_SUSPENDED_EVENT(b, ch) && AlmostFull[b])) {  // bad event (or suspended in an almost full condition) => discard it and try next one
			Stats.EvCnt[b][ch]++;
			Stats.CurrEventTime[b][ch] = EvQ[b][ch][EvRpnt[b][ch]].TimeStamp;
			EvRpnt[b][ch] = (EvRpnt[b][ch] + 1) % EV_QUEUE_SIZE;
			Qnev[b][ch]--;
			if (Qnev[b][ch] <= EV_QUEUE_ALMFULL_LEVEL)
				AlmostFull[b] &= ~(1<<ch);
			if (ENROLOG) fprintf(rolog, "Discarded event from queue in ch %d: WP=%d RP=%d\n", ch, EvWpnt[b][ch], EvRpnt[b][ch]);
			continue;
		} else if (IS_SUSPENDED_EVENT(b, ch))  {  // can't use this event; need to read more data to unblock it
			NeedNewData = 1;
			break;
		} else if (IS_UNPROCESSED_EVENT(b, ch)) {  // can't use this event; just wait for the next EventSelection call
			break;
		}
	}

	// If good event is found, get it from the queue
	if (ret == 1)
		PopEvent(b, ch, evnt);		

	// check for board empty (all queues must be empty)
	if (Qnev[b][ch] == 0) { 
		BoardEmpty[b] = 1;
		AllEmpty = 1;
		if (ENROLOG) fprintf(rolog, "Ch %d is empty\n", ch);
		for(b=0; b<WDcfg.NumBrd; b++) {
			for(ch=0; ch<WDcfg.NumCh; ch++) {
				if (Qnev[b][ch] > 0) {
					BoardEmpty[b] = 0;
					AllEmpty = 0;
					break;
				}
			}
			if (AllEmpty == 0)
				break;
		}
	} else {
		AllEmpty = 0;
	}


	// if the board is empty or need fresh data to unblock suspended events, then read more data and try to refill the queue;
	// Don't do that if the boards has at least one queue almost full (to avoid overflow)
	if (AllEmpty) { 
		uint32_t nb;
		if (!WDcfg.ReplayDataFile) {
			for(b=0; b<WDcfg.NumBrd; b++) {
				ret1 = ReadDataFromBoard(b, &nb);
				if (ret1 < 0) {
					printf("Data Readout Error\n");
					return -1;    
				}
				Stats.TotNb += nb;
			}
		} else {
			ret1 = ReadDataFromFile(&nb);
			Stats.TotNb += nb;
		}
	}

	if (Qnev[b][ch] <= EV_QUEUE_ALMFULL_LEVEL)
		AlmostFull[b] &= ~(1<<ch);
	return ret;
}

// --------------------------------------------------------------------------------------------------------- 
// --------------------------------------------------------------------------------------------------------- 
int IsAnyChannelAlmostFull()
{
	int b, ch;
	for(b=0; b<WDcfg.NumBrd; b++) {
		for(ch=0; ch<WDcfg.NumCh; ch++) {
			if (Qnev[b][ch] > EV_QUEUE_ALMFULL_LEVEL) {
				AlmostFull[b] = 1;
				return 1;
			}
		}
	}
	return 0;
}

int AllChannelsReady()
{
	int b, ch;
	for(b=0; b<WDcfg.NumBrd; b++) {
		for(ch=0; ch<WDcfg.NumCh; ch++) {
			if (WDcfg.EnableInput[b][ch] && (Qnev[b][ch] == 0))
				return 0;
		}
	}
	return 1;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: 
// Inputs:		
//				
// Outputs:		
// Return:		
// --------------------------------------------------------------------------------------------------------- 
int GetEvent_Correlated(int StartCh, int StartBrd, GenericDPPEvent_t evnt[MAX_NBRD][MAX_NCH], int EvRdy[MAX_NBRD][MAX_NCH])
{
	int ret = 0, ret1, NeedNewData = 0;
	int b, ch;
	int nevFound = 0;
	int NaIfound;
	GenericDPPEvent_t StartEvent, *StopEvent, Trash;
	uint64_t StartTime;

	memset(EvRdy, 0, MAX_NBRD*MAX_NCH*sizeof(int));

	while ((!IsAnyChannelAlmostFull()) && (!AllChannelsReady())) {
		uint32_t nb;
		if (!WDcfg.ReplayDataFile) {
			for(b=0; b<WDcfg.NumBrd; b++) {
				ret1 = ReadDataFromBoard(b, &nb);
				if (ret1 < 0) {
					printf("Data Readout Error\n");
					return -1;    
				}
				Stats.TotNb += nb;
			}
		} else {
			ret1 = ReadDataFromFile(&nb);
			Stats.TotNb += nb;
		}
		if (Qnev[StartBrd][StartCh] == 0) {
			for(b=0; b<WDcfg.NumBrd; b++) {
				for(ch=0; ch<WDcfg.NumCh; ch++) {
					while(PopEvent(b,ch,&Trash));
				}
			}
			break;
		}
	}


	if (PopEvent(StartBrd, StartCh, &StartEvent))   
		StartTime = StartEvent.TimeStamp;	
	else
		return 0;
	nevFound = 1;

	NaIfound = 0;
	for(b=0; b<WDcfg.NumBrd; b++) {
		for(ch=0; ch<WDcfg.NumCh; ch++) {
			if ((b == StartBrd) && (ch == StartCh)) continue;
			if (Qnev[b][ch] == 0) 
				continue;

			while(1) {
				if (!PeekEvent(b, ch, &StopEvent))
					break;
				/*if (WDcfg.EnableEnergyFilter && (((float)StopEvent.Energy < WDcfg.EnergyLLD[b][ch]) || (((float)StopEvent.Energy > WDcfg.EnergyULD[b][ch])))) {
					PopEvent(b, ch, &Trash);
					continue;
				}*/
				/*if (WDcfg.EnablePSDFilter && ((StopEvent.psd < WDcfg.PsdLLD[b][ch]) || ((StopEvent.psd > WDcfg.PsdULD[b][ch])))) {
					PopEvent(b, ch, &Trash);
					continue;
				}*/
				if (StopEvent->TimeStamp < StartTime) { // stop is too old; discard
					PopEvent(b, ch, &Trash);
					continue;
				}
				if (StopEvent->TimeStamp < (StartTime + WDcfg.TimeCorrelWindow)) {
					PopEvent(b, ch, &evnt[b][ch]);
					EvRdy[b][ch] = 1;
					nevFound++;
					if (b>0)  // this is a stop from a NaI detector (boards 1, 2 and 3)
						NaIfound = 1;
					memcpy(&evnt[StartBrd][StartCh], &StartEvent, sizeof(GenericDPPEvent_t));
					//evnt[StartBrd][StartCh].TimeStamp -= 300;
					EvRdy[StartBrd][StartCh] = 1;
					Stats.Matched[b][ch]++;
					break;
				}
				break;
			}
		}
	}

	// events without any stop from NaI are discarded (return empty event array)
	if (!NaIfound) {
		nevFound = 0;
		memset(EvRdy, 0, MAX_NBRD*MAX_NCH*sizeof(int));
	}


	// check for board empty (all queues must be empty)
	if (Qnev[b][ch] == 0) { 
		BoardEmpty[b] = 1;
		if (ENROLOG) fprintf(rolog, "Ch %d is empty\n", ch);
		for(b=0; b<WDcfg.NumBrd; b++) {
			for(ch=0; ch<WDcfg.NumCh; ch++) {
				if (Qnev[b][ch] > 0) {
					BoardEmpty[b] = 0;
					break;
				}
			}
		}
	}

	return nevFound;
}


