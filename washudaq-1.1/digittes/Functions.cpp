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
    #include <time.h>
    #include <sys/timeb.h>
    #include <process.h>
#else
    #include <unistd.h>
    #include <stdint.h>   /* C99 compliant compilers: uint64_t */
    #include <ctype.h>    /* toupper() */
    #include <sys/time.h>
#endif

// Global Variables
uint64_t PlotTime=0;

int InitListFiles = 1;
FILE *list[MAX_NBRD][MAX_NCH];   // list output files


// --------------------------------------------------------------------------------------------------------- 
// Description: get time from the computer
// Return:		time in ms
// --------------------------------------------------------------------------------------------------------- 
long get_time()
{
    long time_ms;
#ifdef WIN32
    struct _timeb timebuffer;
    _ftime( &timebuffer );
    time_ms = (long)timebuffer.time * 1000 + (long)timebuffer.millitm;
#else
    struct timeval t1;
    struct timezone tz;
    gettimeofday(&t1, &tz);
    time_ms = (t1.tv_sec) * 1000 + t1.tv_usec / 1000;
#endif
    return time_ms;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: clear the console
// --------------------------------------------------------------------------------------------------------- 
void ClearScreen()
{
#ifdef WIN32
	system("cls");  
#else
	// HACK : how to clear screen in Linux?
#endif
}


// --------------------------------------------------------------------------------------------------------- 
// Description: Set bits of a registers 
// Inputs:		handle = board handle 
//				addr = address of the register
//				start_bit = first bit of the parameter being written 
//				end_bit = last bit of the parameter being written 
//				val: value to write
// Return:		0=OK, negative number = error code
// --------------------------------------------------------------------------------------------------------- 
int RegisterSetBits(int handle, uint16_t addr, int start_bit, int end_bit, int val)
{
	uint32_t mask=0, reg;
	int ret;
	int i;

	ret = CAEN_DGTZ_ReadRegister(handle, addr, &reg);   
	for(i=start_bit; i<=end_bit; i++)
		mask |= 1<<i;
	reg = reg & ~mask | ((val<<start_bit) & mask);
	ret |= CAEN_DGTZ_WriteRegister(handle, addr, reg);   
	return ret;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Save an histogram to output file
// Inputs:		FileName = filename radix (ch and board index will be added)
//				Nbin = number of bins
//				Histo = histogram to save
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int SaveHistogram(char *FileName, int Nbin, uint32_t *Histo)
{
    FILE *fh, *ansi42;
    int i;
    char str[200];

	fh = fopen(FileName, "w");
    if (fh == NULL)
		return -1;
	if (WDcfg.HistoOutputFormat == HISTO_FILE_FORMAT_ANSI42) {
		ansi42 = fopen("ansi42template.txt", "r");
		if (ansi42 != NULL) {
			while(!feof(ansi42)) {
				fgets(str, 200, ansi42);
				if (strstr(str, "*PutChannelDataHere*")) {
					for(i=0; i<Nbin; i++) 
						fprintf(fh, "%d\n", Histo[i]);
				} else {
					fprintf(fh, "%s", str);
				}
			}
			fclose(ansi42);
		}
	} else if (WDcfg.HistoOutputFormat == HISTO_FILE_FORMAT_1COL) {
		for(i=0; i<Nbin; i++) 
			fprintf(fh, "%d\n", Histo[i]);
	} else if (WDcfg.HistoOutputFormat == HISTO_FILE_FORMAT_2COL) {
		for(i=0; i<Nbin; i++) 
			fprintf(fh, "%d %d\n", i, Histo[i]);
	}
    fclose(fh);
    return 0;
}


// --------------------------------------------------------------------------------------------------------- 
// Description: Save all histograms to output file
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int SaveAllHistograms() 
{
	int b, ch, ret=0;
	char prefix[300], fname[300], fext[10];

	if (WDcfg.RunNumberInDataFiles)
		sprintf(prefix, "%sRun%d_", WDcfg.DataFilePath, WDcfg.RunNumber);
	else
		sprintf(prefix, "%s", WDcfg.DataFilePath);

	if (WDcfg.HistoOutputFormat == HISTO_FILE_FORMAT_ANSI42)
		sprintf(fext, "n42");
	else
		sprintf(fext, "txt");

	for(b=0; b<WDcfg.NumBrd; b++) {
		for(ch=0; ch<WDcfg.NumCh; ch++) {
			if (WDcfg.EnableInput[b][ch]) {
				/* Save Histograms to file for each board/channel */
				sprintf(fname, "%sEhisto_%d_%d.%s", prefix, b, ch, fext);
				ret |= SaveHistogram(fname, WDcfg.EHnbin, Stats.EHisto[b][ch]);
				sprintf(fname, "%sThisto_%d_%d.%s", prefix, b, ch, fext);
				ret |= SaveHistogram(fname, WDcfg.THnbin, Stats.THisto[b][ch]);
				/*sprintf(fname, "%sPSDhisto_%d_%d.%s", prefix, b, ch, fext);
				ret |= SaveHistogram(fname, 1024, Stats.PSDHisto[b][ch]);*/
			}
		}
	}
	if (NUM_VIRTUAL_BOARDS) {
		sprintf(fname, "%sEsum.%s", prefix, fext);
		ret |= SaveHistogram(fname, WDcfg.EHnbin, Stats.EHisto[WDcfg.NumBrd][0]);
		sprintf(fname, "%sTsum.%s", prefix, fext);
		ret |= SaveHistogram(fname, WDcfg.THnbin, Stats.THisto[WDcfg.NumBrd][0]);
	}
	return ret;
}



// --------------------------------------------------------------------------------------------------------- 
// Description: Save one event to the List file
// Inputs:		list = output files
//				b = board index
//				ch = channel
//				timestamp = timestamp of the event
//				energy = energy of the event
//				format = OUTFILE_ASCII or OUTFILE_BINARY
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int SaveList(int b, int ch, uint64_t timestamp, uint16_t energy, float psd)
{
	char prefix[300];

	if (InitListFiles) {
		int bb, cc;
		InitListFiles = 0;
		for(bb=0; bb<MAX_NBRD; bb++) 
			for(cc=0; cc<MAX_NCH; cc++)
				list[bb][cc] = NULL;
	}

	if (WDcfg.RunNumberInDataFiles)
		sprintf(prefix, "%sRun%d_", WDcfg.DataFilePath, WDcfg.RunNumber);
	else
		sprintf(prefix, "%s", WDcfg.DataFilePath);


	if (list[b][ch]==NULL) {
		char fname[100];					
		if (WDcfg.OutFileFormat == OUTFILE_ASCII) {
			sprintf(fname, "%sList_%d_%d.txt", prefix, b, ch);
			list[b][ch] = fopen(fname, "w");
		} else {
			sprintf(fname, "%sList_%d_%d.dat", prefix, b, ch);
			list[b][ch] = fopen(fname, "wb");
		}
		if (list[b][ch]==NULL)
			return -1;
	}

	if (WDcfg.OutFileFormat == OUTFILE_ASCII) {
		fprintf(list[b][ch], "%lld %d %f\n", timestamp, energy, psd);
	} else {
		char buff[100];
		memcpy(buff, &timestamp, sizeof(timestamp));
		memcpy(buff+sizeof(timestamp), &energy, sizeof(energy));
		fwrite(buff, 1, sizeof(timestamp)+sizeof(energy), list[b][ch]);
	}
	return 0;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Plot the waveforms of one event
// Inputs:		wplot = gnuplot pipe
//				Wfm = Event to plot
//				title = title of the plot
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int PlotWaveforms(FILE *wplot, Waveform_t Wfm, char *title)
{
	int i, cl=1;
	char comma;
	FILE *wpdata;

	wpdata = fopen("PlotData.txt", "w");
	if (wpdata == NULL) {
		printf("Can't open plot data\n");
		return -1;
	}
								
	/* Save Waveform and Plot it using gnuplot */
	for(i=0; i<(int)Wfm.Ns; i++) {
		if (Wfm.TraceSet[0] != -1) fprintf(wpdata, "%d ", (short)Wfm.AnalogTrace[0][i]);  // analog trace 1
		if (Wfm.TraceSet[1] != -1) fprintf(wpdata, "%d ", (short)Wfm.AnalogTrace[1][i]);  // analog trace 2
		if (Wfm.TraceSet[2] != -1) fprintf(wpdata, "%d ", Wfm.DigitalTrace[0][i]*100 +  50);  // digital trace 0
		if (Wfm.TraceSet[3] != -1) fprintf(wpdata, "%d ", Wfm.DigitalTrace[1][i]*100 + 250);  // digital trace 1
		if (Wfm.TraceSet[4] != -1) fprintf(wpdata, "%d ", Wfm.DigitalTrace[2][i]*100 + 450);  // digital trace 2
		if (Wfm.TraceSet[5] != -1) fprintf(wpdata, "%d ", Wfm.DigitalTrace[3][i]*100 + 650);  // digital trace 3
		fprintf(wpdata, "\n");
	}
	fclose(wpdata);
	fprintf(wplot, "set title '1: %s'\n", title);
	fprintf(wplot, "plot");
	comma = ' '; // first command after "plot" doesn't have comma
	for(i=0; i<MAX_NTRACES; i++) {
		if (Wfm.TraceSet[i] != -1) {
			fprintf(wplot, "%c 'PlotData.txt' u ($0*%f):%d t 'T%d: %s' w step ls %d", comma, (float)WDcfg.Tsampl/1000.0, cl++, i+1, TraceNames[i][Wfm.TraceSet[i]], i+1);
			comma = ',';
		}
	}
	fprintf(wplot, "\n");
	fflush(wplot);
	return 0;
}



// --------------------------------------------------------------------------------------------------------- 
// Description: Plot one histogram
// Inputs:		hplot = gnuplot pipe
//				Histo = histogram to plot
//				Nbin = number of bins of the histogram
//				title = title of the plot
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int PlotHisto(FILE *hplot, uint32_t *Histo , int Nbin, float Xmin, float Xmax, char *title, char *xlabel)
{
	int i;
	FILE *phdata;
	float xa, xb;

	xb = Xmin;
	xa = (Xmax-Xmin)/Nbin;
	phdata = fopen("PlotHistoData.txt", "w");
	for(i=0; i<Nbin; i++)
		fprintf(phdata, "%d  \n", Histo[i]);
	fclose(phdata);
	fprintf(hplot, "set title '%s'\n", title);
    fprintf(hplot, "set xlabel '%s'\n", xlabel);
    fprintf(hplot, "set ylabel 'Counts'\n");
	fprintf(hplot, "plot 'PlotHistoData.txt' using ($0*%f+%f):($1) title 'BinSize = %f' with step\n", xa, xb, xa);
	fflush(hplot);
	return 0;
}

// --------------------------------------------------------------------------------------------------------- 
// Description: Plot 2D histogram
// Inputs:		hplot = gnuplot pipe
//				Histo2D = Histogram to plot
//				nx, ny = number of bins of the histogram (x and y axes)
//				title = title of the plot
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int PlotHisto2D(FILE *hplot, uint32_t *Histo2D, int nx, int ny, char *title)
{
    int i,j;
	FILE *w1;
	float xbinsize=1.0; 
	float ybinsize=1.0/512;
	uint64_t CurrentTime = get_time();

	if ((CurrentTime - PlotTime) < 5000)
		return 0;
    fprintf(hplot, "set xlabel 'Total Energy'\n");
    fprintf(hplot, "set ylabel 'PSD (Total-Fast)/Total'\n");
	fprintf(hplot, "set title '%s'\n", title);
	//fprintf(hplot, "set xrange [0:%d]\n", (int)(nx*xbinsize));
	//fprintf(hplot, "set yrange [0:%d]\n", (int)(ny*ybinsize));
	w1 = fopen("histogram2d.txt", "w");			
	for(i=0; i<nx; i++) {
		for(j=0; j<ny; j++)
			fprintf(w1, "%f %f %d\n", i*xbinsize, j*ybinsize, Histo2D[i*ny+j]);
		fprintf(w1, "\n");
	}
	fclose(w1);
	Sleep(10);
	//fprintf(hplot, "unset grid; set palette model CMY rgbformulae 7,5,15\n");
	fprintf(hplot, "unset grid; set palette model CMY rgbformulae 15,7,3\n");
	fprintf(hplot, "plot 'histogram2d.txt' with image\n");
	fflush(hplot);	
	PlotTime = CurrentTime;
	
    return 0;
}



// --------------------------------------------------------------------------------------------------------- 
// Description: Save all the regsiters of the borad to a file
// Inputs:		handle = handle of the board
// Return:		0=OK, -1=error
// --------------------------------------------------------------------------------------------------------- 
int SaveRegImage(int handle) 
{
	FILE *regs;
	char fname[100];
	int ret;
	uint32_t addr, reg, ch;

	sprintf(fname, "reg_image_%d.txt", handle);
	regs=fopen(fname, "w");
	if (regs==NULL)
		return -1;

	fprintf(regs, "[COMMON REGS]\n");
	for(addr=0x8100; addr <= 0x8200; addr += 4) {
		ret = CAEN_DGTZ_ReadRegister(handle, addr, &reg);
		if (ret==0)
			fprintf(regs, "%04X : %08X\n", addr, reg);
	}
	for(addr=0xEF00; addr <= 0xEF34; addr += 4) {
		ret = CAEN_DGTZ_ReadRegister(handle, addr, &reg);
		if (ret==0)
			fprintf(regs, "%04X : %08X\n", addr, reg);
	}
	for(addr=0xF000; addr <= 0xF088; addr += 4) {
		ret = CAEN_DGTZ_ReadRegister(handle, addr, &reg);
		if (ret==0)
			fprintf(regs, "%04X : %08X\n", addr, reg);
	}

	for(ch=0; ch<8; ch++) {
		fprintf(regs, "[CHANNEL %d]\n", ch);
		for(addr=0x1000+(ch<<8); addr <= (0x10FF+(ch<<8)); addr += 4) {
			ret = CAEN_DGTZ_ReadRegister(handle, addr, &reg);
			if (ret==0)
				fprintf(regs, "%04X : %08X\n", addr, reg);
		}
	}

	fclose(regs);
	return 0;
}



