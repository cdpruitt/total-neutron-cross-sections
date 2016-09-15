#ifndef RAW_H
#define RAW_H

// Raw data is stored as hexadecimal words (16 bits long) in the .evt files
int const BUFFER_WORDS = 1; // number of chars per buffer word
int const BUFFER_BYTES = BUFFER_WORDS*2; // number of bytes per buffer word

// keep track of event statistics
long rawNumberOfEvents = 0;
long rawNumberOfDPPs = 0;
long rawNumberOfWaveforms = 0;
long rawNumberOfCh0Waveforms = 0;
long rawNumberOfCh2Waveforms = 0;
long rawNumberOfCh4Waveforms = 0;

#endif
