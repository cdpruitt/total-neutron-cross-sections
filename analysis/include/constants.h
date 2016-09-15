#ifndef CONSTANTS_H
#define CONSTANTS_H

/* event variables */

// variables for holding tree event data to be histogrammed
unsigned int evtNo, macroNo, targetPos;
unsigned int sgQ, lgQ;
double completeTime, macroTime;

vector<int> *waveform; // for holding one event's waveform data

int numberGoodFits = 0;
int numberBadFits = 0;
int numberOnePeakFits = 0;
int numberOnePeakExpBackFits = 0; // Successfully fit as one peak riding on
// an exponential tail
int numberTwoPeakFits = 0;        // Successfully fit as two peaks

// total number of micropulses processed per target (for performing dead time
// calculation)
vector<long> microsPerTarget(6,0);
