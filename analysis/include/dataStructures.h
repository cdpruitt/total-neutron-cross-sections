#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <vector>

/******************************************************************************/
/* Event structure */

// To populate events from the raw file into a ROOT tree, we provide
// an event structure listing all the variables that comprise an event.
// When a new event is added to the tree, a subset of these variables is
// included in the new event being saved.
struct RawEvent
{
    // Event header variables:
    unsigned int size;    // size of event, in bytes
    unsigned int evtType; // "event type", either 1 (DPP) or 2 (waveform)
    unsigned int chNo;    // "channel number" indicates event origin (i.e.,
    // detector, monitor, target changer)
    unsigned int timetag;       // "coarse timestamp of event" includes 32 bits of time
    // information. Units are the same as the sample period of
    // the digitizer (e.g. 2 ns or 5 ns)

    // Event body variables:
    unsigned int extraSelect; // indicates the meaning of the "extras" words
    unsigned int sgQ;     // "short gate integrated charge" provides the charge
    // integral over an adjustable range of the event's peak 
    unsigned int lgQ;     // "long gate integrated charge" provides the charge
    // integral over an adjustable range of the event's peak 
    unsigned int nSamp; // number of samples in the event's waveform
    std::vector<int>* waveform; // "digital waveform of event" is a series of waveform
                                // samples for each event

    // Variables extracted from "extras"
    unsigned int baseline;
    unsigned int flags;
    unsigned int extTime; // "extended timestamp of event" extends timetag with
    // 16 additional bits for times greater than 2^32 sample periods.
    double fineTime;// "fine timestamp of event" sub-divides timetag with
    // 10 additional bits of time granularity, with units of
    // (sample period)/2^10 units.
    unsigned int PZC; // positive zero-crossing (for manual CFD calculation)
    unsigned int NZC; // negative zero-crossing (for manual CFD calculation)
};

// used to sort events into channel-specific trees, but do no processing
struct SeparatedEvent
{
    unsigned int timetag; // 1 sample granularity, 32 bits
    unsigned int extTime;
    double fineTime; // provide additional bits of granularity 
    unsigned int evtNo;
    unsigned int evtType;
    unsigned int chNo;

    unsigned int sgQ; // integrated charge of event from short gate
    unsigned int lgQ; // integrated charge of event from long gate
    std::vector<int>* waveform; // contains all waveform samples for each event to allow for corrections in analysis
};

// used to store processed events after they have been mated with a macropulse
struct ProcessedEvent
{
    unsigned int macroNo; // label each event by macropulse
    long double macroTime; // provide the macropulse "zero" time
    unsigned int evtNo; // uniquely label each event in a macropulse
    long double completeTime; // the event's 48-bit timestamp
    double fineTime; // the event's fine time
    unsigned int targetPos; // target position
    unsigned int sgQ, lgQ; // the event's short and long integrated charge gates
    std::vector<int> *waveform; // waveform data for this event
};

// used to keep track of the macropulse structure of the sub-run.
struct TargetChangerEvent
{
    unsigned int macroNo; // label each event by macropulse
    unsigned int modeChange; // indicate the first event after a mode change
    long double macroTime; // the event's time-zero reference (the macropulse start)
    double fineTime; // the event's calculated fine time
    unsigned int targetPos; // target position
    unsigned int lgQ; // the event's long integrated charge gate
    std::vector<int> *waveform; // waveform data for this event
};

#endif
