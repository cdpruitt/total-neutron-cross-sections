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
    std::vector<int> waveform; // "digital waveform of event" is a series of waveform samples for each event

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

    double completeTime; // digitizer assigned "raw time" (before corrections)
    unsigned int cycleNumber; // number of times the digitizer has cycled between DPP and waveform mode during this subrun
};

// used to sort events into channel-specific trees, but do no processing
struct SeparatedEvent
{
    unsigned int timetag; // 1 sample granularity, 32 bits
    unsigned int extTime;
    double fineTime; // provide additional bits of granularity 
    unsigned int eventNo;
    unsigned int evtType;
    unsigned int chNo;

    unsigned int sgQ; // integrated charge of event from short gate
    unsigned int lgQ; // integrated charge of event from long gate
    std::vector<int>* waveform = new std::vector<int>; // contains all waveform samples for each event to allow for corrections in analysis
};

// create struct for holding detector event data and link to trees
struct DetectorEvent
{
    double macroTime = 0;
    int macroNo = 0;
    int targetPos = 0;

    int cycleNumber = 0;
    double completeTime = 0;
    unsigned int timetag = 0;
    unsigned int extTime = 0;
    double fineTime = 0;
    int eventNo = 0;
    int sgQ = 0;
    int lgQ = 0;
    unsigned int baseline = 0;

    bool vetoed = false;

    std::vector<int> waveform;
};

struct MacropulseEvent
{
    MacropulseEvent() {}
    MacropulseEvent(
            int mn,
            int ne,
            int tp) :
        macroNo(mn), numberOfEventsInMacro(ne), targetPos(tp) {}

    int cycleNumber = 0;
    int macroNo = 0;
    double macroTime = 0;
    int targetPos = 0;
    int lgQ = 0;
    std::vector<int> waveform;

    int numberOfEventsInMacro = 0;
    int numberOfMonitorsInMacro = 0;

    bool isGoodMacro = 0;
};

// used to store processed events after they have been mated with a macropulse
struct ProcessedEvent
{
    unsigned int macroNo; // label each event by macropulse
    double macroTime; // provide the macropulse "zero" time
    unsigned int eventNo; // uniquely label each event in a macropulse
    double completeTime; // the event's 48-bit timestamp
    double fineTime; // the event's fine time
    int targetPos; // target position
    unsigned int sgQ, lgQ; // the event's short and long integrated charge gates
    std::vector<int> waveform; // waveform data for this event
};

// used to keep track of the macropulse structure of the sub-run.
struct TargetChangerEvent
{
    unsigned int macroNo; // label each event by macropulse
    unsigned int modeChange; // indicate the first event after a mode change
    double macroTime; // the event's time-zero reference (the macropulse start)
    double fineTime; // the event's calculated fine time
    int targetPos; // target position
    unsigned int lgQ; // the event's long integrated charge gate
    std::vector<int> waveform; // waveform data for this event
};

#endif
