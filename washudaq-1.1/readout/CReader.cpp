/**
#******************************************************************************
#
# Via Vetraia, 11 - 55049 - Viareggio ITALY
# +390594388398 - www.caen.it
#
#***************************************************************************//**
# 

##
# @file CReader.cpp
# @brief common utility functions for all readers
# @author Ron Fox (rfoxkendo@gmail.com)

*/
#include "CReader.h"

#include <stdint.h>
#include <stdlib.h>
#include <iostream>

/**
 * setBits
 *   Set a bit field in an existing value
 *
 *  @param value - Value to modify
 *  @param startBit - First bit number of field to modify (from 0)
 *  @param endBit   - Last bit number of field to modify (inclusive).
 *  @param bits     - Bits to stick into that field.
 *  @return modified 'value'.
 */
uint32_t
CReader::setBits(
    uint32_t value, unsigned startBit, unsigned endBit, uint32_t bits
)
{
    // Build a mask of the bits in the field:
    
    uint32_t mask = 0;
    for (int i = startBit; i <= endBit; i++) {
        mask |= 1 << i;
    }
    // Clear that bit field in value:
    
    value &= ~mask;
    
    // Or in the new field:
    
    value |= (bits << startBit) & mask;
    
    return value;
}
/**
 * findEarliestChannel
 *    Determines which channel has the smallest coarse timestamp.  This is
 *    quite stupid, not taking into account the potentials for a wrap around
 *    to muck up ordering.
 *
 *  @param dppBuffer - the decoded DPP events each channel has a  pointer
 *                     to a set of events.
 *  @param nDppEvents - Array of number of events in each channel.
 *  @param indices    - Array of indices into the next event for each channel.
 *  @return int   - If >= 0 a channel number if < 0 there are no more events.
 *                  (all indices are at the end of the event list for all channels).
 */
int
CReader::findEarliestChannel(
    CAEN_DGTZ_DPP_PSD_Event_t** dppBuffer, uint32_t* nDppEvents, int* indices
)
{
    int channel       = -1;
    uint64_t earliest = INT64_MAX;    // coarse stamp is unsigned 32 bits...

    for (int i = 0; i < CAEN_DGTZ_MAX_CHANNEL; i++) {
        if (indices[i] < nDppEvents[i] ) {
            uint32_t time = dppBuffer[i][indices[i]].TimeTag;
            if (time < earliest) {
                earliest = time;
                channel  = i;
            }
        }
    }
    return channel;
}
/**
 * checkStatus
 * 
 *    If status is not CAEN_DGTZ_Success error exit with message:
 *
 * @param status - Status from a CAEN_DGTZ call.
 * @param msg    - msg - message emitted if status != CAEN_DGTZ_Success.
 */
void
CReader::checkStatus(CAEN_DGTZ_ErrorCode code, const char* msg)
{
    if (code != CAEN_DGTZ_Success) {
        std::cerr << "Error : " << code << std::endl;
        std::cerr << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}