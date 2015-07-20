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

/**
*  @author Ron Fox <ron@caentech.com>
*  @brief  Event building.
*/

#include <stdlib.h>

#include "evb.h"

extern DPP_Config_t	WDcfg;



//---------------------------------------------------------------------------------------
// Description: Adds a new fragment to the built event.
// Inputs:  b         - Board the Fragment comes from.
//          c         - Channel the Fragment comes from.
//          fragment - pointer to the fragment to add.
//          nf       - number of fragments already in the array.
//  In/outs pBuiltEvnts - pointer to the array of fragments being built.
// Return nf++ (new number of fragments).
//
int addFragment(int b, int c, pFragment_t pBuiltEvents, GenericDPPEvent_t *fragment, int nf)
{
	pBuiltEvents[nf].board = b;
	pBuiltEvents[nf].channel = c;
	memcpy(&(pBuiltEvents[nf].fragment), fragment, sizeof(GenericDPPEvent_t));

	return ++nf;
}

//---------------------------------------------------------------------------------------
// Description: Build coincident data into a 'full' event.
// Inputs:      mindt    - Minimum time differnce
//              deltaT       - The coincidence window for the build (in ns).
//              bTrigger - Board that has the trigger channel.
//              chTrigger - Channel in the trigger board that has the trigger channel.
// Outputs:    EventData  - A pointer to a dynamically allocated array of fragments that
//                          were put together...note that the actual storage allocated
//                          may be larger than the used storage.
//             nFrags     - Number of elements in EventData that contain fragments. 
// Return:   0  - Success, -1, failure.
//
// NOTE:  The caller, uses this via e.g.:
//        int ret, nFrags;
//        pFragment_t EventData;
//       ....
//        ret =  GetBuiltEvent(100, 1, 5, &EventData, &nFrags);
//        ... // Process the array of fragments in EventData (nFrags worth of them).
//        free(EventData);
//
int GetBuiltEvent(int bTrigger, int chTrigger, int EvRdy[MAX_NBRD][MAX_NCH], GenericDPPEvent_t evnt[MAX_NBRD][MAX_NCH], pFragment_t pBuiltEvents, int* nFrags)
{
	// Allocate the worst case event data:

	int b; 
	int c;
	int nf = 0;                           // Count of fragments.

    memset(pBuiltEvents,0,WDcfg.NumBrd*WDcfg.NumCh*sizeof(Fragment_t));
	// If there are no events in the trigger queue try to fill it

    if (EvRdy[bTrigger][chTrigger] == 1) 
	    nf = addFragment(bTrigger, chTrigger, pBuiltEvents, &(evnt[bTrigger][chTrigger]), nf);   // Add fragment to event.

	// Get data from the other queues...those that correspond to enabled inputs that is.

	for (b = 0; b < WDcfg.NumBrd; b++) { 
		for (c = 0; c < WDcfg.NumCh; c++) {
			if (((b != bTrigger) || (c != chTrigger))) {
                if (EvRdy[b][c] == 1) {
						nf = addFragment(b, c, pBuiltEvents,  &(evnt[b][c]), nf);
				}
			}
		}
	}
    *nFrags = nf;
	return 0;
}