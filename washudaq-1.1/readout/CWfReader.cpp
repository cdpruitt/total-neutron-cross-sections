/**
#******************************************************************************
#
# Via Vetraia, 11 - 55049 - Viareggio ITALY
# +390594388398 - www.caen.it
#
#***************************************************************************//**
# 

##
# @file CWfReader.cpp
# @brief Reader for waveforms.
# @author Ron Fox (rfoxkendo@gmail.com)

*/

#include "CWfReader.h"
#include "DPPConfig.h"
#include "config.h"
#include <stdlib.h>
#include <stdexcept>
#include <iostream>

/**
 * constructor
 *    @param handle - the handle to the open digitizer.
 *    @param config - The new configuration to load when we switch to waveform
 *                    mode.
 *    @param oldconfig - The configuration to restore when we switch back to
 *                       DPP Mode.
 *   @param  fd     - File descriptor open on the file to which we write.
 *   @param triggers - Number of triggers to take in this mode (might be 1).
 */
CWfReader::CWfReader(
    int handle,
    dictionary* config, dictionary* oldconfig, int fd, unsigned triggers
) :
    m_nHandle(handle),
    m_pWfConfig(config),
    m_pOldConfig(oldconfig),
    m_nFd(fd),
    m_nTriggers(triggers)
{
    
}

/**
 * setup
 *   Set up the digitizer.  Only some parameters can be modified from the base
 *   configuration in oldconfig:
 *
 *          Common parameters:
 *          
 *   *   RECORD_LENGTH - (global?) number of samples in the waveform.
 *   *   SELF_TRIGGER  - Self trigger mode (NORMAL/COINCIDENCE)
 *   *   MAX_NUM_EVENTS_BLT - Max number of events transferred per read.
 *   *   NEVT_AGGR       - Number of events required to claim data avail.
 *
 *           Per channel parameters:
 *           
 *   *   ENABLE_INPUT - enables/disables the input.
 *   *   PRE_TRIGGER   - Samples accepted before the trigger.
 *   *   CHANNEL_TRIGGER - DISABLED, ACQUISITION_ONLY, ACQUISITION_AND_TRGOUT
 *   *   TRG_THRESHOLD  - Trigger threshold level.
 *   *   TRG_HOLDOFF
 *   *   CFD_DELAY      - CFD delay parameter.
 *   *   CFD_ATTENUATION - CFD attenuation parameter
 *   *   DISCRIMINATOR_MODE - Discriminator mode: CFD or LED
 *
 *  @note Default values are taken from the old settings except for
 *        NEVT_AGGR(1) and MAX_NUM_EVENTS_BLT(1)
 *  @note All other parameters are not modified.
 *  @note It is assumed that the digitizer is stoppped on entry, on exit however
 *        data taking is active.
 *        
 */
void
CWfReader::setup()
{
    CAEN_DGTZ_ErrorCode status;
    
    
    int reg;

    status = CAEN_DGTZ_SetRecordLength(
            m_nHandle, getIntParam(m_pWfConfig, "RECORD_LENGTH", -1, getIntParam(
            m_pOldConfig, "RECORD_LENGTH", -1, 94)

        )
    );
    checkStatus(status, "Setting record Length [global]");    
    status = CAEN_DGTZ_SetDPPEventAggregation(m_nHandle, 0, 0);
    checkStatus(status, "Setting DPP Event aggregation");
    status = CAEN_DGTZ_SetNumEventsPerAggregate(
        m_nHandle, getIntParam(m_pWfConfig, "NEVT_AGGR", -1, getIntParam(
            m_pOldConfig, "NEVT_AGGR", -1, 0)
        )
    );
    checkStatus(status, "Setting numevents/agg");
    status = CAEN_DGTZ_SetMaxNumAggregatesBLT(
        m_nHandle, getIntParam(
            m_pWfConfig, "MAX_NUM_EVENTS_BLT", -1, getIntParam(
                m_pOldConfig, "MAX_NUM_EVENTS_BLT", -1, 1
            )
        )
    );
    checkStatus(status, "Setting max num agg/blt");
    
    // EXTERNAL trigger mode:
    
    status = CAEN_DGTZ_SetExtTriggerInputMode(
        m_nHandle, getExternalTrigger(m_pWfConfig)
    );
    checkStatus(status, "Setting the external trigger mode");
        
    // Per channel configuration:
    CAEN_DGTZ_BoardInfo_t info;
    status = CAEN_DGTZ_GetInfo(m_nHandle, &info);
    checkStatus(status, "Getting board info");
    CAEN_DGTZ_DPP_PSD_Params_t params;
    
    int chanMask =  0;
    params.trgho = getIntParam(
            m_pWfConfig, "TRG_HOLDOFF", -1, getIntParam(
                m_pOldConfig, "TRG_HOLDOFF", -1, 0
        )
    );
    params.blthr = getIntParam(m_pOldConfig, "PSD_BL_THRESHOLD", -1, 0);
    params.bltmo = getIntParam(m_pOldConfig, "PSD_BL_TIMEOUT", -1, 255);
    params.purh  = getPileupRejectionMode(m_pOldConfig);
    params.purgap= getIntParam(m_pOldConfig, "PSD_PUR_GAP", -1, 0);
    params.trgho = getIntParam(m_pOldConfig, "TRG_HOLDOFF", -1, 0);

    for (int c = 0; c < info.Channels; c++) {
        if(getParamString(m_pWfConfig, "ENABLE_INPUT", c, "NO") == "YES") {
            chanMask |= 1 << c;
         
            status = CAEN_DGTZ_SetDPPPreTriggerSize(
                m_nHandle, c, getIntParam(
                    m_pWfConfig, "PRE_TRIGGER", c, getIntParam(
                        m_pOldConfig, "PRE_TRIGGER", c, 0
                    )
                )
            );
            checkStatus(status, "Setting dpp pretrigger size(channel)");
            // DPP Parameters can only be set so we need to fish out the
            // current values from the DPP configuration in m_pOldConfig.
            
            params.thr[c] = getIntParam(
                m_pWfConfig, "TRG_THRESHOLD", c, getIntParam(
                    m_pOldConfig, "TRG_THRESHOLD", c, 50
                )
            );
            
            if (getParamString(
                m_pWfConfig, "CHANNEL_TRIGGER", c, getParamString(
                    m_pOldConfig, "CHANNEL_TRIGGER", c, "DISABLED"))    ==
                "DISABLED"
            ) {
                params.selft[c] = 0;
            } else {
                params.selft[c] = 1;
            }
            // CFD parameters.
            
            uint32_t trigReg;
            status = CAEN_DGTZ_ReadRegister(
                m_nHandle, CAEN_DGTZ_SAM_TRIGGER_REG_ADD + (c << 8), &trigReg
            );
            checkStatus(status, "gettting samtrigger reg value");
            trigReg = setBits(
                trigReg, 0, 7, getIntParam(
                    m_pWfConfig, "CFD_DELAY", c, getIntParam(
                        m_pOldConfig, "CFD_DELAY", c, 40
                    )
                )
            );
            trigReg = setBits(
                trigReg, 8,9, getIntParam(
                    m_pWfConfig, "CFD_ATTENUATION", c, getIntParam(
                        m_pOldConfig, "CFD_ATTENUATION", c, 0
                    )
                )
            );
            int cfdInterp = getIntParam(
                m_pWfConfig, "CFD_INTERPOLATE", c, getIntParam(
                    m_pOldConfig, "CFD_INTERPOLATE", c, 0
                )
            );
            trigReg = setBits(trigReg, 10,11, cfdInterp);
            CAEN_DGTZ_WriteRegister(
                m_nHandle, CAEN_DGTZ_SAM_TRIGGER_REG_ADD + (c << 8), trigReg
            );
            checkStatus(status, "Writing CFD parameters");
            // Do I need to save the interpolation?
            
            // Discriminator mode (LED or CFD):
            
            uint32_t dppCtlReg;
            int      mode;
            status = CAEN_DGTZ_ReadRegister(m_nHandle, 0x1080 + (c << 8), &dppCtlReg);
            checkStatus(status, "Getting reg with disc. mode");
            if (getParamString(
                m_pWfConfig, "DISCRIMINATOR_MODE", c, getParamString(
                    m_pOldConfig, "DISCRIMINATOR_MODE", c, "CFD"
                )
            ) == "CFD" ) {
                mode = 1;
            } else {
                mode = 0;
            }
            dppCtlReg = setBits(dppCtlReg, 6,6, mode);
            status = CAEN_DGTZ_WriteRegister(m_nHandle, 0x1080 + (c << 8), dppCtlReg);
            checkStatus(status, "Setting disc. mode");
            // DPP Parameters:
            
            params.lgate[c] = getIntParam(m_pOldConfig, "PSD_LONG_GATE", c, 60);
            params.sgate[c] = getIntParam(m_pOldConfig, "PSD_SHORT_GATE", c, 20);
            params.pgate[c] = getIntParam(m_pOldConfig, "PSD_PRE_GATE", c, 16);
            params.nsbl[c]  = getIntParam(m_pOldConfig, "PSD_BL_SAMPLES", c, 3);
            params.csens[c] = getIntParam(m_pOldConfig, "PSD_SEL_CHARGE_SENSE", c, 0);
            
        }
        status = CAEN_DGTZ_SetDPPParameters(m_nHandle, chanMask, &params);
        checkStatus(status, "Setting dpp params");
        status = CAEN_DGTZ_SetChannelEnableMask(m_nHandle, chanMask);
        checkStatus(status, "Setting channel enables");
       
    }
     // Enable data acquisition.
        
    status = CAEN_DGTZ_SetDPPAcquisitionMode(
        m_nHandle, CAEN_DGTZ_DPP_ACQ_MODE_Mixed,
        CAEN_DGTZ_DPP_SAVE_PARAM_TimeOnly
    );
    checkStatus(status, "Setting dpp acq mode");
    // checkStatus(CAEN_DGTZ_ClearData(m_nHandle), "Clearing the data");
    status = CAEN_DGTZ_SWStartAcquisition(m_nHandle);
    checkStatus(status, "Starting acq");
    std::cout << "WFmode acquisition started\n";
}
/**
 * ReadEvents
 *   Read at most the specified number of triggers-- writing them to file.
 */
int CWfReader::ReadEvents()
{
    // Allocate the event storage and  waveform storage:
    
    uint32_t rawBufferSize;
    char*    rawBuffer(0);
    CAEN_DGTZ_MallocReadoutBuffer(m_nHandle, &rawBuffer, &rawBufferSize);
    
    uint32_t dppBufferSize;
    CAEN_DGTZ_DPP_PSD_Event_t  *dppBuffer[CAEN_DGTZ_MAX_CHANNEL];
    int nDppEvents[CAEN_DGTZ_MAX_CHANNEL];
    CAEN_DGTZ_MallocDPPEvents(m_nHandle, (void**)dppBuffer, &dppBufferSize);
    
    CAEN_DGTZ_DPP_PSD_Waveforms_t  *pWaveform(0);
    uint32_t waveformSize;
    CAEN_DGTZ_MallocDPPWaveforms(m_nHandle, (void**)&pWaveform, &waveformSize);
    
    // The number of samples in a waveform:
    
    unsigned samples = getIntParam(
        m_pWfConfig, "RECORD_LENGTH", -1, getIntParam(
            m_pOldConfig, "RECORD_LENGTH", -1, 94
        )
    );
    
    // Read data until we have at least the requested number of triggers.
    
    CAEN_DGTZ_ErrorCode status;
    int trigsRead = 0;
    while (trigsRead < m_nTriggers) {
        uint32_t rawReadSize(0);
        uint32_t eventsRead[CAEN_DGTZ_MAX_CHANNEL];
        status = CAEN_DGTZ_ReadData(
            m_nHandle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT, rawBuffer,
            &rawReadSize
        );
        if (status != CAEN_DGTZ_Success) {
            break;
        }
        status = CAEN_DGTZ_GetDPPEvents(
            m_nHandle, rawBuffer, rawReadSize, (void**)dppBuffer, eventsRead
        );
        if (status != CAEN_DGTZ_Success) {
            break;
        }
        // Write the events:
        
        
        trigsRead += writeEvents(dppBuffer, eventsRead, pWaveform, samples);
        
        
    }
    // Free allocations:
    
    CAEN_DGTZ_FreeDPPWaveforms(m_nHandle, pWaveform);
    CAEN_DGTZ_FreeDPPEvents(m_nHandle, (void**)dppBuffer);
    CAEN_DGTZ_FreeReadoutBuffer(&rawBuffer);
    
    return trigsRead;
}
/**
 * finalize
 *   Restore the bits of the DPP setup that we butchered:
 */
void
CWfReader::finalize()
{
    
    CAEN_DGTZ_SWStopAcquisition(m_nHandle);
}    

/*----------------------------------------------------------------------------
 * Private utilities
 */

/**
 * writeEvents
 *    Writes waveform events to file.
 *    Set writeEvent for information about how this works.
 * @param dppBuffer - buffer of events from the decode.  This is an array of
 *                    pointers to the events from each channe.
 * @param eventsRead - An array, one element per channel, of counts of events in
 *                    each channel.
 * @param pWaveform - Space allocated for the waveform as decoded from each event.
 * @param samples   - Number of samples in each waveform.
 * @return int - number of triggers processed.
 */
int
CWfReader::writeEvents(
    CAEN_DGTZ_DPP_PSD_Event_t** dppBuffer, uint32_t* nRead,
    CAEN_DGTZ_DPP_PSD_Waveforms_t* pWaveform, unsigned int samples
)
{
    int nTriggers(0);
    int indices[CAEN_DGTZ_MAX_CHANNEL];
    memset(indices, 0, sizeof(indices));
    
    while(1) {
        int ch = findEarliestChannel(dppBuffer, nRead, indices);
        if (ch == -1) break;
        
        CAEN_DGTZ_DPP_PSD_Event_t* pEvent = & (dppBuffer[ch][indices[ch]]);
        indices[ch]++;
        CAEN_DGTZ_DecodeDPPWaveforms(m_nHandle, pEvent, pWaveform);
        
        writeEvent(ch, pEvent, pWaveform, samples);
        nTriggers++;
    }
    
    return nTriggers;
}
/**
 * writeEvent
 *   Write a single event to file. An event has the format:
 *     -  uint32_t size - Number of bytes in the event.
 *     -  uint32_t type - Type of event (2 for waveform).
 *     -  uint32_t channel - Channel number.
 *     -  uint32_t timeTag - Time tag at which the event was taken.
 *     -  uint32_t nSamples - Number of waveform samples.
 *     -  uint16_t samples[] - The signal waveform.
 * @param ch - The channel being written.
 * @param pEvent - Pointer to the DPP Event.
 * @param pWaveform - Pointrer to the waveform struct.
 * @param samples - Number of samples in the waveform.
 */
void
CWfReader::writeEvent(
    int ch, CAEN_DGTZ_DPP_PSD_Event_t* pEvent,
    CAEN_DGTZ_DPP_PSD_Waveforms_t* pWaveform, int nSamples
)
{
    // Compute the event size (self inclusive) and create the output bufer.
    
    uint32_t size = sizeof(uint32_t)*5 + sizeof(uint16_t)*nSamples;
    uint32_t type = 2;
    uint8_t  buffer[size];
    
    // Fill the event buffer 
    
    uint8_t* p = buffer;
    memcpy(p, &size, sizeof(uint32_t));            p += sizeof(uint32_t);
    memcpy(p, &type, sizeof(uint32_t));            p += sizeof(uint32_t);
    memcpy(p, &ch, sizeof(uint32_t));              p += sizeof(uint32_t);
    memcpy(p, &(pEvent->TimeTag), sizeof(uint32_t));  p += sizeof(uint32_t);
    memcpy(p, &pWaveform->Ns, sizeof(uint32_t));         p += sizeof(uint32_t);
    memcpy(p, pWaveform->Trace1, nSamples*sizeof(uint16_t));
    
    
    // Write the event:
    
    
    int n = write(m_nFd, buffer, size);
    if (n == -1) {
        throw std::runtime_error("Unable to write waveform event");
    }
    
}