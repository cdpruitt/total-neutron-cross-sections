/******************************************************************************
*
* CAEN SpA - System integration division
* Via Vetraia, 11 - 55049 - Viareggio ITALY
* +390594388398 - www.caen.it
*
***************************************************************************/
/**
* \note TERMS OF USE:
* This program is free software; you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the Free Software
* Foundation. This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The user relies on the
* software, documentation and results solely at his own risk.
******************************************************************************/

/**
 * @file CDPPReader.cpp
 * @brief Implements the class that reads data in DPP mode.
 */


#include "CDPPReader.h"
#include <CAENDigitizer.h>
#include "config.h"
#include "DPPConfig.h"
#include <iostream>
#include <cstdlib>
#include <stdexcept>



/**
 * constructor
 *
 * @param handle - Handle to the digitizer.
 * @param config - Dictionary pointer to the parsed configuration.
 * @param fd     - File descriptor to which data are written.
 * @param triggers - Number of triggers to read before returning on ReadEvents().
 */

CDPPReader::CDPPReader(int handle, dictionary* config, int fd, unsigned triggers) :
  m_nHandle(handle),
  m_pConfig(config),
  m_nFd(fd),
  m_nTriggers(triggers) {}

/**
 * setup
 *   Prepares the digitizer for readout.  This is where/when the configuration is
 *   processed and the digitizer set up in accordance with that configuration.
 *
 * @note - we are assuming the digitizer is a V1730.
 */
void
CDPPReader::setup()
{
    std::cout << "Setting up DPP Mode\n";
        
    int reg;
    CAEN_DGTZ_ErrorCode status;
    
    CAEN_DGTZ_BoardInfo_t info;
    status = CAEN_DGTZ_GetInfo(m_nHandle, &info);
    checkStatus(status, "Getting board info");
    unsigned nChans = info.Channels;
    std::cerr << "Board with " << nChans << std::endl;
    
    // Desk top digitizers have a FAN_SPEED param (global)...not named in header
    
    status = CAEN_DGTZ_WriteRegister(
        m_nHandle, 0x8168, (getIntParam(m_pConfig, "FAN_SPEED", -1, 0) << 3) & 8
    );
    checkStatus(status, "Setting fan speed");
    // ACQUISITION_MODE
    
    status = CAEN_DGTZ_SetDPPAcquisitionMode(
        m_nHandle, getAcquisitionMode(m_pConfig),
        CAEN_DGTZ_DPP_SAVE_PARAM_EnergyAndTime
    );
    checkStatus(status, "Setting acquisition mode.");
    
    // Baseline mode on by default
    
    if (getIntParam(m_pConfig, "PSD_SEL_BASELINE", -1, 1)) {
        reg = CAEN_DGTZ_BROAD_CH_CONFIGBIT_SET_ADD;
    } else {
        reg = CAEN_DGTZ_BROAD_CH_CLEAR_CTRL_ADD;
    }
    status = CAEN_DGTZ_WriteRegister(m_nHandle, reg, 1 << 17);
    checkStatus(status, "Setting PSD_SEL_BASELINE");
    // Record length:
    
    status = CAEN_DGTZ_SetRecordLength(
        m_nHandle, getIntParam(m_pConfig, "RECORD_LENGTH", -1, 96), -1
    );
    checkStatus(status, "Setting record length");
    // Figure out what to do with the TRIGGER_MODE(NORMAL, COINCIDENCE).
    // If coincidence mode, the trigger validation direction sigs must be reversed.
    
    if (getParamString(m_pConfig, "TRIGGER_MODE", -1, "NORMAL") == "NORMAL") {
        reg = CAEN_DGTZ_BROAD_CH_CLEAR_CTRL_ADD;            // BIt clear.
    } else {
        reg = CAEN_DGTZ_BROAD_CH_CONFIGBIT_SET_ADD;        // BIt set
    }
    status = CAEN_DGTZ_WriteRegister(m_nHandle, reg, 1 << 2 );    // 
    checkStatus(status, "Setting self trigger mode");
    
    // SELF_TRIGGER itself??
    
    // Front panel signal levels:
    
    status = CAEN_DGTZ_SetIOLevel(m_nHandle, getFpLevel(m_pConfig));
    checkStatus(status, "Setting I/O level");
    
    // GATED_START means start via front panel..otherwise start is via software:
    
    if (getParamString(m_pConfig, "GATED_START", -1, "DISABLED") == "DISABLED") {
        status = CAEN_DGTZ_SetAcquisitionMode(m_nHandle, CAEN_DGTZ_SW_CONTROLLED);
        checkStatus(status, "Setting acquisition mode SW controlled");
        status = CAEN_DGTZ_SetRunSynchronizationMode(
            m_nHandle, CAEN_DGTZ_RUN_SYNC_Disabled
        );
        checkStatus(status, "Disabling synchronized start");
        
    } else {
        status = CAEN_DGTZ_SetAcquisitionMode(m_nHandle, CAEN_DGTZ_S_IN_CONTROLLED);
        checkStatus(status, "Setting Acq mode to S_IN controlled");
        status = CAEN_DGTZ_WriteRegister(m_nHandle, 0x8170, 0);   // Only one board so no delay.
        checkStatus(status, "Setting delay 0");
    }
    // Is the external trigger enabled or only internal?  We don't allow
    // software triggers.
    
    status = CAEN_DGTZ_SetExtTriggerInputMode(m_nHandle, getExternalTrigger(m_pConfig));
    checkStatus(status, "Setting the external trigger mode");
    
    // Event aggregation (NEVT_AGGR and MAX_NUM_AGGREGATES_BLT):
    // My opinion - the second parameter is the one that really affects performance
    // not the first.
    
    status = CAEN_DGTZ_SetDPPEventAggregation(m_nHandle, 0, 0);
    checkStatus(status, "Setting DPPEvent aggregation");
    
    status = CAEN_DGTZ_SetNumEventsPerAggregate(
        m_nHandle, getIntParam(m_pConfig, "NEVT_AGGR", -1, 0)
    );
    checkStatus(status, "Setting num events/aggregate");
    
    status = CAEN_DGTZ_SetMaxNumAggregatesBLT(
        m_nHandle, getIntParam(m_pConfig, "MAX_NUM_EVENTS_BLT", -1, 0)
    );
    checkStatus(status, "Setting max number of events/blt");
    
    // Per channel settings;
    

    
    uint32_t chanMask(0);
    CAEN_DGTZ_DPP_PSD_Params_t params;
    params.trgho = getIntParam(m_pConfig, "TRG_HOLDOFF", -1, 0);
    params.purh  = getPileupRejectionMode(m_pConfig);
    params.purgap= getIntParam(m_pConfig, "PSD_PUR_GAP", -1, 0);
    params.blthr = getIntParam(m_pConfig, "PSD_BL_THRESHOLD", -1, 0);
    params.bltmo = getIntParam(m_pConfig, "PSD_BL_TIMEOUT", -1, 255);

    for (int  c = 0; c < nChans; c++) {
        // DC Offset and Pre Trigger:
        
        if (getParamString(m_pConfig, "ENABLE_INPUT", c, "NO") == "YES") {
            chanMask |= 1 << c;
            status = CAEN_DGTZ_SetChannelDCOffset(
                m_nHandle, c, getIntParam(m_pConfig, "DC_OFFSET", c, 32767)
            );
            checkStatus(status, "Setting channel DCOffset");
            
            status = CAEN_DGTZ_SetDPPPreTriggerSize(
                m_nHandle, c, getIntParam(m_pConfig, "PRE_TRIGGER", c, 0)
            );
            checkStatus(status, "Setting channel pre trigger");
            
            // Get the parameters for setDPPParameters
            

            params.thr[c] = getIntParam(m_pConfig, "TRG_THRESHOLD", c, 50);
            
            if (getParamString(m_pConfig, "CHANNEL_TRIGGER", c, "DISABLED") ==
                "DISABLED"
            ) {
                params.selft[c] = 0;
            } else {
                params.selft[c] = 1;
                params.trgc[c]  = CAEN_DGTZ_DPP_TriggerConfig_Threshold;
            }
            params.lgate[c] = getIntParam(m_pConfig, "PSD_LONG_GATE", c, 60);
            params.sgate[c] = getIntParam(m_pConfig, "PSD_SHORT_GATE", c, 20);
            params.pgate[c] = getIntParam(m_pConfig, "PSD_PRE_GATE", c, 16);
            params.nsbl[c]  = getIntParam(m_pConfig, "PSD_BL_SAMPLES", c, 3);
            params.csens[c] = getIntParam(m_pConfig, "PSD_SEL_CHARGE_SENSE", c, 0);
            params.tvaw[c]  = getIntParam(m_pConfig, "TRIGGER_VALIDATION_WINDOW", c, 50);
            
            // Set the CFD parameters; bit fields in the CAEN_DGTZ_SAM_TRIGGER_REG.
            
            uint32_t trigReg;
            status = CAEN_DGTZ_ReadRegister(
                m_nHandle, CAEN_DGTZ_SAM_TRIGGER_REG_ADD + (c << 8), &trigReg
            );
            checkStatus(status, "Reading SAM_TRIGGER_REG - to fold in CFD parameters");
            
            trigReg = setBits(
                trigReg, 0, 7, getIntParam(m_pConfig, "CFD_DELAY", c, 40)/2
            );
            trigReg = setBits(
                trigReg, 8, 9, getIntParam(m_pConfig, "CFD_ATTENUATION", c, 0)
            );
            int cfdInterp = getIntParam(m_pConfig, "CFD_INTERPOLATE", c, 0);
            trigReg = setBits(trigReg, 10, 11, cfdInterp);
            m_nCFDDt = (1 + cfdInterp*2) * 1000 *2;      // Dt in ps.
            status = CAEN_DGTZ_WriteRegister(
                m_nHandle, CAEN_DGTZ_SAM_TRIGGER_REG_ADD + (c << 8), trigReg
            );
            checkStatus(status, "Writing SM_TRIGGER_REG_ADD with CFD Parameters");
            

            // LED or CFD?  DIS
            
            uint32_t discMode;
            status = CAEN_DGTZ_ReadRegister(
                m_nHandle, 0x1080 + (c << 8), &discMode
            );
            int usecfd;
            checkStatus(status, "Reading the 0x1080 - to fold in disc mode");
            if (getParamString(m_pConfig, "DISC_MODE", -1, "LED") == "LED") {
                usecfd = 0;
            } else {
                usecfd = 1;
            }
            discMode = setBits(discMode, 6,6, usecfd);
            status = CAEN_DGTZ_WriteRegister(
                m_nHandle, 0x1080 + (c << 8), discMode
            );
            checkStatus(status, "Writing 0x1080 with discriminator mode set");
            
            // Enable interpolation:
            
            uint32_t interpReg;
            status = CAEN_DGTZ_ReadRegister(
                m_nHandle, 0x1084 + (c << 8), &interpReg
            );
            checkStatus(status, "Reading reg 0x1084 to fold request to supply interp");
            interpReg = setBits(interpReg, 8,9, 1);
            interpReg = setBits(interpReg, 10, 11, 1);
            status = CAEN_DGTZ_WriteRegister(
                m_nHandle, 0x1084  + (c << 8), interpReg
            );
            checkStatus(status, "Setting CFD Interp register");
            
            

            // Gain:
            
            std::string dynRange = getParamString(
                m_pConfig, "DYNAMIC_RANGE", c, ".5"
            );
            uint32_t rangeValue = (dynRange == ".5") ? 1 : 0;
            status = CAEN_DGTZ_WriteRegister(
                m_nHandle, 0x1028  + (c << 8), rangeValue
            );
            
            // Resolution and sampling frequency:
            
            std::string resolution = getParamString(
                m_pConfig, "RESOLUTION", c, "14"
            );
            uint32_t resValue = 0;
            if (resolution == "14") {
                resValue = 0;
            } else if (resolution == "13") {
                resValue = 1;
            } else if (resolution == "12") {
                resValue = 2;
            } else {
                resValue = 3;             // 10 bits.
            }
            
            std::string frequency = getParamString(
                m_pConfig, "FREQUENCY", c, "500"
            );
            uint32_t freqValue = 0;
            if (frequency == "250") freqValue = 1;
            
            uint32_t freqResReg;
            status = CAEN_DGTZ_ReadRegister(
                m_nHandle, 0x1080 + (c << 8), &freqResReg
            );
            checkStatus(status, "Reading 1x80 to set freq/res");
            freqResReg = setBits(freqResReg, 12, 13, resValue);
            freqResReg = setBits(freqResReg, 14, 14, freqValue);
            
            status = CAEN_DGTZ_WriteRegister(
                m_nHandle, 0x1080 + (c << 8), freqResReg
            );
            checkStatus(status, "Writing 1n80 to set freq/res.");
            
            // pulse Polarity. POLARITY
            
            
            status = CAEN_DGTZ_SetChannelPulsePolarity(
                m_nHandle, c, getChannelPulsePolarity(m_pConfig, c)
            );
            
            // PSD cut:  PSD_CUT - DISABLED, GAMMA, or NEUTRON
            //           PSD_CUT_LEVEL - float.
            
            std::string psdCutType = getParamString(
                m_pConfig, "PSD_CUT", c, "DISABLED"
            );
            if (psdCutType != "DISABLED") {
                uint32_t dppCutLevel =
                    1024*(getDoubleParam(m_pConfig, "PSD_CUT_LEVEL", c, 0.5));
                status = CAEN_DGTZ_WriteRegister(
                    m_nHandle, 0x1078 + (c << 8), dppCutLevel
                );
                checkStatus(status, "Setting PSD cut level");
                uint32_t dppControl;
                status = CAEN_DGTZ_ReadRegister(
                    m_nHandle, 0x1080 + (c << 8), &dppControl
                );
                checkStatus(status, "Reading 0x1n80 prior to settting psd cut type");
                uint32_t cutType = (psdCutType == "GAMMA") ? 1 : 2;
                dppControl = setBits(dppControl, 27,28, cutType);
                status = CAEN_DGTZ_WriteRegister(
                    m_nHandle, 0x1080 + (c << 8), dppControl
                );
                checkStatus(status, "Writing PSD Cut type");
            }
            
        }
        
    }
    
    status = CAEN_DGTZ_SetDPPParameters(m_nHandle, chanMask, &params);
    checkStatus(status, "Writing DPP parameters");
    
    std::cerr << "Channels enabled: " << std::hex << chanMask << std::dec << std::endl;
    
    status = CAEN_DGTZ_SetChannelEnableMask(m_nHandle, chanMask);
    checkStatus(status, "Setting channel enable mask");
    
    // Enable channels and set coincidencde masks(?)
    
    status = CAEN_DGTZ_WriteRegister(m_nHandle, 0x8180, 0x103);
    checkStatus(status, "Writing register 0x8180 - coinc masks?");
    
    status = CAEN_DGTZ_WriteRegister(m_nHandle, 0x8184, 0x103);
    checkStatus(status, "Writing 0x8184 - coinc masks??");
    
    // Enable data taking:
    // checkStatus(CAEN_DGTZ_ClearData(m_nHandle), "Clearing the data");
    status = CAEN_DGTZ_SWStartAcquisition(m_nHandle);
    checkStatus(status, "Starting acquisition");
    std::cout << "DPP mode acquisition started\n";
}

/**
 * finalize
 *    Stop taking data
 */
void
CDPPReader::finalize()
{
    CAEN_DGTZ_SWStopAcquisition(m_nHandle);
}
/**
 * readEvents
 *    Read at most m_nTriggers events.
 *    If a read fails or produces no events, then we return the number of events
 *    we have so far.
 * @return int - Number of events actually read.
 */
int
CDPPReader::ReadEvents()
{
    // Allocate the buffers
    
    CAEN_DGTZ_ErrorCode status;
    uint32_t rawBufferSize;
    char*    rawBuffer(0);
    CAEN_DGTZ_DPP_PSD_Event_t*  dppBuffer[CAEN_DGTZ_MAX_CHANNEL];
    uint32_t                    dppBufferSize;
    
    status = CAEN_DGTZ_MallocReadoutBuffer(m_nHandle, &rawBuffer, &rawBufferSize);
    checkStatus(status, "Allocating raw readout buffer");
    
    status = CAEN_DGTZ_MallocDPPEvents(m_nHandle, (void**)dppBuffer, &dppBufferSize);
    checkStatus(status, "Allocating the DPPEvents array");
    
    CAEN_DGTZ_DPP_PSD_Waveforms_t* pWaveforms(0);
    uint32_t waveformSize(0);
    status = CAEN_DGTZ_MallocDPPWaveforms(m_nHandle, (void**)&pWaveforms, &waveformSize);
    
    // We don't read or care about waveforms so don't malloc them.
    
    // Read the data:

    
    int nRead = 0;
    while (nRead < m_nTriggers) {
        uint32_t rawReadSize;
        uint32_t nDppEvents[CAEN_DGTZ_MAX_CHANNEL];
        
        // Do the raw read:
        
        CAEN_DGTZ_ErrorCode readStatus = CAEN_DGTZ_ReadData(
            m_nHandle, CAEN_DGTZ_SLAVE_TERMINATED_READOUT_MBLT,
            rawBuffer, &rawReadSize
        );
        if (readStatus != CAEN_DGTZ_Success) {
            break;
        }
        // Decode the DPP data:
        
        CAEN_DGTZ_ErrorCode decodeStatus = CAEN_DGTZ_GetDPPEvents(
            m_nHandle, rawBuffer, rawReadSize, (void**)dppBuffer, nDppEvents
        );
        if (decodeStatus != CAEN_DGTZ_Success) {
            break;
        }
        nRead += writeEvents(dppBuffer, nDppEvents, pWaveforms);
    }
    
    // Free the buffers
    
    CAEN_DGTZ_FreeDPPWaveforms(m_nHandle, pWaveforms);
    CAEN_DGTZ_FreeDPPEvents(m_nHandle, (void**)dppBuffer);
    CAEN_DGTZ_FreeReadoutBuffer(&rawBuffer);

    // Return the number we actually read:
    
    return nRead;
}

/*------------------------------------------------------------------------------
 * Utilities:
 */



/**
 * writeEvents
 *    Write the events to m_nFd  Note that since there are
 *    very few assurances about struct padding (C++/C even worse) we're writing
 *    an item at a time and computing the exact size based on that.
 *
 * @param dppBuffer - Pointers to the DPP data.
 * @param sizes     - Array of sizes.
 * @param pWaveform - Pointer to allocated waveforms struct.
 * @return unsigned - Number of triggered items actually written.'
 *
 * @note data are written in order of their coarse timestamp.
 * @note the fine timestamp is computed.
 */
unsigned
CDPPReader::writeEvents(
    CAEN_DGTZ_DPP_PSD_Event_t**  dppBuffer, uint32_t* nDppEvents,
    CAEN_DGTZ_DPP_PSD_Waveforms_t* pWaveform
)
{
    // we need to keep track of where we are in each of the channels.  There
    // are CAEN_DGTZ_MAX_CHANNEL ints.
    
    int indices[CAEN_DGTZ_MAX_CHANNEL];
    memset(indices, 0, CAEN_DGTZ_MAX_CHANNEL*sizeof(int));
    int nWritten = 0;
    
    while(1) {
        int ch = findEarliestChannel(dppBuffer, nDppEvents, indices);
        if (ch == -1) {
            return nWritten;              // Wrote all.
        }
        
        CAEN_DGTZ_DPP_PSD_Event_t* pEvent = &(dppBuffer[ch][indices[ch]]);
        indices[ch]++;
        CAEN_DGTZ_DecodeDPPWaveforms(m_nHandle, pEvent, pWaveform);
        
        writeEvent(ch, pEvent, pWaveform);
        nWritten++;
    }
}

/**
 * writeEvent
 *   Write a single event to file.
 *   This code writes the following to the file:
 * \verbatim
 *    uint32_t   size   - size of the event.
 *    uint32_t   type   - Type of event (1 for DPP event).
 *    uint32_t   chan   - Channel number.
 *    uint32_t   timeTag- The coarse time tag.
 *    uint16_t   interpTime - interpolated time (integer ps).
 *    int16_t   qShort
 *    int16_t   qLong
 *    int16_t   baseline
 *    int16_t   Pile up rejection.
 *    uint32_t  Number of waveform samples (could be zero)
 *    uint16_t[] wave form samples
 *  \endverbatim
 *
 *  @param ch    - The channel.
 *  @param pEvent - Pointer to a single event.
 *  @param pWaveform - Decoded waveform.
 *  
 */
void
CDPPReader::writeEvent(
    uint32_t ch, CAEN_DGTZ_DPP_PSD_Event_t* pEvent,
    CAEN_DGTZ_DPP_PSD_Waveforms_t* pWaveform
)
{
    // Compute the interpolated time
    // from zp, zn in the extra word;
    
    uint16_t zn = pEvent->Extras & 0xffff;
    uint16_t zp = (pEvent->Extras >> 16) & 0xffff;
    uint16_t zc = 0;
    if ((zn < 8192) && (zp >= 8192)) { // only compute if there was a zc.
        zc = m_nCFDDt * (8192 - zn) / (zp - zn);
    }
    
    
    // Compute the size of the event.
    
    uint32_t eventSize = sizeof(uint32_t)*5 + sizeof(uint16_t)*(5 + pWaveform->Ns);
    uint32_t eventType = 1;
    
    // Buffer all the data so that we can do a single write.
    
    uint8_t event[eventSize];
    uint8_t* p  = event;
    
    memcpy(p, &eventSize, sizeof(uint32_t));             p += sizeof(uint32_t);
    memcpy(p, &eventType, sizeof(uint32_t));             p += sizeof(uint32_t);
    memcpy(p, &ch,        sizeof(uint32_t));             p += sizeof(uint32_t);
    memcpy(p, &(pEvent->TimeTag), sizeof(uint32_t));     p += sizeof(uint32_t);
    memcpy(p, &zc, sizeof(uint16_t));                    p += sizeof(uint16_t);
    memcpy(p, &(pEvent->ChargeShort), sizeof(uint16_t)); p += sizeof(uint16_t);
    memcpy(p, &(pEvent->ChargeLong), sizeof(uint16_t));  p += sizeof(uint16_t);
    memcpy(p, &(pEvent->Baseline), sizeof(int16_t));     p += sizeof(int16_t);
    memcpy(p, &(pEvent->Pur), sizeof(int16_t));          p += sizeof(int16_t);
    memcpy(p, &(pWaveform->Ns), sizeof(uint32_t));       p += sizeof(uint32_t);
    memcpy(p, pWaveform->Trace1, pWaveform->Ns*sizeof(uint16_t));
    
    // Write the data:
    
    int n = write(m_nFd, event, eventSize);
    if (n == -1) {
        throw std::runtime_error("Unable to write a DPP event!!");
     
    }
    
    
}