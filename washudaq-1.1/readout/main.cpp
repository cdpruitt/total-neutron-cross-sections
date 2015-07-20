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
 * @file main.cpp
 * @brief Main program unit.
 */


#include <config.h>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <CAENDigitizer.h>

#include "DPPConfig.h"
#include "CDPPReader.h"
#include "CWfReader.h"
#include "keyb.h"

struct Configuration {
  std::string dppConfig;
  std::string wfConfig;
  int  dppTriggers;
  int  wfTriggers;
};


/**
 * Usage 
 *   
 * Output proper program usage and then exits with bad status.
 */
static void Usage()
{
  std::cerr << "Usage:\n";
  std::cerr << "  readout mainconfig outputfile\n";
  std::cerr << "Where: \n";
  std::cerr << "  mainconfig - is the top level configuration file.\n";

  std::exit(EXIT_FAILURE);
}

/**
 * checkConfigAccess 
 *  Ensure that a configuration file is non-empty and also accessible
 *
 * @param which - Desdcribes the config file (part of exception msg).
 * @param file  - Filename path.
 * @throw std::runtime_error - if file is empty or not accessible.
 */
static void
checkConfigAccess(std::string which, std::string file)
{
  std::string msg = which;
  if (file == "") {
    
    msg += " configuration file needs to be defined";
    throw std::runtime_error(msg);
  }
  if (access(file.c_str(), R_OK)) {
    msg += " configuration file is not accessible";
    throw std::runtime_error(msg);
  }
}
/**
 * checkPositive
 *   Require an integer be greater than 0.
 *
 * @param which - describes the integer.
 * @param value - The integer value.
 * @throw std::runtime_error if the tests fail.
 */
static void
checkPositive(std::string which, int value)
{
  if (value <= 0) {
    std::string msg = which;
    msg += " must be strictly greater than 0 and was not";
    throw std::runtime_error(msg);
  }
}

/**
 * unpackConfig
 *   Unpacks a configuration struct from the parsed configuration file:
 *   -   It's an error for a required element to not be there.
 *   -   It's an error to have one of the configuration files missing.
 * 
 * @param dict - The parsed toplevel config file.
 */

static Configuration 
unpackConfig(dictionary* dict)
{
  Configuration result;

  // Get the parameters with good defaults for the ints, but unreasonable names
  // for the child config files.

  result.dppConfig   = getParamString(dict, "dppconfig", -1, "");
  result.wfConfig    = getParamString(dict, "waveformconfig", -1, "");
  result.dppTriggers = getIntParam(dict, "dpptriggers", -1, 10000);
  result.wfTriggers  = getIntParam(dict, "waveformtriggers", -1, 1);

  // Ensure that the configs are nonempty and accessible:

  checkConfigAccess("DPP Config file", result.dppConfig);
  checkConfigAccess("Waveform config file", result.wfConfig);

  // The trigger counts must be positive:

  checkPositive("DPP Trigger count", result.dppTriggers);
  checkPositive("Waveform trigger count", result.wfTriggers);

  return result;
}



/**
 * Main program.
 *  Usage:
 *     readout master-config-file
 **/
int main(int argc, char** argv)
{
  if (argc != 3) {
    Usage();
  }
  // Load the main configuration file:


  dictionary* mainConfig;


  try {
    mainConfig = loadConfig(argv[1]);
  }
  catch(std::exception& e) {
    std::cerr << "Failed to process main configuration file:\n";
    std::cerr << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Configuration cfg;
  try {
    cfg = unpackConfig(mainConfig);
  }
  catch (std::exception& e) {
    std::cerr << "Failed to get what we need from the main config\n";
    std::cerr << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Open the output file:
 
  int fd = open(argv[2], O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fd == -1) {
    perror("Unable to open the output file: ");
    std::exit(EXIT_FAILURE);
  }


  // Load the dpp config file.. we need to see enough info to open the handle.
  // If we can we open the digitizer here:

  dictionary* dppConfig;
  int handle;

  try {
    dppConfig = loadConfig(cfg.dppConfig.c_str());
    ConnectionParams conSpec = getConnectionParams(dppConfig);
    CAEN_DGTZ_ErrorCode stat = CAEN_DGTZ_OpenDigitizer(
        conSpec.linkType, conSpec.linkNum, conSpec.ConetNode, conSpec.VMEBaseAddress, &handle
    );
    if (stat != CAEN_DGTZ_Success) {
      throw std::runtime_error("Failed to open the digitizer specified by the DPP config file");
    }

  }
  catch (std::exception& e ) {
    std::cerr << "Failed to process the dpp configuration file \n";
    std::cerr << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Process the waveform mode config file:

  dictionary* wfConfig;
  try {
    wfConfig = loadConfig(cfg.wfConfig.c_str());
  }
  catch (std::exception& e) {
    std::cerr << "Failed to process the wave form configuration file\n";
    std::cerr << e.what() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  
  CDPPReader dppReader(handle, dppConfig, fd, cfg.dppTriggers);
  CWfReader  wfReader(handle, wfConfig, dppConfig, fd, cfg.wfTriggers);
  
  // Alternate between dppReader and wfReader checking for keyboard hits in the
  // middle:
  
  
  while (1) {
    int nevents = 0;
    dppReader.setup();

    while (nevents < cfg.dppTriggers) {
      nevents += dppReader.ReadEvents();
      if (kbhit() && (getch() == 'x')) {
        goto done;
      }
  
    }

    dppReader.finalize();
    
    nevents = 0;
    wfReader.setup();
    while (nevents < cfg.wfTriggers) {
      nevents += wfReader.ReadEvents();
      if (kbhit() && (getch() == 'x')) {
       goto done;
      }

    }
    wfReader.finalize();            // Restore the changed bits of dppreader.
    
  }
done:
  close(fd);
  std::exit(EXIT_SUCCESS);
  std::cerr << "Exit requested\n";
}
