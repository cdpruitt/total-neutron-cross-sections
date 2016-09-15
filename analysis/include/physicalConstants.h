#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

/******************************************************************************/
/* Physical constants used in the experiment */
/******************************************************************************/

// Facility parameters

const double MACRO_FREQUENCY = 120;   // frequency of beam macropulses, in Hz
const double MACRO_PERIOD = pow(10,9)/MACRO_FREQUENCY;
                                      // macropulse period, in ns
const double MACRO_LENGTH = 650000;   // macropulse duration, in ns
const double MICRO_LENGTH = 1788.814; // micropulse duration, in ns
const double FLIGHT_DISTANCE = 2672;  // detector distance from neutron
                                      // source, in cm

// Digitizer constants

const double SAMPLE_PERIOD = 2;       // digitizer sample rate, in ns
const double MACROPULSE_OFFSET = 842; // timing delay between the digitizer (channel 0)
                                      // and the facility's RF clock (due to cable delay,
                                      // NIM logic, etc.), in ns
const double SCAVENGER_OFFSET = -12.51; // timing delay between channels 0 and 6
                                        // on the digitizer (due to cable delay,
                                        // NIM logic, etc.), in ns

// Nuclear physics constants

const double C = 299792458; // speed of light in m/s
const double NEUTRON_MASS = 939.56536; // in MeV/c^2
const double AVOGADROS_NUMBER = 6.02214*pow(10.,23.); // in atoms/mol

#endif
