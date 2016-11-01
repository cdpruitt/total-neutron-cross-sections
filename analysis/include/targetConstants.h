#ifndef TARGET_CONSTANTS_H
#define TARGET_CONSTANTS_H

const int NUMBER_OF_TARGETS = 6; // includes the 'blank' position

// physical target data, listed in order given in 'targetNames':

/*const double targetLength[NUMBER_OF_TARGETS] =    {0,     1.366,  2.737,  1.369,  1.373,  1.371};  // in cm
const double targetDiameter[NUMBER_OF_TARGETS] =  {0.827, 0.826,  0.827,  0.825,  0.827,  0.827};  // in cm
const double targetMass[NUMBER_OF_TARGETS] =      {0,     1.2363, 2.4680, 4.9749, 5.3286, 5.5505}; // in grams
const double targetMolMass[NUMBER_OF_TARGETS] =   {0,     12.01,  12.01,  112,    118.7,  124};    // in g/mol

*/const std::vector<std::string> targetNames = {"blank", "shortCarbon", "longCarbon", "Sn112", "NatSn", "Sn124"}; 

#endif
