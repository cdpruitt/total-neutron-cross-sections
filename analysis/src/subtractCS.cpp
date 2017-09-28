/* How to use:
 *
 * ./subtractCS [raw cross section filename] [raw cross section graph name]
 *              [subtrahend filename] [subtrahend graph name]
 *              [subtrahend factor] [result divisor] [output graph name]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"
#include "../include/experimentalConfig.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>

using namespace std;

ExperimentalConfig experimentalConfig;

int main(int, char* argv[])
{
    subtractCS(argv[1], argv[2], argv[3], argv[4], stod(argv[5]), stod(argv[6]), argv[7]);
    return 0;
}
