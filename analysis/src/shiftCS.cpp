/* How to use:
 *
 * ./shiftCS [raw cross section filename] [raw cross section graph name]
 *              [subtrahend filename] [subtrahend graph name]
 *              [subtrahend factor] [result divisor] [output graph name]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"
#include "../include/config.h"
#include "../include/CSUtilities.h"

#include "TFile.h"
#include "TGraphAsymmErrors.h"

#include <string>
#include <iostream>

using namespace std;

Config config;

int main(int, char* argv[])
{
    shiftCS(argv[1], argv[2], stod(argv[3]), argv[4], argv[5]);
    return 0;
}
