/* How to use:
 *
 * ./subtractCS [raw cross section filename] [raw cross section graph name]
 *              [subtrahend filename] [subtrahend graph name]
 *              [subtrahend multiplier]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>

using namespace std;

int main(int, char* argv[])
{
    subtractCS(argv[1], argv[2], argv[3], argv[4], stod(argv[5]), stod(argv[6]));
    return 0;
}
