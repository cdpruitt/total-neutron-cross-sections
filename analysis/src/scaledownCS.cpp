/* How to use:
 *
 * ./scaledown [CS to-be-corrected file name] [CS to-be-corrected graph name]
 *             [scaledown]
 *             [output file name] [output graph name]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"
#include "../include/config.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>

using namespace std;

Config config;

int main(int, char* argv[])
{
    scaledownCS(argv[1], argv[2], stoi(argv[3]), argv[4], argv[5]);

    return 0;
}
