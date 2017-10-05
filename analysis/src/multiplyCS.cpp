/* How to use:
 *
 * ./multiplyCS [raw cross section filename] [raw cross section graph name]
 *              [multiplier] [output graph name]
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
    multiplyCS(argv[1], argv[2], stod(argv[3]), argv[4]);
    return 0;
}
