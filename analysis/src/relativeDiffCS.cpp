/* How to use:
 *
 * ./relativeDiffCS [first cross section filename] [first cross section graph name]
 *              [second cross section filename] [second cross section graph name]
 *              [output graph name]
 */

#include "../include/dataSet.h"
#include "../include/crossSection.h"
#include "../include/plots.h"
#include "../include/config.h"
#include "../include/CSUtilities.h"

#include "TFile.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>

using namespace std;

Config config;

int main(int, char* argv[])
{
    relativeDiffCS(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); 

    return 0;
}
