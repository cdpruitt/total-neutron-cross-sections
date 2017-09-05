/* How to use:
 *
 * ./applyCSCorrectionFactor [CS correction file name] [CS correction graph name]
 *                        [CS to-be-corrected file name] [CS to-be-corrected graph name]
 *                        [output file name] [output graph name]
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
    applyCSCorrectionFactor(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]); 

    return 0;
}
