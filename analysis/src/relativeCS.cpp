/* How to use:
 *
 * ./relativeCS [first cross section filename] [first cross section graph name]
 *              [second cross section filename] [second cross section graph name]
 *              [output graph name]
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
    relativeCS(argv[1], argv[2], argv[3], argv[4], argv[5]); 

    return 0;
}
