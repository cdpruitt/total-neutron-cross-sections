#ifndef TARGET_H
#define TARGET_H

#include <string>

#include "TH1I.h"
#include "TGraphErrors.h"

#include "crossSection.h"

const std::string TARGET_DATA_FILE_PATH = "targetData/";
const std::string TARGET_DATA_FILE_EXTENSION = ".txt";

class Target
{
    public:
        // creates a new target, based on a text file in ./targetData/
        Target(std::string targetDataLocation);
        std::string getName();
        double getLength();
        double getDiameter();
        double getMass();
        double getMolMass();

    private:
        // physical target parameters
        std::string name;
        double length;
        double diameter;
        double mass;
        double molMass;
};

#endif
