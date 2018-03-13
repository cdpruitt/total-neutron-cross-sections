#ifndef TARGET_H
#define TARGET_H

#include <string>

#include "TH1I.h"
#include "TGraphAsymmErrors.h"

class Target
{
    public:
        // creates a new target, based on a text file in ./targetData/
        Target();
        Target(std::string targetDataLocation);
        std::string getName() const;
        double getLength() const;
        double getDiameter() const;
        double getMass() const;
        double getMolarMass() const;

        double getDiameterUncertainty() const;
        double getMassUncertainty() const;
        double getMolarMassUncertainty() const;

        void setName(std::string n);
        void setLength(double l);
        void setDiameter(double d);
        void setMass(double m);
        void setMolarMass(double mm);

        void setDiameterUncertainty(double d);
        void setMassUncertainty(double m);
        void setMolarMassUncertainty(double mm);

    private:
        // physical target parameters
        std::string name;
        double length;
        double diameter;
        double mass;
        double molMass;

        double diameterUncertainty;
        double massUncertainty;
        double molMassUncertainty;
};

#endif
