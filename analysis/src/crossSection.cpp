#include <iostream>

#include "../include/crossSection.h"
#include "../include/DataPoint.h"

using namespace std;

CrossSection::CrossSection()
{
}

void CrossSection::addDataPoint(DataPoint dataPoint)
{
    data.push_back(dataPoint);
}

DataPoint CrossSection::getDataPoint(int i)
{
    if((size_t)i>data.size())
    {
        cout <<
        "Error: tried to retrieve a non-existent cross section data point" <<
        endl;

        exit(1);
    }

    return data[i];
}

int CrossSection::getNumberOfPoints()
{
    return data.size();
}

vector<double> CrossSection::getEnergyValues()
{
    vector<double> energyValues;
    for(DataPoint d : data)
    {
        energyValues.push_back(d.getXValue());
    }
    return energyValues;
}

vector<double> CrossSection::getEnergyErrors()
{
    vector<double> energyErrors;
    for(DataPoint d : data)
    {
        energyErrors.push_back(d.getXError());
    }
    return energyErrors;
}

vector<double> CrossSection::getCrossSectionValues()
{
    vector<double> crossSectionValues;
    for(DataPoint d : data)
    {
        crossSectionValues.push_back(d.getYValue());
    }
    return crossSectionValues;
}

vector<double> CrossSection::getCrossSectionErrors()
{
    vector<double> crossSectionErrors;
    for(DataPoint d : data)
    {
        crossSectionErrors.push_back(d.getYError());
    }
    return crossSectionErrors;
}
