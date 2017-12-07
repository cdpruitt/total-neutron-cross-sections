#include "../include/crossSection.h"
#include "../include/dataSet.h"
#include "../include/dataPoint.h"

CrossSection calculateRelativeDiff(CrossSection a, CrossSection b)
{
    CrossSection relative; 

    DataSet aData = a.getDataSet();
    DataSet bData = b.getDataSet();

    if(aData.getNumberOfPoints()!=bData.getNumberOfPoints())
    {
        cerr << "Error: can't calculate relative cross section from "
             << "data sets of different sizes. Returning empty cross section."
             << endl;
        return relative;
    }

    DataSet relativeDataSet;

    // for each point, calculate the relative cross section, including error
    for(int i=0; i<aData.getNumberOfPoints(); i++)
    {
        DataPoint aPoint = aData.getPoint(i);
        DataPoint bPoint = bData.getPoint(i);

        double aXValue = aPoint.getXValue();
        double bXValue = bPoint.getXValue();
        
        if(aXValue != bXValue)
        {
            cerr << "Error: can't calculate relative cross section from "
                 << "data points with different x values. Returning cross section."
                 << endl;
            return relative;
        }

        double aYValue = aPoint.getYValue();
        double bYValue = bPoint.getYValue();

        double yValue = (aYValue-bYValue)/(aYValue+bYValue);

        // calculate cross section error
        double aError = getPartialError(aPoint, bPoint, a.getArealDensity());
        double bError = getPartialError(bPoint, aPoint, b.getArealDensity());
        double totalError = pow(pow(aError,2)+pow(bError,2),0.5);

        relativeDataSet.addPoint(
                DataPoint(aXValue, aPoint.getXError(), 100*yValue, 100*totalError, /* convert to % */
                          aPoint.getBlankMonitorCounts()+bPoint.getBlankMonitorCounts(),
                          aPoint.getTargetMonitorCounts()+bPoint.getTargetMonitorCounts(),
                          aPoint.getBlankDetCounts()+bPoint.getBlankDetCounts(),
                          aPoint.getTargetDetCounts()+bPoint.getTargetDetCounts()));
    }

    relative.addDataSet(relativeDataSet);
    return relative;
}
