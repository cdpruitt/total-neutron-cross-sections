#ifndef DATA_POINT_H
#define DATA_POINT_H

class DataPoint
{
    public:
        DataPoint(double xValue, double xError,
                  double yValue, double yError);
        double getXValue();
        double getXError();
        double getYValue();
        double getYError();

    private:
        double xValue;
        double xError;
        double yValue;
        double yError;
};

#endif
