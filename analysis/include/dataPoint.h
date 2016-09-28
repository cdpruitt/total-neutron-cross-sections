#ifndef DATA_POINT_H
#define DATA_POINT_H

class DataPoint
{
    public:
        DataPoint(double xValue, double xError,
                  double yValue, double yError);
        double getXValue() const;
        double getXError() const;
        double getYValue() const;
        double getYError() const;

        friend DataPoint operator-(const DataPoint& minuend, const DataPoint& subtrahend);

    private:
        double xValue;
        double xError;
        double yValue;
        double yError;
};

#endif
