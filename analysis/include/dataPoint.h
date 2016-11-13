#ifndef DATA_POINT_H
#define DATA_POINT_H

class DataPoint
{
    public:
        DataPoint();
        DataPoint(double xValue, double xError,
                  double yValue, double yError);
        double getXValue() const;
        double getXError() const;
        double getYValue() const;
        double getYError() const;

        void setXValue(double xv);
        void setXError(double xe);
        void setYValue(double yv);
        void setYError(double ye);

        DataPoint mergePoints(DataPoint point2);

        friend DataPoint operator-(const DataPoint& minuend, const DataPoint& subtrahend);
        friend DataPoint operator/(const DataPoint& dividend, const DataPoint& divisor);

        friend DataPoint operator*(const DataPoint& multiplicand, const double multiplier);
        friend DataPoint operator/(const DataPoint& dividend, const double divisor);

    private:
        double xValue;
        double xError;
        double yValue;
        double yError;
};

#endif
