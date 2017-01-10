#ifndef DATA_POINT_H
#define DATA_POINT_H

class DataPoint
{
    public:
        DataPoint();
        DataPoint(double xValue, double xError,
                  double yValue, double yError);
        DataPoint(double xValue, double xError,
                  double yValue, double yError,
                  long blankMonitorCounts,
                  long targetMonitorCounts,
                  long blankDetCounts,
                  long targetDetCounts);

        double getXValue() const;
        double getXError() const;
        double getYValue() const;
        double getYError() const;

        void setXValue(double xv);
        void setXError(double xe);
        void setYValue(double yv);
        void setYError(double ye);

        long getBlankMonitorCounts() const;
        long getTargetMonitorCounts() const;
        long getBlankDetCounts() const;
        long getTargetDetCounts() const;

        void setBlankMonitorCounts(long bMC);
        void setTargetMonitorCounts(long tMC);
        void setBlankDetCounts(long bDC);
        void setTargetDetCounts(long tDC);

        DataPoint mergePoints(DataPoint point2);

        friend DataPoint operator+(const DataPoint& augend, const DataPoint& addend);
        friend DataPoint operator-(const DataPoint& minuend, const DataPoint& subtrahend);
        friend DataPoint operator/(const DataPoint& dividend, const DataPoint& divisor);

        friend DataPoint operator*(const DataPoint& multiplicand, const double multiplier);
        friend DataPoint operator/(const DataPoint& dividend, const double divisor);

    private:
        // cross section data
        double xValue;
        double xError;
        double yValue;
        double yError;

        // data for propagating error
        double blankMonitorCounts;
        double targetMonitorCounts;
        double blankDetCounts;
        double targetDetCounts;
};

#endif
