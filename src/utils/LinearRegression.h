#ifndef LINEAR_REGRESSION_H
#define LINEAR_REGRESSION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>

/** 
  * Linear Regression calculation class adapted from http://david.swaim.com/cpp/linreg.htm
  *
  * This class implements a standard linear regression on
  * experimental data using a least squares fit to a straight
  * line graph.  Calculates coefficients a and b of the equation:
  *
  * y = a + b * x
  *
  * for data points of x and y.  Also calculates the coefficient of
  * determination, the coefficient of correlation, and standard
  * error of estimate.
  *
  * The value n (number of points) must be greater than 2 to
  * calculate the regression.  This is primarily because the
  * standard error has a (N-2) in the denominator.
  *
  * Check haveData() to see if there is enough data in
  * LinearRegression to get values.
  *
  **/
namespace AMP {

class LinearRegression {
  friend std::ostream& operator<<(std::ostream& os, LinearRegression& linReg) {
    if (linReg.haveData()) 
      os<<"f(x) = "<<linReg.getA()<<" + ( "<<linReg.getB()<<" * x )";
    return os;
  }

public:
  LinearRegression(const std::vector<double>& x, const std::vector<double>& y) 
    : nPoints(0), a(0.0), b(0.0), sumX(0.0), sumY(0.0), sumXsquared(0.0), sumYsquared(0.0), sumXY(0.0), coefD(-1.0), coefC(-1.0), stdError(-1.0) {
    unsigned int sizeData = x.size();
    if(sizeData != y.size()) std::cout<<"throw an error\n";
    for(unsigned int i = 0; i < sizeData; ++i) addXY(x[i], y[i]);
  }
	
  void addXY(const double& x, const double& y) {
    ++nPoints; sumX += x; sumY += y; sumXsquared += x*x; sumYsquared += y*y; sumXY += x*y;
    Calculate();
  }

  // Must have at least 3 points to calculate
  // standard error of estimate.  Do we have enough data?
  bool haveData() const { return (nPoints > 2 ? true : false); }
  unsigned int items() const { return nPoints; }
	
  double getA() const { return a; }
  double getB() const { return b; }
	
  double getCoefDeterm() const  { return coefD; }
  double getCoefCorrel() const { return coefC; }
  double getStdErrorEst() const { return stdError; }
  double estimateY(double x) const { return (a + b * x); }
	
protected:
  unsigned int nPoints;     // number of data points input so far
  double a, b;        // coefficients of f(x) = a + b*x
  double sumX, sumY;  // sums of x and y
  double sumXsquared, // sum of x squares
         sumYsquared; // sum y squares
  double sumXY;       // sum of x*y
	
  double coefD,       // coefficient of determination
         coefC,       // coefficient of correlation
         stdError;    // standard error of estimate
	
  void Calculate() {
    if (haveData()) {
      if (std::fabs( double(nPoints) * sumXsquared - sumX * sumX) > DBL_EPSILON) {
        b = ( double(nPoints) * sumXY - sumY * sumX) / ( double(nPoints) * sumXsquared - sumX * sumX);
        a = (sumY - b * sumX) / double(nPoints);

        double sx = b * ( sumXY - sumX * sumY / double(nPoints) );
        double sy2 = sumYsquared - sumY * sumY / double(nPoints);
        double sy = sy2 - sx;

        coefD = sx / sy2;
        coefC = sqrt(coefD);
        stdError = sqrt(sy / double(nPoints - 1));
      }
      else {
        a = b = coefD = coefC = stdError = 0.0;
      }
    }
  }
};

} // end namespace AMP

#endif // LINEAR_REGRESSION_H
