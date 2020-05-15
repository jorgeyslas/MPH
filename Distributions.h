// Densities and CDF's of different discrete and continuous distributions


#ifndef Distributions_h
#define Distributions_h

#include <vector>
#include <cmath>

//Continuos densities
double dUniform(double x, const std::vector<double> & parameters);
double dNormal(double x, const std::vector<double> & parameters);
double dLogNormal(double x, const std::vector<double> & parameters);
double dWeibull(double x,const std::vector<double>  & parameters);
double dInverseGaussian(double x, const std::vector<double> & parameters);
double dParetoTypeI(double x, const std::vector<double> & parameters);
double dParetoTypeII(double x, const std::vector<double> & parameters);
double dExponential(double x, const std::vector<double> & parameters);
double dGumbel(double x, const std::vector<double> & parameters);
double dGEVD(double x, const std::vector<double> & parameters);
double dGompertz(double x, const std::vector<double> & parameters);

//Continuous distributions
double cdfStandardNormal(double x);
double cdfNormal(double x, const std::vector<double> & parameters);
double cdfLogNormal(double x, const std::vector<double> & parameters);
double cdfWeibull(double x, const std::vector<double> & parameters);
double cdfParetoTypeI(double x, const std::vector<double> & parameters);
//  Pareto Type I on 1
double cdfPareto(double x, const std::vector<double> & parameters);
double cdfParetoTypeII(double x, const std::vector<double> & parameters);
double cdfExponential(double x, const std::vector<double> & parameters);
double cdfGumbel(double x, const std::vector<double> & parameters);
double cdfGEVD(double x, const std::vector<double> & parameters);
double cdfGompertz(double x, const std::vector<double> & parameters);

//Tails of continuous distributions
double tailNormal(double x, const std::vector<double> & parameters);
double tailLogNormal(double x, const std::vector<double> & parameters);
double tailWeibull(double x, const std::vector<double> & parameters);
double tailParetoTypeI(double x, const std::vector<double> & parameters);
double tailPareto(double x, const std::vector<double> & parameters);
double tailParetoTypeII(double x, const std::vector<double> & parameters);
double tailExponential(double x, const std::vector<double> & parameters);
double tailGumbel(double x, const std::vector<double> & parameters);
double tailGEVD(double x, const std::vector<double> & parameters);
double tailGompertz(double x, const std::vector<double> & parameters);


//Discrete densities
double dPoisson( double k, double lambda);

//Discrete distributions
double cdfPoisson( double k, double lambda);
#endif /* Distributions_h */
