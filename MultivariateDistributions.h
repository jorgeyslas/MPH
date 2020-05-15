//
//  MultivariateDistributions.h
//  PhaseType2.2
//
//  Created by Jorge Yslas on 05/12/19.
//  Copyright © 2019 Jorge Yslas. All rights reserved.
//

#ifndef MultivariateDistributions_h
#define MultivariateDistributions_h

#include <cmath>
#include "Matrix.h"
#include "MatrixOperations.h"

class MultivariateDistribution {
public:
    virtual double jointDensity (double x1, double x2) = 0;
    virtual double densityMarginal1(double x) = 0;
    virtual double densityMarginal2(double x) = 0;
    virtual double densitySum(double x) = 0;
    virtual ~MultivariateDistribution(){}
};


class BivariatePareto: public MultivariateDistribution {
    double m_alpha;
    double m_theta1;
    double m_theta2;
public:
    BivariatePareto(double alpha, double theta1, double theta2);
    virtual double jointDensity (double x1, double x2);
    virtual double densityMarginal1(double x);
    virtual double densityMarginal2(double x);
    virtual double densitySum(double x);
    virtual ~BivariatePareto(){}
};


class BivariateNormal: public MultivariateDistribution {
    double m_mu1;
    double m_mu2;
    double m_sigma1;
    double m_sigma2;
    double m_rho;
public:
    BivariateNormal(double mu1, double mu2, double sigma1, double sigma2, double rho);
    virtual double jointDensity (double x1, double x2);
    virtual double densityMarginal1(double x);
    virtual double densityMarginal2(double x);
    virtual double densitySum(double x);
    virtual ~BivariateNormal(){}
};


class MarshallOlkinBivariateExponential: public MultivariateDistribution {
    double m_lambda1;
    double m_lambda2;
    double m_lambda12;
public:
    MarshallOlkinBivariateExponential(double lambda1, double lambda2, double lambda12);
    virtual double jointDensity (double x1, double x2);
    virtual double densityMarginal1(double x);
    virtual double densityMarginal2(double x);
    virtual double densitySum(double x);
    virtual ~MarshallOlkinBivariateExponential(){}
};



class BivariateWeibull: public MultivariateDistribution {
    double m_lambda1;
    double m_lambda2;
    double m_beta;
    double m_gamma;
public:
    BivariateWeibull(double lambda1, double lambda2, double beta, double gamma);
    virtual double jointDensity (double x1, double x2);
    virtual double densityMarginal1(double x);
    virtual double densityMarginal2(double x);
    virtual double densitySum(double x);
    virtual ~BivariateWeibull(){}
};

#endif /* MultivariateDistributions_h */
