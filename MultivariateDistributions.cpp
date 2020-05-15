//
//  MultivariateDistributions.cpp
//  PhaseType2.2
//
//  Created by Jorge Yslas on 05/12/19.
//  Copyright © 2019 Jorge Yslas. All rights reserved.
//

#include "MultivariateDistributions.h"


BivariatePareto::BivariatePareto(double alpha, double theta1, double theta2 ) : m_alpha(alpha), m_theta1(theta1), m_theta2(theta2) {
}

double BivariatePareto::jointDensity (double x1, double x2) {
    return ( m_alpha * (m_alpha + 1) / (m_theta1 * m_theta2) * pow(x1 / m_theta1 + x2 / m_theta2 + 1, -m_alpha - 2));
}

double BivariatePareto::densityMarginal1(double x) {
    return ((m_alpha / m_theta1) * pow((x + m_theta1) / m_theta1, -m_alpha - 1));
}

double BivariatePareto::densityMarginal2(double x) {
    return ((m_alpha / m_theta2) * pow((x + m_theta2) / m_theta2, -m_alpha - 1));
}

double BivariatePareto::densitySum(double x) {
    if (m_theta1 == m_theta2) {
        return ((x * m_alpha * (m_alpha + 1.0) / (m_theta1 * m_theta1)) * pow(x / m_theta1 + 1.0, -m_alpha - 2.0));
    }
    else {
        return ((m_alpha / (m_theta1 - m_theta2)) * pow(x / m_theta1 + 1.0, -m_alpha - 1.0) - (m_alpha / (m_theta1 - m_theta2)) * pow(x / m_theta2 + 1.0, -m_alpha - 1.0));
    }
}



BivariateNormal::BivariateNormal(double mu1, double mu2, double sigma1, double sigma2, double rho ) : m_mu1(mu1), m_mu2(mu2), m_sigma1(sigma1), m_sigma2(sigma2), m_rho(rho) {
}
double BivariateNormal::jointDensity (double x1, double x2) {
    return 1.0 / (2 * std::atan(1) * 4 * m_sigma1 * m_sigma2 * sqrt(1 - m_rho * m_rho)) * exp(-1.0 / (2 * (1 - m_rho * m_rho)) * (((x1 - m_mu1) * (x1 - m_mu1) / (m_sigma1 * m_sigma1)) + ((x2 - m_mu2) * (x2 - m_mu2) / (m_sigma2 * m_sigma2)) - 2 * m_rho * (x1 - m_mu1) * (x2 - m_mu2) / (m_sigma1 * m_sigma2)));
}
double BivariateNormal::densityMarginal1(double x) {
    return (1 / (sqrt(2 * std::atan(1) * 4) * m_sigma1 * exp((x - m_mu1) * (x - m_mu1) / (2 * m_sigma1 * m_sigma1))));
}
double BivariateNormal::densityMarginal2(double x) {
    return (1 / (sqrt(2 * std::atan(1) * 4) * m_sigma2 * exp((x - m_mu2) * (x - m_mu2) / (2 * m_sigma2 * m_sigma2))));
}
double BivariateNormal::densitySum(double x) {
    double mu{m_mu1 + m_mu2};
    double sigma{sqrt(m_sigma1 * m_sigma1 + m_sigma2 * m_sigma2 + 2 * m_rho * m_sigma1 * m_sigma2)};
    return (1 / (sqrt(2 * std::atan(1) * 4) * sigma * exp((x - mu) * (x - mu) / (2 * sigma * sigma))));
}




MarshallOlkinBivariateExponential::MarshallOlkinBivariateExponential(double lambda1, double lambda2, double lambda12 ) : m_lambda1(lambda1), m_lambda2(lambda2), m_lambda12(lambda12) {
}
double MarshallOlkinBivariateExponential::jointDensity (double x1, double x2) {
    if (x2 < x1) {
        return m_lambda2 * (m_lambda1 + m_lambda12) * exp(-(m_lambda1 + m_lambda12) * x1 - m_lambda2 * x2);
    } else if (x2 > x1) {
        return m_lambda1 * (m_lambda2 + m_lambda12) * exp(-(m_lambda2 + m_lambda12) * x2 - m_lambda1 * x1);
    } else {
        return m_lambda12 * exp(-m_lambda12 * x1);
    }
    return 0;
}
double MarshallOlkinBivariateExponential::densityMarginal1(double x) {
    return (m_lambda1 + m_lambda12) * exp(-(m_lambda1 + m_lambda12) * x);
}
double MarshallOlkinBivariateExponential::densityMarginal2(double x) {
    return (m_lambda2 + m_lambda12) * exp(-(m_lambda2 + m_lambda12) * x);
}
double MarshallOlkinBivariateExponential::densitySum(double x) {
    Numeric_lib::Matrix<double,2> pi(1,3);
    pi(0,0) = 1;
    Numeric_lib::Matrix<double,2> T(3,3);
    T(0,0) = -0.5 * (m_lambda1 + m_lambda2 + m_lambda12);
    T(0,1) = 0.5 * m_lambda2;
    T(0,2) = 0.5 * m_lambda1;
    T(1,1) = -(m_lambda1  + m_lambda12);
    T(2,2) = -(m_lambda2  + m_lambda12);
    Numeric_lib::Matrix<double,2> e(3,1);
    e += 1.0;
    Numeric_lib::Matrix<double,2> t(3,1);
    matrixProduct(T * (-1.0), e, t);
    
    if (x == 0) {
        return (1.0 - matrixProduct(pi, e)(0,0));
    }
    return (matrixProduct(pi, matrixProduct(matrixExponential(T * x), t))(0,0));
}




BivariateWeibull::BivariateWeibull(double lambda1, double lambda2, double beta, double gamma ) : m_lambda1(lambda1), m_lambda2(lambda2), m_beta(beta), m_gamma(gamma) {
}
double BivariateWeibull::jointDensity (double x1, double x2) {
   return m_gamma * m_beta * m_beta * m_lambda1 * m_lambda2 * pow(x1 * x2, m_beta - 1) * pow(m_lambda1 * pow(x1, m_beta) + m_lambda2 * pow(x2, m_beta), m_gamma - 2) * exp(- pow(m_lambda1 * pow(x1, m_beta) + m_lambda2 * pow(x2, m_beta), m_gamma)) * ( m_gamma * pow(m_lambda1 * pow(x1, m_beta) + m_lambda2 * pow(x2, m_beta), m_gamma) - m_gamma + 1);
}
double BivariateWeibull::densityMarginal1(double x) {
    return m_gamma * m_beta * pow(m_lambda1, m_gamma) * pow(x, m_gamma * m_beta - 1) * exp(-pow(m_lambda1 * pow(x, m_beta), m_gamma));
}
double BivariateWeibull::densityMarginal2(double x) {
    return m_gamma * m_beta * pow(m_lambda2, m_gamma) * pow(x, m_gamma * m_beta - 1) * exp(-pow(m_lambda2 * pow(x, m_beta), m_gamma));
}
double BivariateWeibull::densitySum(double x) {
    return 0;
}

