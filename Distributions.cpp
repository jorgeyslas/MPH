// Densities and CDF's of different discrete and continuous distributions

#include "Distributions.h"

//Univariate continuous distributions
//  Densities

//      Uniform (a,b) a=parameters[0], b=parameters[1]
double dUniform(double x, const std::vector<double> & parameters) {
    if ((x >= parameters[0]) && (x <= parameters[1])) {
        return (1.0 / (parameters[1] - parameters[0]));
    }
    else {
        return 0.0;
    }
}
//      Normal (mu, sigma) mu=parameters[0], sigma=parameters[1]
double dNormal(double x, const std::vector<double> & parameters) {
    return (1 / (sqrt(2 * std::atan(1) * 4) * parameters[1] * exp((x - parameters[0]) * (x - parameters[0]) / (2 * parameters[1] * parameters[1]))));
}
//      LogNormal (mu, sigma) mu=parameters[0], sigma=parameters[1]
double dLogNormal(double x, const std::vector<double> & parameters) {
    if (x > 0) {
        return (1 / (sqrt(2 * std::atan(1) * 4) * x * parameters[1] * exp((log(x) - parameters[0]) * (log(x) - parameters[0]) / (2 * parameters[1] * parameters[1]))));
    }
    else {
        return 0.0;
    }
    
}
//    Returns a "big" real number.
double bigNum() {
    return 1.0E+30;
}
//      Weibull (tau, beta) tau=parameters[0] beta=parameters[1], tau is the one that determines the tail
double dWeibull(double x, const std::vector<double> & parameters) {
    if (x == 0 && parameters[0] < 1) {
        return (bigNum());
    }
    return (parameters[0] * pow(parameters[1], parameters[0]) * pow(x, parameters[0] - 1) / exp(pow((parameters[1] * x), parameters[0])));
}
//      Inverse Gaussian (lambda, mu) lambda=parameters[0] mu=parameters[1]
double dInverseGaussian(double x, const std::vector<double> & parameters) {
    return (sqrt(parameters[0] / (2 * std::atan(1) * 4)) * pow(x, -1.5) * exp(-0.5 * parameters[0] * (x - parameters[1]) * (x - parameters[1]) / (parameters[1] * parameters[1] * x)));
}
//      Pareto type I (alpha, kappa) alpha determines the tail, only values larger than kappa
double dParetoTypeI(double x, const std::vector<double> & parameters) {
    if (x >= parameters[1]) {
        return (parameters[0] * pow(parameters[1] / x, parameters[0]) / x);
    }
    else {
        return 0.0;
    }
}
//      Pareto Lomax or type II (alpha, kappa) alpha determines the tail, values larger than 0
double dParetoTypeII(double x, const std::vector<double> & parameters) {
    return (parameters[0] * pow(parameters[1] / (x + parameters[1]), parameters[0]) / (x + parameters[1]));
}
//      Exponential (lambda)
double dExponential(double x, const std::vector<double> & parameters) {
    return (parameters[0] * exp(-parameters[0] * x));
}
//      Gumbel (mu, sigma). Suport on (-Inf, Inf)
double dGumbel(double x, const std::vector<double> & parameters) {
    double z{(x - parameters[0]) / parameters[1]};
    return ( exp(-(z + exp(-z))) / parameters[1] );
}
//      GEVD (mu, sigma, xi). Suport on (-Inf, Inf) for xi = 0, (mu - sigma / xi, Inf) for xi > 0 and (-Inf, mu - sigma / xi) for xi < 0
double dGEVD(double x, const std::vector<double> & parameters) {
    if (parameters[2] == 0) {
        return dGumbel(x, parameters);
    }
    else if (parameters[2] > 0) {
        if (x > parameters[0] - parameters[1] / parameters[2]) {
            double z{pow(1 + parameters[2] * (x - parameters[0]) / parameters[1], -1 / parameters[2])};
            return ( pow(z, parameters[2] + 1) * exp(-z) / parameters[1] );
        }
        else {
            return 0;
        }
    }
    else {
        if (x < parameters[0] - parameters[1] / parameters[2]) {
            double z{pow(1 + parameters[2] * (x - parameters[0]) / parameters[1], -1 / parameters[2])};
            return ( pow(z, parameters[2] + 1) * exp(-z) / parameters[1] );
        }
        else {
            return 0;
        }
    }
}
//      Gompertz with scale parameter alpha and shape beta
double dGompertz(double x, const std::vector<double> & parameters) {
    return parameters[0] * exp(parameters[1] * x) * exp(- parameters[0] * (exp(parameters[1] * x) - 1) / parameters[1] );
}


//  Distribution functions
//      Standard normal
double cdfStandardNormal(double x)
{
    double a1{0.398942280444};
    double a2{0.399903438504};
    double a3{5.75885480458};
    double a4{29.8213557808};
    double a5{2.62433121679};
    double a6{48.6959930692};
    double a7{5.92885724438};
    double b0{0.398942280385};
    double b1{3.8052E-08};
    double b2{1.00000615302};
    double b3{3.98064794E-04};
    double b4{1.98615381364};
    double b5{0.151679116635};
    double b6{5.29330324926};
    double b7{4.8385912808};
    double b8{15.1508972451};
    double b9{0.742380924027};
    double b10{30.789933034};
    double b11{3.99019417011};
    double cdf{};
    double q{};
    double y{};

    if (fabs(x) <= 1.28) {
        y = 0.5 * x * x;
        q = 0.5 - fabs(x) * (a1 - a2*y / (y + a3 - a4 / (y + a5 + a6 / (y + a7))));
    }
    else if (fabs(x) <= 12.7) {
        y = 0.5 * x * x;
        q = exp(-y) * b0 / (fabs(x) - b1 + b2 / (fabs(x) + b3 + b4 / (fabs(x) - b5 + b6 / (fabs(x) + b7 - b8 / (fabs(x) + b9 + b10 / (fabs(x) + b11))))));
    }
    else {
        q = 0.0;
    }
    if (x < 0.0) {
        cdf = q;
    }
    else {
        cdf = 1.0 - q;
    }
    return cdf;
}
//      Normal(mu,sigma)
double cdfNormal(double x, const std::vector<double> & parameters) {
    return cdfStandardNormal((x - parameters[0]) / parameters[1]);
}
//      LogNormal
double cdfLogNormal(double x, const std::vector<double> & parameters) {
    return cdfNormal(log(x), parameters);
}
//      Weibull
double cdfWeibull(double x, const std::vector<double> & parameters) {
    return (1.0 - exp(pow(parameters[1] * x, -parameters[0])));
}
//      Pareto Type I
double cdfParetoTypeI(double x, const std::vector<double> & parameters) {
    if (x > parameters[1]) {
        return (1.0 - pow(parameters[1] / x, parameters[0]));
    }
    else {
        return 0.0;
    }
}
//      Pareto Type I on 1
double cdfPareto(double x, const std::vector<double> & parameters) {
    if (x > 1.0) {
        return (1.0 - pow(x, -parameters[0]));
    }
    else {
        return 0.0;
    }
}
//      Pareto Lomax or type II (alpha, kappa) alpha determines the tail, values larger than 0
double cdfParetoTypeII(double x, const std::vector<double> & parameters) {
    return (1.0 - pow(parameters[1] / (x + parameters[1]), parameters[0]));
}
//      Exponential (lambda)
double cdfExponential(double x, const std::vector<double> & parameters) {
    return (1.0 - exp(-parameters[0] * x));
}
//      Gumbel (mu, sigma)
double cdfGumbel(double x, const std::vector<double> & parameters) {
    return (exp(-exp(-(x - parameters[0]) / parameters[1])));
}
//      GEVD (mu, sigma, xi). Suport on (-Inf, Inf) for xi = 0, (mu - sigma / xi, Inf) for xi > 0 and (-Inf, mu - sigma / xi) for xi < 0
double cdfGEVD(double x, const std::vector<double> & parameters) {
    if (parameters[2] == 0) {
           return cdfGumbel(x, parameters);
    }
    else if (parameters[2] > 0) {
        if (x > parameters[0] - parameters[1] / parameters[2]) {
               return (exp(-pow(1 + parameters[2] * (x - parameters[0]) / parameters[1], -1 / parameters[2]) )  );
        }
        else {
            return 0;
        }
    }
    else {
        if (x < parameters[0] - parameters[1] / parameters[2]) {
            return (exp(-pow(1 + parameters[2] * (x - parameters[0]) / parameters[1], -1 / parameters[2]) )  );
        }
        else {
            return 1;
        }
    }
}
//      Gompertz with scale parameter alpha and shape beta
double cdfGompertz(double x, const std::vector<double> & parameters) {
    return 1 - exp(- parameters[0] * (exp(parameters[1] * x) - 1) / parameters[1] );
}


//  Survival functions or tails
//      Normal(mu,sigma)
double tailNormal(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfNormal(x, parameters));
}
//      LogNormal
double tailLogNormal(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfLogNormal(x, parameters));
}
//      Weibull
double tailWeibull(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfWeibull(x, parameters));
}
//      Pareto Type I
double tailParetoTypeI(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfParetoTypeI(x, parameters));
}
//      Pareto Type I on 1
double tailPareto(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfPareto(x, parameters));
}
//      Pareto Lomax or type II (alpha, kappa) alpha determines the tail, values larger than 0
double tailParetoTypeII(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfParetoTypeII(x, parameters));
}
//      Exponential (lambda)
double tailExponential(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfExponential(x, parameters));
}
double tailGumbel(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfGumbel(x, parameters));
}
double tailGEVD(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfGEVD(x, parameters));
}
double tailGompertz(double x, const std::vector<double> & parameters) {
    return (1.0 - cdfGompertz(x, parameters));
}

//Discrete Distributions
//  Densities
//      Poisson
double dPoisson( double k, double lambda) {
  return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}

//  Distribution functions
//      Poisson
double cdfPoisson( double k, double lambda) {
    double cumProb{0.0};
    for (int n{0}; n <= k; ++n) {
        cumProb += dPoisson(n, lambda);
    }
  return cumProb;
}
