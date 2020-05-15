// Simulation of phase-type distributions

#include "Simulation.h"

//General functions

Numeric_lib::Matrix<double,2> EmbeddedMC(const Numeric_lib::Matrix<double,2> & T, const Numeric_lib::Matrix<double,2> & t) {
//  Returns the transition probabilities of the embeded markov chain determined by the transition matrix
    long p{T.dim1()};
    Numeric_lib::Matrix<double,2> Q(p + 1,p + 1);
    
    for (int i{0}; i < p; ++i) {
        for (int j{0}; j < p + 1; ++j) {
            if (j != i && j < p) {
                 Q(i,j) = -1.0 * T(i,j) / T(i,i);
            }
            else if(j == p) {
                Q(i,j) = -1.0 * t(i,0) / T(i,i);
            }
        }
    }
    Q(p,p) = 1;
    
    return Q;
}

long initialState(const Numeric_lib::Matrix<double,2> & cumulatedPi, double u) {
//Given the accumulated values of the initial probabilities (Pi) and a uniform value (u) returns the states that satisfies cumPi_(k-1)<u<cumPi_(k)
    
    if (u <= cumulatedPi(0,0)) {
        return 0;
    }
    
    for( int i{1}; i < cumulatedPi.dim2(); ++i) {
        if (cumulatedPi(0,i - 1) < u && u <= cumulatedPi(0,i)) {
            return i;
        }
    }
    return 0;
}

long newState(long previousState, const Numeric_lib::Matrix<double,2> & cumulatedEmbeddedMC, double u) {
//  Given the values of a probabilities matrix (Q), a uniform value (u) and a previous state,returns the new states of the Markov chain
    
    if (u <= cumulatedEmbeddedMC(previousState,0)) {
        return 0;
    }
    
    for (int i{1}; i < cumulatedEmbeddedMC.dim2(); ++i) {
        if (cumulatedEmbeddedMC(previousState,i - 1) < u && u <= cumulatedEmbeddedMC(previousState,i)) {
            return i;
        }
    }
    
    return 0;
}



//Simulation of univariate PH
void simulatePH(PhaseType & ph, int sampleSize, AbsRNG & rUniform, std::ofstream & outputFile, int transformation, const std::vector<double> & beta) {
//  Simulates a univariate phase-type distribution
    Numeric_lib::Matrix<double,2> T(ph.getT());
    Numeric_lib::Matrix<double,2> cumulatedEmbeddedMC(cumulateMatrix(EmbeddedMC(ph.getT(), ph.gett())));
    Numeric_lib::Matrix<double,2> cumulatedPi(cumulateMatrix(ph.getPi()));
    
    
    long p{ph.getp()};
    long state{0};
    for (int i{0}; i < sampleSize; ++i) {
        double time{0.0};
        state = initialState(cumulatedPi, rUniform());
        while (state != p) {
            time += log(1.0 - rUniform()) / T(state,state);
            state = newState(state, cumulatedEmbeddedMC, rUniform());
        }
        
        if (transformation == 1) {
            time = exp(time) - 1;
        }
        else if (transformation == 2) {
            time = beta[0] * (exp(time) - 1);
        }
        else if (transformation == 3) {
            time = pow(time, 1.0 / beta[0]);
        }
        else if (transformation == 4) {
            time = log(beta[0] * time + 1) / beta[0];
        }
        else if (transformation == 5) {
            time = beta[0] - beta[1] * log(time);
        }
        else if (transformation == 6) {
            time = beta[0] + beta[1]  * (pow(time, - beta[2]) - 1) / beta[2] ;
        }
        
        if (i < sampleSize - 1) {
            outputFile << time << '\n';
        }
        else {
            outputFile << time;
        }
    }
}

//Simulation of Multivariate PH
void simulateMPH(MPH & mph, int sampleSize, AbsRNG & rUniform, std::ofstream & outputFile, int transformation, std::vector<double> beta) {
    Numeric_lib::Matrix<double,2> T(mph.getT());
    Numeric_lib::Matrix<double,2> R(mph.getR());
    Numeric_lib::Matrix<double,2> cumulatedEmbeddedMC(cumulateMatrix(EmbeddedMC(mph.getT(), mph.gett())));
    Numeric_lib::Matrix<double,2> cumulatedPi(cumulateMatrix(mph.getPi()));
    
    long p{mph.getp()};
    long dim{mph.getd()};
    long state{0};
    for (int i{0}; i < sampleSize; ++i) {
        std::vector<double> sample(dim);
        double time{0.0};
        state = initialState(cumulatedPi, rUniform());
        while (state != p) {
            time = log(1.0 - rUniform()) / T(state,state);
            for (int j{0}; j < dim; ++j) {
                sample[j] += R(state,j) * time;
            }
            state = newState(state, cumulatedEmbeddedMC, rUniform());
        }
        for (int j{0}; j < dim; ++j) {
            if (transformation == 1) {
                sample[j] = (exp(sample[j]) - 1);
            }
            else if (transformation == 2) {
                sample[j] = beta[j] * (exp(sample[j]) - 1);
            }
            else if (transformation == 3) {
                sample[j] =  pow(sample[j], 1.0 / beta[j]);
            }
            else if (transformation == 4) {
                sample[j] =  log(sample[j] * beta[j * dim] + 1) / beta[j * dim];
            }
            else if (transformation == 5) {
                sample[j] = beta[j * dim] - beta[j * dim + 1] * log(sample[j]);
            }
            else if (transformation == 6) {
                sample[j] = beta[j * dim] + beta[j * dim + 1] * (pow(sample[j], - beta[j * dim + 2]) - 1) / beta[j * dim + 2];
            }
            
            if (j < dim - 1) {
                outputFile << sample[j] << '\t';
            }
            else {
                outputFile << sample[j];
            }
        }
        if (i < sampleSize - 1) {
            outputFile << '\n';
        }
    }
}



double upperTailDependence(MPH & mph, int sampleSize, double u, double epsilon, AbsRNG & rUniform) {
    Numeric_lib::Matrix<double,2> T(mph.getT());
    Numeric_lib::Matrix<double,2> R(mph.getR());
    Numeric_lib::Matrix<double,2> cumulatedEmbeddedMC(cumulateMatrix(EmbeddedMC(mph.getT(), mph.gett())));
    Numeric_lib::Matrix<double,2> cumulatedPi(cumulateMatrix(mph.getPi()));
    
    long p{mph.getp()};
    long dim{mph.getd()};
    long state{0};
    
    double q1{mph.quantile(0, u, epsilon)};
    double q2{mph.quantile(1, u, epsilon)};
    
    
    std::cout << "Quantiles:\n" << q1 << '\t' << q2 << '\n';
    
    int prob_indicator{0};
    
    int i{0};
    
    while (i < sampleSize) {
        std::vector<double> sample(dim);
        double time{0.0};
        state = initialState(cumulatedPi, rUniform());
        while (state != p) {
            time = log(1.0 - rUniform()) / T(state,state);
            for (int j{0}; j < dim; ++j) {
                sample[j] += R(state,j) * time;
            }
            state = newState(state, cumulatedEmbeddedMC, rUniform());
        }
        if (sample[1] > q2) {
            ++i;
            if (sample[0] > q1) {
                ++prob_indicator;
            }
        }
    }
    return static_cast<double>(prob_indicator) / sampleSize;
}



std::vector<double> upperTailDependence3d(MPH & mph, int sampleSize, double u, double epsilon, AbsRNG & rUniform) {
    Numeric_lib::Matrix<double,2> T(mph.getT());
    Numeric_lib::Matrix<double,2> R(mph.getR());
    Numeric_lib::Matrix<double,2> cumulatedEmbeddedMC(cumulateMatrix(EmbeddedMC(mph.getT(), mph.gett())));
    Numeric_lib::Matrix<double,2> cumulatedPi(cumulateMatrix(mph.getPi()));
    
    long p{mph.getp()};
    long dim{mph.getd()};
    long state{0};
    
    double q1{mph.quantile(0, u, epsilon)};
    double q2{mph.quantile(1, u, epsilon)};
    double q3{mph.quantile(2, u, epsilon)};
    
    std::vector<double> lambdaU{};
    
    std::cout << "Quantiles:\n" << q1 << '\t' << q2 << '\t' << q3 << '\n';
    
    int prob_indicator1{0};
    int prob_indicator2{0};
    int prob_indicator3{0};
    
    int i{0};
    int i1{0};
    int i2{0};
    int i3{0};
    
    while (i < sampleSize) {
        std::vector<double> sample(dim);
        double time{0.0};
        state = initialState(cumulatedPi, rUniform());
        while (state != p) {
            time = log(1.0 - rUniform()) / T(state,state);
            for (int j{0}; j < dim; ++j) {
                sample[j] += R(state,j) * time;
            }
            state = newState(state, cumulatedEmbeddedMC, rUniform());
        }
        if (sample[1] > q2) { //1-2
            ++i1;
            if (sample[0] > q1) {
                ++prob_indicator1;
            }
        }
        if (sample[2] > q3) { // 1 - 3
            ++i2;
            if (sample[0] > q1) {
                ++prob_indicator2;
            }
        }
        if (sample[2] > q3) { // 2 - 3
            ++i3;
            if (sample[1] > q2) {
                ++prob_indicator3;
            }
        }
        i = std::min(i1, std::min(i2,i3));
    }
    
    lambdaU.push_back(static_cast<double>(prob_indicator1) / (i1 + 1)) ;
    lambdaU.push_back(static_cast<double>(prob_indicator2) / (i2 + 1)) ;
    lambdaU.push_back(static_cast<double>(prob_indicator3) / (i3 + 1)) ;
    
    return lambdaU;
   // return static_cast<double>(prob_indicator) / sampleSize;
}



double lowerTailDependence(MPH & mph, int sampleSize, double u, double epsilon, AbsRNG & rUniform) {
    Numeric_lib::Matrix<double,2> T(mph.getT());
    Numeric_lib::Matrix<double,2> R(mph.getR());
    Numeric_lib::Matrix<double,2> cumulatedEmbeddedMC(cumulateMatrix(EmbeddedMC(mph.getT(), mph.gett())));
    Numeric_lib::Matrix<double,2> cumulatedPi(cumulateMatrix(mph.getPi()));
    
    long p{mph.getp()};
    long dim{mph.getd()};
    long state{0};
    
    double q1{mph.quantile(0, u, epsilon)};
    double q2{mph.quantile(1, u, epsilon)};
    
    std::cout << "Quantiles:\n" << q1 << '\t' << q2 << '\n';
    
    int prob_indicator{0};
    
    int i{0};
    
    while (i < sampleSize) {
        std::vector<double> sample(dim);
        double time{0.0};
        state = initialState(cumulatedPi, rUniform());
        while (state != p) {
            time = log(1.0 - rUniform()) / T(state,state);
            for (int j{0}; j<dim; ++j) {
                sample[j] += R(state,j) * time;
            }
            state = newState(state, cumulatedEmbeddedMC, rUniform());
        }
        if (sample[1] <= q2) {
            ++i;
            if (sample[0] <= q1) {
                ++prob_indicator;
            }
        }
    }
    return static_cast<double>(prob_indicator) / sampleSize;
}


// Simulation for NPH
double getLevel(std::vector<DiscreteDistribution> & theCdf, ScaleDistribution *N, double u){

    DiscreteDistribution auxiliar;
    double cum{};
    int k{};

    if (theCdf.size() > 0) {
        cum = theCdf[theCdf.size() - 1].cdf;
        k = static_cast<int>(theCdf.size());
    }

    while (cum < u) {
        ++k;
        auxiliar.x = N->si(k);
        cum += N->density(k);
        auxiliar.cdf = cum;
        theCdf.push_back(auxiliar);
    }


    if (u <= theCdf[0].cdf) {
        return theCdf[0].x;
    }

    for (int i{1}; i < theCdf.size(); ++i) {
        if(theCdf[i - 1].cdf < u && u <= theCdf[i].cdf) {
            return theCdf[i].x;
        }
    }
    return (theCdf[theCdf.size() - 1].x); // Creo que no lo necesito
}


void simulateNPH(PhaseType & ph, ScaleDistribution *N, int sampleSize, AbsRNG & rUniform, std::ofstream & outputFile) {
    
    long p{ph.getp()};
    Numeric_lib::Matrix<double,2> T(p,p);
    Numeric_lib::Matrix<double,2> Torg(ph.getT());
    Numeric_lib::Matrix<double,2> t(p, 1);
    Numeric_lib::Matrix<double,2> cumulatedEmbeddedMC(p + 1, p + 1);
    Numeric_lib::Matrix<double,2> cumulatedPi(cumulateMatrix(ph.getPi()));
    Numeric_lib::Matrix<double,2> e(ph.gete());
    
    long state{0};
    double level{};
    
    std::vector<DiscreteDistribution> theCdf;
    
    for (int i{0}; i < sampleSize; ++i) {
        level = getLevel(theCdf, N, rUniform());
        T = Torg / level;
        t = matrixProduct( T * (-1.0), e);
        cumulatedEmbeddedMC = cumulateMatrix(EmbeddedMC(T, t));
        double time{0.0};
        state = initialState(cumulatedPi, rUniform());
        while (state != p) {
            time += log(1.0 - rUniform()) / T(state,state);
            state = newState(state, cumulatedEmbeddedMC, rUniform());
        }
        
        if (i < sampleSize - 1) {
            outputFile << time << '\n';
        }
        else {
            outputFile << time;
        }
    }
}
