//Implementation of the EM algorithm

#include "EM.h"

double norm(double x1, double x2) {
    return sqrt(x1 * x1 + x2 * x2);
}

double norm(double x1, double x2, double x3) {
    return sqrt(x1 * x1 + x2 * x2 + x3 * x3);
}

int askTypeOfStepRK() {
    int typeOfStep{};
    do {
        std::cout << "\nChoose step-length for the Runge-Kutta procedure:\n";
        std::cout << "      1. Default value\n";
        std::cout << "      2. Your own choice of value\n";
        std::cout << "Select 1 or 2 :";
        std::cin >> typeOfStep;
        if (typeOfStep != 1 && typeOfStep != 2) {
            std::cout << "Please enter a valid option\n";
        }
    } while (typeOfStep != 1 && typeOfStep != 2);
    return typeOfStep;
}


double askStepLengthRK() {
    double RKstep{};
    do {
        std::cout << "\nStep-length: ";
        std::cin >> RKstep;
        if (RKstep <= 0) {
            std::cout << "The Step-length has to be postive.";
        }
    } while (RKstep <= 0);
    return RKstep;
}


// Provides a way to compare the elements of the structure "Sample"
bool compare(Sample &a, Sample &b) {
    return a.obs < b.obs;
}

// Sort a vector of the type sample.
void sortObservations(std::vector<Sample> & vec) {
    long newsize{0};
    newsize = vec.size();
    
    sort(vec.begin(), vec.end(), compare);
    
    for (long i = vec.size() - 2; i >= 0; --i) {
        if (vec[i].obs == vec[i + 1].obs) {
            vec[i].weight += vec[i + 1].weight;
            --newsize;
            for (long j = i + 1; j < newsize; ++j) {
                vec[j] = vec[j + 1] ;
            }
            vec.erase(vec.begin() + newsize);
        }
    }
}

//Prints a vector of doubles 
void printVector(const std::vector<double> & a) {
    for (int i{0}; i < a.size(); ++i) {
        std::cout << a[i] << " ";
    }
    std::cout << '\n';
}

void printVector(const std::vector<int> & a) {
    for (int i{0}; i < a.size(); ++i) {
        std::cout << a[i] << " ";
    }
    std::cout << '\n';
}


// Read univariate data from file and put it in a vector of type "Sample"
void readDataForPH(std::ifstream & inputFile, std::vector<Sample> & observations, std::vector<Sample> & censored, int indWeighted, int indCensored) {
    Sample aux;
    int typeObs{0};
    
    while (!inputFile.eof()) {
        if (indCensored == 2) {
            inputFile >> typeObs;
        }
        
        inputFile >> aux.obs;
        
        if (indWeighted == 2) {
            inputFile >> aux.weight;
        }
        else {
            aux.weight = 1.0;
        }
        
        if (typeObs == 0) {
            observations.push_back(aux);
        }
        else {
            censored.push_back(aux);
        }
    }
}


//Given a choice of density (See Distributions.h) put a sample in a vector of type "Sample" for approximation of a continous density.
void inputDensity(std::function<double(double, std::vector<double>)> density, const std::vector<double> & parameters, double truncationPoint, double maxProbability, double maxDeltat, std::vector<Sample> & observations) {
    Sample auxiliar;
    
    double deltat{};
    double t{0.0};
    
    while (t < truncationPoint) {
        if (density(t, parameters) < maxProbability / maxDeltat) {
            deltat = maxDeltat;
        }
        else {
            deltat = maxProbability / density(t, parameters);
        }
        auxiliar.weight = deltat / 6 * (density(t, parameters) + 4 * density(t + deltat / 2, parameters) + density(t + deltat, parameters));
        while (auxiliar.weight > maxProbability) {
            deltat = deltat * 0.9;
            auxiliar.weight = deltat / 6 * (density(t ,parameters) + 4 * density(t + deltat / 2, parameters) + density(t + deltat, parameters));
        }
        if (auxiliar.weight > 0) {
            auxiliar.obs = (t * density(t, parameters) + 4 * (t + deltat / 2) * density(t + deltat / 2, parameters) + (t + deltat) * density(t + deltat, parameters)) / (density(t, parameters) + 4 * density(t + deltat / 2, parameters) + density(t + deltat, parameters));
            observations.push_back(auxiliar);
        }
        t += deltat;
    }
}

// Gives basic information (size, sum of weigths, mean and sd) of a vector of type "Sample"
void sampleStatisticsPH(const std::vector<Sample> & observations) {
    double sumOfx{0.0};
    double sumOfxSquared{0.0};
    double sumOfWeigths{0.0};
    
    for (int i{0}; i < observations.size(); ++i) {
        sumOfx += observations[i].obs * observations[i].weight;
        sumOfxSquared += observations[i].obs * observations[i].obs * observations[i].weight;
        sumOfWeigths += observations[i].weight;
    }
    std::cout << "Sample size: " << observations.size() << '\n';
    std::cout << "Sum of weigths: " << sumOfWeigths << '\n';
    std::cout << "Sample mean: " << sumOfx / (sumOfWeigths) << '\n';
    std::cout << "Sample standard deviation: " << sqrt(sumOfxSquared / (sumOfWeigths) - pow(sumOfx / (sumOfWeigths), 2.0)) << '\n';
}

// Gives basic information (size, sum of weigths, mean and sd) of a vector of type "Sample"
void infoDataPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored) {
    double sumOfWeigths{0.0};
    double sumOfWeigthsRC{0.0};
    
    for (int i{0}; i < observations.size(); ++i) {
        sumOfWeigths += observations[i].weight;
    }
    for (int i{0}; i < censored.size(); ++i) {
        sumOfWeigthsRC += censored[i].weight;
    }

    std::cout << "Uncensored data.\n";
    std::cout << "    Sample size: " << observations.size() << '\n';
    std::cout << "    Sum of weigths: " << sumOfWeigths << '\n';
    std::cout << "Right-censored data.\n";
    std::cout << "    Sample size: " << censored.size() << '\n';
    std::cout << "    Sum of weigths: " << sumOfWeigthsRC << '\n';
    
}


// Given a file with the structure of a MPH (0's and 1's) returns the dimension of the MPH and the size p of the PH component - For a file for a PH returns dim = 0.
DimensionsMPH dimensionsOfFile(std::ifstream & inputFile) {
    DimensionsMPH sizeFile;
    
    sizeFile.m_p = std::count(std::istreambuf_iterator<char>(inputFile), std::istreambuf_iterator<char>(), '\n') + 1;
    
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);
    
    
    std::string firstline;
    
    getline(inputFile, firstline);
    int count{0};
    long length = firstline.length();
    for (int i {0}; i < length; ++i) {
        if (isspace(firstline[i]))
            ++count;
    }
    
    // If a R matrix is not given then it will be zero
    sizeFile.m_dim = count - sizeFile.m_p;
    
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);
    
    return sizeFile;
}

// Given a file with the structure of a MPH (0's and 1's)  puts the corresponding values in PiLegal and TLegal
void readStructureForPH(std::ifstream & inputFile, Numeric_lib::Matrix<double,2> & piLegal ,Numeric_lib::Matrix<double,2> & TLegal) {
    long p{dimensionsOfFile(inputFile).m_p};
    
    for (int i{0}; i < p; ++i) {
        inputFile >> piLegal(0,i);
        for (int j{0}; j < p; ++j) {
            inputFile >> TLegal(i,j);
        }
    }
}


int askForSampleType() {
    int input{};
    do {
        std::cout << "\nType of input:\n";
        std::cout << "      1. Sample file\n";
        std::cout << "      2. Density\n";
        std::cout << "Select 1 or 2: ";
        std::cin >> input;
        if (input != 1 && input != 2) {
            std::cout << "Please enter a valid option\n";
        }
    } while (input != 1 && input != 2);
    
    return input;
}

int askForStructure() {
    int structure{};
    do {
        std::cout << "\nSelect structure type:\n";
        std::cout << "      1. General phase-type with random initialization\n";
        std::cout << "      2. User specified structure from file and random initialization\n";
        std::cout << "      3. User specified starting values\n";
        std::cout << "Select(1-3): ";
        std::cin >> structure;
        if (structure != 1 && structure != 2 && structure != 3) {
            std::cout << "Please enter a valid option\n";
        }
    } while (structure != 1 && structure != 2 && structure != 3);
    
    return structure;
}

void askContinousDist(std::vector<Sample> & observations) {
    int typeContDist{};
    do {
        std::cout << "\nType of density:\n";
        std::cout << "      1. Rectangle\n";
        std::cout << "      2. Normal\n";
        std::cout << "      3. Lognormal\n";
        std::cout << "      4. Weibull\n";
        std::cout << "      5. Inverse Gaussian\n";
        std::cout << "      6. Pareto\n";
        std::cout << "      7. Gumbel\n";
        std::cout << "      8. Generalized extreme value distribution\n";
        std::cout << "      9. Gompertz\n";
        std::cout << "Select(1-8): ";
        std::cin >> typeContDist;
        if (typeContDist < 1 || typeContDist > 9) {
            std::cout << "Please enter a valid option\n";
        }
    } while (typeContDist < 1 || typeContDist > 9);
    
    
    std::vector<double>  parameters;
    std::function<double(double, std::vector<double>)> density;
    double auxValue{};
    switch(typeContDist) {
        case 1:
            std::cout << "\nUniform (rectangle) distribution between a and b:\n";
            std::cout << "a: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "b: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dUniform;
            break;
        case 2:
            std::cout << "\nNormal distribution with mean mu and standard deviation sigma:\n";
            std::cout << "mu: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "sigma: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dNormal;
            break;
        case 3:
            std::cout << "\nlognormal distribution with parameters mu and sigma:\n";
            std::cout << "mu: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "sigma: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dLogNormal;
            break;
        case 4:
            std::cout << "\nWeibull distribution with shape parameter tau and scale parameter lambda:\n";
            std::cout << "tau: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "lambda: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dWeibull;
            break;
        case 5:
            std::cout << "\nInverse Gaussian distribution with parameters lambda and mu:\n";
            std::cout << "lambda: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "mu:";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dInverseGaussian;
            break;
        case 6:
            std::cout << "\nPareto distribution with shape alpha and scale parameter kappa:\n";
            std::cout << "alpha: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "kappa: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dParetoTypeII;
            break;
        case 7:
            std::cout << "\nGumbel distribution with location parameter mu and scale parameter sigma:\n";
            std::cout << "mu: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "sigma: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dGumbel;
            break;
        case 8:
            std::cout << "\nGEVD distribution with location parameter mu, scale parameter sigma and shape parameter xi:\n";
            std::cout << "mu: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "sigma: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "xi: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dGEVD;
            break;
        case 9:
            std::cout << "\nGompertz distribution with scale parameter alpha and shape parameter beta:\n";
            std::cout << "alpha: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            std::cout << "beta: ";
            std::cin >> auxValue;
            parameters.push_back(auxValue);
            density = dGumbel;
            break;
    }
    
    std::cout << "Time at which density can be truncated: ";
    double truncationPoint{};
    std::cin >> truncationPoint;

    std::cout << "Maximum acceptable probability in one point: ";
    double maxProbability{};
    std::cin >> maxProbability;
    
    std::cout << "Maximum time interval corresponding to one point: ";
    double maxDeltat{};
    std::cin >> maxDeltat;
    
    inputDensity(density, parameters, truncationPoint, maxProbability, maxDeltat, observations);
    
}


double derivativeMatrixWeibull(const std::vector<Sample> & observations, PhaseType & PH, const std::vector<double> &  beta) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh +=  observations[k].weight * ( ( matrixProduct(PH.getPi(), matrixProduct(matrixExponential(PH.getT() * pow(observations[k].obs, beta[0])), matrixProduct(PH.getT(), PH.gett())))(0,0) * log(observations[k].obs) * pow(observations[k].obs, beta[0]) ) / PH.density(pow(observations[k].obs, beta[0])) + 1 / beta[0] + log(observations[k].obs) );
    }
    return logLh;
}

double derivativeMatrixPareto(const std::vector<Sample> & observations, PhaseType & PH, const std::vector<double> &  beta) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight *  ( ( matrixProduct(PH.getPi(), matrixProduct(matrixExponential(PH.getT() * log(observations[k].obs / beta[0] + 1)), matrixProduct(PH.getT(), PH.gett())))(0,0) * (- 1.0 * observations[k].obs) / ( beta[0] * (observations[k].obs + beta[0])) ) / PH.density(log(observations[k].obs / beta[0] + 1)) - 1 / (beta[0] + observations[k].obs) );
    }
    return logLh;
}

double derivativeMatrixGompertz(const std::vector<Sample> & observations, PhaseType & PH, const std::vector<double> & beta) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( ( matrixProduct(PH.getPi(), matrixProduct(matrixExponential(PH.getT() * (exp(observations[k].obs * beta[0]) - 1) / beta[0] ), matrixProduct(PH.getT(), PH.gett())))(0,0) * ( exp(observations[k].obs * beta[0]) * (observations[k].obs / beta[0] - 1 / (beta[0] * beta[0])) + 1 / (beta[0] * beta[0]) )) / PH.density((exp(observations[k].obs * beta[0]) - 1) / beta[0])  + observations[k].obs );
    }
    return logLh;
}

double firstDerivativeGumbel(const std::vector<Sample> & observations, PhaseType & PH, const std::vector<double> & beta) { //wrt mu
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( ( matrixProduct(PH.getPi(), matrixProduct(matrixExponential(PH.getT() * exp( -(observations[k].obs - beta[0]) / beta[1]) ), matrixProduct(PH.getT(), PH.gett())))(0,0) * ( exp( -(observations[k].obs - beta[0]) / beta[1]) / beta[1] ) ) / PH.density(exp( -(observations[k].obs - beta[0]) / beta[1]))  + 1 / beta[1] );
    }
    return logLh;
}
double secondDerivativeGumbel(const std::vector<Sample> & observations, PhaseType & PH, const std::vector<double> & beta) { //wrt sigma
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( ( matrixProduct(PH.getPi(), matrixProduct(matrixExponential(PH.getT() * exp( -(observations[k].obs - beta[0]) / beta[1]) ), matrixProduct(PH.getT(), PH.gett())))(0,0) * ( exp( -(observations[k].obs - beta[0]) / beta[1]) * (observations[k].obs - beta[0]) / (beta[1] * beta[1]) ) ) / PH.density(exp( -(observations[k].obs - beta[0]) / beta[1]))  - 1 / beta[1] + (observations[k].obs - beta[0]) / (beta[1] * beta[1]) );
    }
    return logLh;
}

std::vector<double> DerivativeGumbel(const std::vector<Sample> & observations, PhaseType & PH, const std::vector<double> & beta) {
    std::vector<double> derivatives{};
    derivatives.push_back(0);
    derivatives.push_back(0);
    
    double factor{};
    
    for (int k{0}; k < observations.size(); ++k) {
        factor = ( matrixProduct(PH.getPi(), matrixProduct(matrixExponential(PH.getT() * exp( -(observations[k].obs - beta[0]) / beta[1]) ), matrixProduct(PH.getT(), PH.gett())))(0,0) * exp( -(observations[k].obs - beta[0]) / beta[1]) ) / PH.density(exp( -(observations[k].obs - beta[0]) / beta[1])) ;
        
        derivatives[0] += observations[k].weight * ( factor / beta[1]   + 1 / beta[1] );
        derivatives[1] += observations[k].weight * ( factor * (observations[k].obs - beta[0]) / (beta[1] * beta[1])  - 1 / beta[1] + (observations[k].obs - beta[0]) / (beta[1] * beta[1]) );
    }
    return derivatives;
}

std::vector<double> DerivativeGEVD(const std::vector<Sample> & observations, PhaseType & PH, const std::vector<double> & beta) {
    std::vector<double> derivatives{};
    derivatives.push_back(0);
    derivatives.push_back(0);
    derivatives.push_back(0);
    
    double factor{};
    double y{};
    
    for (int k{0}; k < observations.size(); ++k) {
        y = 1 + beta[2] * (observations[k].obs - beta[0]) / beta[1];
        factor = ( matrixProduct(PH.getPi(), matrixProduct(matrixExponential(PH.getT() * pow(y, -1.0 / beta[2]) ), matrixProduct(PH.getT(), PH.gett())))(0,0) ) / PH.density(pow(y, -1.0 / beta[2]) ) ;
        
        derivatives[0] += observations[k].weight * ( factor * pow(y, -(1 + 1 / beta[2])) / (beta[1]) + (beta[2] + 1) / (beta[1] * y)  ); //Mu
        derivatives[1] += observations[k].weight * ( factor * pow(y, -(1 + 1 / beta[2])) * (observations[k].obs - beta[0]) / (beta[1] * beta[1]) - 1 / beta[1] + (observations[k].obs - beta[0]) * (beta[2] + 1) / (beta[1] * beta[1] * y) ); //Sigma
        derivatives[2] += observations[k].weight * ( factor * pow(y, -1 / beta[2]) * ( log(y) / (beta[2] * beta[2]) - (observations[k].obs - beta[0]) / (beta[2] * beta[1] * y) )  + log(y) / (beta[2] * beta[2]) - (1 + 1 / beta[2]) * (observations[k].obs - beta[0]) / (beta[1] * y) ); //Xi
    }
    return derivatives;
}


// Transform data for a IPH fit
void transformData(const std::vector<Sample> & observations, std::vector<Sample> & TransformObs, int transformation, const std::vector<double> & beta) {
    
    if (transformation == 1) { // Pareto
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].obs = log(observations[i].obs + 1);
        }
    }
    else if (transformation == 2) { // Pareto parameter-dependent
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].obs = log(observations[i].obs / beta[0] + 1);
        }
    }
    else if (transformation == 3) { // Weibull
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].obs = pow(observations[i].obs, beta[0]);
        }
    }
    else if (transformation == 4) { // Gompertz
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].obs = (exp(observations[i].obs * beta[0]) - 1) / beta[0];
        }
    }
    else if (transformation == 5) { // Gumbel
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].obs = exp( -(observations[i].obs - beta[0]) / beta[1]) ;
        }
    }
    else if (transformation == 6) { // GEVD
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].obs = pow( 1 + beta[2] * (observations[i].obs - beta[0]) / beta[1] , -1 / beta[2]) ;
        }
    }
    
}

void reversTransformData(const std::vector<Sample> & observations, std::vector<Sample> & TransformObs, int transformation, const std::vector<double> & beta) {
    size_t N{observations.size()};
    if (transformation == 5) { // Gumbel
        for (int i{0}; i < N; ++i) {
            TransformObs[i].obs = exp( -(observations[N - i - 1].obs - beta[0]) / beta[1]) ;
            TransformObs[i].weight = observations[N - i - 1].weight;
        }
    }
    else if (transformation == 6) { // GEVD
        for (int i{0}; i < N; ++i) {
            TransformObs[i].obs = pow( 1 + beta[2] * (observations[N - i - 1].obs - beta[0]) / beta[1] , -1 / beta[2]);
            TransformObs[i].weight = observations[N - i - 1].weight;
        }
    }
    
}



// Runge Kutta for the calculation of the a,b and c vectors in the EM step
void rungeKutta(Numeric_lib::Matrix<double,2> & avector, Numeric_lib::Matrix<double,2> & bvector, Numeric_lib::Matrix<double,2> & cmatrix, double dt, double h, const Numeric_lib::Matrix<double,2> & T, const Numeric_lib::Matrix<double,2> & t) {
    long p{T.dim1()};
    int j{};
    int m{};
    double eps{};
    double sum{};
    
    int i{};
    i = dt / h;
    
    double h2{};
    h2 = dt / (i + 1);
    
    Numeric_lib::Matrix<double,2> ka(4,p);
    Numeric_lib::Matrix<double,2> kb(4,p);
    Numeric_lib::Matrix<double,3> kc(4,p,p);
    
    for (eps = 0; eps <= dt - h2 / 2; eps += h2) {
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(j,i) * avector(0,j);
            }
            ka(0,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(j,i) * (avector(0,j) + ka(0,j) / 2);
            }
            ka(1,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(j,i) * (avector(0,j) + ka(1,j) / 2);
            }
            ka(2,i) = h2 * sum;
        }
        for (i=0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j)
            {
                sum += T(j,i) * (avector(0,j) + ka(2,j));
            }
            ka(3,i) = h2 * sum;
        }
        
        
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(i,j) * bvector(j,0);
            }
            kb(0,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(i,j) * (bvector(j,0) + kb(0,j) / 2);
            }
            kb(1,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(i,j) * (bvector(j,0) + kb(1,j) / 2);
            }
            kb(2,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(i,j) * (bvector(j,0) + kb(2,j));
            }
            kb(3,i) = h2 * sum;
        }
        
        for (m = 0; m < p; ++m) {
            for (i = 0; i < p; ++i) {
                sum = t(m,0) * avector(0,i);
                for (j = 0; j < p; ++j) {
                    sum += T(m,j) * cmatrix(j,i);
                }
                kc(0,m,i) = h2 * sum;
            }
        }
        for (m = 0; m < p; ++m) {
            for (i = 0; i < p; ++i) {
                sum = t(m,0) * (avector(0,i) + ka(0,i) / 2);
                for (j = 0; j < p; ++j) {
                    sum += T(m,j) * (cmatrix(j,i) + kc(0,j,i) / 2);
                }
                kc(1,m,i) = h2 * sum;
            }
        }
        for (m = 0; m < p; ++m) {
            for (i = 0; i < p; ++i) {
                sum = t(m,0) * (avector(0,i) + ka(1,i) / 2);
                for (j = 0; j < p; ++j) {
                    sum += T(m,j) * (cmatrix(j,i) + kc(1,j,i) / 2);
                }
                kc(2,m,i) = h2 * sum;
            }
        }
        for (m = 0; m < p; ++m) {
            for (i = 0; i < p; ++i) {
                sum = t(m,0) * (avector(0,i) + ka(2,i));
                for (j = 0; j < p; ++j) {
                    sum += T(m,j) * (cmatrix(j,i) + kc(2,j,i));
                }
                kc(3,m,i) = h2 * sum;
            }
        }
        
        for (i = 0; i < p; ++i) {
            avector(0,i) += (ka(0,i) + 2 * ka(1,i) + 2 * ka(2,i) + ka(3,i)) / 6;
            bvector(i,0) += (kb(0,i) + 2 * kb(1,i) + 2 * kb(2,i) + kb(3,i)) / 6;
            for (j = 0; j < p; ++j) {
                cmatrix(i,j) +=(kc(0,i,j) + 2 * kc(1,i,j) + 2 * kc(2,i,j) + kc(3,i,j)) / 6;
            }
        }
    }
}

// Default size of the steps in the RK
double defaultStepLength(const Numeric_lib::Matrix<double,2> & T) {
    double h{-0.1 / T(0,0)};

    for (int i{1}; i < T.dim1(); ++i) {
        if (h > -0.1 / T(i,i)) {
            h = -0.1 / T(i,i);
        }
    }
    return h;
}



// EM using Runge Kutta
void EMstepRungeKutta(double h, PhaseType & PH, const std::vector<Sample> & observations, const std::vector<Sample> & censored) {
    long p{PH.getp()};
    Numeric_lib::Matrix<double,2> pi(PH.getPi());
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    Numeric_lib::Matrix<double,2> e(PH.gete());
    
    Numeric_lib::Matrix<double,2> Bmean(p,1);
    Numeric_lib::Matrix<double,2> Zmean(p,1);
    Numeric_lib::Matrix<double,2> Nmean(p,p + 1);
    
    Numeric_lib::Matrix<double,2> avector(1,p);
    Numeric_lib::Matrix<double,2> bvector(p,1);
    Numeric_lib::Matrix<double,2> cmatrix(p,p);
    
    // initial conditions
    avector = pi;
    bvector = t;
    
    double dt{0.0};
    if (observations.size()>0) {
        dt = observations[0].obs;
    }
    
    double SumOfWeights{0.0};
    double density{0.0};
    
    // E step
    //  Unccensored data
    for (int k{0}; k < observations.size(); ++k) {
            
        SumOfWeights += observations[k].weight;
            
        rungeKutta(avector, bvector, cmatrix, dt, h, T, t);
        density = matrixProduct(pi, bvector)(0,0);
            
        // E-step
        for (int i{0}; i < p; ++i) {
            Bmean(i,0) += pi(0,i) * bvector(i,0) * observations[k].weight / density;
            Nmean(i,p) += avector(0,i) * t(i,0) * observations[k].weight / density;
            Zmean(i,0) += cmatrix(i,i) * observations[k].weight / density;
            for (int j{0}; j < p; ++j) {
                Nmean(i,j) += T(i,j) * cmatrix(j,i) * observations[k].weight / density;
            }
        }
        
        dt = observations[k + 1].obs - observations[k].obs;
    }
    
    //  Right-Censored Data
    double SumOfCensored{0.0};
    if (censored.size() > 0) {
        dt = censored[0].obs;
        cmatrix *= 0;
        avector = pi;
        bvector = e;
    }
    for (int k{0}; k < censored.size(); ++k) {
            
        SumOfCensored += censored[k].weight;
            
        rungeKutta(avector, bvector, cmatrix, dt, h, T, e);
        density = matrixProduct(pi, bvector)(0,0);
            
        //E-step
        for (int i{0}; i < p; ++i) {
            Bmean(i,0) += pi(0,i) * bvector(i,0) * censored[k].weight / density;
            Zmean(i,0) += cmatrix(i,i) * censored[k].weight / density;
            for (int j{0}; j < p; ++j) {
                Nmean(i,j) += T(i,j) * cmatrix(j,i) * censored[k].weight / density;
            }
        }
        
        dt = censored[k + 1].obs - censored[k].obs;
    }
    
    // M step
    for (int i{0}; i < p; ++i) {
        pi(0,i) = Bmean(i,0) / (SumOfWeights + SumOfCensored);
        if (pi(0,i) < 0) {
            pi(0,i) = 0;
        }
        t(i,0) = Nmean(i,p) / Zmean(i,0);
        if (t(i,0) < 0) {
            t(i,0) = 0;
        }
        T(i,i) = -t(i,0);
        for (int j{0}; j < p; ++j) {
            if (i != j) {
                T(i,j) = Nmean(i,j) / Zmean(i,0);
                if (T(i,j) < 0) {
                    T(i,j) = 0;
                }
                T(i,i) -= T(i,j);
            }
        }
    }
    PH.changeAll(pi, T);
}



void a_rungekutta(Numeric_lib::Matrix<double,2> & avector, double dt, double h, const Numeric_lib::Matrix<double,2> & T) {
    long p{T.dim1()};
    int j{};
    double eps{};
    double sum{};

    int i{};
    i = dt / h;
    
    double h2{};
    h2 = dt / (i + 1);
    
    Numeric_lib::Matrix<double,2> ka(4,p);

    for (eps = 0; eps <= dt - h2 / 2; eps += h2) {
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(j,i) * avector(0,j);
            }
            ka(0,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(j,i) * (avector(0,j) + ka(0,j) / 2);
            }
            ka(1,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(j,i) * (avector(0,j) + ka(1,j) / 2);
            }
            ka(2,i) = h2 * sum;
        }
        for (i = 0; i < p; ++i) {
            sum = 0;
            for (j = 0; j < p; ++j) {
                sum += T(j,i) * (avector(0,j) + ka(2,j));
            }
            ka(3,i) = h2 * sum;
        }

        for (i = 0; i < p; ++i) {
            avector(0,i) += (ka(0,i) + 2 * ka(1,i) + 2 * ka(2,i) + ka(3,i)) / 6;
        }
    }
}


// EM using Matlab algorithm
void EMstep(const std::vector<Sample> & observations, const std::vector<Sample> & censored,  PhaseType & PH) {
    long p{PH.getp()};
    Numeric_lib::Matrix<double,2> pi(PH.getPi());
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    Numeric_lib::Matrix<double,2> e(PH.gete());
    
    Numeric_lib::Matrix<double,2> Bmean(p,1);
    Numeric_lib::Matrix<double,2> Zmean(p,1);
    Numeric_lib::Matrix<double,2> Nmean(p,p + 1);
    
    Numeric_lib::Matrix<double,2> avector(1,p);
    Numeric_lib::Matrix<double,2> bvector(p,1);
    Numeric_lib::Matrix<double,2> cmatrix(p,p);
    Numeric_lib::Matrix<double,2> aux_exp(p,p);
    
    Numeric_lib::Matrix<double,2> J(2 * p,2 * p);
    Numeric_lib::Matrix<double,2> tProductPi(p,p);
    matrixProduct(t, pi, tProductPi);
    
    double SumOfWeights{0.0};
    double density{0.0};
    
    //E-step
    //  Unccensored data
    for (int k{0}; k < observations.size(); ++k) {
        SumOfWeights += observations[k].weight;
        
        VanLoanForIntegrals(observations[k].obs, T, T, tProductPi, J);
        
        for (int i{0}; i < p; ++i) {
            for (int j{0}; j < p; ++j) {
                aux_exp(i,j) = J(i,j);
                cmatrix(i,j) = J(i,j + p);
            }
        }
        
        matrixProduct(pi, aux_exp, avector);
        matrixProduct(aux_exp, t, bvector);
        density = matrixProduct(pi, bvector)(0,0);
        
        //E-step
        for (int i{0}; i < p; ++i) {
            Bmean(i,0) += pi(0,i) * bvector(i,0) * observations[k].weight / density;
            Nmean(i,p) += avector(0,i) * t(i,0) * observations[k].weight / density;
            Zmean(i,0) += cmatrix(i,i) * observations[k].weight / density;
            for (int j{0}; j < p; ++j) {
                Nmean(i,j) += T(i,j) * cmatrix(j,i) * observations[k].weight / density;
            }
        }
    }
    //  Right-Censored Data
    double SumOfCensored{0.0};
    if (censored.size() > 0) {
        matrixProduct(e, pi, tProductPi);
    }
    for (int k{0}; k < censored.size(); ++k) {
        SumOfCensored += censored[k].weight;
        
        VanLoanForIntegrals(censored[k].obs, T, T, tProductPi, J);
        
        for (int i{0}; i < p; ++i) {
            for (int j{0}; j < p; ++j) {
                aux_exp(i,j) = J(i,j);
                cmatrix(i,j) = J(i,j + p);
            }
        }
        
        matrixProduct(aux_exp, e, bvector);
        density = matrixProduct(pi, bvector)(0,0);
        
        //E-step
        for (int i{0}; i < p; ++i) {
            Bmean(i,0) += pi(0,i) * bvector(i,0) * censored[k].weight / density;
            Zmean(i,0) += cmatrix(i,i) * censored[k].weight / density;
            for (int j{0}; j < p; ++j) {
                Nmean(i,j) += T(i,j) * cmatrix(j,i) * censored[k].weight / density;
            }
        }
    }
    
    // M step
    for (int i{0}; i < p; ++i) {
        pi(0,i) = Bmean(i,0) / (SumOfWeights + SumOfCensored);
        if (pi(0,i) < 0) {
            pi(0,i) = 0;
        }
        t(i,0) = Nmean(i,p) / Zmean(i,0);
        if (t(i,0) < 0) {
            t(i,0) = 0;
        }
        T(i,i) = -t(i,0);
        for (int j{0}; j < p; ++j) {
            if (i != j) {
                T(i,j) = Nmean(i,j) / Zmean(i,0);
                if (T(i,j) < 0) {
                    T(i,j) = 0;
                }
                T(i,i) -= T(i,j);
            }
        }
    }
    PH.changeAll(pi, T);
}



// Loglikelihood using Matlab matrix exponential
double logLikelihoodPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, int transformation, std::vector<double> beta) {
    double logLh{0.0};
    if (transformation == 0) { // No transformation
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * log(PH.density(observations[k].obs));
        }
        for (int k{0}; k < censored.size(); ++k) {
            logLh += censored[k].weight * log(PH.tail(censored[k].obs));
        }
    }
    else if (transformation == 1) { // Pareto
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * (log(PH.density(log(observations[k].obs + 1))) - log(observations[k].obs + 1) );
        }
        for (int k{0}; k < censored.size(); ++k) {
            logLh += censored[k].weight * log(PH.tail(log(observations[k].obs + 1)));
        }
    }
    else if (transformation == 2) { // Pareto parameter-dependent
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * (log(PH.density(log(observations[k].obs / beta[0] + 1))) - log(observations[k].obs / beta[0] + 1) - log(beta[0]) );
        }
        for (int k{0}; k < censored.size(); ++k) {
            logLh += censored[k].weight * log(PH.tail(log(observations[k].obs / beta[0] + 1)));
        }
    }
    else if (transformation == 3) { // Weibull
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * (log(PH.density(pow(observations[k].obs, beta[0]))) + log(beta[0]) + (beta[0] - 1) * log(observations[k].obs ));
        }
        for (int k{0}; k < censored.size(); ++k) {
            logLh += censored[k].weight * log(PH.tail(pow(censored[k].obs, beta[0])));
        }
    }
    else if (transformation == 4) { // Gompertz
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * (log(PH.density((exp(beta[0] * observations[k].obs) - 1) / beta[0])) + beta[0] * observations[k].obs);
        }
        for (int k{0}; k < censored.size(); ++k) {
            logLh += censored[k].weight * log(PH.tail((exp(beta[0] * censored[k].obs) - 1) / beta[0]));
        }
    }
    else if (transformation == 5) { // Gumbel
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * (log(PH.density(exp(-(observations[k].obs - beta[0]) / beta[1]))) - log(beta[1]) -  (observations[k].obs - beta[0]) / beta[1]);
        }
        for (int k{0}; k < censored.size(); ++k) {
            logLh += censored[k].weight * log(PH.tail(exp(-(censored[k].obs - beta[0]) / beta[1])));
        }
    }
    else if (transformation == 6) { // GEVD
        double z{};
        for (int k{0}; k < observations.size(); ++k) {
            z = pow(1 + (beta[2] / beta[1]) * (observations[k].obs - beta[0]) , - 1 / beta[2]);
            logLh += observations[k].weight * (log(PH.density(z)) - log(beta[1]) + (beta[2] + 1) * log(z));
        }
        for (int k{0}; k < censored.size(); ++k) {
            logLh += censored[k].weight * log(PH.tail(pow(1 + (beta[2] / beta[1]) * (censored[k].obs - beta[0]) , - 1 / beta[2])));
        }
    }

    return logLh;
}


// Loglikelihood using RK
double logLikelihoodPH_RK(double h, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, int transformation, std::vector<double> beta) {
    long p{PH.getp()};
    Numeric_lib::Matrix<double,2> avector(1,p);
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    Numeric_lib::Matrix<double,2> e(PH.gete());
    
    
    // Uncensored data
    //   initial condition
    avector = PH.getPi();
    
    double dt{0.0};

    double density{0.0};
    
    double logLh{0.0};
    
    if (transformation == 0) { // No transformation
        // Non censored data
        if (observations.size() > 0) {
               dt = observations[0].obs;
        }
        for (int k{0}; k < observations.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, t)(0,0);
            logLh += observations[k].weight * log(density);
            dt = observations[k + 1].obs - observations[k].obs;
        }
        //Right censored data
        if (censored.size() > 0) {
            dt = censored[0].obs;
            avector = PH.getPi();
        }
        for (int k{0}; k < censored.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, e)(0,0);
            logLh += censored[k].weight * log(density);
            dt = censored[k + 1].obs - censored[k].obs;
        }
    }
    else if (transformation == 1 || transformation == 2) {  // Pareto
        if (transformation==1) {
            beta.push_back(0);
            beta[0] = 1;
        }
        
        // Non censored data
        if (observations.size() > 0) {
               dt = log(observations[0].obs / beta[0] + 1);
        }
        for (int k{0}; k < observations.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, t)(0,0);
            logLh += observations[k].weight * (log(density)  - log(observations[k].obs / beta[0] + 1) - log(beta[0]));
            dt = log(observations[k + 1].obs / beta[0] + 1) - log(observations[k].obs / beta[0] + 1);
        }
        //Right censored data
        if (censored.size() > 0) {
            dt = log(censored[0].obs / beta[0] + 1);
            avector = PH.getPi();
        }
        for (int k{0}; k < censored.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, e)(0,0);
            logLh += censored[k].weight * log(density);
            dt = log(censored[k + 1].obs / beta[0] + 1) - log(censored[k].obs / beta[0] + 1);
        }
    }
    else if (transformation == 3) { // Weibull
        // Non censored data
        if (observations.size() > 0) {
               dt = pow(observations[0].obs, beta[0]);
        }
        for (int k{0}; k < observations.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, t)(0,0);
            logLh += observations[k].weight * (log(density) + log(beta[0]) + (beta[0] - 1) * log(observations[k].obs));
            dt = pow(observations[k + 1].obs, beta[0]) - pow(observations[k].obs, beta[0]);
        }
        //Right censored data
        if (censored.size() > 0) {
            dt = pow(censored[0].obs, beta[0]);
            avector = PH.getPi();
        }
        for (int k{0}; k < censored.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, e)(0,0);
            logLh += censored[k].weight * log(density);
            dt = pow(censored[k + 1].obs, beta[0]) - pow(censored[k].obs, beta[0]);
        }
    }
    else if (transformation == 4) { // Gompertz
        // Non censored data
        if (observations.size() > 0) {
               dt = (exp(beta[0] * observations[0].obs) - 1) / beta[0];
        }
        for (int k{0}; k < observations.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, t)(0,0);
            logLh += observations[k].weight * (log(density) +  beta[0] * observations[k].obs);
            dt = (exp(beta[0] * observations[k + 1].obs) - 1) / beta[0] - (exp(beta[0] * observations[k].obs) - 1) / beta[0];
        }
        //Right censored data
        if (censored.size() > 0) {
            dt = (exp(beta[0] * censored[0].obs) - 1) / beta[0];
            avector = PH.getPi();
        }
        for (int k{0}; k < censored.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, e)(0,0);
            logLh += censored[k].weight * log(density);
            dt = (exp(beta[0] * censored[k + 1].obs) - 1) / beta[0] - (exp(beta[0] * censored[k].obs) - 1) / beta[0];
        }
    }
    else if (transformation == 5) { // Gumbel
        // Non censored data
        size_t N{observations.size()};
        if (N > 0) {
               dt = exp(-(observations[N - 1].obs - beta[0]) / beta[1]);
        }
        for (size_t k{1}; k <= N; ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, t)(0,0);
            logLh += observations[N - k].weight * (log(density) - log(beta[1]) -  (observations[N - k].obs - beta[0]) / beta[1]);
            dt = exp(-(observations[N - k - 1].obs - beta[0]) / beta[1]) - exp(-(observations[N - k].obs - beta[0]) / beta[1]);
        }
        //Right censored data
        N = censored.size();
        if (N > 0) {
            dt = exp(-(censored[N - 1].obs - beta[0]) / beta[1]);
            avector = PH.getPi();
        }
        for (int k{0}; k < censored.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, e)(0,0);
            logLh += censored[N - k].weight * log(density);
            dt = exp(-(censored[N - k - 1].obs - beta[0]) / beta[1]) - exp(-(censored[N - k].obs - beta[0]) / beta[1]);
        }
    }
    else if (transformation == 6) { // GEVD
        // Non censored data
        size_t N{observations.size()};
        if (N > 0) {
               dt = pow(1 + (beta[2] / beta[1]) * (observations[N - 1].obs - beta[0]) , - 1 / beta[2]);
        }
        for (size_t k{1}; k <= N; ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, t)(0,0);
            logLh += observations[N - k].weight * (log(density)  - log(beta[1]) - (1 + 1 / beta[2]) * log(1 + (beta[2] / beta[1]) * (observations[N - k].obs - beta[0])) );
            dt = pow(1 + (beta[2] / beta[1]) * (observations[N - k - 1].obs - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (observations[N - k].obs - beta[0]) , - 1 / beta[2]);
        }
        //Right censored data
        N = censored.size();
        if (N > 0) {
            dt = pow(1 + (beta[2] / beta[1]) * (censored[N - 1].obs - beta[0]) , - 1 / beta[2]);
            avector = PH.getPi();
        }
        for (int k{0}; k < censored.size(); ++k) {
            a_rungekutta(avector, dt, h, T);
            density = matrixProduct(avector, e)(0,0);
            logLh += censored[N - k].weight * log(density);
            dt = pow(1 + (beta[2] / beta[1]) * (censored[N - k - 1].obs - beta[0]) , - 1 / beta[2]) - pow(1 + (beta[2] / beta[1]) * (censored[N - k].obs - beta[0]) , - 1 / beta[2]);
        }
    }
    
    
    return logLh;
}


// Loglikelihood using Uniformization
double logLikelihoodPH_Uni(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, double epsilon, int transformation, std::vector<double> beta) {
    Numeric_lib::Matrix<double,2> pi(PH.getPi());
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    Numeric_lib::Matrix<double,2> e(PH.gete());
    
    //Uncensored data

    double density{0.0};
    double logLh{0.0};
    
    double alpha{0.0};
    alpha = matrixMaxDiagonal(T * (-1.0));
    
    int N{findN(epsilon, 1)};
    
    std::vector<Numeric_lib::Matrix<double,2>> theVector;
    
    vectorOfMatrices(theVector, T, alpha, N);
    
    if (transformation == 0) { // No transformation
        //Non-censored
        for (int k{0}; k < observations.size(); ++k) {
            double x{observations[k].obs};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), t)(0,0);
            logLh += observations[k].weight * log(density);
        }
        //Right censored data
        for (int k{0}; k < censored.size(); ++k) {
            double x{censored[k].obs};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), e)(0,0);
            logLh += censored[k].weight * log(density);
        }
    }
    
    else if (transformation == 1 || transformation == 2) { // Pareto
        if (transformation==1) {
            beta.push_back(0);
            beta[0] = 1;
        }
        //Non-censored
        for (int k{0}; k < observations.size(); ++k) {
            double x{log(observations[k].obs / beta[0] + 1)};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), t)(0,0);
            logLh += observations[k].weight * (log(density)  - log(observations[k].obs / beta[0] + 1) - log(beta[0]));
        }
        //Right censored data
        for (int k{0}; k < censored.size(); ++k) {
            double x{log(censored[k].obs / beta[0] + 1)};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), e)(0,0);
            logLh += censored[k].weight * log(density);
        }
    }
    else if (transformation == 3) { // Weibull
        //Non-censored
        for (int k{0}; k < observations.size(); ++k) {
            double x{pow(observations[k].obs, beta[0])};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), t)(0,0);
            logLh += observations[k].weight * (log(density) + log(beta[0]) + (beta[0] - 1) * log(observations[k].obs));
        }
        //Right censored data
        for (int k{0}; k < censored.size(); ++k) {
            double x{pow(censored[k].obs, beta[0])};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), e)(0,0);
            logLh += censored[k].weight * log(density);
        }
    }
    else if (transformation == 4) { // Gompertz
        //Non-censored
        for (int k{0}; k < observations.size(); ++k) {
            double x{(exp(beta[0] * observations[k].obs) - 1) / beta[0]};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), t)(0,0);
            logLh += observations[k].weight * (log(density) + beta[0] * observations[k].obs);
        }
        //Right censored data
        for (int k{0}; k < censored.size(); ++k) {
            double x{(exp(beta[0] * censored[k].obs) - 1) / beta[0]};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), e)(0,0);
            logLh += censored[k].weight * log(density);
        }
    }
    else if (transformation == 5) { // Gumbel
        //Non-censored
        for (int k{0}; k < observations.size(); ++k) {
            double x{exp(-(observations[k].obs - beta[0]) / beta[1])};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), t)(0,0);
            logLh += observations[k].weight * (log(density) - log(beta[1]) -  (observations[k].obs - beta[0]) / beta[1]);
        }
        //Right censored data
        for (int k{0}; k < censored.size(); ++k) {
            double x{exp(-(censored[k].obs - beta[0]) / beta[1])};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), e)(0,0);
            logLh += censored[k].weight * log(density);
        }
    }
    else if (transformation == 6) { // GEVD
        //Non-censored
        for (int k{0}; k < observations.size(); ++k) {
            double x{pow(1 + (beta[2] / beta[1]) * (observations[k].obs - beta[0]) , - 1 / beta[2])};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), t)(0,0);
            logLh += observations[k].weight * (log(density) - log(beta[1]) - (1 + 1 / beta[2]) * log(1 + (beta[2] / beta[1]) * (observations[k].obs - beta[0])) );
        }
        //Right censored data
        for (int k{0}; k < censored.size(); ++k) {
            double x{pow(1 + (beta[2] / beta[1]) * (censored[k].obs - beta[0]) , - 1 / beta[2])};
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, T);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, T);
                pow2Matrix(n, T);
            }
            density = matrixProduct(matrixProduct(pi, T), e)(0,0);
            logLh += censored[k].weight * log(density);
        }
    }
    
    
    return logLh;
}


// EM using Matlab matrix exponential
void EMIterate(int stepsEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, int printLikehood, int transformation, std::vector<double> & beta, double lambda, double epsilon) {
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    //double epsilon{0.001}; //Stopping criteria for gradient ascent
    
    logLh = logLikelihoodPH(observations, censored, phaseTypeEM, transformation, beta);
    logLh_ant = logLh;
    
    if (transformation == 0) {
        for (int i{1}; i <= stepsEM; ++i) {
            EMstep(observations, censored, phaseTypeEM);
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH(observations, censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    
    else if (transformation == 1) {
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
        
        transformData(observations, TransformObs, transformation, beta);
        transformData(censored, TransformCensored, transformation, beta);
        
        for (int i{1}; i <= stepsEM; ++i) {
            EMstep(TransformObs, TransformCensored, phaseTypeEM);
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH(observations, censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    
    else if (transformation == 2 || transformation == 3 || transformation == 4) {
        double fxnprime{};
        
        std::function<double(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
        
        if (transformation == 2) {
            derivative = derivativeMatrixPareto;
        }
        else if (transformation == 3)  {
            derivative = derivativeMatrixWeibull;
        }
        else if (transformation == 4)  {
            derivative = derivativeMatrixGompertz;
        }
        
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
        
        for (int i{1}; i <= stepsEM; ++i) {
            transformData(observations, TransformObs, transformation, beta);
            transformData(censored, TransformCensored, transformation, beta);
            
            EMstep(TransformObs, TransformCensored, phaseTypeEM);
            
            //Gradient ascent
            // Stoping condition relative error lest than epsilon
//            double beta_new{};
//            while (relError > epsilon) {
//                fxnprime = FnForNR(observations, phaseTypeEM, beta);
//                beta_new = beta + lambda * fxnprime;
//                relError = abs( (beta_new - beta) / beta );
//                beta = beta_new;
//                //std::cout<< fxnprime << '\t' << beta << '\n';
//            }
            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
            
            while (abs(fxnprime) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime;
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
            
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH(observations, censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". Beta: " << beta[0] << '\n';
            }
        }
    }
    
//    else if (transformation == 5) { // Gumbel
//        double fxnprime1{};
//        double fxnprime2{};
//
//        std::function<double(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> partial1;
//        std::function<double(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> partial2;
//
//        partial1 = firstDerivativeGumbel;
//        partial2 = secondDerivativeGumbel;
//
//        std::vector<Sample> TransformObs{observations};
//        std::vector<Sample> TransformCensored{censored};
//
//        for (int i{1}; i <= stepsEM; ++i) {
//            transformData(observations, TransformObs, transformation, beta);
//            transformData(censored, TransformCensored, transformation, beta);
//
//            EMstep(TransformObs, TransformCensored, phaseTypeEM);
//
//
//            // Stoping condition derivarite close to zero
//            fxnprime1 = partial1(observations, phaseTypeEM, beta);
//            fxnprime2 = partial2(observations, phaseTypeEM, beta);
//
//            while (norm(fxnprime1, fxnprime2) > epsilon) {
//                beta[0] = beta[0] + lambda * fxnprime1;
//                beta[1] = beta[1] + lambda * fxnprime2;
//                fxnprime1 = partial1(observations, phaseTypeEM, beta);
//                fxnprime2 = partial2(observations, phaseTypeEM, beta);
//            }
//
//            if(i % printLikehood == 0) {
//                logLh = logLikelihoodPH(observations, censored, phaseTypeEM, transformation, beta);
//                relError = abs((logLh - logLh_ant) / logLh_ant);
//                logLh_ant = logLh;
//                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". mu: " << beta[0] << ". sigma: " << beta[1] << '\n';
//            }
//        }
//    }
    else if (transformation == 5) { // Gumbel
        std::vector<double> fxnprime;
            
        std::function<std::vector<double>(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
            
        derivative = DerivativeGumbel;
            
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
            
        for (int i{1}; i <= stepsEM; ++i) {
            transformData(observations, TransformObs, transformation, beta);
            transformData(censored, TransformCensored, transformation, beta);
                
            EMstep(TransformObs, TransformCensored, phaseTypeEM);
                

            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
                
            while (norm(fxnprime[0], fxnprime[1]) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime[0];
                beta[1] = beta[1] + lambda * fxnprime[1];
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
            
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH(observations, censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". mu: " << beta[0] << ". sigma: " << beta[1] << '\n';
            }
        }
    }
    else if (transformation == 6) { // GEVD
        std::vector<double> fxnprime;
            
        std::function<std::vector<double>(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
            
        derivative = DerivativeGEVD;
            
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
            
        for (int i{1}; i <= stepsEM; ++i) {
            transformData(observations, TransformObs, transformation, beta);
            transformData(censored, TransformCensored, transformation, beta);
                
            EMstep(TransformObs, TransformCensored, phaseTypeEM);
                

            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
            
            while (norm(fxnprime[0], fxnprime[1], fxnprime[2]) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime[0];
                beta[1] = beta[1] + lambda * fxnprime[1];
                beta[2] = beta[2] + lambda * fxnprime[2];
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
            
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH(observations, censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". mu: " << beta[0] << ". sigma: " << beta[1]  << ". xi: " << beta[2] << ". mu - sigma / xi: " << beta[0] - beta[1] / beta[2] << '\n';
            }
        }
    }
    
    
}


void EMIterate_RK(int stepsEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, int printLikehood, int transformation, std::vector<double> & beta, double lambda, double epsilon) {
    
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    int typeOfStep{askTypeOfStepRK()};
          
    double RKstep{0.0};
    RKstep = (typeOfStep == 2) ? askStepLengthRK() : defaultStepLength(phaseTypeEM.getT());
          
    logLh = logLikelihoodPH_RK(RKstep, observations, censored, phaseTypeEM, transformation, beta);
    logLh_ant = logLh;
    
    
    if (transformation == 0) {
        for (int i{1}; i <= stepsEM; ++i) {
            if (typeOfStep == 1) {
                RKstep = defaultStepLength(phaseTypeEM.getT());
            }
            EMstepRungeKutta(RKstep, phaseTypeEM, observations, censored);
                  
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH_RK(RKstep, observations , censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    else if (transformation == 1) {
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
        
        transformData(observations, TransformObs, transformation, beta);
        transformData(censored, TransformCensored, transformation, beta);
        
        for (int i{1}; i <= stepsEM; ++i) {
            if (typeOfStep == 1) {
                RKstep = defaultStepLength(phaseTypeEM.getT());
            }
            EMstepRungeKutta(RKstep, phaseTypeEM, TransformObs, TransformCensored);
                  
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH_RK(RKstep, observations , censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    else if (transformation == 2 || transformation == 3 || transformation == 4) { //Pareto, Weibull and Gompertz
        double fxnprime{};
        
        std::function<double(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
        
        if (transformation == 2) {
            derivative = derivativeMatrixPareto;
        }
        else if (transformation == 3)  {
            derivative = derivativeMatrixWeibull;
        }
        else if (transformation == 4)  {
            derivative = derivativeMatrixGompertz;
        }
        
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
        
        for (int i{1}; i <= stepsEM; ++i) {
            if (typeOfStep == 1) {
                RKstep = defaultStepLength(phaseTypeEM.getT());
            }
            
            transformData(observations, TransformObs, transformation, beta);
            transformData(censored, TransformCensored, transformation, beta);
            
            EMstepRungeKutta(RKstep, phaseTypeEM, TransformObs, TransformCensored);
                
            //Gradient ascent
            // Stoping condition relative error lest than epsilon
//            double beta_new{};
//            while (relError > epsilon) {
//                fxnprime = FnForNR(observations, phaseTypeEM, beta);
//                beta_new = beta + lambda * fxnprime;
//                relError = abs( (beta_new - beta) / beta );
//                beta = beta_new;
//                //std::cout<< fxnprime << '\t' << beta << '\n';
//            }
            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
            while (abs(fxnprime) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime;
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
            
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH_RK(RKstep, observations , censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". Beta: " << beta[0] << '\n';
            }
        }
    }
    else if (transformation == 5) { // Gumbel
        std::vector<double> fxnprime;
                
        std::function<std::vector<double>(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
                
        derivative = DerivativeGumbel;
            
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
        
        for (int i{1}; i <= stepsEM; ++i) {
            if (typeOfStep == 1) {
                RKstep = defaultStepLength(phaseTypeEM.getT());
            }
            
            reversTransformData(observations, TransformObs, transformation, beta);
            reversTransformData(censored, TransformCensored, transformation, beta);
                
            EMstepRungeKutta(RKstep, phaseTypeEM, TransformObs, TransformCensored);
            

            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
                    
            while (norm(fxnprime[0], fxnprime[1]) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime[0];
                beta[1] = beta[1] + lambda * fxnprime[1];
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
            
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH_RK(RKstep, observations , censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". mu: " << beta[0] << ". sigma: " << beta[1] << '\n';
            }
        }
    }
    else if (transformation == 6) { // GEVD
        std::vector<double> fxnprime;
                
        std::function<std::vector<double>(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
                
        derivative = DerivativeGEVD;
            
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
        
        for (int i{1}; i <= stepsEM; ++i) {
            if (typeOfStep == 1) {
                RKstep = defaultStepLength(phaseTypeEM.getT());
            }
            
            reversTransformData(observations, TransformObs, transformation, beta);
            reversTransformData(censored, TransformCensored, transformation, beta);
                
            EMstepRungeKutta(RKstep, phaseTypeEM, TransformObs, TransformCensored);
            

            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
            
            while (norm(fxnprime[0], fxnprime[1], fxnprime[2]) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime[0];
                beta[1] = beta[1] + lambda * fxnprime[1];
                beta[2] = beta[2] + lambda * fxnprime[2];
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
            if(i % printLikehood == 0) {
                logLh = logLikelihoodPH_RK(RKstep, observations , censored, phaseTypeEM, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". mu: " << beta[0] << ". sigma: " << beta[1]  << ". xi: " << beta[2] << ". mu - sigma / xi: " << beta[0] - beta[1] / beta[2] << '\n';
            }
        }
    }
}


void EMIterate_Uni(int stepsEM, double epsilon, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, int printLikehood, int transformation, std::vector<double> & beta, double lambda, double epsilonGD) {
    
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    logLh = logLikelihoodPH_Uni(observations, censored, phaseTypeEM, epsilon, transformation, beta);
    logLh_ant = logLh;
    
    if (transformation == 0) {
        for (int i{1}; i <= stepsEM; ++i) {
            EMstepUniformization(observations, censored, phaseTypeEM, epsilon);
            if (i % printLikehood == 0) {
                logLh = logLikelihoodPH_Uni(observations, censored, phaseTypeEM, epsilon, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    else if (transformation == 1) {
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
    
        transformData(observations, TransformObs, transformation, beta);
        transformData(censored, TransformCensored, transformation, beta);
        
        for (int i{1}; i <= stepsEM; ++i) {
            EMstepUniformization(TransformObs, TransformCensored, phaseTypeEM, epsilon);
            if (i % printLikehood == 0) {
                logLh = logLikelihoodPH_Uni(observations, censored, phaseTypeEM, epsilon, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    else if (transformation == 2 || transformation == 3 || transformation == 4) {
        double fxnprime{};
        
        std::function<double(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
        
        if (transformation == 2) {
            derivative = derivativeMatrixPareto;
        }
        else if (transformation == 3)  {
            derivative = derivativeMatrixWeibull;
        }
        else if (transformation == 4)  {
            derivative = derivativeMatrixGompertz;
        }
        
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
        
        for (int i{1}; i <= stepsEM; ++i) {
            
            transformData(observations, TransformObs, transformation, beta);
            transformData(censored, TransformCensored, transformation, beta);
            
            EMstepUniformization(TransformObs, TransformCensored, phaseTypeEM, epsilon);
            
            //Gradient ascent
            // Stoping condition relative error lest than epsilon
//            double beta_new{};
//            while (relError > epsilon) {
//                fxnprime = FnForNR(observations, phaseTypeEM, beta);
//                beta_new = beta + lambda * fxnprime;
//                relError = abs( (beta_new - beta) / beta );
//                beta = beta_new;
//                //std::cout<< fxnprime << '\t' << beta << '\n';
//            }
            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
            while (abs(fxnprime) > epsilonGD) {
                beta[0] = beta[0] + lambda * fxnprime;
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
            
            
            if (i % printLikehood == 0) {
                logLh = logLikelihoodPH_Uni(observations, censored, phaseTypeEM, epsilon, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". Beta: " << beta[0] << '\n';
            }
        }
    }
    else if (transformation == 5) { // Gumbel
        std::vector<double> fxnprime;
                    
        std::function<std::vector<double>(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
            
        derivative = DerivativeGumbel;
            
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
            
        for (int i{1}; i <= stepsEM; ++i) {
                
            transformData(observations, TransformObs, transformation, beta);
            transformData(censored, TransformCensored, transformation, beta);
                
            EMstepUniformization(TransformObs, TransformCensored, phaseTypeEM, epsilon);
                
            //Gradient ascent
            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
                        
            while (norm(fxnprime[0], fxnprime[1]) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime[0];
                beta[1] = beta[1] + lambda * fxnprime[1];
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
                
            if (i % printLikehood == 0) {
                logLh = logLikelihoodPH_Uni(observations, censored, phaseTypeEM, epsilon, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". mu: " << beta[0] << ". sigma: " << beta[1] << '\n';
            }
        }
    }
    else if (transformation == 6) { // GEVD
        std::vector<double> fxnprime;
                    
        std::function<std::vector<double>(const std::vector<Sample> &, PhaseType &, const std::vector<double> &)> derivative;
        
        derivative = DerivativeGEVD;
        
        std::vector<Sample> TransformObs{observations};
        std::vector<Sample> TransformCensored{censored};
            
        for (int i{1}; i <= stepsEM; ++i) {
                
            transformData(observations, TransformObs, transformation, beta);
            transformData(censored, TransformCensored, transformation, beta);
                
            EMstepUniformization(TransformObs, TransformCensored, phaseTypeEM, epsilon);
                
            //Gradient ascent
            // Stoping condition derivarite close to zero
            fxnprime = derivative(observations, phaseTypeEM, beta);
            
            while (norm(fxnprime[0], fxnprime[1], fxnprime[2]) > epsilon) {
                beta[0] = beta[0] + lambda * fxnprime[0];
                beta[1] = beta[1] + lambda * fxnprime[1];
                beta[2] = beta[2] + lambda * fxnprime[2];
                fxnprime = derivative(observations, phaseTypeEM, beta);
            }
                
            if (i % printLikehood == 0) {
                logLh = logLikelihoodPH_Uni(observations, censored, phaseTypeEM, epsilon, transformation, beta);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". mu: " << beta[0] << ". sigma: " << beta[1]  << ". xi: " << beta[2] << ". mu - sigma / xi: " << beta[0] - beta[1] / beta[2] << '\n';
            }
        }
    }
}


// EM using uniformization

//  Find n such that P(N > n) = epsilon with N Poisson distributed
int findN(double epsilon, double lambda) {
    int n{0};
    double cumProb{0.0};
    
    do {
        cumProb += dPoisson(n, lambda);
        ++n;
    } while (cumProb < 1.0 - epsilon);
    
    return (n - 1);
}

// Computes the elements T^n / n! till the value size
void vectorOfMatrices(std::vector<Numeric_lib::Matrix<double,2>> & theVector, const Numeric_lib::Matrix<double,2> & T, double alpha, int size) {
    
    Numeric_lib::Matrix<double,2> P(addMatrices(identityMatrix(static_cast<int>(T.dim1())), T / alpha));
    
    
    theVector.push_back(identityMatrix(static_cast<int>(T.dim1())));
    
    for (int k{1}; k <= size; ++k) {
        theVector.push_back( matrixProduct(P * (1.0 / k), theVector[k - 1]));
    }
    
}

// Computes (T)^(2^n)
void pow2Matrix(int n , Numeric_lib::Matrix<double,2> & T) {
    Numeric_lib::Matrix<double,2> auxMat(T.dim1(),T.dim2());
    
    for (int i{1}; i <= n; ++i) {
        auxMat = matrixProduct(T, T);
        T = auxMat;
    }
}

// Computes e^(Tx) base on the values on powerVector
void matrixExpSum(double x, int n, const std::vector<Numeric_lib::Matrix<double,2>> & powerVector, double alpha, Numeric_lib::Matrix<double,2> & resultmatrix) {
    
    resultmatrix = powerVector[0];
    
    for (int i{1}; i <= n; ++i) {
        resultmatrix = addMatrices(resultmatrix, powerVector[i] * exp(i * log(alpha * x)));
    }
    resultmatrix *= exp(-alpha * x);
}


void EMstepUniformization(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, double epsilon) {
    long p{PH.getp()};
    Numeric_lib::Matrix<double,2> pi(PH.getPi());
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    Numeric_lib::Matrix<double,2> e(PH.gete());
    
    Numeric_lib::Matrix<double,2> Bmean(p,1);
    Numeric_lib::Matrix<double,2> Zmean(p,1);
    Numeric_lib::Matrix<double,2> Nmean(p,p + 1);
    
    Numeric_lib::Matrix<double,2> avector(1,p);
    Numeric_lib::Matrix<double,2> bvector(p,1);
    Numeric_lib::Matrix<double,2> cmatrix(p,p);
    Numeric_lib::Matrix<double,2> aux_exp(p,p);
    
    Numeric_lib::Matrix<double,2> J(2 * p,2 * p);
    
    Numeric_lib::Matrix<double,2> tProductPi(p,p);
    matrixProduct(t, pi, tProductPi);
    
    matrixForVanLoan(T, T, tProductPi, J);
    
    
    double alpha{0.0};
    alpha = matrixMaxDiagonal(J * (-1.0));
    
    
    int N{findN(epsilon, 1)};
    
    std::vector<Numeric_lib::Matrix<double,2>> theVector;
    
    vectorOfMatrices(theVector, J, alpha, N);
    
    double SumOfWeights{0.0};
    double density{0.0};
    
    // E-step
    //  Unccensored data
    for (int k{0}; k < observations.size(); ++k) {
        SumOfWeights += observations[k].weight;
        
        double x{observations[k].obs};
        
        if (x * alpha <= 1.0) {
            matrixExpSum(x, N, theVector, alpha, J);
        }
        else {
            int n{};
            n = log(alpha * x) / log(2);
            ++n;
            
            matrixExpSum(x / pow(2.0, n), N, theVector, alpha, J);
            
            pow2Matrix(n, J);
        }
        
        
        for (int i{0}; i < p; ++i) {
            for (int j{0}; j < p; ++j) {
                aux_exp(i,j) = J(i,j);
                cmatrix(i,j) = J(i,j + p);
            }
        }
        
        matrixProduct(pi, aux_exp, avector);
        matrixProduct(aux_exp, t, bvector);
        density=matrixProduct(pi, bvector)(0,0);
        
        //E-step
        for (int i{0}; i < p; ++i) {
            Bmean(i,0) += pi(0,i) * bvector(i,0) * observations[k].weight / density;
            Nmean(i,p) += avector(0,i) * t(i,0) * observations[k].weight / density;
            Zmean(i,0) += cmatrix(i,i) * observations[k].weight / density;
            for (int j{0}; j < p; ++j) {
                Nmean(i,j) += T(i,j) * cmatrix(j,i) * observations[k].weight / density;
            }
        }
    }
    //  Right-Censored Data - Aun no lo he cambiado 
    double SumOfCensored{0.0};
    if (censored.size() > 0) {
        matrixProduct(e, pi, tProductPi);
        matrixForVanLoan(T, T, tProductPi, J);
        theVector.clear();
        vectorOfMatrices(theVector, J, alpha, N);
    }
    for (int k{0}; k < censored.size(); ++k) {
        SumOfCensored += censored[k].weight;
        
        double x{censored[k].obs};
        
        if (x * alpha <= 1.0) {
            matrixExpSum(x, N, theVector, alpha, J);
        }
        else {
            int n{};
            n = log(alpha * x) / log(2);
            ++n;
            
            matrixExpSum(x / pow(2.0, n), N, theVector, alpha, J);
            
            pow2Matrix(n, J);
        }
        
        
        for (int i{0}; i < p; ++i) {
            for (int j{0}; j < p; ++j) {
                aux_exp(i,j) = J(i,j);
                cmatrix(i,j) = J(i,j + p);
            }
        }
        
        matrixProduct(aux_exp, e, bvector);
        density = matrixProduct(pi, bvector)(0,0);
        
        //E-step
        for (int i{0}; i < p; ++i) {
            Bmean(i,0) += pi(0,i) * bvector(i,0) * censored[k].weight / density;
            Zmean(i,0) += cmatrix(i,i) * censored[k].weight / density;
            for (int j{0}; j < p; ++j) {
                Nmean(i,j) += T(i,j) * cmatrix(j,i) * censored[k].weight / density;
            }
        }
    }
    
    // M step
    for (int i{0}; i < p; ++i) {
        pi(0,i) = Bmean(i,0) / (SumOfWeights + SumOfCensored);
        if (pi(0,i) < 0) {
            pi(0,i) = 0;
        }
        t(i,0) = Nmean(i,p) / Zmean(i,0);
        if (t(i,0) < 0) {
            t(i,0) = 0;
        }
        T(i,i) = -t(i,0);
        for (int j{0}; j < p; ++j) {
            if (i != j) {
                T(i,j) = Nmean(i,j) / Zmean(i,0);
                if (T(i,j) < 0) {
                    T(i,j) = 0;
                }
                T(i,i) -= T(i,j);
            }
        }
    }
    PH.changeAll(pi, T);
}




/* ***************
        MPH
*************** */

// Read univariate data from file and put it in a vector of type "Sample"
long readDataForMPH(std::ifstream & inputFile, std::vector<std::vector<Sample>> & observations,std::vector<Sample> & sumObservations, std::vector<std::vector<Sample>> & censored, std::vector<Sample> & sumCensored, int indCensored, int transformation) {
    
    
    // Determines how many colums the file has
    std::string firstline;
    
    getline(inputFile, firstline);
    int count{1};
    long length = firstline.length();
    for (int i{0}; i < length; ++i) {
        if (isspace(firstline[i]))
            ++count;
    }
    
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);
    
    // The file has more colums for censored data then we adjust appropiately
    long divisor{1};
    if (indCensored == 2) {
        ++divisor;
    }
    long dim{count / divisor};
    
    std::vector<Sample> auxVector{};
    // Initializes the vector of obsevations according to the dimension
    for (int i{0}; i < dim; ++i) {
        observations.push_back(auxVector);
        censored.push_back(auxVector);
    }
    
    
    Sample aux;
    int typeObs{0};
    
    // Read data
    while (!inputFile.eof()) {
        double sumObs{0.0};
        int typeObsSum{0};
        
        for (int i{0}; i < dim; ++i) {
            if (indCensored == 2) {
                inputFile >> typeObs;
                if (typeObs == 1) {
                    typeObsSum = 1;
                }
            }
            
            inputFile >> aux.obs;
            
            if (transformation == 1) {
                aux.obs = log(aux.obs + 1); // Matrix-Pareto
            }
            
            sumObs += aux.obs;
            aux.weight = 1.0;
            
            if (typeObs == 0) {
                observations[i].push_back(aux);
            }
            else {
                censored[i].push_back(aux);
            }
        }
        
        aux.obs = sumObs;
        aux.weight = 1.0;
        
        if (typeObsSum == 0) {
            sumObservations.push_back(aux);
        }
        else {
            sumCensored.push_back(aux);
        }
    }
    return dim;
}


// Sign of x
int signFn(double x) {
    return (x < 0) ? -1 : (x > 0);
}


// Gives basic information (size, sum of weigths, mean and sd) of a vector of type "Sample" - This only makes sense for unweithed file from sample (for a density you may check the means an SD but correlation does not makes sense becase of the weights (one would be assuming independence if you multiply the weights
void sampleStatisticsMPH(const std::vector<std::vector<Sample>> & observations) {
    int dim{static_cast<int>(observations.size())};
    
    Numeric_lib::Matrix<double,2> meanVec(1,dim);
    Numeric_lib::Matrix<double,2> sdVec(1,dim);
    
    for (int j{0}; j < dim; ++j) {
        double sumOfx{0.0};
        double sumOfxSquared{0.0};
        
        for (int i{0}; i < observations[j].size(); ++i) {
            sumOfx += observations[j][i].obs;
            sumOfxSquared += observations[j][i].obs * observations[j][i].obs;
        }
        
        meanVec(0,j) = sumOfx / observations[j].size();
        sdVec(0,j) = sqrt(sumOfxSquared / (observations[j].size()) - pow(sumOfx / (observations[j].size()), 2.0));
    }
    
    std::cout << "Dimension: " << observations.size() << '\n';
    std::cout << "Sample size: " << observations[0].size() << '\n';
    std::cout << "Sample mean: " << '\n';
    printMatrix(meanVec);
    std::cout << "Sample standard deviation: " << '\n';
    printMatrix(sdVec);
    
    Numeric_lib::Matrix<double,2> correlations(dim,dim);
    Numeric_lib::Matrix<double,2> kendallsTau(dim,dim);
    
    for (int j{0}; j < dim; ++j) {
        correlations(j,j) = 1;
        kendallsTau(j,j) = 1;
        for (int k{j+1}; k < dim; ++k) {
            double sumOfxy{0.0};
            double sumOfSigns{0.0};
            for (int i{0}; i < observations[j].size(); ++i) {
                sumOfxy += observations[j][i].obs * observations[k][i].obs;
                for (int l{i+1}; l < observations[j].size(); ++l) {
                    sumOfSigns += signFn((observations[j][i].obs - observations[j][l].obs) * (observations[k][i].obs - observations[k][l].obs));
                }
            }
            sumOfSigns = (2.0/(observations[j].size() * (observations[j].size() - 1))) * sumOfSigns;
            kendallsTau(j,k) = sumOfSigns;
            kendallsTau(k,j) = kendallsTau(j,k);
            
            sumOfxy = sumOfxy / observations[j].size();
            correlations(j,k) = (sumOfxy - meanVec(0,j) * meanVec(0,k)) / (sdVec(0,j) * sdVec(0,k));
            correlations(k,j) = correlations(j,k);
            
        }
    }
    std::cout << "Correlations: " << '\n';
    printMatrix(correlations);
    
    std::cout << "Kendall's Tau: " << '\n';
    printMatrix(kendallsTau);
    
}


void infoDataMPH(const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored) {
    std::cout << "Dimension: " << observations.size() << '\n';
    std::cout << "Uncensored data.\n";
    std::cout << "    Sample size: " << '\n';
    for (int j{0}; j < observations.size(); ++j) {
        std::cout << '\t' << observations[j].size();
    }
    std::cout << "\nRight-censored data.\n";
    std::cout << "    Sample size: " << '\n';
    for (int j{0}; j < censored.size(); ++j) {
        std::cout << '\t' << censored[j].size();
    }
    std::cout << '\n';
}


// Given a file with the structure of a MPH (0's and 1's)  puts the corresponding values in PiLegal, TLegal and RLegal
void readStructureForMPH(std::ifstream & inputFile, Numeric_lib::Matrix<double,2> & piLegal ,Numeric_lib::Matrix<double,2> & TLegal, Numeric_lib::Matrix<double,2> & RLegal) {
    long p{dimensionsOfFile(inputFile).m_p};
    long d{dimensionsOfFile(inputFile).m_dim};
    
    for (int i{0}; i < p; ++i) {
        inputFile >> piLegal(0,i);
        for (int j{0}; j < p; ++j) {
            inputFile >> TLegal(i,j);
        }
        for (int j{0}; j < d; ++j) {
            inputFile >> RLegal(i,j);
        }
    }
}


void secondEMstep(const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored, MPH & mph) {
    double lowerbound{1.0E-15}; //Makes a reward zero if it is below of this value - Otherwise it will never be zero
    long p = mph.getp();
    long dim = mph.getd();
    
    Numeric_lib::Matrix<double,2> Zmean(p,dim);
    Numeric_lib::Matrix<double,2> Zsum(p,1);
    
    Numeric_lib::Matrix<double,2> R(p,dim);
    
    //E-step
    for (int j{0}; j < dim; ++j) {
        
        std::vector<int> newStates;
        
        PhaseType PHmarginal(mph.marginal(j, newStates));
        
        long pmarginal = PHmarginal.getp();
        
        Numeric_lib::Matrix<double,2> pi(PHmarginal.getPi());
        Numeric_lib::Matrix<double,2> T(PHmarginal.getT());
        Numeric_lib::Matrix<double,2> t(PHmarginal.gett());
        Numeric_lib::Matrix<double,2> e(PHmarginal.gete());
        
        Numeric_lib::Matrix<double,2> bvector(pmarginal, 1);
        Numeric_lib::Matrix<double,2> cmatrix(pmarginal, pmarginal);
        Numeric_lib::Matrix<double,2> aux_exp(pmarginal, pmarginal);
        
        Numeric_lib::Matrix<double,2> J(2 * pmarginal, 2 * pmarginal);
        Numeric_lib::Matrix<double,2> tProductPi(pmarginal, pmarginal);
        matrixProduct(t, pi, tProductPi);
        
        double density{0.0};
        
        //  Unccensored data
        for (int k{0}; k < observations[j].size(); ++k) {
            VanLoanForIntegrals(observations[j][k].obs, T, T, tProductPi, J);
            
            for (int i{0}; i < pmarginal; ++i) {
                for (int j{0}; j < pmarginal; ++j) {
                    aux_exp(i,j) = J(i,j);
                    cmatrix(i,j) = J(i,j + pmarginal);
                }
            }
            
            matrixProduct(aux_exp, t, bvector);
            density = matrixProduct(pi, bvector)(0,0);
            
            //E-step
            for (int i{0}; i < pmarginal; ++i) {
                Zmean(newStates[i],j) += cmatrix(i,i) * observations[j][k].weight / density;
            }
        }
        

        //  Right-Censored Data
        if (censored[j].size() > 0) {
            matrixProduct(e, pi, tProductPi);
        }
        for (int k{0}; k < censored[j].size(); ++k) {
    
            VanLoanForIntegrals(censored[j][k].obs, T, T, tProductPi, J);
    
            for (int i{0}; i < pmarginal; ++i) {
                for (int j{0}; j < pmarginal; ++j) {
                    aux_exp(i,j) = J(i,j);
                    cmatrix(i,j) = J(i,j + pmarginal);
                }
            }
    
            matrixProduct(aux_exp, e, bvector);
            density = matrixProduct(pi, bvector)(0, 0);
    
            //E-step
            for (int i{0}; i < pmarginal; ++i) {
                Zmean(newStates[i],j) += cmatrix(i,i) * censored[j][k].weight / density;
            }
        }
    }
    
    //M-step
    for (int i{0}; i < p; ++i) {
        for (int j{0}; j < dim; ++j) {
            Zsum(i,0) += Zmean(i,j);
        }
    }
    
    for (int j{0}; j < dim; ++j) {
        for (int i{0}; i < p; ++i){
            R(i,j) = Zmean(i,j) / Zsum(i,0);
            if (R(i,j) < lowerbound) {
                R(i,j) = 0;
            }
//            if(std::isnan(R(i,j))){
//                mph.printR();
//                double xa;
//                std::cin >> xa;
//            }
        }
    }
    
    mph.changeR(R);
}



void EMIterateMPH(int stepsFirstEM, int stepsSecondEM, const std::vector<std::vector<Sample>> & observations, const std::vector<Sample> & sumObservations, const std::vector<std::vector<Sample>> & censored, const std::vector<Sample> & sumCensored, MPH & mph, PhaseType & PHsumOfMarginals, int printLikehood) {
    
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    std::vector<double> dummyVector{};
    
    //First EM - Fit the sum
    logLh = logLikelihoodPH(sumObservations, sumCensored, PHsumOfMarginals, 0, dummyVector);
    logLh_ant = logLh;
    for (int i{1}; i <= stepsFirstEM; ++i) {
        EMstep(sumObservations, sumCensored, PHsumOfMarginals);
        if (i % printLikehood == 0) {
            logLh = logLikelihoodPH(sumObservations, sumCensored, PHsumOfMarginals, 0, dummyVector);
            relError = abs((logLh - logLh_ant) / logLh_ant);
            logLh_ant = logLh;
            std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
        }
    }
    mph.changePi(PHsumOfMarginals.getPi());
    mph.changeT(PHsumOfMarginals.getT());
    
    //Second EM
    for (int i{1}; i <= stepsSecondEM; ++i) {
        secondEMstep(observations, censored, mph);
        if (i % printLikehood == 0) {
            std::cout << "Second EM iteration " << i << '\n';
        }
    }
}


void secondEMstepUniformization(const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored, MPH & mph, double epsilon) {
    double lowerbound{1.0E-10}; //Makes a reward zero if it is below of this value - Otherwise it will never be zero
    long p = mph.getp();
    long dim = mph.getd();
    
    Numeric_lib::Matrix<double,2> Zmean(p,dim);
    Numeric_lib::Matrix<double,2> Zsum(p,1);
    
    Numeric_lib::Matrix<double,2> R(p,dim);
    
    int N{findN(epsilon, 1)};
    std::vector<Numeric_lib::Matrix<double,2>> theVector;
    
    //E-step
    //  Unccensored data
    for (int j{0}; j < dim; ++j) {
        
        std::vector<int> newStates;
        
        PhaseType PHmarginal(mph.marginal(j, newStates));
        
        long pmarginal = PHmarginal.getp();
        
        Numeric_lib::Matrix<double,2> pi(PHmarginal.getPi());
        Numeric_lib::Matrix<double,2> T(PHmarginal.getT());
        Numeric_lib::Matrix<double,2> t(PHmarginal.gett());
        Numeric_lib::Matrix<double,2> e(PHmarginal.gete());
        
        Numeric_lib::Matrix<double,2> bvector(pmarginal,1);
        Numeric_lib::Matrix<double,2> cmatrix(pmarginal,pmarginal);
        Numeric_lib::Matrix<double,2> aux_exp(pmarginal,pmarginal);
        
        Numeric_lib::Matrix<double,2> J(2 * pmarginal,2 * pmarginal);
        Numeric_lib::Matrix<double,2> tProductPi(pmarginal,pmarginal);
        matrixProduct(t, pi, tProductPi);
        
        matrixForVanLoan(T, T, tProductPi, J);
        
        double alpha{0.0};
        alpha = matrixMaxDiagonal(J * (-1.0));
        
        theVector.clear();
        
        vectorOfMatrices(theVector, J, alpha, N);
        
        
        double density{0.0};
        
        for (int k{0}; k < observations[j].size(); ++k) {
            
            double x{observations[j][k].obs};
                   
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, J);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                       
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, J);
                       
                pow2Matrix(n, J);
            }
            
            
            for (int i{0}; i < pmarginal; ++i) {
                for (int j{0}; j < pmarginal; ++j) {
                    aux_exp(i,j) = J(i,j);
                    cmatrix(i,j) = J(i,j + pmarginal);
                }
            }
            
            matrixProduct(aux_exp, t, bvector);
            density = matrixProduct(pi, bvector)(0,0);
            
            //E-step
            for (int i{0}; i < pmarginal; ++i) {
                Zmean(newStates[i],j) += cmatrix(i,i) * observations[j][k].weight / density;
            }
        }
        

        //  Right-Censored Data
        if (censored[j].size() > 0) {
            matrixProduct(e, pi, tProductPi);
            matrixForVanLoan(T, T, tProductPi, J);
            theVector.clear();
            vectorOfMatrices(theVector, J, alpha, N);
        }
        for (int k{0}; k < censored[j].size(); ++k) {
    
            double x{censored[j][k].obs};
            
            if (x * alpha <= 1.0) {
                matrixExpSum(x, N, theVector, alpha, J);
            }
            else {
                int n{};
                n = log(alpha * x) / log(2);
                ++n;
                
                matrixExpSum(x / pow(2.0, n), N, theVector, alpha, J);
                
                pow2Matrix(n, J);
            }
    
            for (int i{0}; i < pmarginal; ++i) {
                for (int j{0}; j < pmarginal; ++j) {
                    aux_exp(i,j) = J(i,j);
                    cmatrix(i,j) = J(i,j + pmarginal);
                }
            }
    
            matrixProduct(aux_exp, e, bvector);
            density = matrixProduct(pi, bvector)(0,0);
    
            //E-step
            for (int i{0}; i < pmarginal; ++i) {
                Zmean(newStates[i],j) += cmatrix(i,i) * censored[j][k].weight / density;
            }
        }
    }
    
    //M-step
    for (int i{0}; i < p; ++i) {
        for (int j{0}; j < dim; ++j) {
            Zsum(i,0) += Zmean(i,j);
        }
    }
    
    for (int j{0}; j < dim; ++j) {
        for (int i{0}; i < p; ++i) {
            R(i,j) = Zmean(i,j) / Zsum(i,0);
            if (R(i,j) < lowerbound) {
                R(i,j) = 0;
            }
        }
    }
    
    mph.changeR(R);
}

void EMIterateMPH_Uni(int stepsFirstEM, int stepsSecondEM, double epsilon, const std::vector<std::vector<Sample>> & observations, const std::vector<Sample> & sumObservations, const std::vector<std::vector<Sample>> & censored, const std::vector<Sample> & sumCensored, MPH & mph, PhaseType & PHsumOfMarginals, int printLikehood) {
    
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    std::vector<double> dummyVector{};
    //First EM - Fit the sum
    logLh = logLikelihoodPH_Uni(sumObservations, sumCensored, PHsumOfMarginals, epsilon, 0, dummyVector);
    logLh_ant = logLh;
    for (int i{1}; i <= stepsFirstEM; ++i) {
        EMstepUniformization(sumObservations, sumCensored, PHsumOfMarginals, epsilon);
        if (i % printLikehood == 0) {
            logLh = logLikelihoodPH_Uni(sumObservations, sumCensored, PHsumOfMarginals, epsilon, 0, dummyVector);
            relError = abs((logLh - logLh_ant) / logLh_ant);
            logLh_ant = logLh;
            std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
        }
    }
    mph.changePi(PHsumOfMarginals.getPi());
    mph.changeT(PHsumOfMarginals.getT());
    
    
    //Second EM
    for (int i{1}; i <= stepsSecondEM; ++i) {
        secondEMstepUniformization(observations, censored, mph, epsilon);
        if (i % printLikehood == 0) {
            std::cout << "Second EM iteration " << i << '\n';
        }
    }
}




void secondEMstepRungeKutta(int typeOfStep, double h, const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored, MPH & mph) {
    double lowerbound{1.0E-10}; //Makes a reward zero if it is below of this value - Otherwise it will never be zero
    long p = mph.getp();
    long dim = mph.getd();
    
    Numeric_lib::Matrix<double,2> Zmean(p,dim);
    Numeric_lib::Matrix<double,2> Zsum(p,1);
    
    Numeric_lib::Matrix<double,2> R(p,dim);
    
    //E-step
    for (int j{0}; j < dim; ++j) {
        
        std::vector<int> newStates;
        
        PhaseType PHmarginal(mph.marginal(j, newStates));
        
        long pmarginal = PHmarginal.getp();
        
        Numeric_lib::Matrix<double,2> pi(PHmarginal.getPi());
        Numeric_lib::Matrix<double,2> T(PHmarginal.getT());
        Numeric_lib::Matrix<double,2> t(PHmarginal.gett());
        Numeric_lib::Matrix<double,2> e(PHmarginal.gete());
        
        Numeric_lib::Matrix<double,2> avector(1,pmarginal);
        Numeric_lib::Matrix<double,2> bvector(pmarginal,1);
        Numeric_lib::Matrix<double,2> cmatrix(pmarginal,pmarginal);
        
        if (typeOfStep == 1) {
            h = defaultStepLength(T);
        }
        
        //  Unccensored data
        
        // initial conditions
        avector = pi;
        bvector = t;
        
        double dt{0.0};
        if (observations[j].size() > 0) {
            dt = observations[j][0].obs;
        }
        
        double density{0.0};
        
        for (int k{0}; k < observations[j].size(); ++k) {
            
            rungeKutta(avector, bvector, cmatrix, dt, h, T, t);
            
            density = matrixProduct(pi, bvector)(0,0);
            
            //E-step
            for (int i{0}; i < pmarginal; ++i) {
                Zmean(newStates[i],j) += cmatrix(i,i) * observations[j][k].weight / density;
            }
            
            dt = observations[j][k + 1].obs - observations[j][k].obs;
        }
        

        //  Right-Censored Data
        if (censored[j].size() > 0) {
            dt = censored[j][0].obs;
            cmatrix *= 0;
            avector = pi;
            bvector = e;
        }
        for (int k{0}; k < censored[j].size(); ++k) {
    
            rungeKutta(avector, bvector, cmatrix, dt, h, T, e);
            
            density = matrixProduct(pi, bvector)(0,0);
    
            //E-step
            for (int i{0}; i < pmarginal; ++i) {
                Zmean(newStates[i],j) += cmatrix(i,i) * censored[j][k].weight / density;
            }
            
            dt = censored[j][k + 1].obs - censored[j][k].obs;
        }
    }
    
    //M-step
    for (int i{0}; i < p; ++i) {
        for (int j{0}; j < dim; ++j) {
            Zsum(i,0) += Zmean(i,j);
        }
    }
    
    for (int j{0}; j < dim; ++j) {
        for (int i{0}; i < p; ++i) {
            R(i,j) = Zmean(i,j) / Zsum(i,0);
            if (R(i,j) < lowerbound) {
                R(i,j) = 0;
            }
        }
    }
    
    mph.changeR(R);
}



void EMIterateMPH_RK(int stepsFirstEM, int stepsSecondEM, const std::vector<std::vector<Sample>> & observations, const std::vector<Sample> & sumObservations, const std::vector<std::vector<Sample>> & censored, const std::vector<Sample> & sumCensored, MPH & mph, PhaseType & PHsumOfMarginals, int printLikehood) {
    
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    int typeOfStep{askTypeOfStepRK()};

    double RKstep{0.0};
    
    RKstep = (typeOfStep == 2) ? askStepLengthRK() : defaultStepLength(PHsumOfMarginals.getT());
    
    std::vector<double> dummyVector{};
    //First EM - Fit the sum
    logLh = logLikelihoodPH_RK(RKstep, sumObservations, sumCensored, PHsumOfMarginals, 0, dummyVector);
    logLh_ant = logLh;
    for (int i{1}; i <= stepsFirstEM; ++i) {
        if (typeOfStep == 1) {
            RKstep = defaultStepLength(PHsumOfMarginals.getT());
        }
        EMstepRungeKutta(RKstep, PHsumOfMarginals, sumObservations, sumCensored);

        if (i % printLikehood == 0) {
            logLh = logLikelihoodPH_RK(RKstep, sumObservations, sumCensored, PHsumOfMarginals, 0, dummyVector);
            relError = abs((logLh - logLh_ant) / logLh_ant);
            logLh_ant = logLh;
            std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
        }
    }
    mph.changePi(PHsumOfMarginals.getPi());
    mph.changeT(PHsumOfMarginals.getT());
    
    
    //Second EM
    for (int i{1}; i <= stepsSecondEM; ++i) {
        secondEMstepRungeKutta(typeOfStep, RKstep, observations, censored, mph);
        if (i % printLikehood == 0) {
            std::cout << "Second EM iteration " << i << '\n';
        }
    }
}



/*
 Bivariate
 */


// Read univariate data from file and put it in a vector of type "Sample"
void readDataForBivPH(std::ifstream & inputFile, std::vector<BivariateSample> & observations, int transformation) {
    BivariateSample auxSample;
    
    while (!inputFile.eof()) {
        inputFile >> auxSample.x1;
        inputFile >> auxSample.x2;
        
//        if (transformation == 1) {  // Matrix-Pareto
//            auxSample.x1 = log(auxSample.x1 + 1);
//            auxSample.x2 = log(auxSample.x2 + 1);
//        }
        
        auxSample.weight = 1.0;
        observations.push_back(auxSample);
    }
}


void sampleStatisticsBivPH(const std::vector<BivariateSample> & observations) {
    
    Numeric_lib::Matrix<double,2> meanVec(1,2);
    Numeric_lib::Matrix<double,2> sdVec(1,2);
    Numeric_lib::Matrix<double,2> correlations(2,2);
    
    double sumOfx1{0.0};
    double sumOfx2{0.0};
    double sumOfx1Squared{0.0};
    double sumOfx2Squared{0.0};
    double sumOfxy{0.0};
    for (int i{0}; i < observations.size(); ++i) {
        sumOfx1 += observations[i].x1;
        sumOfx2 += observations[i].x2;
        sumOfx1Squared += observations[i].x1 * observations[i].x1;
        sumOfx2Squared += observations[i].x2 * observations[i].x2;
        sumOfxy += observations[i].x1 * observations[i].x2;
    }
        
    meanVec(0,0) = sumOfx1 / observations.size();
    meanVec(0,1) = sumOfx2 / observations.size();
    sdVec(0,0) = sqrt(sumOfx1Squared / (observations.size()) - pow(sumOfx1 / (observations.size()), 2.0));
    sdVec(0,1) = sqrt(sumOfx2Squared / (observations.size()) - pow(sumOfx2 / (observations.size()), 2.0));
    
    std::cout << "Dimension: " << 2 << '\n';
    std::cout << "Sample size: " << observations.size() << '\n';
    std::cout << "Sample mean: " << '\n';
    printMatrix(meanVec);
    std::cout << "Sample standard deviation: " << '\n';
    printMatrix(sdVec);
    
    correlations(0,0) = 1;
    correlations(1,1) = 1;
        
    sumOfxy = sumOfxy / observations.size();
    correlations(0,1) = (sumOfxy - meanVec(0,0) * meanVec(0,1)) / (sdVec(0,0) * sdVec(0,1));
    correlations(1,0) = correlations(0,1);
    
    std::cout << "Correlations: " << '\n';
    printMatrix(correlations);
    
}




DimensionsMPH dimensionsOfBivFile(std::ifstream & inputFile) {
    DimensionsMPH sizeFile;
    
    sizeFile.m_p = std::count(std::istreambuf_iterator<char>(inputFile),
                std::istreambuf_iterator<char>(), '\n') + 1;
    
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);
    
    
    double aux{0.0};
    long p2{0};
    for (int i{0}; i < sizeFile.m_p; ++i) {
        for (int j{0}; j < 1 + sizeFile.m_p + 2; ++j) {
            if (j < sizeFile.m_p + 2) {
                inputFile >> aux;
            }
            else {
                inputFile >> aux;
                p2 += static_cast<long>(aux);
            }
        }
    }
    long p1{sizeFile.m_p - p2};
    
    sizeFile.m_dim = p2;
    sizeFile.m_p = p1;
    
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);
    
    return sizeFile;
}


void readStructureForBivPH(std::ifstream & inputFile, Numeric_lib::Matrix<double,2> & alphaLegal ,Numeric_lib::Matrix<double,2> & T11Legal, Numeric_lib::Matrix<double,2> & T12Legal, Numeric_lib::Matrix<double,2> & T22Legal) {
    long p1{T11Legal.dim1()};
    long p2{T22Legal.dim1()};
    
    double helpvar{};
    
    for (int i{0}; i < p1; i++) {
        inputFile >> alphaLegal(0,i);
        for (int j{0}; j < p1; j++) {
            inputFile >> T11Legal(i,j);
        }
        for (int j{0}; j < p2; j++) {
            inputFile >> T12Legal(i,j);
        }
        inputFile >> helpvar; // Avoids reading R
        inputFile >> helpvar;
    }
    for (int i{0}; i < p2; i++) {
        inputFile >> helpvar;
        if (helpvar != 0) {
            std::cout << "Warning: the extructure of alpha is not correct";
        }
        for (int j{0}; j < p1; j++) {
            inputFile >> helpvar;
            if(helpvar != 0){
                std::cout << "Warning: the extructure of T is not correct";
            }
        }
        for (int j{0}; j < p2; j++) {
            inputFile >> T22Legal(i,j);
        }
        inputFile >> helpvar;
        inputFile >> helpvar;
    }
}



//  EM for a MPH fit
void EMstepForBiv(const std::vector<BivariateSample> & observations, BivariatePH & bivPH) {
    long p1{bivPH.getp1()};
    long p2{bivPH.getp2()};
    long p{p1 + p2};
    
    double density{0};
    double aux;
    
    Numeric_lib::Matrix<double,2> alpha(bivPH.getAlpha());
    Numeric_lib::Matrix<double,2> T11(bivPH.getT11());
    Numeric_lib::Matrix<double,2> T12(bivPH.getT12());
    Numeric_lib::Matrix<double,2> T22(bivPH.getT22());
    
    Numeric_lib::Matrix<double,2> Bmean(p1,1);
    Numeric_lib::Matrix<double,2> Zmean(p,1);
    Numeric_lib::Matrix<double,2> Nmean(p,p + 1);
    
    Numeric_lib::Matrix<double,2> bmatrix1(p1,p1);
    Numeric_lib::Matrix<double,2> bmatrix2(p2,p2);
    Numeric_lib::Matrix<double,2> cmatrix1(p1,p1);
    Numeric_lib::Matrix<double,2> cmatrix2(p2,p2);
    Numeric_lib::Matrix<double,2> aux_exp1(p1,p1);
    Numeric_lib::Matrix<double,2> aux_exp2(p2,p2);
    
    Numeric_lib::Matrix<double,2> J1(2 * p1,2 * p1);
    Numeric_lib::Matrix<double,2> J2(2 * p2,2 * p2);
    
    Numeric_lib::Matrix<double,2> exitvec(bivPH.gett());
    
    Numeric_lib::Matrix<double,2> auxMatrix1(p1,1);
    Numeric_lib::Matrix<double,2> auxMatrix2(p2,p1);
    Numeric_lib::Matrix<double,2> auxMatrix3(1,p2);
    
    double sumOfWeights{0.0};
    //E step
    for (int k{0}; k < observations.size(); ++k) {
        
        sumOfWeights += observations[k].weight;
        
        matrixExponential(T22 * observations[k].x2, aux_exp2);
        matrixProduct(T12, matrixProduct(aux_exp2, matrixProduct(exitvec, alpha)), bmatrix1);
        
        VanLoanForIntegrals(observations[k].x1, T11, T11, bmatrix1, J1);
        
        for (int i{0}; i < p1; ++i) {
            for (int j{0}; j < p1; ++j) {
                aux_exp1(i,j) = J1(i,j);
                cmatrix1(i,j) = J1(i,j + p1);
            }
        }
        
        matrixProduct(exitvec, matrixProduct(alpha, matrixProduct(aux_exp1, T12)), bmatrix2);
        
        VanLoanForIntegrals(observations[k].x2, T22, T22, bmatrix2, J2);
        
        for (int i{0}; i < p2; ++i) {
            for (int j{0}; j < p2; ++j) {
                cmatrix2(i,j) = J2(i,j + p2);
            }
        }
        density = matrixProduct(alpha, matrixProduct(aux_exp1, matrixProduct(T12, matrixProduct(aux_exp2, exitvec))))(0,0);
        
        
        //E-step
        matrixProduct(aux_exp1, matrixProduct(T12, matrixProduct(aux_exp2, exitvec)), auxMatrix1);
        matrixProduct(aux_exp2, matrixProduct(exitvec, matrixProduct(alpha, aux_exp1)), auxMatrix2);
        for (int i{0}; i < p1; ++i) {
            aux = auxMatrix1(i,0);
            Bmean(i,0) += alpha(0,i) * aux * observations[k].weight / density;
            Zmean(i,0) += cmatrix1(i,i) * observations[k].weight / density;
            for (int j{0}; j < p1; ++j) {
                Nmean(i,j) += T11(i,j) * cmatrix1(j,i) * observations[k].weight / density;
            }
            for (int j{0}; j < p2; ++j) {
                aux = auxMatrix2(j,i);
                Nmean(i,j + p1) += T12(i,j) * aux * observations[k].weight / density;
            }
        }
        
        matrixProduct(alpha, matrixProduct(aux_exp1, matrixProduct(T12, aux_exp2)), auxMatrix3);
        for (int i{0}; i < p2; ++i) {
            Zmean(i + p1,0) += cmatrix2(i,i) * observations[k].weight / density;
            aux = auxMatrix3(0,i);
            Nmean(i + p1,p) += aux * exitvec(i,0) * observations[k].weight / density;
            for (int j{0}; j < p2; ++j){
                Nmean(i + p1,j + p1) += T22(i,j) * cmatrix2(j,i) * observations[k].weight / density;
            }
        }
    }
    
    // M step
    for (int i{0}; i < p1; ++i) {
        alpha(0,i) = Bmean(i,0) / sumOfWeights;
        if (alpha(0,i) < 0) {
            alpha(0,i) = 0;
        }
        T11(i,i) = 0;
        for (int j{0}; j < p1; ++j) {
            if (i != j) {
                T11(i,j) = Nmean(i,j) / Zmean(i,0);
                if (T11(i,j) < 0) {
                    T11(i,j) = 0;
                }
                T11(i,i) -= T11(i,j);
            }
        }
        for (int j{0}; j < p2; ++j) {
            T12(i,j) = Nmean(i,j + p1) / Zmean(i,0);
            if (T12(i,j) < 0){
                T12(i,j) = 0;
            }
            T11(i,i) -= T12(i,j);
        }
    }
    for (int i{0}; i < p2; ++i) {
        exitvec(i,0) = Nmean(i + p1,p) / Zmean(i + p1,0);
        if (exitvec(i,0) < 0) {
            exitvec(i,0) = 0;
        }
        T22(i,i) = -exitvec(i,0);
        for (int j{0}; j < p2; ++j) {
            if (i != j) {
                T22(i,j) = Nmean(i + p1,j + p1) / Zmean(i + p1,0);
                if (T22(i,j) < 0){
                    T22(i,j) = 0;
                }
                T22(i,i) -= T22(i,j);
            }
        }
    }
    
    bivPH.changeAll(alpha, T11, T12, T22);
}

double logLikelihoodBivPH(const std::vector<BivariateSample> & observations, BivariatePH & bivPH, int transformation, double beta1, double beta2) {
    double logLh{0.0};
    if (transformation == 0) {
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * log(bivPH.density(observations[k].x1, observations[k].x2));
        }
    }
    else if(transformation == 1 || transformation==2) {
        if (transformation == 1) {
            beta1 = 1;
            beta2 = 1;
        }
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * ( log(bivPH.density(log(observations[k].x1 / beta1 + 1), log(observations[k].x2 / beta2 + 1)))  - log(observations[k].x1 + beta1)  - log(observations[k].x2 + beta2));
        }
    }
    else if (transformation == 3) {
        for (int k{0}; k < observations.size(); ++k) {
            logLh += observations[k].weight * ( log(bivPH.density(pow(observations[k].x1, beta1), pow(observations[k].x2, beta2))) + log(beta1) + (beta1 - 1) * log(observations[k].x1) + log(beta2) + (beta2 - 1) * log(observations[k].x2) );
        }
    }
    
    return logLh;
}

double partial1MW(const std::vector<BivariateSample> & observations, BivariatePH & PH, double beta1, double beta2) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( ( matrixProduct(PH.getAlpha(), matrixProduct(matrixExponential(PH.getT11() * pow(observations[k].x1, beta1)), matrixProduct(PH.getT11(), matrixProduct(PH.getT12(), matrixProduct(matrixExponential(PH.getT22() * pow(observations[k].x2, beta2)), matrixProduct(PH.getT22() * (-1.0), PH.gete()))))))(0,0) * pow(observations[k].x1, beta1) * log(observations[k].x1) ) / PH.density(pow(observations[k].x1, beta1), pow(observations[k].x2, beta2)) + 1 / beta1 + log(observations[k].x1) );
    }
    return logLh;
}

double partial2MW(const std::vector<BivariateSample> & observations, BivariatePH & PH, double beta1, double beta2) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( ( matrixProduct(PH.getAlpha(), matrixProduct(matrixExponential(PH.getT11() * pow(observations[k].x1, beta1)), matrixProduct(PH.getT12(), matrixProduct(matrixExponential(PH.getT22() * pow(observations[k].x2, beta2)), matrixProduct(PH.getT22(), matrixProduct(PH.getT22() * (-1.0), PH.gete()))))))(0,0) * pow(observations[k].x2, beta2) * log(observations[k].x2) ) / PH.density(pow(observations[k].x1, beta1), pow(observations[k].x2, beta2)) + 1 / beta2 + log(observations[k].x2) );
    }
    return logLh;
}

double partial1MP(const std::vector<BivariateSample> & observations, BivariatePH & PH, double beta1, double beta2) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( ( matrixProduct(PH.getAlpha(), matrixProduct(matrixExponential(PH.getT11() * log(observations[k].x1 / beta1 + 1)), matrixProduct(PH.getT11(), matrixProduct(PH.getT12(), matrixProduct(matrixExponential(PH.getT22() * log(observations[k].x2 / beta2 + 1)), matrixProduct(PH.getT22() * (-1.0), PH.gete()))))))(0,0) * (- observations[k].x1) / (beta1 * (observations[k].x1 + beta1)) ) / PH.density(log(observations[k].x1 / beta1 + 1), log(observations[k].x2 / beta2 + 1)) - 1 / (observations[k].x1 + beta1) );
    }
    return logLh;
}

double partial2MP(const std::vector<BivariateSample> & observations, BivariatePH & PH, double beta1, double beta2) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( ( matrixProduct(PH.getAlpha(), matrixProduct(matrixExponential(PH.getT11() * log(observations[k].x1 / beta1 + 1)), matrixProduct(PH.getT12(), matrixProduct(matrixExponential(PH.getT22() * log(observations[k].x2 / beta2 + 1)), matrixProduct(PH.getT22(), matrixProduct(PH.getT22() * (-1.0), PH.gete()))))))(0,0)  * (- observations[k].x2) / (beta2 * (observations[k].x2 + beta2)) ) / PH.density(log(observations[k].x1 / beta1 + 1), log(observations[k].x2 / beta2 + 1)) - 1 / (observations[k].x2 + beta2) );
    }
    return logLh;
}


void transformDataBiv(const std::vector<BivariateSample> & observations, std::vector<BivariateSample> & TransformObs, int transformation, double beta1, double beta2) {
    if (transformation == 1 || transformation == 2) {
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].x1 = log(observations[i].x1 / beta1 + 1);
            TransformObs[i].x2 = log(observations[i].x2 / beta2 + 1);
        }
    }
    else if (transformation == 3) {
        for (int i{0}; i < observations.size(); ++i) {
            TransformObs[i].x1 = pow(observations[i].x1, beta1);
            TransformObs[i].x2 = pow(observations[i].x2, beta2);
        }
    }
    
}

void EMIterateBivPH(int stepsEM, const std::vector<BivariateSample> & observations, BivariatePH & bivPH, int printLikehood, int transformation, double & beta1, double  & beta2, double lambda, double epsilon) {
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    logLh = logLikelihoodBivPH(observations, bivPH, transformation, beta1, beta2);
    logLh_ant = logLh;
    if (transformation == 0) {
        for (int i{1}; i <= stepsEM; ++i) {
            EMstepForBiv(observations, bivPH);
            if (i % printLikehood == 0) {
                logLh = logLikelihoodBivPH(observations, bivPH, transformation, beta1, beta2);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    else if (transformation == 1) {
        std::vector<BivariateSample> transformObs{observations};
        transformDataBiv(observations, transformObs, transformation, 1, 1);
        
        for (int i{1}; i <= stepsEM; ++i) {
            EMstepForBiv(transformObs, bivPH);
            if (i % printLikehood == 0) {
                logLh = logLikelihoodBivPH(observations, bivPH, transformation, beta1, beta2);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << '\n';
            }
        }
    }
    else if (transformation == 2 || transformation == 3) {
        double fxnprime1{};
        double fxnprime2{};
        
        std::function<double(const std::vector<BivariateSample> &, BivariatePH &, double, double)> partial1;
        std::function<double(const std::vector<BivariateSample> &, BivariatePH &, double, double)> partial2;
        
        if (transformation == 2) {
            partial1 = partial1MP;
            partial2 = partial2MP;
        }
        else {
            partial1 = partial1MW;
            partial2 = partial2MW;
        }
        
        
        std::vector<BivariateSample> transformObs{observations};
        
        for (int i{1}; i <= stepsEM; ++i) {
            transformDataBiv(observations, transformObs, transformation, beta1, beta2);
            EMstepForBiv(transformObs, bivPH);
            
            // Gradient ascent
            fxnprime1 = partial1(observations, bivPH, beta1, beta2);
            fxnprime2 = partial2(observations, bivPH, beta1, beta2);
            
            while (norm(fxnprime1, fxnprime2) > epsilon) {
                beta1 = beta1 + lambda * fxnprime1;
                beta2 = beta2 + lambda * fxnprime2;
                fxnprime1 = partial1(observations, bivPH, beta1, beta2);
                fxnprime2 = partial2(observations, bivPH, beta1, beta2);
            }
            
            if (i % printLikehood == 0) {
                logLh = logLikelihoodBivPH(observations, bivPH, transformation, beta1, beta2);
                relError = abs((logLh - logLh_ant) / logLh_ant);
                logLh_ant = logLh;
                std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << ". Beta1: " << beta1 << ". Beta2: " << beta2 << '\n';
            }
        }
    }
}


/* Test - Bivariate Weibull */



double logLikelihoodBivMW(const std::vector<BivariateSample> & observations, BivariatePH & bivPH, double beta1, double beta2) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * ( log(bivPH.density(pow(observations[k].x1, beta1), pow(observations[k].x2, beta2))) + log(beta1) + (beta1 - 1) * log(observations[k].x1) + log(beta2) + (beta2 - 1) * log(observations[k].x2) );
    }
    return logLh;
}




void EMIterateBivMW(int stepsEM, const std::vector<BivariateSample> & observations, BivariatePH & bivPH, int printLikehood, double & beta1, double & beta2, double lambda) {
   
    const std::string LogL{"LLBivMW.txt"};
    std::ofstream outFile{LogL};
    
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    
    double fxnprime1{};
    double fxnprime2{};
    double epsilon{0.001};
    
    logLh = logLikelihoodBivMW(observations, bivPH, beta1, beta2);
    logLh_ant = logLh;
    
    std::vector<BivariateSample> TransformObs{observations};
    
    for (int i{1}; i <= stepsEM; ++i) {
        
        transformDataBiv(observations, TransformObs,3, beta1, beta2);
        
        EMstepForBiv(TransformObs, bivPH);
        
        
        fxnprime1 = partial1MW(observations, bivPH, beta1, beta2);
        fxnprime2 = partial2MW(observations, bivPH, beta1, beta2);
        //std::cout<< fxnprime1 << '\t' << beta1 << '\n' << fxnprime2 << '\t' << beta2 << '\n';
        while (norm(fxnprime1, fxnprime2) > epsilon) {
            beta1 = beta1 + lambda * fxnprime1;
            beta2 = beta2 + lambda * fxnprime2;
            fxnprime1 = partial1MW(observations, bivPH, beta1, beta2);
            fxnprime2 = partial2MW(observations, bivPH, beta1, beta2);
            //std::cout<< fxnprime1 << '\t' << beta1 << '\n' << fxnprime2 << '\t' << beta2 << '\n';
        }
        
        if (i % printLikehood == 0) {
            logLh = logLikelihoodBivMW(observations, bivPH, beta1, beta2);
            relError = abs((logLh - logLh_ant) / logLh_ant);
            logLh_ant = logLh;
            std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << " beta 1 " << beta1 << " beta 2 " << beta2 << '\n';
            outFile << logLh << '\t' << beta1 << '\t' << beta2 << '\n';
        }
    }
    outFile.close();
}



/* End test */


void askMultivariateDist(std::vector<std::vector<Sample>> & observations,  std::vector<Sample> & sumObservations) {
    int typeContDist{};
    do {
        std::cout << "\nType of density:\n";
        std::cout << "      1. Bivariate Pareto\n";
        std::cout << "      2. Bivariate M-O Exponential\n";
        std::cout << "      3. Bivariate Normal\n";
        std::cout << "      4. Bivariate Weibull\n";
        std::cout << "Select(1-4): ";
        std::cin >> typeContDist;
        if (typeContDist < 1 || typeContDist > 4) {
            std::cout << "Please enter a valid option\n";
        }
    } while (typeContDist < 1 || typeContDist > 4);
    
    
    MultivariateDistribution *theMulDist_prt;
    
    switch(typeContDist) {
        case 1: {
            std::cout << "\nBivariate Pareto with shape parameter alpha and scale parameters theta1 and theta2:\n";
            std::cout << "alpha:";
            double alpha{};
            std::cin >> alpha;
            std::cout << "theta1:";
            double theta1{};
            std::cin >> theta1;
            std::cout << "theta2:";
            double theta2{};
            std::cin >> theta2;
            theMulDist_prt = new BivariatePareto(alpha, theta1, theta2);
            break;
        }
        case 2: {
            std::cout << "\nMarshall Olkin’s Bivariate Exponential with parameters lambda1, lambda2 and lambda12.\n";
            std::cout << "lambda1:";
            double lambda1{};
            std::cin >> lambda1;
            std::cout << "lambda2:";
            double lambda2{};
            std::cin >> lambda2;
            std::cout << "lambda12:";
            double lambda12{};
            std::cin >> lambda12;
            theMulDist_prt = new MarshallOlkinBivariateExponential(lambda1, lambda2, lambda12);
            break;
        }
        case 3: {
            std::cout << "\n Bivariate Normal with parameters mu1, mu2, sigma1, sigma2 and rho.\n";
            std::cout << "mu1:";
            double mu1{};
            std::cin >> mu1;
            std::cout << "mu2:";
            double mu2{};
            std::cin >> mu2;
            std::cout << "sigma1:";
            double sigma1{};
            std::cin >> sigma1;
            std::cout << "sigma2:";
            double sigma2{};
            std::cin >> sigma2;
            std::cout << "rho:";
            double rho{};
            std::cin >> rho;
            theMulDist_prt = new BivariateNormal(mu1, mu2, sigma1, sigma2, rho);
            break;
        }
        default:
            theMulDist_prt = nullptr;
            break;
    }
    
    std::cout << "Time at which the density of the first marginal can be truncated:";
    double truncationPoint1{};
    std::cin >> truncationPoint1;

    
    std::cout << "Time at which the density of the second marginal can be truncated:";
    double truncationPoint2{};
    std::cin >> truncationPoint2;
    
    std::cout << "Maximum acceptable probability in one point: ";
    double maxProbability{};
    std::cin >> maxProbability;
    
    std::cout << "Maximum time interval corresponding to one point: ";
    double maxDeltat{};
    std::cin >> maxDeltat;
    
    inputMultivariateDensity(theMulDist_prt, truncationPoint1, truncationPoint2, maxProbability, maxDeltat, observations, sumObservations);
    
    delete theMulDist_prt;
}


void inputMultivariateDensity( MultivariateDistribution *theMulDist_prt, double truncationPoint1, double truncationPoint2, double maxProbability, double maxDeltat, std::vector<std::vector<Sample>> & observations, std::vector<Sample> & sumObservations) {
    
    Sample auxiliar;
    
    //Marginal 1
    double deltat{};
    double t{0.0};
    
    while (t < truncationPoint1) {
        if (theMulDist_prt->densityMarginal1(t) < maxProbability / maxDeltat) {
            deltat = maxDeltat;
        }
        else {
            deltat = maxProbability / theMulDist_prt->densityMarginal1(t);
        }
        auxiliar.weight = deltat / 6 * (theMulDist_prt->densityMarginal1(t) + 4 * theMulDist_prt->densityMarginal1(t + deltat / 2) + theMulDist_prt->densityMarginal1(t + deltat));
        while (auxiliar.weight > maxProbability) {
            deltat = deltat * 0.9;
            auxiliar.weight = deltat / 6 * (theMulDist_prt->densityMarginal1(t) + 4 * theMulDist_prt->densityMarginal1(t + deltat / 2) + theMulDist_prt->densityMarginal1(t + deltat));
        }
        if (auxiliar.weight > 0) {
            auxiliar.obs = (t * theMulDist_prt->densityMarginal1(t) + 4 * (t + deltat / 2) * theMulDist_prt->densityMarginal1(t + deltat / 2) + (t + deltat) * theMulDist_prt->densityMarginal1(t + deltat / 2)) / (theMulDist_prt->densityMarginal1(t) + 4 * theMulDist_prt->densityMarginal1(t + deltat / 2) + theMulDist_prt->densityMarginal1(t + deltat));
            observations[0].push_back(auxiliar);
        }
        t += deltat;
    }
    
    
    //Marginal 2
    deltat = 0.0;
    t = 0.0;
    
    while (t < truncationPoint2) {
        if (theMulDist_prt->densityMarginal2(t) < maxProbability / maxDeltat) {
            deltat = maxDeltat;
        }
        else {
            deltat = maxProbability / theMulDist_prt->densityMarginal2(t);
        }
        auxiliar.weight = deltat / 6 * (theMulDist_prt->densityMarginal2(t) + 4 * theMulDist_prt->densityMarginal2(t + deltat / 2) + theMulDist_prt->densityMarginal2(t + deltat));
        while (auxiliar.weight > maxProbability) {
            deltat = deltat * 0.9;
            auxiliar.weight = deltat / 6 * (theMulDist_prt->densityMarginal2(t) + 4 * theMulDist_prt->densityMarginal2(t + deltat / 2) + theMulDist_prt->densityMarginal2(t + deltat));
        }
        if (auxiliar.weight > 0) {
            auxiliar.obs = (t * theMulDist_prt->densityMarginal2(t) + 4 * (t + deltat / 2) * theMulDist_prt->densityMarginal2(t + deltat / 2) + (t + deltat) * theMulDist_prt->densityMarginal2(t + deltat / 2)) / (theMulDist_prt->densityMarginal2(t) + 4 * theMulDist_prt->densityMarginal2(t + deltat / 2) + theMulDist_prt->densityMarginal2(t + deltat));
            observations[1].push_back(auxiliar);
        }
        t += deltat;
    }
    
    
    //Sum
    deltat = 0.0;
    t = 0.0;
    
    while (t < truncationPoint1 + truncationPoint2) {
        if (theMulDist_prt->densitySum(t) < maxProbability / maxDeltat) {
            deltat = maxDeltat;
        }
        else {
            deltat = maxProbability / theMulDist_prt->densitySum(t);
        }
        auxiliar.weight = deltat / 6 * (theMulDist_prt->densitySum(t) + 4 * theMulDist_prt->densitySum(t + deltat / 2) + theMulDist_prt->densitySum(t + deltat));
        while (auxiliar.weight > maxProbability) {
            deltat = deltat * 0.9;
            auxiliar.weight = deltat / 6 * (theMulDist_prt->densitySum(t) + 4 * theMulDist_prt->densitySum(t + deltat / 2) + theMulDist_prt->densitySum(t + deltat));
        }
        if (auxiliar.weight > 0) {
            auxiliar.obs = (t * theMulDist_prt->densitySum(t) + 4 * (t + deltat / 2) * theMulDist_prt->densitySum(t + deltat / 2) + (t + deltat) * theMulDist_prt->densitySum(t + deltat / 2)) / (theMulDist_prt->densitySum(t) + 4 * theMulDist_prt->densitySum(t + deltat / 2) + theMulDist_prt->densitySum(t + deltat));
            sumObservations.push_back(auxiliar);
        }
        t += deltat;
    }
}



void sampleStatisticsMultDist(const std::vector<std::vector<Sample>> & observations) {
    int dim{};
    dim = static_cast<int>(observations.size());
    
    Numeric_lib::Matrix<double,2> meanVec(1,dim);
    Numeric_lib::Matrix<double,2> sumWeigths(1,dim);
    Numeric_lib::Matrix<double,2> sdVec(1,dim);
    
    for (int j{0}; j < dim; ++j) {
        double sumOfx{0.0};
        double sumOfxSquared{0.0};
        
        for (int i{0}; i < observations[j].size(); ++i) {
            sumOfx += observations[j][i].obs * observations[j][i].weight;
            sumOfxSquared += observations[j][i].obs * observations[j][i].obs * observations[j][i].weight;
            sumWeigths(0,j) += observations[j][i].weight;
        }
        
        meanVec(0,j) = sumOfx / sumWeigths(0,j);
        sdVec(0,j) = sqrt(sumOfxSquared / (sumWeigths(0,j)) - pow(sumOfx / (sumWeigths(0,j)), 2.0));
    }
    
    std::cout << "Dimension: " << observations.size() << '\n';
    std::cout << "Sample sizes: " << observations[1].size() << '\t' << observations[0].size() << '\n';
    std::cout << "Sum of weights: " << '\n';
    printMatrix(sumWeigths);
    std::cout << "Sample mean: " << '\n';
    printMatrix(meanVec);
    std::cout << "Sample standard deviation: " << '\n';
    printMatrix(sdVec);
    
}



void askBivariateDist(std::vector<BivariateSample> & observations) {
    int typeContDist{};
    do {
        std::cout << "\nType of density:\n";
        std::cout << "      1. Bivariate Pareto\n";
        std::cout << "      2. Bivariate M-O Exponential\n";
        std::cout << "      3. Bivariate Normal\n";
        std::cout << "      4. Bivariate Weibull\n";
        std::cout << "Select(1-4): ";
        std::cin >> typeContDist;
        if (typeContDist < 1 || typeContDist > 4) {
            std::cout << "Please enter a valid option\n";
        }
    } while (typeContDist < 1 || typeContDist > 4);
    
    
    MultivariateDistribution *theMulDist_prt;
    
    switch(typeContDist) {
        case 1: {
            std::cout << "\nBivariate Pareto with shape parameter alpha and scale parameters theta1 and theta2:\n";
            std::cout << "alpha:";
            double alpha{};
            std::cin >> alpha;
            std::cout << "theta1:";
            double theta1{};
            std::cin >> theta1;
            std::cout << "theta2:";
            double theta2{};
            std::cin >> theta2;
            theMulDist_prt = new BivariatePareto(alpha, theta1, theta2);
            break;
        }
        case 2: {
            std::cout << "\nMarshall Olkin’s Bivariate Exponential with parameters lambda1, lambda2 and lambda12.\n";
            std::cout << "lambda1:";
            double lambda1{};
            std::cin >> lambda1;
            std::cout << "lambda2:";
            double lambda2{};
            std::cin >> lambda2;
            std::cout << "lambda12:";
            double lambda12{};
            std::cin >> lambda12;
            theMulDist_prt = new MarshallOlkinBivariateExponential(lambda1, lambda2, lambda12);
            break;
        }
        case 3: {
            std::cout << "\n Bivariate Normal with parameters mu1, mu2, sigma1, sigma2 and rho.\n";
            std::cout << "mu1:";
            double mu1{};
            std::cin >> mu1;
            std::cout << "mu2:";
            double mu2{};
            std::cin >> mu2;
            std::cout << "sigma1:";
            double sigma1{};
            std::cin >> sigma1;
            std::cout << "sigma2:";
            double sigma2{};
            std::cin >> sigma2;
            std::cout << "rho:";
            double rho{};
            std::cin >> rho;
            theMulDist_prt = new BivariateNormal(mu1, mu2, sigma1, sigma2, rho);
            break;
        }
        default:
            theMulDist_prt = nullptr;
            break;
    }
    
    std::cout << "Time at which the density of the first marginal can be truncated:";
    double truncationPoint1{};
    std::cin >> truncationPoint1;

    
    std::cout << "Time at which the density of the second marginal can be truncated:";
    double truncationPoint2{};
    std::cin >> truncationPoint2;
    
    std::cout << "Maximum acceptable probability in one bivariate point: ";
    double maxProbability{};
    std::cin >> maxProbability;
    
    std::cout << "Maximum time interval corresponding to one point of the first component: ";
    double maxDeltat1{};
    std::cin >> maxDeltat1;
    
    std::cout << "Maximum time interval corresponding to one point of the second component: ";
    double maxDeltat2{};
    std::cin >> maxDeltat2;
    
    inputBivariateDensity(theMulDist_prt, truncationPoint1, truncationPoint2, maxProbability, maxDeltat1, maxDeltat2, observations);
    
    delete theMulDist_prt;
}



void inputBivariateDensity( MultivariateDistribution *theMulDist_prt, double truncationPoint1, double truncationPoint2, double maxProbability, double maxDeltat1, double maxDeltat2, std::vector<BivariateSample> & observations) {
    
    BivariateSample auxiliar;
    
    double deltat1{};
    double deltat2{};
    double t1{0.0};
    double t2{0.0};
    
    while (t2 < truncationPoint2) {
        
        t1 = 0.0;
        while (t1 < truncationPoint1) {
            deltat1 = maxDeltat1;
            deltat2 = maxDeltat2;
            
            auxiliar.weight = deltat1 * deltat2 / 36 * ( theMulDist_prt->jointDensity(t1,t2) + 4 * theMulDist_prt->jointDensity(t1, t2 + deltat2 / 2) + theMulDist_prt->jointDensity(t1, t2 + deltat2) + 4 * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2) + 16 * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2 + deltat2 / 2) + 4 * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2 + deltat2) + theMulDist_prt->jointDensity(t1 + deltat1, t2) + 4 * theMulDist_prt->jointDensity(t1 + deltat1, t2 + deltat2 / 2) +
                    theMulDist_prt->jointDensity(t1 + deltat1, t2 + deltat2));
            while (auxiliar.weight > maxProbability) {
                deltat1 = deltat1 * 0.9;
                deltat2 = deltat2 * 0.9;
                auxiliar.weight = deltat1 * deltat2 / 36 * ( theMulDist_prt->jointDensity(t1, t2) + 4 * theMulDist_prt->jointDensity(t1, t2 + deltat2 / 2) + theMulDist_prt->jointDensity(t1, t2 + deltat2) + 4 * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2) + 16 * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2 + deltat2 / 2) + 4 * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2 + deltat2) + theMulDist_prt->jointDensity(t1 + deltat1, t2) + 4 * theMulDist_prt->jointDensity(t1 + deltat1, t2 + deltat2 / 2) +
                theMulDist_prt->jointDensity(t1 + deltat1, t2 + deltat2));
            }
            
            if (auxiliar.weight > 0) {
                auxiliar.x1 = (t1 * theMulDist_prt->jointDensity(t1, t2) + 4 * (t1 + deltat1 / 2) * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2) + (t1 + deltat1) * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2)) / (theMulDist_prt->jointDensity(t1, t2) + 4 * theMulDist_prt->jointDensity(t1 + deltat1 / 2, t2) + theMulDist_prt->jointDensity(t1 + deltat1, t2));
                
                auxiliar.x2 = (t2 * theMulDist_prt->jointDensity(t1, t2) + 4 * (t2 + deltat2 / 2) * theMulDist_prt->jointDensity(t1, t2 + deltat2 / 2) + (t2 + deltat2) * theMulDist_prt->jointDensity(t1, t2 + deltat2 / 2)) / (theMulDist_prt->jointDensity(t1, t2) + 4 * theMulDist_prt->jointDensity(t1, t2 + deltat2 / 2) + theMulDist_prt->jointDensity(t1, t2 + deltat2));
                
                observations.push_back(auxiliar);
            }
            t1 += deltat1;
        }
        
        t2 += deltat2;
    }

    
}


// EM for NPH

double NPHdensity(double x, PhaseType & PH, ScaleDistribution *N) {
    double epsilon{1e-6};
    
    double density{};
    double previous{};
    double relError{1};
    int i{1};
    
    while (relError > epsilon){
        density += (N->density(i)) * PH.densityLevel(x, N->si(i));
        
        if (i > 1) {
            relError = abs((density - previous) / previous);
        }
        
        previous = density;
        ++i;
    }
    return density;
}

double NPHtail(double x, PhaseType & PH, ScaleDistribution *N) {
    double epsilon{1e-6};
    
    double density{};
    double previous{};
    double relError{1};
    int i{1};
    
    while (relError > epsilon){
        density += (N->density(i)) * PH.tail(x / N->si(i));
        
        if (i > 1) {
            relError = abs((density - previous) / previous);
        }
        
        previous = density;
        ++i;
    }
    return density;
}

void relativeErrorMatrices(const Numeric_lib::Matrix<double,2> & M1, const Numeric_lib::Matrix<double,2> & M2, Numeric_lib::Matrix<double,2> & E) {
    
    for (int i{0}; i < M1.dim1(); ++i) {
        for (int j{0}; j < M1.dim2(); ++j) {
            if (M1(i,j) == M2(i,j)) {
                E(i,j) = 0;
            } else {
                E(i,j) = abs((M1(i,j) - M2(i,j)) / M2(i,j));
            }
        }
    }
}


std::vector<double> thetaFn(std::vector<double> & weigths, double epsilon, PhaseType & PH, ScaleDistribution *N, double theta_org, const std::vector<Sample> & observations, const std::vector<Sample> & censored) {
    
    Numeric_lib::Matrix<double,2> pi(PH.getPi());
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    
    std::vector<double> result;
    double density{};
    double fn=0, fnp=0, prev1=0, prev2=0;
    
    double theta{N->getParameter()};
    
    double relativeError{1.0};
    double phdensity{};
    int level{1};
    
    while (relativeError > epsilon) {
        
        
        if (weigths.size() < level) {
            
            N->updateParameter(theta_org);
            
            double si{N->si(level)};
            double densityScale{N->density(level)};
            
            double factor{};
            //Uncensored
            for (int k{0}; k < observations.size(); ++k) {
                phdensity = PH.density(observations[k].obs / si) / si;
                density = NPHdensity(observations[k].obs, PH, N);
                factor += densityScale * phdensity * observations[k].weight  / density;
            }
            //Censored
            for (int k{0}; k < censored.size(); ++k) {
                phdensity = PH.tail(censored[k].obs / si);
                density = NPHtail(censored[k].obs, PH, N);
                factor += densityScale * phdensity * censored[k].weight  / density;
            }
            
            weigths.push_back(factor);
            
            N->updateParameter(theta);
        }
       
        fn += weigths[level - 1] * (N->firstDerivative(level)) / (N->density(level));
        fnp += weigths[level - 1] * (((N->density(level)) * (N->secondDerivative(level)) - pow(N->firstDerivative(level), 2)) / pow(N->density(level), 2));
        
        
        if (level > 1) {
            relativeError = fmax(abs((fn - prev1) / prev1), abs((fnp - prev2) / prev2));
        }
        prev1 = fn;
        prev2 = fnp;
        ++level;
    }
    result.push_back(fn);
    result.push_back(fnp);
    return (result);
}


double thetaFnGA(std::vector<double> & weigths, double epsilon, PhaseType & PH, ScaleDistribution *N, double theta_org, const std::vector<Sample> & observations, const std::vector<Sample> & censored, const std::vector<double> & densityVector, const std::vector<double> & densityVectorCensored) {
    
    Numeric_lib::Matrix<double,2> pi(PH.getPi());
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    
    double density{};
    double fn=0, prev1=0;
    
    double theta{N->getParameter()};
    
    double relativeError{1.0};
    double phdensity{};
    int level{1};
    
    while (relativeError > epsilon) {
        
        
        if (weigths.size() < level) {
            
            N->updateParameter(theta_org);
            
            double si{N->si(level)};
            double densityScale{N->density(level)};
            
            double factor{};
            //Uncensored
            for (int k{0}; k < observations.size(); ++k) {
                phdensity = PH.density(observations[k].obs / si) / si;
                density = densityVector[k];
                factor += densityScale * phdensity * observations[k].weight  / density;
            }
            //Censored
            for (int k{0}; k < censored.size(); ++k) {
                phdensity = PH.tail(censored[k].obs / si);
                density = densityVectorCensored[k];
                factor += densityScale * phdensity * censored[k].weight  / density;
            }
            
            weigths.push_back(factor);
            
            N->updateParameter(theta);
        }
       
        fn += weigths[level - 1] * (N->firstDerivative(level)) / (N->density(level));
        
        
        if (level > 1) {
            relativeError = abs((fn - prev1) / prev1);
        }
        prev1 = fn;
        ++level;
    }
    return (fn);
}

void EMstepNPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored,  PhaseType & PH, ScaleDistribution *N, int fixParameter, int maxMethod, double lambda) {
    long p{PH.getp()};
    Numeric_lib::Matrix<double,2> pi(PH.getPi());
    Numeric_lib::Matrix<double,2> T(PH.getT());
    Numeric_lib::Matrix<double,2> t(PH.gett());
    Numeric_lib::Matrix<double,2> e(PH.gete());
    
    Numeric_lib::Matrix<double,2> Bmean(p,1);
    Numeric_lib::Matrix<double,2> Zmean(p,1);
    Numeric_lib::Matrix<double,2> Nmean(p,p + 1);
    
    Numeric_lib::Matrix<double,2> Bmean_prev(p,1);
    Numeric_lib::Matrix<double,2> Zmean_prev(p,1);
    Numeric_lib::Matrix<double,2> Nmean_prev(p,p + 1);
    
    Numeric_lib::Matrix<double,2> relativeErrorB(p,1);
    Numeric_lib::Matrix<double,2> relativeErrorZ(p,1);
    Numeric_lib::Matrix<double,2> relativeErrorN(p,p + 1);
    
    Numeric_lib::Matrix<double,2> avector(1,p);
    Numeric_lib::Matrix<double,2> bvector(p,1);
    Numeric_lib::Matrix<double,2> cmatrix(p,p);
    Numeric_lib::Matrix<double,2> aux_exp(p,p);
    
    Numeric_lib::Matrix<double,2> J(2 * p,2 * p);
    Numeric_lib::Matrix<double,2> tProductPi(p,p);
    
    std::vector<double> densityVector; //Vector that contains the densities to avoid multiple computations
    std::vector<double> LVector; // Vector that contains the expectec value of L - to avoid multiple computations
    
    std::vector<double> densityVectorCensored; //Vector that contains the tail of censored observations to avoid multiple computations
    
    int level{1};
    double SumOfWeights{0.0};
    double SumOfCensored{0.0};
    double density{0.0};
    double maxError{1.0};
    
    //Constants
    double epsilonMatrices{1e-3}; // For the maximun error in the computation of B, Z, N
    double epsilonConv{1e-3}; // For N-R
    
    double phdensity{};;
    
    //E-step
    while (maxError > epsilonMatrices) {
        double si{N->si(level)};
        double densityScale{N->density(level)};
        double Lweight{};
        
        //  Unccensored data
        for (int k{0}; k < observations.size(); ++k) {
            
            matrixProduct(t / si, pi, tProductPi);
            
            VanLoanForIntegrals(observations[k].obs, T / si, T / si, tProductPi, J);
            
            for (int i{0}; i < p; ++i) {
                for (int j{0}; j < p; ++j) {
                    aux_exp(i,j) = J(i,j);
                    cmatrix(i,j) = J(i,j + p);
                }
            }
            
            matrixProduct(pi, aux_exp, avector);
            matrixProduct(aux_exp, t / si, bvector);
            
            
            if (level == 1) { // Computations required only once
                //Sum of Weights
                SumOfWeights += observations[k].weight;
                
                //The densitity in each observation
                density = NPHdensity(observations[k].obs, PH, N);
                densityVector.push_back(density);
            }
            else {
                //Recovers the value of the density for further levels
                density = densityVector[k];
            }
            
            if (fixParameter == 1) { // In case we need to estimate the parameter of the scale component, we need to compute E(L)
                phdensity = PH.density(observations[k].obs / si) / si;
                Lweight += densityScale * phdensity * observations[k].weight  / density;
            }
            
            
            //E-step
            for (int i{0}; i < p; ++i) {
                Bmean(i,0) += densityScale * pi(0,i) * bvector(i,0) * observations[k].weight / density;
                Nmean(i,p) += densityScale * avector(0,i) * t(i,0) * observations[k].weight / (density * si);
                Zmean(i,0) += densityScale * cmatrix(i,i) * observations[k].weight / (density * si);
                for (int j{0}; j < p; ++j) {
                    Nmean(i,j) += densityScale * T(i,j) * cmatrix(j,i) * observations[k].weight / (density * si);
                }
            }
        }
        
        //  Right-Censored Data
        if (censored.size() > 0) {
            matrixProduct(e, pi, tProductPi);
        }
        for (int k{0}; k < censored.size(); ++k) {
            
            VanLoanForIntegrals(censored[k].obs, T / si, T / si, tProductPi, J);
            
            for (int i{0}; i < p; ++i) {
                for (int j{0}; j < p; ++j) {
                    aux_exp(i,j) = J(i,j);
                    cmatrix(i,j) = J(i,j + p);
                }
            }
            
            matrixProduct(aux_exp, e, bvector);
            
            if (level == 1) { // Sum only once
                SumOfCensored += censored[k].weight;
                
                density = NPHtail(censored[k].obs, PH, N);
                densityVectorCensored.push_back(density);
            }
            else {
                density = densityVectorCensored[k];
            }
            
            
            if (fixParameter == 1) { // In case we need to estimate the parameter of the scale component, we need to compute E(L)
                phdensity = PH.tail(observations[k].obs / si);
                Lweight += densityScale * phdensity * observations[k].weight  / density;
            }
            
            //E-step
            for (int i{0}; i < p; ++i) {
                Bmean(i,0) += density * pi(0,i) * bvector(i,0) * censored[k].weight / density;
                Zmean(i,0) += density * cmatrix(i,i) * censored[k].weight / (density * si);
                for (int j{0}; j < p; ++j) {
                    Nmean(i,j) += density * T(i,j) * cmatrix(j,i) * censored[k].weight / (density * si);
                }
            }
        }
        
        //Compute errors
        if (level>1){
            relativeErrorMatrices(Bmean, Bmean_prev, relativeErrorB);
            relativeErrorMatrices(Nmean, Nmean_prev, relativeErrorN);
            relativeErrorMatrices(Zmean, Zmean_prev, relativeErrorZ);
            
            maxError = fmax(matrixMax(relativeErrorB), fmax(matrixMax(relativeErrorN), matrixMax(relativeErrorZ)));
        }
        
        Bmean_prev = Bmean;
        Zmean_prev = Zmean;
        Nmean_prev = Nmean;
        
        if (fixParameter == 1) { // Save value
            LVector.push_back(Lweight);
        }
        
        ++level;
        
    }
    
    
    // M step
    //Scale component
    if (fixParameter == 1) { // Case we need to estimate parameter
        int idN{N->getid()};
        
        if (idN == 1) { // Riemann Zeta - NR but simplify
            double epsilonTheta{1e-2};
            double theta{N->getParameter()};
            double theta_prev{theta};
            double relativeError{1};
            double relErrorTheta{1};
            double Fn{};
            double FnPrime{};
            double factor{0};
            double prev{};
            
            level = 1;
            while(relativeError > epsilonConv){
                double densityScale{N->density(level)};
                
                if (LVector.size() < level) {
                    //Uncensored
                    for (int k{0}; k < observations.size(); ++k) {
                        phdensity = PH.density(observations[k].obs / level) / level;
                        density = densityVector[k];
                        factor += log(level) * densityScale * phdensity * observations[k].weight  / density;
                    }
                    //Censored
                    for (int k{0}; k < censored.size(); ++k) {
                        phdensity = PH.tail(censored[k].obs / level);
                        density = densityVectorCensored[k];
                        factor += log(level) * densityScale * phdensity * censored[k].weight  / density;
                    }
                }
                else {
                    factor += log(level) * LVector[level - 1];
                }
                
                if (level > 1) {
                    relativeError = abs((factor - prev) / prev);
                }
                
                prev = factor;
                ++level;
            }
            factor = factor / (SumOfWeights + SumOfCensored);
            
            while (relErrorTheta > epsilonTheta) {
                Fn = (N->firstDerivative(1)) / (N->auxiliarFn()) + factor;
                FnPrime = ((N->secondDerivative(1)) * (N->auxiliarFn()) - (N->firstDerivative(1)) * (N->firstDerivative(1))) / ((N->auxiliarFn()) * (N->auxiliarFn()));
                theta = theta_prev - Fn / FnPrime;
                relErrorTheta = abs((theta - theta_prev) / theta_prev); // Se podria cambiar a convergencia cerca de cero de Fn - Igual que en los otros
                theta_prev = theta;
                N->updateParameter(theta);
            }
        }
        else if (idN == 2) { // Pareto Geometric - Explicit solution
            double relativeError{1};
            double c{};
            double factor{0};
            double prev{};
                   
            c = log(N->auxiliarFn());
            
            level = 1;
                   
            while(relativeError > epsilonConv){
                
                double si{N->si(level)};
                double densityScale{N->density(level)};
                
                if (LVector.size() < level) {
                    //Uncensored
                    for (int k{0}; k < observations.size(); ++k) {
                        phdensity = PH.density(observations[k].obs / si) / si;
                        density = densityVector[k];
                        factor += level * densityScale * phdensity * observations[k].weight  / density;
                    }
                    //Censored
                    for (int k{0}; k < censored.size(); ++k) {
                        phdensity = PH.tail(censored[k].obs / si);
                        density = densityVectorCensored[k];
                        factor += level * densityScale * phdensity * censored[k].weight  / density;
                    }
                    
                }
                else {
                    factor += level * LVector[level - 1];
                }
                
                if(level>1) {
                    relativeError = abs((factor - prev) / prev);
                }
                prev = factor;
                ++level;
            }
            
            N->updateParameter(-log(1 - (SumOfWeights + SumOfCensored) / factor) / c);
        }
        else { // NR
            double epsilonFn{1e-3};
            double theta_org{N->getParameter()};
            double theta{theta_org};
            
            if (maxMethod == 1) { //Gradient ascent
                double fnPrime{};
                fnPrime = thetaFnGA(LVector, epsilonFn, PH, N, theta_org, observations, censored, densityVector, densityVectorCensored);
                while (abs(fnPrime) > epsilonConv) {
                    theta = theta + lambda * fnPrime;
                    N->updateParameter(theta);
                    fnPrime = thetaFnGA(LVector, epsilonFn, PH, N, theta_org, observations, censored, densityVector, densityVectorCensored);
                }
            }
            else {
                std::vector<double> functionval;
                
                double theta_prev{theta};
                double relativeError{1};
                
                while (relativeError > epsilonConv) {
                    functionval = thetaFn(LVector, epsilonFn, PH, N, theta_org, observations, censored);
                    theta = theta_prev - functionval[0] / functionval[1];
                    relativeError = abs((theta - theta_prev) / theta_prev); // Creo q deberia de usar functionval[0] cerca a cero - si no podria hacer mas de la cuenta
                    theta_prev = theta;
                    N->updateParameter(theta);
                }
            }
        }
    
    }
    
    
    //Parameters of the PH
    for (int i{0}; i < p; ++i) {
        pi(0,i) = Bmean(i,0) / (SumOfWeights + SumOfCensored);
        if (pi(0,i) < 0) {
            pi(0,i) = 0;
        }
        t(i,0) = Nmean(i,p) / Zmean(i,0);
        if (t(i,0) < 0) {
            t(i,0) = 0;
        }
        T(i,i) = -t(i,0);
        for (int j{0}; j < p; ++j) {
            if (i != j) {
                T(i,j) = Nmean(i,j) / Zmean(i,0);
                if (T(i,j) < 0) {
                    T(i,j) = 0;
                }
                T(i,i) -= T(i,j);
            }
        }
    }
    PH.changeAll(pi, T);
}


// Loglikelihood using Matlab matrix exponential
double logLikelihoodNPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, ScaleDistribution *N) {
    double logLh{0.0};
    for (int k{0}; k < observations.size(); ++k) {
        logLh += observations[k].weight * log(NPHdensity(observations[k].obs, PH, N));
    }
    for (int k{0}; k < censored.size(); ++k) {
        logLh += censored[k].weight * log(NPHtail(censored[k].obs, PH, N));
    }
    return logLh;
}


void EMIterateNPH(int stepsEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, ScaleDistribution *N, int fixParameter, int printLikehood, int maxMethod, double lambda) {
    double relError{1.0};
    double logLh{0.0};
    double logLh_ant{0.0};
    
    logLh = logLikelihoodNPH(observations, censored, phaseTypeEM, N);
    logLh_ant = logLh;
    for (int i{1}; i <= stepsEM; ++i) {
        EMstepNPH(observations, censored, phaseTypeEM, N, fixParameter, maxMethod, lambda);
        if(i % printLikehood == 0) {
            logLh = logLikelihoodNPH(observations, censored, phaseTypeEM, N);
            relError = abs((logLh - logLh_ant) / logLh_ant);
            logLh_ant = logLh;
            std::cout << "Iteration " << i << ". Loglikehood: " << logLh << ". Relative difference: " << relError << " Parameter: " << N->getParameter() << '\n';
        }
    }
}
