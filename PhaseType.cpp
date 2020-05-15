// Implementation of the functions of the phaseType class 

#include "PhaseType.h"


/* *********************
        PH class
********************* */

//  Constructors
//      Constructor using a pi and a T as input
PhaseType::PhaseType(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T) : m_pi(pi), m_T(T), m_p(T.dim1()), m_e(T.dim1(),1), m_t(T.dim1(),1) {
    m_e += 1;
    matrixProduct(m_T * (-1.0), m_e, m_t);
}
//      Initializated in zeros with dimension p
PhaseType::PhaseType(long p) : m_pi(1,p), m_T(p,p), m_p(p), m_e(p,1), m_t(p,1) {
    m_e += 1;
}

void PhaseType::readParametersFromFile(std::ifstream & inputFile) {
    double aux{0.0};
    for (int i{0}; i < m_p; ++i) {
        for (int j{0}; j < 1 + m_p; ++j) {
            if (j == 0) {
                inputFile >> aux;
                m_pi(0,i) = aux;
            }
            else {
                inputFile >> aux;
                m_T(i,j - 1) = aux;
            }
        }
    }
    matrixProduct(m_T * (-1.0), m_e, m_t);
}


//  Printing Data
//      Prints Pi - Initial distribution
void PhaseType::printPi() {
     printMatrix(m_pi);
 }
//      Prints T - intensity matrix
 void PhaseType::printT() {
     printMatrix(m_T);
 }
//      Prints t - Exit rates
 void PhaseType::printt() {
    printMatrix(m_t);
}
//      Prints t - vector of ones
 void PhaseType::printe() {
    printMatrix(m_e);
}
// Prints pi,T and t
void PhaseType::print() {
    std::cout << "Initial probabilities (Pi):\n";
    printPi();
    std::cout << "Transition matrix (T):\n";
    printT();
    std::cout << "Exit rates (t):\n";
    printt();
}


//  Get data
Numeric_lib::Matrix<double,2> PhaseType::getPi() {
    return m_pi;
}
Numeric_lib::Matrix<double,2> PhaseType::getT() {
    return m_T;
}
Numeric_lib::Matrix<double,2> PhaseType::gett() {
    return m_t;
}
Numeric_lib::Matrix<double,2> PhaseType::gete() {
    return m_e;
}
long PhaseType::getp() {
    return m_p;
}


//  Change Data
void PhaseType::changePi(const Numeric_lib::Matrix<double,2> & pi) {
    m_pi = pi;
}
void PhaseType::changeT(const Numeric_lib::Matrix<double,2> & T) {
    m_T = T;
    matrixProduct(m_T * (-1.0), m_e, m_t);
}
void PhaseType::changeAll(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T) {
    m_pi = pi;
    m_T = T;
    matrixProduct(m_T * (-1.0), m_e, m_t);
}

//  Properties
double PhaseType::density(double x) {
    if (x == 0) {
        return (1.0 - matrixProduct(m_pi, m_e)(0,0));
    }
    return (matrixProduct(m_pi, matrixProduct(matrixExponential(m_T * x), m_t))(0,0));
}
double PhaseType::tail(double x) {
    return (matrixProduct(m_pi, matrixProduct(matrixExponential(m_T * x), m_e))(0,0));
}
double PhaseType::distribution(double x) {
    return (1.0 - tail(x));
}

double PhaseType::LaplaceTransform(double s) {
    return matrixProduct(m_pi, matrixProduct(matrixInverse(addMatrices(identityMatrix(static_cast<int>(m_p)) * s, m_T * (-1.0))), m_t))(0,0);
}

double PhaseType::quantile(double u, double epsilon) {
    int n{10000};
    double delta{mean() / n};
    double x0{0};
     
    while (distribution(x0) - u < 0) {
        x0 += delta;
    } 
    double x1{fmax(x0 - delta, 0)};
    double xm{x0};
    double fx{distribution(x0) - u};
    
    if (fx > 0) {
        xm = 0.5 * (x0 + x1);
        double fxm{distribution(xm) - u};
            
        while (abs((fxm - fx) / fx) > epsilon) {
            fx = fxm;
            if (fxm < 0) {
                x1 = xm;
            }
            else {
                x0 = xm;
            }
            xm = 0.5 * (x0 + x1);
            fxm = distribution(xm) - u;
        }
    }
    return xm;
}

double PhaseType::densityLevel(double x, double si) {
    if (x == 0) {
        return (1.0 - matrixProduct(m_pi, m_e)(0,0));
    }
    return (matrixProduct(m_pi, matrixProduct(matrixExponential(m_T * (x / si)), m_t * (1 / si)))(0,0));
}



//  Functions for the EM
//      Generates random parameters for pi, T and t based on the marks in piLegal and TLegal
void PhaseType::randomPhase(const Numeric_lib::Matrix<double,2> & piLegal, const Numeric_lib::Matrix<double,2> & TLegal, double scaleFactor, AbsRNG & rUniform) {
    long p{m_T.dim1()};
    double sum{0.0};
        
    for (int i{0}; i < p; ++i) {
        if (piLegal(0,i) == 1) {
            m_pi(0,i) = rUniform();
            sum += m_pi(0,i);
        }
    }
        
    for (int i{0}; i < p; ++i) {
        m_pi(0,i) = m_pi(0,i) / sum;
    }
        
    for (int i{0}; i < p; ++i) {
        for (int j{0}; j < p; ++j) {
            if ((i != j) && (TLegal(i,j) == 1)) {
                m_T(i,j) = rUniform();
                m_T(i,i) -= m_T(i,j);
            }
        }
    }
    
    for (int i{0}; i < p; ++i) {
        if (TLegal(i,i) == 1) {
            m_t(i,0) = rUniform();
            m_T(i,i) -= m_t(i,0);
        }
    }
    m_T *= (p / scaleFactor);
    m_t *= (p / scaleFactor);
        
}

double PhaseType::moment(int k) {
    Numeric_lib::Matrix<double,2> U(matrixInverse(m_T * (-1.0)));
    
    U = matrixPower(k, U);
    
    return factorial(k) * matrixProduct(m_pi, matrixProduct(U, m_e))(0,0);
}

double PhaseType::mean() {
    return moment(1);
}

double PhaseType::var() {
    return moment(2) - pow(moment(1), 2.0);
}

double PhaseType::sd() {
    return sqrt(var());
}



void PhaseType::exportToFile(std::ofstream & outputFile) {
    for (int i{0}; i < m_p; ++i) {
        outputFile << m_pi(0,i) << '\t';
        for (int j{0}; j < m_p; ++j) {
            if (j == m_p - 1) {
                outputFile << m_T(i,j);
            }
            else {
                outputFile << m_T(i,j) << '\t';
            }
        }
        if(i < m_p-1) {
            outputFile << '\n';
        }
    }
}


/* *********************
        MPH class
********************* */

//  Constructors
//      Constructor using a pi, T and R as input
MPH::MPH(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T, const Numeric_lib::Matrix<double,2> & R) : m_pi(pi), m_T(T), m_R(R), m_p(T.dim1()), m_d(R.dim2()), m_e(T.dim1(),1), m_t(T.dim1(),1) {
    m_e += 1;
    matrixProduct(m_T * (-1.0), m_e, m_t); // Lo puedo usar para la distribucion de la suma
}
//      Initializated in zeros with dimension p
MPH::MPH(long p, long d) : m_pi(1,p), m_T(p,p), m_R(p,d), m_p(p),m_d(d), m_e(p,1), m_t(p,1) {
    m_e += 1;
}


void MPH::readParametersFromFile(std::ifstream & inputFile) {
    double aux{0.0};
    for (int i{0}; i < m_p; ++i) {
        for (int j{0}; j < 1 + m_p + m_d; ++j) {
            if (j == 0) {
                inputFile >> aux;
                m_pi(0,i) = aux;
            }
            else if (j <= m_p) {
                inputFile >> aux;
                m_T(i,j - 1) = aux;
            }
            else {
                inputFile >> aux;
                m_R(i,j - m_p - 1) = aux;
            }
        }
    }
    matrixProduct(m_T * (-1.0), m_e, m_t);
}



//  Printing Data
//      Prints Pi - Initial distribution
void MPH::printPi() {
     printMatrix(m_pi);
 }
//      Prints T - intensity matrix
 void MPH::printT() {
     printMatrix(m_T);
 }
//
void MPH::printR() {
    printMatrix(m_R);
}
//      Prints t - Exit rates
 void MPH::printt() {
    printMatrix(m_t);
}
//      Prints t - vector of ones
 void MPH::printe() {
    printMatrix(m_e);
}
// Prints pi,T and t
void MPH::print() {
    std::cout << "Initial probabilities (Pi):\n";
    printPi();
    std::cout << "Transition matrix (T):\n";
    printT();
    std::cout << "Rewards matrix (R):\n";
    printR();
}

//  Change Data
void MPH::changePi(const Numeric_lib::Matrix<double,2> & pi) {
    m_pi = pi;
}
void MPH::changeT(const Numeric_lib::Matrix<double,2> & T) {
    m_T = T;
    matrixProduct(m_T * (-1.0), m_e, m_t);
}
void MPH::changeR(const Numeric_lib::Matrix<double,2> & R) {
    m_R = R;
}
void MPH::changeAll(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T, const Numeric_lib::Matrix<double,2> & R) {
    m_pi = pi;
    m_T = T;
    m_R = R;
    matrixProduct(m_T * (-1.0), m_e, m_t);
}


//  Get data
Numeric_lib::Matrix<double,2> MPH::getPi() {
    return m_pi;
}
Numeric_lib::Matrix<double,2> MPH::getT() {
    return m_T;
}
Numeric_lib::Matrix<double,2> MPH::getR() {
    return m_R;
}
Numeric_lib::Matrix<double,2> MPH::gett() {
    return m_t;
}
Numeric_lib::Matrix<double,2> MPH::gete() {
    return m_e;
}
long MPH::getp() {
    return m_p;
}
long MPH::getd() {
    return m_d;
}



void MPH::exportToFile( std::ofstream & outputFile ) {
    for (int i{0}; i < m_p; ++i) {
        outputFile << m_pi(0,i) << '\t';
        for (int j{0}; j < m_p; ++j) {
            outputFile << m_T(i,j) << '\t';
        }
        for (int j{0}; j < m_d; ++j) {
            if (j == m_d - 1) {
                outputFile << m_R(i,j);
            }
            else {
                outputFile << m_R(i,j) << '\t';
            }
        }
        if (i < m_p - 1) {
            outputFile << '\n';
        }
    }
}


PhaseType MPH::linearCombination(const std::vector<double> & w, std::vector<int> & newStates ) {
    long p{m_p};
    
    int NumZeros{0};
    std::vector<int> deleteRows; //states to be deleted
    std::vector<int> keepRows; //states to keep
    
    std::vector<double> Rw(matrixVectorProduct(m_R, w)); //overload fuction
    
    for (int j{0}; j < p; ++j) {
        if (Rw[j] == 0) {
            deleteRows.push_back(j);
            ++NumZeros;
        }
        else {
            keepRows.push_back(j);
        }
    }
    
    newStates = keepRows;
    
    PhaseType wPH(p - NumZeros);
    
    if (NumZeros == 0) {
        Numeric_lib::Matrix<double,2> diagonal(p,p);
        
        for (int i{0}; i < p; ++i) {
            diagonal(i,i) = 1.0 / Rw[i];
        }
        diagonal = matrixProduct(diagonal, m_T);
        
        wPH.changeAll(m_pi, diagonal);
        //return diagonal;
    }
    else {
        long n1 = deleteRows.size();
        long n2 = keepRows.size();
        
        Numeric_lib::Matrix<double,2> Spp(n2,n2);
        Numeric_lib::Matrix<double,2> Sp0(n2,n1);
        Numeric_lib::Matrix<double,2> S0p(n1,n2);
        Numeric_lib::Matrix<double,2> S00(n1,n1);
        
        Numeric_lib::Matrix<double,2> Taux(n2,n2);
        Numeric_lib::Matrix<double,2> diagonal(n2,n2);
        
        Numeric_lib::Matrix<double,2> pi0(1,n1);
        Numeric_lib::Matrix<double,2> pip(1,n2);
        
        Numeric_lib::Matrix<double,2> piaux(1,n2);
        
        for (int i{0}; i < n2; i++) {
            for (int j = 0; j < n2; j++) {
                Spp(i,j) = m_T(keepRows[i],keepRows[j]);
            }
            for (int j{0}; j < n1; j++) {
                Sp0(i,j) = m_T(keepRows[i],deleteRows[j]);
            }
            pip(0,i) = m_pi(0,keepRows[i]);
        }
        for (int i{0}; i < n1; i++) {
            for (int j{0}; j < n2; j++) {
                S0p(i,j) = m_T(deleteRows[i],keepRows[j]);
            }
            for (int j{0}; j < n1; j++){
                S00(i,j) = m_T(deleteRows[i],deleteRows[j]);
            }
            pi0(0,i) = m_pi(0,deleteRows[i]);
        }
        
        piaux = addMatrices(pip, matrixProduct(pi0, matrixProduct(matrixInverse(S00 * (-1.0)), S0p)));
        
        Taux = addMatrices(Spp, matrixProduct(Sp0, matrixProduct(matrixInverse(S00 * (-1.0)), S0p)));
        
        for (int i{0}; i < n2; ++i) {
            diagonal(i,i) = 1.0 / Rw[keepRows[i]];
        }
        diagonal = matrixProduct(diagonal, Taux);
        
        wPH.changeAll(piaux, diagonal);
        
        //return diagonal;
    }
    
    return wPH;
}


PhaseType MPH::marginal(int j, std::vector<int> & newStates) {
    std::vector<double> e(m_p);
    
    e[j] = 1.0;
    
    return linearCombination(e, newStates);
}

double MPH::jointMGF(const std::vector<double> & w) {
    
    std::vector<double> Rw(matrixVectorProduct(m_R, w));
    
    Numeric_lib::Matrix<double,2> diagonal(m_p,m_p);
    
    for (int i{0}; i < m_p; ++i) {
        diagonal(i,i) = -Rw[i];
    }
    return matrixProduct(m_pi, matrixProduct(matrixInverse(addMatrices(diagonal, m_T * (-1.0))), m_t))(0,0);
}


//      Generates random parameters for pi, T, R and t based on the marks in piLegal, TLegal and RLegal
void MPH::randomPhase(const Numeric_lib::Matrix<double,2> & piLegal, const Numeric_lib::Matrix<double,2> & TLegal, const Numeric_lib::Matrix<double,2> & RLegal, double scaleFactor, AbsRNG & rUniform) {
    long p{m_T.dim1()};
    long dim{m_R.dim2()};
    double sum{0.0};
        
    for (int i{0}; i < p; ++i) {
        if (piLegal(0,i) == 1) {
            m_pi(0,i) = rUniform();
            sum += m_pi(0,i);
        }
    }
        
    for (int i{0}; i < p; ++i) {
        m_pi(0,i) = m_pi(0,i) / sum;
    }
        
    for (int i{0}; i < p; ++i) {
        for (int j{0}; j < p; ++j) {
            if ((i != j) && (TLegal(i,j) == 1)) {
                m_T(i,j) = rUniform();
                m_T(i,i) -= m_T(i,j);
            }
        }
    }
        
    for (int i{0}; i < p; ++i) {
        if (TLegal(i,i) == 1) {
            m_t(i,0) = rUniform();
            m_T(i,i) -= m_t(i,0);
        }
    }
    
    for (int i{0}; i < p; ++i) {
        sum = 0;
        for (int j{0}; j < dim; ++j) {
            if (RLegal(i,j) == 1) {
                m_R(i,j) = rUniform();
                sum += m_R(i,j);
            }
        }
        for (int j{0}; j < dim; ++j) {
            m_R(i,j) = m_R(i,j) / sum;
        }
    }
    m_T *= (p / scaleFactor);
    m_t *= (p / scaleFactor);
        
}


double MPH::moment(int k, int j) {
    std::vector<int> states;
    
    return (marginal(j, states).moment(k));
    
}

double MPH::mean(int j) {
    return moment(1, j);
}

std::vector<double> MPH::mean() {
    std::vector<double> theMean;
    for (int j{0}; j < m_d; ++j) {
        theMean.push_back(moment(1, j));
    }
    return theMean;
}

double MPH::var(int j) {
    return moment(2, j) - pow(moment(1, j), 2.0);
}

std::vector<double> MPH::var() {
    std::vector<double> theVar;
    for (int j{0}; j < m_d; ++j) {
        theVar.push_back(moment(2, j) - pow(moment(1, j), 2.0));
    }
    return theVar;
}

double MPH::sd(int j) {
    return sqrt(var(j));
}

std::vector<double> MPH::sd() {
    std::vector<double> theSd{var()};
    for (int j{0}; j < m_d; ++j) {
        theSd[j] = sqrt(theSd[j]);
    }
    return theSd;
}

double MPH::crossMoment(int i, int j) {
    Numeric_lib::Matrix<double,2> U(matrixInverse(m_T * (-1.0)));
    double theCrossMoment{0.0};
    
    theCrossMoment=matrixProduct(m_pi, matrixProduct(U, matrixProduct(diagonalVector(columVector(m_R, i)), matrixProduct(U, columVector(m_R, j)))))(0,0) + matrixProduct(m_pi, matrixProduct(U, matrixProduct(diagonalVector(columVector(m_R, j)), matrixProduct(U, columVector(m_R, i)))))(0,0);
    
    return theCrossMoment;
}

Numeric_lib::Matrix<double,2> MPH::allCrossMoments() {
    Numeric_lib::Matrix<double,2> moments(m_d,m_d);
    for (int i{0}; i < m_d; ++i) {
        for (int j{0}; j < m_d; ++j) {
            moments(i,j) = crossMoment(i, j);
        }
    }
    return moments;
}


double MPH::covariance(int i, int j) {
    return (crossMoment(i, j) - mean(i) * mean(j));
}

Numeric_lib::Matrix<double,2> MPH::covarianceMatrix() {
    Numeric_lib::Matrix<double,2> cov(m_d,m_d);
    for (int i{0}; i < m_d; ++i) {
        for (int j{0}; j < m_d; ++j) {
            cov(i,j) = covariance(i, j);
        }
    }
    return cov;
}

double MPH::correlation(int i, int j) {
    return (covariance(i, j) / (sd(i) * sd(j)));
}


Numeric_lib::Matrix<double,2> MPH::correlationMatrix() {
    Numeric_lib::Matrix<double,2> correlations(m_d,m_d);
    for (int i{0}; i < m_d; ++i) {
        for (int j{0}; j < m_d; ++j) {
            correlations(i,j) = correlation(i, j);
        }
    }
    return correlations;
}

double MPH::quantile(int j, double u, double epsilon){
    std::vector<int> states;
    
    return (marginal(j, states).quantile(u, epsilon));
}





/* ******************************
       Bivariate PH class
****************************** */

BivariatePH::BivariatePH(const Numeric_lib::Matrix<double,2> & alpha, const Numeric_lib::Matrix<double,2> & T11, const Numeric_lib::Matrix<double,2> & T12, const Numeric_lib::Matrix<double,2> & T22) : m_alpha(alpha), m_T11(T11), m_T12(T12), m_T22(T22), m_p1(T11.dim1()), m_p2(T22.dim1()), m_e(T22.dim1(),1), m_t(T22.dim1(),1), m_p(T11.dim1() + T22.dim1()), m_R(T11.dim1() + T22.dim1(),2), m_T(T11.dim1() + T22.dim1(),T11.dim1() + T22.dim1()), m_pi(1,T11.dim1() + T22.dim1()) {
    m_e += 1;
    matrixProduct(m_T22 * (-1.0), m_e, m_t);
    
    //Put proper values in T
    for (int i{0}; i < m_p; ++i) {
        for (int j{0}; j < m_p; ++j) {
            if( i < m_p1 && j < m_p1) {
                m_T(i,j) = m_T11(i,j);
            }
            else if (i >= m_p1 && j < m_p1) {
                m_T(i,j) = 0;
            }
            else if (i < m_p1 && j >= m_p1) {
                m_T(i,j) = m_T12(i,j - m_p1);
            }
            else {
                m_T(i,j) = m_T22(i - m_p1,j - m_p1);
            }
        }
        //Values of R
        if (i < m_p1) {
            m_R(i,0) = 1;
        }
        else {
            m_R(i,1) = 1;
        }
    }
    
    //Put proper values in pi
    for (int i{0}; i < m_p1; ++i) {
        m_pi(0,i) = m_alpha(0,i);
    }
    
}


BivariatePH::BivariatePH(long p1, long p2) : m_pi(1,p1 + p2), m_T(p1 + p2,p1 + p2), m_R(p1 + p2,2), m_p(p1 + p2), m_p1(p1), m_p2(p2), m_e(p2,1), m_t(p2,1), m_alpha(1,p1), m_T12(p1,p2), m_T22(p2,p2), m_T11(p1,p1) {
    m_e += 1;
    
    for (int i{0}; i < m_p; ++i) {
        //Values of R
        if (i < m_p1) {
            m_R(i,0) = 1;
        }
        else {
            m_R(i,1) = 1;
        }
    }
}


void BivariatePH::updatePi() {
    //Put proper values in pi
    for (int i{0}; i < m_p1; ++i) {
        m_pi(0,i) = m_alpha(0,i);
    }
}

void BivariatePH::updateT() {
    for (int i{0}; i < m_p; ++i) {
        for (int j{0}; j < m_p; ++j) {
            if( i < m_p1 && j < m_p1) {
                m_T(i,j) = m_T11(i,j);
            }
            else if (i >= m_p1 && j < m_p1) {
                m_T(i,j) = 0;
            }
            else if (i < m_p1 && j >= m_p1) {
                m_T(i,j) = m_T12(i,j - m_p1);
            }
            else {
                m_T(i,j) = m_T22(i - m_p1,j - m_p1);
            }
        }
    }

}

void BivariatePH::updateR() {
    for (int i{0}; i < m_p; ++i) {
        //Values of R
        if (i < m_p1) {
            m_R(i,0) = 1;
        }
        else {
            m_R(i,1) = 1;
        }
    }
}

void BivariatePH::readParametersFromFile(std::ifstream & inputFile) {
    double helpvar{};
    
    for (int i{0}; i < m_p1; ++i) {
        inputFile >> m_alpha(0, i);
        for (int j{0}; j < m_p1; ++j) {
            inputFile >> m_T11(i,j);
        }
        for (int j{0}; j < m_p2; ++j) {
            inputFile >> m_T12(i,j);
        }
        //El archivo implicitamente tiene R para sacar las dimensiones p1 y p2, por eso hay que leer R. Esto para mantener el mismo formato con el MPH
        inputFile >> helpvar;
        inputFile >> helpvar;
    }
    
    for (int i{0}; i < m_p2; ++i) {
        inputFile >> helpvar;
        if (helpvar != 0) {
            std::cout << "Warning: the extructure of alpha is not correct";
        }
        for (int j{0}; j < m_p1; ++j) {
            inputFile >> helpvar;
            if (helpvar != 0) {
                std::cout << "Warning: the extructure of T is not correct";
            }
        }
        for (int j{0}; j < m_p2; ++j) {
            inputFile >> m_T22(i,j);
        }
        //For R
        inputFile >> helpvar;
        inputFile >> helpvar;
    }
    
    matrixProduct(m_T22 * (-1.0), m_e, m_t);
    updatePi();
    updateT();
}

//  Printing Data
//      Prints Pi - Initial distribution
void BivariatePH::printPi() {
     printMatrix(m_pi);
 }
//      Prints Pi - Initial distribution
void BivariatePH::printAlpha() {
     printMatrix(m_alpha);
 }
//      Prints T - intensity matrix
 void BivariatePH::printT11() {
     printMatrix(m_T11);
 }
//      Prints T - intensity matrix
void BivariatePH::printT12() {
    printMatrix(m_T12);
}
//      Prints T - intensity matrix
void BivariatePH::printT22() {
    printMatrix(m_T22);
}
//      Prints T - intensity matrix
void BivariatePH::printT() {
    printMatrix(m_T);
}
//
void BivariatePH::printR() {
    printMatrix(m_R);
}
//      Prints t - Exit rates
 void BivariatePH::printt() {
    printMatrix(m_t);
}
//      Prints t - vector of ones
 void BivariatePH::printe() {
    printMatrix(m_e);
}
// Prints pi,T and t
void BivariatePH::print() {
    std::cout << "Initial probabilities (Pi):\n";
    printPi();
    std::cout << "Transition matrix (T):\n";
    printT();
    std::cout << "Rewards matrix (R):\n";
    printR();
}


//  Change Data
void BivariatePH::changeAlpha(const Numeric_lib::Matrix<double,2> & alpha) {
    m_alpha = alpha;
    updatePi();
}
void BivariatePH::changeT(const Numeric_lib::Matrix<double,2> & T11, const Numeric_lib::Matrix<double,2> & T12, const Numeric_lib::Matrix<double,2> & T22) {
    
    m_T11 = T11;
    m_T12 = T12;
    m_T22 = T22;
    
    matrixProduct(m_T22 * (-1.0), m_e, m_t);
    
    updateT();
    
}
void BivariatePH::changeAll(const Numeric_lib::Matrix<double,2> & alpha, const Numeric_lib::Matrix<double,2> & T11, const Numeric_lib::Matrix<double,2> & T12, const Numeric_lib::Matrix<double,2> & T22) {
    m_alpha = alpha;
    
    m_T11 = T11;
    m_T12 = T12;
    m_T22 = T22;
    
    matrixProduct(m_T22 * (-1.0), m_e, m_t);
    
    updatePi();
    
    updateT();
}


//  Get data
Numeric_lib::Matrix<double,2> BivariatePH::getAlpha() {
    return m_alpha;
}
Numeric_lib::Matrix<double,2> BivariatePH::getPi() {
    return m_pi;
}
Numeric_lib::Matrix<double,2> BivariatePH::getT11() {
    return m_T11;
}
Numeric_lib::Matrix<double,2> BivariatePH::getT12() {
    return m_T12;
}
Numeric_lib::Matrix<double,2> BivariatePH::getT22() {
    return m_T22;
}
Numeric_lib::Matrix<double,2> BivariatePH::getT() {
    return m_T;
}
Numeric_lib::Matrix<double,2> BivariatePH::getR() {
    return m_R;
}
Numeric_lib::Matrix<double,2> BivariatePH::gett() {
    return m_t;
}
Numeric_lib::Matrix<double,2> BivariatePH::gete() {
    return m_e;
}
long BivariatePH::getp1() {
    return m_p1;
}
long BivariatePH::getp2() {
    return m_p2;
}
long BivariatePH::getp() {
    return m_p;
}


PhaseType BivariatePH::marginal1() {
    PhaseType theMarginal(m_alpha, m_T11);
    return theMarginal;
}

PhaseType BivariatePH::marginal2() {
    PhaseType theMarginal(matrixProduct(m_alpha, matrixProduct(matrixInverse(m_T11 * (-1.0)), m_T12)), m_T22);
    return theMarginal;
}

double BivariatePH::density(double x1,double x2) {
    double densityVal{};
    densityVal = matrixProduct(m_alpha, matrixProduct(matrixExponential(m_T11 * x1), matrixProduct(m_T12, matrixProduct(matrixExponential(m_T22 * x2), m_t))))(0,0);
    return densityVal;
}
double BivariatePH::jointTail(double x1,double x2) {
    double tailVal{};
    tailVal = matrixProduct(m_alpha, matrixProduct(matrixInverse(m_T11 * (-1.0)), matrixProduct(matrixExponential(m_T11 * x1), matrixProduct(m_T12, matrixProduct( matrixExponential(m_T22 * x2), m_e)))))(0,0);
    return tailVal;
}
double BivariatePH::moment(int k, int j) {
    if (j == 0) {
        return marginal1().moment(k);
    }
    else {
        return marginal2().moment(k);
    }
}
double BivariatePH::crossMoment() {
    Numeric_lib::Matrix<double,2> U(matrixInverse(m_T * (-1.0)));
    double theCrossMoment{0.0};
    
    theCrossMoment = matrixProduct(m_pi, matrixProduct(U, matrixProduct(diagonalVector(columVector(m_R, 0)), matrixProduct(U, columVector(m_R, 1)))))(0,0) + matrixProduct(m_pi, matrixProduct(U, matrixProduct(diagonalVector(columVector(m_R, 1)), matrixProduct(U, columVector(m_R, 0)))))(0,0);
    
    return theCrossMoment;
}
double BivariatePH::covariance() {
    return crossMoment() - mean1() * mean2();
}
double BivariatePH::correlation() {
    return covariance() / (sd1() * sd2());
}
double BivariatePH::mean1() {
    return moment(1, 0);
}
double BivariatePH::mean2() {
    return moment(1, 1);
}
std::vector<double> BivariatePH::mean() {
    std::vector<double> theMean{};
    theMean.push_back(mean1());
    theMean.push_back(mean2());
    return theMean;
}
double BivariatePH::var1() {
    return moment(2, 0) - pow(moment(1, 0), 2.0);
}
double BivariatePH::var2() {
    return moment(2, 1) - pow(moment(1, 1), 2.0);
}
std::vector<double> BivariatePH::var() {
    std::vector<double> theVar{};
    theVar.push_back(var1());
    theVar.push_back(var2());
    return theVar;
}
double BivariatePH::sd1() {
    return sqrt(var1());
}
double BivariatePH::sd2() {
    return sqrt(var2());
}
std::vector<double> BivariatePH::sd() {
    std::vector<double> theSd{var()};
    for (int j{0}; j < 2; ++j) {
        theSd[j] = sqrt(theSd[j]);
    }
    return theSd;
}


double BivariatePH::jointMGF(double theta1, double theta2) {
    std::vector<double>  w;
    w.push_back(theta1);
    w.push_back(theta2);
    
    std::vector<double> Rw(matrixVectorProduct(m_R, w));
    
    Numeric_lib::Matrix<double,2> diagonal(m_p,m_p);
    
    for(int i{0}; i < m_p; ++i) {
        diagonal(i,i) = -Rw[i];
    }
    return matrixProduct(m_pi, matrixProduct(matrixInverse(addMatrices(diagonal, m_T * (-1.0))), m_t))(0,0);
}


void BivariatePH::randomPhase(const Numeric_lib::Matrix<double,2> & piLegal, const Numeric_lib::Matrix<double,2> & T11Legal, const Numeric_lib::Matrix<double,2> & T12Legal, const Numeric_lib::Matrix<double,2> & T22Legal, double scale, AbsRNG & rUniform) {
    
    double sum{0.0};
    
    for (int i{0}; i < m_p1; ++i) {
        if (piLegal(0,i) == 1) {
            m_alpha(0,i) = rUniform();
            sum += m_alpha(0, i);
        }
    }
    m_alpha = m_alpha / sum;
    
    for (int i{0}; i < m_p1; ++i) {
        for (int j{0}; j < m_p1; ++j) {
            if ((i != j) && (T11Legal(i,j) == 1)) {
                m_T11(i,j) = rUniform();
                m_T11(i,i) -= m_T11(i,j);
            }
        }
        for (int j{0}; j < m_p2; ++j) {
            if (T12Legal(i,j) == 1) {
                m_T12(i,j)=rUniform();
                m_T11(i,i) -= m_T12(i,j);
            }
        }
    }
    
    for (int i{0}; i < m_p2; ++i) {
        for (int j{0}; j < m_p2; ++j)
            if ((i != j) && (T22Legal(i,j) == 1)) {
                m_T22(i,j)=rUniform();
                m_T22(i,i) -= m_T22(i,j);
            }
    }
    
    double r{};
    for (int i{0}; i < m_p2; ++i) {
        if (T22Legal(i,i) == 1) {
            r = rUniform();
            m_T22(i,i) -= r;
            m_t(i,0) = r;
        }
    }
    
    m_T11 *= scale;
    m_T12 *= scale;
    m_T22 *= scale;
    m_t *= scale;
    
    updatePi();
    updateT();
    
}


void BivariatePH::exportToFile(std::ofstream & outputFile) {
    for (int i{0}; i < m_p; ++i) {
        outputFile << m_pi(0,i) << '\t';
        
        for (int j{0}; j < m_p; ++j) {
            outputFile << m_T(i,j) << '\t';
        }
        
        outputFile << m_R(i,0) << '\t';
        outputFile << m_R(i,1);
        
        if (i < m_p - 1) {
            outputFile << '\n';
        }
    }
}
