// Class for phase-type distributions

#ifndef PhaseType_h
#define PhaseType_h


#include "Matrix.h"
#include "MatrixOperations.h"
#include "RandomNumbers.h"

#include <fstream>
#include <cmath>


class PhaseType {
    long m_p; //dimension of the phase-type
    Numeric_lib::Matrix<double,2> m_pi; //Initial distribution
    Numeric_lib::Matrix<double,2> m_T; //Intensity matrix
    Numeric_lib::Matrix<double,2> m_e; //Vector of ones
    Numeric_lib::Matrix<double,2> m_t; //Exit rates
    
public:
    //Constructors
    //  With pi and T given
    PhaseType(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T);
    //  Initalizated on zeros with dimension p
    PhaseType(long p);
    
    
    //Read parameters from file
    void readParametersFromFile(std::ifstream & inputFile);
    
    //Print data
    void printPi();
    void printT();
    void printt();
    void printe();
    void print();
    
    //Get data
    Numeric_lib::Matrix<double,2> getPi();
    Numeric_lib::Matrix<double,2> getT();
    Numeric_lib::Matrix<double,2> gett();
    Numeric_lib::Matrix<double,2> gete();
    long getp();
    
    
    //Modify data
    void changePi(const Numeric_lib::Matrix<double,2> & pi);
    void changeT(const Numeric_lib::Matrix<double,2> & T);
    void changeAll(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T);
    
    //Properties
    double density(double x);
    double tail(double x);
    double distribution(double x);
    double moment(int k);
    double mean();
    double var();
    double sd();
    double LaplaceTransform(double s);
    double quantile(double u, double epsilon);
    double densityLevel(double x, double si); //For NPH
    
    //Functions for the EM
    //  Random initialitation based on the structure of Pilegal and TLegal
    void randomPhase(const Numeric_lib::Matrix<double,2> & piLegal, const Numeric_lib::Matrix<double,2> & TLegal, double scale, AbsRNG & rUniform);
    
    // Writes the parameter in a file
    void exportToFile(std::ofstream & outputFile);

    
};




class MPH {
    long m_p; //Number of transition
    long m_d; //dimension of the vector
    Numeric_lib::Matrix<double,2> m_pi; //Initial distribution
    Numeric_lib::Matrix<double,2> m_T; //Intensity matrix
    Numeric_lib::Matrix<double,2> m_R; //Rewards matrix
    Numeric_lib::Matrix<double,2> m_e; //Vector of ones
    Numeric_lib::Matrix<double,2> m_t; //Exit rates
    
public:
    //Constructors
    //  With pi, T and R given
    MPH(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T,  const Numeric_lib::Matrix<double,2> & R);
    //  Initalizated on zeros with dimensiona p and d
    MPH(long p, long d);
    
    //Read parameters from file
    void readParametersFromFile(std::ifstream & inputFile);
    
    //Print data
    void printPi();
    void printT();
    void printR();
    void printt();
    void printe();
    void print();
    
    //Get data
    Numeric_lib::Matrix<double,2> getPi();
    Numeric_lib::Matrix<double,2> getT();
    Numeric_lib::Matrix<double,2> getR();
    Numeric_lib::Matrix<double,2> gett();
    Numeric_lib::Matrix<double,2> gete();
    long getp();
    long getd();
    
    
    //Modify data
    void changePi(const Numeric_lib::Matrix<double,2> & pi);
    void changeT(const Numeric_lib::Matrix<double,2> & T);
    void changeR(const Numeric_lib::Matrix<double,2> & R);
    void changeAll(const Numeric_lib::Matrix<double,2> & pi, const Numeric_lib::Matrix<double,2> & T, const Numeric_lib::Matrix<double,2> & R);
    
    //Properties
    PhaseType linearCombination(const std::vector<double> & w, std::vector<int> & newStates );
    PhaseType marginal(int j, std::vector<int> & newStates);
    
    //double density(double x1,double x2 ); //Maybe only for the Bivariate
    double moment(int k, int j); //k-moment of the j marginal
    double crossMoment(int i, int j); // E(XiXj)
    Numeric_lib::Matrix<double,2> allCrossMoments(); // E(XiXj)
    double covariance(int i, int j);
    Numeric_lib::Matrix<double,2> covarianceMatrix();
    double correlation(int i, int j);
    Numeric_lib::Matrix<double,2> correlationMatrix();
    double mean(int j);
    std::vector<double> mean(); //Mean vector
    double var(int j);
    std::vector<double> var(); //vector with variances
    double sd(int j);
    std::vector<double> sd(); //vector with sd's
    double jointMGF(const std::vector<double> & w);
    double quantile(int j, double u, double epsilon);
    
    //Functions for the EM
    //  Random initialitation based on the structure of Pilegal and TLegal
    void randomPhase(const Numeric_lib::Matrix<double,2> & piLegal, const Numeric_lib::Matrix<double,2> & TLegal, const Numeric_lib::Matrix<double,2> & RLegal, double scale, AbsRNG & rUniform);
    
    // Writes the parameter in a file
    void exportToFile(std::ofstream & outputFile);

};





class BivariatePH {
    long m_p1; //Number of transition
    long m_p2; //Number of transition
    long m_p; //Number of transition
    Numeric_lib::Matrix<double,2> m_alpha; //Initial distribution - dimension p1
    Numeric_lib::Matrix<double,2> m_pi; //Initial distribution - dimensio p (alpha, 0)
    Numeric_lib::Matrix<double,2> m_T11; //Intensity matrix
    Numeric_lib::Matrix<double,2> m_T12; //Intensity matrix
    Numeric_lib::Matrix<double,2> m_T22; //Intensity matrix
    Numeric_lib::Matrix<double,2> m_T; //Intensity matrix (T11 , T12 ; 0 , T22)
    Numeric_lib::Matrix<double,2> m_R; //Rewards matrix - fix form
    Numeric_lib::Matrix<double,2> m_e; //Vector of ones - of dimension p2 for the exit
    Numeric_lib::Matrix<double,2> m_t; //Exit rates - dimension p2
    
public:
    //Constructors
    //  With pi, T and R given
    BivariatePH(const Numeric_lib::Matrix<double,2> & alpha, const Numeric_lib::Matrix<double,2> & T11, const Numeric_lib::Matrix<double,2> & T12, const Numeric_lib::Matrix<double,2> & T22);
    //  Initalizated on zeros with dimensiona p and d
    BivariatePH(long p1, long p2);
    
    
    //Update pi and T given alpha, T11, T12 and T22
    void updatePi();
    void updateT();
    void updateR();
    
    //Read parameters from file
    void readParametersFromFile(std::ifstream & inputFile);
    
    //Print data
    void printAlpha();
    void printPi();
    void printT11();
    void printT12();
    void printT22();
    void printT();
    void printR();
    void printt();
    void printe();
    void print();
    
    //Get data
    Numeric_lib::Matrix<double,2> getAlpha();
    Numeric_lib::Matrix<double,2> getPi();
    Numeric_lib::Matrix<double,2> getT11();
    Numeric_lib::Matrix<double,2> getT12();
    Numeric_lib::Matrix<double,2> getT22();
    Numeric_lib::Matrix<double,2> getT();
    Numeric_lib::Matrix<double,2> getR();
    Numeric_lib::Matrix<double,2> gett();
    Numeric_lib::Matrix<double,2> gete();
    long getp();
    long getp1();
    long getp2();
    
    
    //Modify data
    void changeAlpha(const Numeric_lib::Matrix<double,2> & alpha);
    void changeT11(const Numeric_lib::Matrix<double,2> & T11);
    void changeT12(const Numeric_lib::Matrix<double,2> & T12);
    void changeT22(const Numeric_lib::Matrix<double,2> & T22);
    void changeT(const Numeric_lib::Matrix<double,2> & T11, const Numeric_lib::Matrix<double,2> & T12, const Numeric_lib::Matrix<double,2> & T22);
    void changeAll(const Numeric_lib::Matrix<double,2> & alpha, const Numeric_lib::Matrix<double,2> & T11, const Numeric_lib::Matrix<double,2> & T12, const Numeric_lib::Matrix<double,2> & T22);
    
    //Properties
    PhaseType marginal1();
    PhaseType marginal2();
    
    double density(double x1,double x2);
    double jointTail(double x1,double x2);
    double moment(int k, int j); //k-moment of the j marginal
    double crossMoment(); // E(XiXj)
    double covariance();
    double correlation();
    double mean1();
    double mean2();
    std::vector<double> mean(); //Mean vector
    double var1();
    double var2();
    std::vector<double> var(); //vector with variances
    double sd1();
    double sd2();
    std::vector<double> sd(); //vector with sd's
    double jointMGF(double theta1, double theta2);
    
    //Functions for the EM
    //  Random initialitation based on the structure of Pilegal and TLegal
    void randomPhase(const Numeric_lib::Matrix<double,2> & piLegal, const Numeric_lib::Matrix<double,2> & T11Legal, const Numeric_lib::Matrix<double,2> & T12Legal, const Numeric_lib::Matrix<double,2> & T22Legal, double scale, AbsRNG & rUniform);
    
    // Writes the parameter in a file
    void exportToFile(std::ofstream & outputFile);

};





#endif /* PhaseType_h */
