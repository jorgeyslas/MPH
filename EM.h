//Implementation of the EM algorithm

#ifndef EM_h
#define EM_h

#include "Matrix.h"
#include "PhaseType.h"
#include "Distributions.h"
#include "MultivariateDistributions.h"
#include "ScaleComponent.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>

//Structure that keeps the dimensions of a multivariate PH - This is use on file reading, to return these values
struct DimensionsMPH {
    long m_p;
    long m_dim;
};

// Structure of a univariate sample: The value of the observation (obs) and the weight of it (weight)
struct Sample {
    double obs;
    double weight;
};


struct BivariateSample {
    double x1;
    double x2;
    double weight;
};


// Provides a way to compare the elements of the structure "Sample"
bool compare(Sample &a, Sample &b);
// Sort a vector of the type "Sample".
void sortObservations(std::vector<Sample> & vec);

//Prints a vector of doubles 
void printVector(const std::vector<double> & a);
void printVector(const std::vector<int> & a);


//EM for a univariate PH
//  Data reading
//      Read univariate data without weigths from file and put it in a vector of type "Sample"
void readDataForPH(std::ifstream & inputFile, std::vector<Sample> & observations, std::vector<Sample> & censored, int indWeighted, int indCensored);
//      Given a choice of density (See Distributions.h) put a sample in a vector of type "Sample" for approximation of a continous density.
void inputDensity(std::function<double(double, std::vector<double>)> density, const std::vector<double> & parameters, double truncationPoint, double maxProbability, double maxDeltat, std::vector<Sample> & observations);
//      Gives basic information (size, sum of weigths, mean and sd) of a vector of type "Sample"
void sampleStatisticsPH(const std::vector<Sample> & observations);
void infoDataPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored);

//      Given a file with the structure of a MPH (0's and 1's) returns the dimension of the MPH and the size p of the PH component - For a file for a PH returns dim=0.
DimensionsMPH dimensionsOfFile(std::ifstream & inputFile);
//      Given a file with the structure of a MPH (0's and 1's)  puts the corresponding values in PiLegal and TLegal
void readStructureForPH(std::ifstream & inputFile, Numeric_lib::Matrix<double,2> & piLegal ,Numeric_lib::Matrix<double,2> & TLegal); //por q no use const ref?


//  Interface - Option selection
//      Type of sample - From file or a density
int askForSampleType();
//      In case of a fit to a density provides the list of options of densities
void askContinousDist(std::vector<Sample> & observations);
//      Initial structure of the PH - General and random, Given by file and random or given by file
int askForStructure();
//      Ask if step lenght is given by the user or computed automatically
int askTypeOfStepRK();
//      Ask for step lenght in the Runge-Kutta procedure
double askStepLengthRK();


//  EM algorithm for a univariate PH
//      Runge-Kutta for the calculation of the vectors a,b and c of the paper of Asmussen
void rungeKutta(Numeric_lib::Matrix<double,2> & avector, Numeric_lib::Matrix<double,2> & bvector, Numeric_lib::Matrix<double,2> & cmatrix, double dt, double h, const Numeric_lib::Matrix<double,2> & T, const Numeric_lib::Matrix<double,2> & t);
//      Default step length for R-K
double defaultStepLength(const Numeric_lib::Matrix<double,2> & T);
//      EM step using Runge-Kutta
void EMstepRungeKutta(double h, PhaseType & phaseTypeEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored);
void a_rungekutta(Numeric_lib::Matrix<double,2> & avector, double dt, double h, const Numeric_lib::Matrix<double,2> & T);

double logLikelihoodPH_RK(double h, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, int transformation, std::vector<double> beta);

void EMIterate_RK(int stepsEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, int printLikehood, int transformation, std::vector<double> & beta, double lambda, double epsilon);

//      EM step using the matrix exponential algorithm of Matlab
void EMstep(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM);

// Loglikelihood using Matlab matrix exponential
double logLikelihoodPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, int transformation, std::vector<double> beta);

void EMIterate(int stepsEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, int printLikehood, int transformation, std::vector<double> & beta, double lambda, double epsilon);

//Uniformization
int findN(double epsilon, double lambda);
void vectorOfMatrices(std::vector<Numeric_lib::Matrix<double,2>> & theVector, const Numeric_lib::Matrix<double,2> & T, double alpha, int size);
void pow2Matrix(int n, Numeric_lib::Matrix<double,2> & T);
void matrixExpSum(double x, int n, const std::vector<Numeric_lib::Matrix<double,2>> & powerVector, double alpha, Numeric_lib::Matrix<double,2> & resultmatrix);
void EMstepUniformization(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, double epsilon);
double logLikelihoodPH_Uni( const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, double epsilon, int transformation, std::vector<double> beta);

void EMIterate_Uni(int stepsEM, double epsilon, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, int printLikehood, int transformation, std::vector<double> & beta, double lambda, double epsilonGD);


// NPH
void EMstepNPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored,  PhaseType & PH, ScaleDistribution *N, int fixParameter, int maxMethod, double lambda);
double logLikelihoodNPH(const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & PH, ScaleDistribution *N);
void EMIterateNPH(int stepsEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, ScaleDistribution *N, int fixParameter, int printLikehood, int maxMethod, double lambda);

//MAtrix Weibull test
void EMIterateMW(int stepsEM, const std::vector<Sample> & observations, const std::vector<Sample> & censored, PhaseType & phaseTypeEM, int printLikehood, double & beta, int maxMethod, double lambda);

/* ********
    MPH
******** */

void readStructureForMPH(std::ifstream & inputFile, Numeric_lib::Matrix<double,2> & piLegal , Numeric_lib::Matrix<double,2> & TLegal, Numeric_lib::Matrix<double,2> & RLegal);
long readDataForMPH(std::ifstream & inputFile, std::vector<std::vector<Sample>> & observations,std::vector<Sample> & sumObservations,  std::vector<std::vector<Sample>>  & censored, std::vector<Sample> & sumCensored, int indCensored, int transformation);
void sampleStatisticsMPH(const std::vector<std::vector<Sample>> & observations);
void infoDataMPH(const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored);

void secondEMstep(const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored, MPH & mph);
void secondEMstepRungeKutta(int typeOfStep, double h, const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored, MPH & mph);
void secondEMstepUniformization(const std::vector<std::vector<Sample>> & observations, const std::vector<std::vector<Sample>> & censored, MPH & mph, double epsilon);

void EMIterateMPH(int stepsFirstEM, int stepsSecondEM, const std::vector<std::vector<Sample>> & observations, const std::vector<Sample> & sumObservations, const std::vector<std::vector<Sample>> & censored, const std::vector<Sample> & sumCensored, MPH & mph, PhaseType & PHsumOfMarginals, int printLikehood);


void EMIterateMPH_RK(int stepsFirstEM, int stepsSecondEM, const std::vector<std::vector<Sample>> & observations,const std::vector<Sample> & sumObservations , const std::vector<std::vector<Sample>> & censored, const std::vector<Sample> & sumCensored, MPH & mph,  PhaseType & PHsumOfMarginals, int printLikehood);


void EMIterateMPH_Uni(int stepsFirstEM, int stepsSecondEM, double epsilon, const std::vector<std::vector<Sample>> & observations,const std::vector<Sample> & sumObservations , const std::vector<std::vector<Sample>> & censored, const std::vector<Sample> & sumCensored, MPH & mph,  PhaseType & PHsumOfMarginals, int printLikehood);


void inputMultivariateDensity(MultivariateDistribution *theMulDist_prt, double truncationPoint1, double truncationPoint2, double maxProbability, double maxDeltat, std::vector<std::vector<Sample>> & observations,  std::vector<Sample> & sumObservations);

void askMultivariateDist(std::vector<std::vector<Sample>> & observations,  std::vector<Sample> & sumObservations);

void sampleStatisticsMultDist(const std::vector<std::vector<Sample>> & observations);


// Bivariate case

void readDataForBivPH(std::ifstream & inputFile, std::vector<BivariateSample> & observations, int transformation);

void sampleStatisticsBivPH(const std::vector<BivariateSample> & observations);

void readStructureForBivPH(std::ifstream & inputFile, Numeric_lib::Matrix<double,2> & alphaLegal, Numeric_lib::Matrix<double,2> & T11Legal, Numeric_lib::Matrix<double,2> & T12Legal, Numeric_lib::Matrix<double,2> & T22Legal);

DimensionsMPH dimensionsOfBivFile(std::ifstream & inputFile);


void EMstepForBiv(const std::vector<BivariateSample> & observations, BivariatePH & bivPH);

void EMIterateBivPH(int stepsEM, const std::vector<BivariateSample> & observations, BivariatePH & bivPH , int printLikehood, int transformation, double & beta1, double  & beta2, double lambda, double epsilon);

void inputBivariateDensity(MultivariateDistribution *theMulDist_prt, double truncationPoint1, double truncationPoint2, double maxProbability, double maxDeltat1, double maxDeltat2, std::vector<BivariateSample> & observations);

void askBivariateDist(std::vector<BivariateSample> & observations);


// Weibull test

void EMIterateBivMW(int stepsEM, const std::vector<BivariateSample> & observations, BivariatePH & bivPH, int printLikehood, double & beta1, double & beta2, double lambda);

#endif /* EM_h */
