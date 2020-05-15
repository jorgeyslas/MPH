// Simulation of phase-type distributions

#ifndef Simulation_h
#define Simulation_h

#include "PhaseType.h"
#include "Matrix.h"
#include "RandomNumbers.h"
#include "ScaleComponent.h"
#include "MatrixOperations.h"

#include<fstream>

struct DiscreteDistribution {
    double x;
    double cdf;
};

//General functions
//  Returns the transition probabilities of the embeded markov chain determined by the transition matrix
Numeric_lib::Matrix<double,2> EmbeddedMC(const Numeric_lib::Matrix<double,2> & T, const Numeric_lib::Matrix<double,2> & t);
//  Given the accumulated values of the initial probabilities (Pi) and a uniform value (u) returns the states that satisfies cumPi_(k-1)<u<cumPi_(k)
long initialState(const Numeric_lib::Matrix<double,2> & cumulatedPi, double u);
//  Given the values of a probabilities matrix (Q), a uniform value (u) and a previous state,returns the new states of the Markov chain
long newState(long previousState, const Numeric_lib::Matrix<double,2> & cumulatedEmbeddedMC, double u);

//Simulates an univariate phase-type distribution - Tranformation: 0. No 1.e^x-1 2.beta(e^x-1) 3. x^(1/beta)
void simulatePH(PhaseType & ph, int SampleSize, AbsRNG & rUniform, std::ofstream & outputFile, int transformation, const std::vector<double> & beta);

//Simulates a multivariate phase-type distribution - Tranformation: 0. No 1.e^x-1 2.beta(e^x-1) 3. x^(1/beta)
void simulateMPH(MPH & mph, int sampleSize, AbsRNG & rUniform, std::ofstream & outputFile, int transformation, std::vector<double> beta);


//Computes upper and lower tail dependence coefficient by simulation
double upperTailDependence(MPH & mph, int sampleSize, double u, double epsilon, AbsRNG & rUniform);
double lowerTailDependence(MPH & mph, int sampleSize, double u, double epsilon, AbsRNG & rUniform);

std::vector<double> upperTailDependence3d(MPH & mph, int sampleSize, double u, double epsilon, AbsRNG & rUniform);

//Simulates a NPH distribution
void simulateNPH(PhaseType & ph, ScaleDistribution *N, int sampleSize, AbsRNG & rUniform, std::ofstream & outputFile);

#endif /* Simulation_h */
