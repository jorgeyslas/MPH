//  Basic operations for objects of the Matrix library

#ifndef MatrixOperations_h
#define MatrixOperations_h

#include "Matrix.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

//General functions
//  Returns a big number
double bigNumber();
//  logarithm base 2
double log2(double x);
//  n!
int factorial(int n);

//Functions for matrices
//  Print a matrix
void printMatrix(const Numeric_lib::Matrix<double,2> & A);
//  Computes the L inf norm of a matrix
double LInfNorm(const Numeric_lib::Matrix<double,2> & A);
//  Returns the maximum entry of a matrix
double matrixMax(const  Numeric_lib::Matrix<double,2> & A);
//  Returns the maximum entry in the diagonal of a matrix
double matrixMaxDiagonal(const  Numeric_lib::Matrix<double,2> & A);

//  Makes C an identity matrix
void identityMatrix(Numeric_lib::Matrix<double,2> & A);
//  Creates an identity matrix of size n
Numeric_lib::Matrix<double,2> identityMatrix (int n);

// Computes A+B saving the result in C
void addMatrices(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & C);
// Computes A+B and returns the result
Numeric_lib::Matrix<double,2> addMatrices(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B);
// Computes A*B saving the result in C
void matrixProduct(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & C);
// Computes A*B and returns the result
Numeric_lib::Matrix<double,2> matrixProduct(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B);

// Finds X in the linear system AX=B and saves the result in X
void solveLinearSystem(Numeric_lib::Matrix<double,2> A, const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & X);
// Finds X in the linear system AX=B and returns the solution
Numeric_lib::Matrix<double,2> solveLinearSystem(Numeric_lib::Matrix<double,2> A, const Numeric_lib::Matrix<double,2> & B);

//  Computes inverse(A)*B and saves the result in X  - This is basically the solution of the linear system AX=B
void matrixInverseByRef(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & X);
//  Computes inverse(A)*B and returns the result
Numeric_lib::Matrix<double,2>  matrixInverse(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B);
//  Computes inverse(A) and saves the result in AInv  - This is the solution of the linear system AX=I, where I is an identity matrix
void  matrixInverseByRef(const Numeric_lib::Matrix<double,2> & A, Numeric_lib::Matrix<double,2> & Ainv);
//  Computes inverse(A) and returns the result
Numeric_lib::Matrix<double,2>  matrixInverse(const Numeric_lib::Matrix<double,2> & A);

//  Computes exp(A) and saves the result in ExpM
void matrixExponential(const Numeric_lib::Matrix<double,2> & A, Numeric_lib::Matrix<double,2> & ExpM);
//  Computes exp(A) and returns the result
Numeric_lib::Matrix<double,2> matrixExponential(const Numeric_lib::Matrix<double,2> & A);

//  Creates a new matrix with entries the cummulated rows of A and saves it in cumulated
void cumulateMatrix(const Numeric_lib::Matrix<double,2> & A, Numeric_lib::Matrix<double,2> & cumulated);
//  Returns a new matrix with entries the cummulated rows of A
Numeric_lib::Matrix<double,2> cumulateMatrix(const Numeric_lib::Matrix<double,2> & A);
//  Computes A^n
Numeric_lib::Matrix<double,2> matrixPower(int n, const Numeric_lib::Matrix<double,2> & A);

// Put the matrix  (A1,B1 ; 0, A2 ) in auxiliarMatrix
void matrixForVanLoan(const Numeric_lib::Matrix<double,2> & A1, const Numeric_lib::Matrix<double,2> & A2, const Numeric_lib::Matrix<double,2> & B1, Numeric_lib::Matrix<double,2> & auxiliarMatrix);
// Computes  exp(t (A1,B1 ; 0, A2 )) and save the result by reference
void VanLoanForIntegrals(double t, const Numeric_lib::Matrix<double,2> & A1, const Numeric_lib::Matrix<double,2> & A2, const Numeric_lib::Matrix<double,2> & B1, Numeric_lib::Matrix<double,2> & auxiliarMatrix);
// Computes  exp(t (A1,B1 ; 0, A2 )) and returns the result
Numeric_lib::Matrix<double,2> VanLoanForIntegrals(double t, const Numeric_lib::Matrix<double,2> & A1, const Numeric_lib::Matrix<double,2> & A2, const Numeric_lib::Matrix<double,2> & B1);


std::vector<double> matrixVectorProduct(const Numeric_lib::Matrix<double,2> & A, const std::vector<double> & B);
Numeric_lib::Matrix<double,2> columVector(const Numeric_lib::Matrix<double,2> & A, int j);
Numeric_lib::Matrix<double,2> diagonalVector(const Numeric_lib::Matrix<double,2> & vec);

#endif /* MatrixOperations_h */
