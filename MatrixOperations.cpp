//  Basic operations for objects of the Matrix library

#include "MatrixOperations.h"


double bigNumber() {
//    Returns a "big" real number.
    return 1.0E+30;
}

double log2(double x) {
//    Returns the logarithm base 2 of the absolute value of a real number.
//    value = Log ( |X| ) / Log ( 2.0 )
    double value{0.0};

    if (x == 0.0) {
        value = - bigNumber();
    }
    else {
        value = log(abs(x)) / log(2.0);
    }

    return value;
}

int factorial(int n) {
// Computes n!
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void printMatrix(const Numeric_lib::Matrix<double,2> & A) {
//Prints a matrix of dimension 2
    for (int i{0}; i < A.dim1(); ++i) {
        for (int j{0}; j < A.dim2(); ++j) {
            std::cout << A(i,j) << '\t';
        }
        std::cout << '\n';
    }
}

double LInfNorm(const Numeric_lib::Matrix<double,2> & A) {
//    Returns the  L-oo norm of a matrix.
//    The matrix L-oo norm is defined as:
//      L-oo A =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
    double value{0.0};

    for (int i{0}; i < A.dim1(); ++i) {
        double row_sum{0.0};
        for (int j{0}; j < A.dim2(); ++j) {
            row_sum += abs(A(i,j));
        }
        value = std::max(value, row_sum);
  }
  return value;
}


void identityMatrix(Numeric_lib::Matrix<double,2> & A) {
//    Modifies a matrix and makes it an identity matrix.
    for (int i{0}; i < A.dim1(); ++i) {
        for (int j{0}; j < A.dim2(); ++j) {
            if (i == j) {
                A(i,j) = 1.0;
            }
            else {
                A(i,j) = 0.0;
            }
        }
    }
}

Numeric_lib::Matrix<double,2> identityMatrix(int n) {
//    Creates an identity matrix of size n.
    Numeric_lib::Matrix<double,2> A(n,n);
    
    identityMatrix(A);
    
    return A;
}


void addMatrices(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & C) {
//  Computes C =  A + B and modifies C by reference
    for (int i{0}; i < A.dim1(); ++i) {
        for (int j{0}; j < A.dim2(); ++j) {
            C(i,j) = A(i,j) + B(i,j);
        }
    }
}

Numeric_lib::Matrix<double,2> addMatrices(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B) {
//  Computes C =  A + B .
    Numeric_lib::Matrix<double,2> C(A.dim1(),A.dim2());
    
    addMatrices(A, B, C);

    return C;
}

//void matrixProduct (const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & C)
////    Multiplies two matrices AxB. Modifying C by reference
//{
//
//    for (int i = 0; i < A.dim1(); ++i )
//    {
//        for (int j = 0; j < B.dim2(); ++j )
//        {
//            C(i,j) = 0.0;
//            for (int k = 0; k < A.dim2(); ++k )
//            {
//                C(i,j) +=  A(i,k) * B(k,j);
//            }
//        }
//    }
//
//}


std::vector<double> matrixVectorProduct(const Numeric_lib::Matrix<double,2> & A, const std::vector<double> & B) {
//    Multiplies two matrices AxB. Modifying C by reference
    std::vector<double> C;
    for (int i{0}; i < A.dim1(); ++i) {
        C.push_back(0.0);
        for (int j{0}; j < A.dim2(); ++j) {
            C[i] +=  A(i,j) * B[j];
        }
    }
    
    return C;
}

//Better version
void matrixProduct(const Numeric_lib::Matrix<double,2>& a, const Numeric_lib::Matrix<double,2>& b, Numeric_lib::Matrix<double,2>& c) {
    
    const size_t rows = a.dim1(), cols = b.dim2(), n = a.dim2();

    for (size_t i = 0; i < rows; ++i) {
        auto ci = c[i];
        for (size_t j = 0; j < cols; ++j) {
            ci[j] = 0;
        }
    }

    for (size_t i = 0; i < rows; ++i) {
        const auto ai = a[i];
        auto ci = c[i];

        for (size_t k = 0; k < n; ++k) {
            const auto bk = b[k];
            const auto aik = ai[k];

            for (size_t j = 0; j < cols; ++j) {
                ci[j] += aik * bk[j];
            }
        }
    }
}

Numeric_lib::Matrix<double,2>  matrixProduct(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B) {
//    Multiplies two matrices AxB.
    Numeric_lib::Matrix<double,2> C(A.dim1(),B.dim2());

    matrixProduct(A, B, C);

    return C;
}


void solveLinearSystem(Numeric_lib::Matrix<double,2> A, const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & X ) {
//    Solves a system with multiple right hand sides. Returns the result by reference
//      AX=B which can be decompose as LUX=B and finds X
//      When B is the identity matrix the solution is the inverse of A
    long i{};
    long ipiv{};
    long j{};
    long jcol{};
    double piv{};

    X = B;

    for (jcol = 1; jcol <= A.dim1(); ++jcol) {
    //  Find the maximum element in column I.
        piv = abs(A(jcol - 1,jcol - 1));
        ipiv = jcol;
        for (i = jcol + 1; i <= A.dim1(); ++i) {
            if (piv < abs(A(i - 1,jcol - 1))) {
                piv = abs(A(i - 1,jcol - 1));
                ipiv = i;
            }
        }

        if (piv == 0.0) {
            std::cerr << "\n";
            std::cerr << "R8MAT_FSS_NEW - Fatal error!\n";
            std::cerr << "  Zero pivot on step " << jcol << "\n";
            exit (1);
        }
    //  Switch rows JCOL and IPIV, and X.
        if (jcol != ipiv) {
            A.swap_rows(jcol - 1,ipiv - 1);
            X.swap_rows(jcol - 1,ipiv - 1);
        }
        
    //  Scale the pivot row.
        double t{};
        t = A(jcol - 1,jcol - 1);
        A(jcol - 1,jcol - 1) = 1.0;
        for (j = jcol+1; j <= A.dim1(); ++j) {
            A(jcol - 1,j - 1) /= t;
        }
        for (j = 0; j < B.dim2(); ++j) {
            X(jcol - 1,j) /= t;
        }
    //  Use the pivot row to eliminate lower entries in that column.
        for (i = jcol + 1; i <= A.dim1(); ++i) {
            if (A(i - 1,jcol - 1) != 0.0) {
                t = - A(i - 1,jcol - 1);
                A(i - 1,jcol - 1) = 0.0;
                for (j = jcol + 1; j <= A.dim1(); ++j) {
                    A(i - 1,j - 1) +=  t * A(jcol - 1,j - 1);
                }
                for (j = 0; j < B.dim2(); ++j) {
                    X(i - 1,j) +=  t * X(jcol - 1,j);
                }
            }
        }
    }

    //  Back solve.
    for (jcol = A.dim1(); 2 <= jcol; --jcol) {
        for (i = 1; i < jcol; ++i) {
          for (j = 0; j < B.dim2(); ++j) {
            X(i - 1,j) -= A(i - 1,jcol - 1) * X(jcol - 1,j);
          }
        }
    }

}

Numeric_lib::Matrix<double,2> solveLinearSystem(Numeric_lib::Matrix<double,2> a, const Numeric_lib::Matrix<double,2> & b) {
//    Solves a system with multiple right hand sides.
//      AX=B which can be decompose as LUX=B and finds X
//      When B is the identity matrix the solution is the inverse of A
    Numeric_lib::Matrix<double,2> x(b.dim1(), b.dim2());
    
    solveLinearSystem(a, b, x);

    return x;
}


Numeric_lib::Matrix<double,2>  matrixInverse(const Numeric_lib::Matrix<double,2> & A, const Numeric_lib::Matrix<double,2> & B) {
//  Computes inverse(A) * B.
//  Parameters:
//    Input, double A[N1*N1], B[N1*N2], the matrices.
//    Output, double C[N1*N2], the result, C = inverse(A) * B.
  return solveLinearSystem(A, B);
}

void  matrixInverseByRef(const Numeric_lib::Matrix<double,2> & A,const Numeric_lib::Matrix<double,2> & B, Numeric_lib::Matrix<double,2> & X) {
//    Computes inverse(A) * B and put the result in X
    solveLinearSystem(A, B, X);
}

Numeric_lib::Matrix<double,2>  matrixInverse(const Numeric_lib::Matrix<double,2> & A) {
//    Computes inverse(A)
    return solveLinearSystem(A, identityMatrix(static_cast<int>(A.dim1())));
}

void  matrixInverseByRef (const Numeric_lib::Matrix<double,2> & A, Numeric_lib::Matrix<double,2> & AInv) {
//    Computes inverse(A) and put it in C
    solveLinearSystem(A, identityMatrix(static_cast<int>(A.dim1())), AInv);
}

void matrixExponential (const Numeric_lib::Matrix<double,2> & A, Numeric_lib::Matrix<double,2> & ExpM ) {
//   MATLAB's built-in matrix exponential algorithm saving the object by reference
    const int q{6};
    
    Numeric_lib::Matrix<double,2> matrixAuxiliar(A.dim1(),A.dim2());
    
    double aNorm{LInfNorm(A)};

    int ee{static_cast<int>(log2(aNorm)) + 1};
  
    int s{std::max(0, ee + 1)};

    double t{1.0 / pow(2.0, s)};
    
    Numeric_lib::Matrix<double,2> x(A.dim1(),A.dim2());
    Numeric_lib::Matrix<double,2> a2(A.dim1(),A.dim2());
    
    a2 = A * t;
    x = a2;

    double c{0.5};
    
    identityMatrix(ExpM);
    
    addMatrices(ExpM, a2 * c, matrixAuxiliar);
    ExpM = matrixAuxiliar;

    Numeric_lib::Matrix<double,2> d(A.dim1(),A.dim2());
    
    identityMatrix(d);

    addMatrices(d, a2 * (-c), matrixAuxiliar);
    d = matrixAuxiliar;


    int p{1};

    for (int k{2}; k <= q; ++k) {
        c = c * static_cast<double>(q - k + 1) / static_cast<double>(k * (2 * q - k + 1));
        
        matrixProduct(a2, x, matrixAuxiliar);
        x = matrixAuxiliar;
        
        addMatrices(x * c, ExpM, matrixAuxiliar);
        ExpM = matrixAuxiliar;

        if (p) {
            addMatrices(x * c, d, matrixAuxiliar);
            d = matrixAuxiliar;
        }
        else {
            addMatrices(x * (-c), d, matrixAuxiliar);
            d = matrixAuxiliar;
        }
        p = !p;
    }
    //  E -> inverse(D) * E
    matrixInverseByRef(d, ExpM, matrixAuxiliar);
    ExpM = matrixAuxiliar;
    //  E -> E^(2*S)
      for (int k = 1; k <= s; ++k) {
          matrixProduct(ExpM, ExpM, matrixAuxiliar);
          ExpM = matrixAuxiliar;
      }

}

Numeric_lib::Matrix<double,2> matrixExponential(const Numeric_lib::Matrix<double,2> & A) {
//   MATLAB's built-in matrix exponential algorithm.
//  Parameters:
//    Input, double A[N*N], the matrix.
//    Output, the estimate for exp(A).
    Numeric_lib::Matrix<double,2> e(A.dim1(),A.dim2());

    matrixExponential(A, e);
    
    return e;
}

void cumulateMatrix(const Numeric_lib::Matrix<double,2> & A, Numeric_lib::Matrix<double,2> & cumulated) {
//  Creates a new matrix with entries the cummulated rows of A - Saving the result by reference
    long p1{A.dim1()};
    long p2{A.dim2()};

    for (int i{0}; i < p1; ++i) {
        for (int j{0}; j < p2; ++j) {
            if (j == 0){
                cumulated(i,j) = A(i,j);
            }
            else {
                cumulated(i,j) = cumulated(i,j-1) + A(i,j);
            }
        }
    }
}

Numeric_lib::Matrix<double,2> cumulateMatrix(const Numeric_lib::Matrix<double,2> & A) {
//  Creates a new matrix with entries the cummulated rows of A
    Numeric_lib::Matrix<double,2> cumulated(A.dim1(),A.dim2());
    
    cumulateMatrix(A, cumulated);
    
    return cumulated;
}

Numeric_lib::Matrix<double,2> matrixPower(int n, const Numeric_lib::Matrix<double,2> & A) {
//  Computes A^n
    if (n == 1) {
        return A;
    }
    else if (n == 2) {
        return matrixProduct(A, A);
    }
    Numeric_lib::Matrix<double,2> previousMatrix(matrixProduct(A, A));
    Numeric_lib::Matrix<double,2> newMatrix(matrixProduct(A, previousMatrix));
    for (int i{4}; i <= n; ++i) {
        previousMatrix = newMatrix;
        newMatrix = matrixProduct(A, previousMatrix);
    }
    return newMatrix;
}

void matrixForVanLoan(const Numeric_lib::Matrix<double,2> & A1, const Numeric_lib::Matrix<double,2> & A2, const Numeric_lib::Matrix<double,2> & B1, Numeric_lib::Matrix<double,2> & auxiliarMatrix) {
// Put the matrix  (A1,B1 ; 0, A2 ) in auxiliarMatrix
    long p1{A1.dim1()};
    long p2{A2.dim1()};
    long p{p1 + p2};
    
    for (int i{0}; i < p; ++i) {
        for (int j{0}; j < p; ++j) {
            if ( i < p1 && j < p1) {
                auxiliarMatrix(i,j) = A1(i,j);
            }
            else if (i >= p1 && j < p1) {
                auxiliarMatrix(i,j) = 0;
            }
            else if (i < p1 && j >= p1) {
                auxiliarMatrix(i,j) = B1(i,j - p1);
            }
            else {
                auxiliarMatrix(i,j) = A2(i - p1,j - p1);
            }
        }
    }
}


void VanLoanForIntegrals(double t, const Numeric_lib::Matrix<double,2> & A1, const Numeric_lib::Matrix<double,2> & A2, const Numeric_lib::Matrix<double,2> & B1, Numeric_lib::Matrix<double,2> & auxiliarMatrix) {
// Computes  exp(t (A1,B1 ; 0, A2 )) and save the result by reference
    matrixForVanLoan(A1, A2, B1, auxiliarMatrix);
    
    auxiliarMatrix = matrixExponential(auxiliarMatrix * t);
}

Numeric_lib::Matrix<double,2> VanLoanForIntegrals(double t, const Numeric_lib::Matrix<double,2> & A1, const Numeric_lib::Matrix<double,2> & A2, const Numeric_lib::Matrix<double,2> & B1) {
// Computes  exp(t (A1,B1 ; 0, A2 ))
    long p{A1.dim1() + A2.dim1()};
    
    Numeric_lib::Matrix<double,2> auxiliarMatrix(p, p);
    
    VanLoanForIntegrals(t, A1, A2, B1, auxiliarMatrix);
    
    return auxiliarMatrix;
}


double matrixMax(const Numeric_lib::Matrix<double,2> & A) {
// Returns the maximum value in a matrix
    double maximum{A(0,0)};
    for (int i{0}; i < A.dim1(); ++i) {
        for (int j{0}; j < A.dim2(); ++j) {
            if (A(i,j) > maximum) {
                maximum = A(i,j);
            }
        }
    }
    return maximum;
}


double matrixMaxDiagonal(const Numeric_lib::Matrix<double,2>& A) {
// Returns the maximum value in the diagonal of a matrix
    double maximum{A(0,0)};
    for (int i{0}; i < A.dim1(); ++i) {
        if (A(i,i) > maximum) {
            maximum = A(i,i);
        }
    }
    return maximum;
}


Numeric_lib::Matrix<double,2> columVector(const Numeric_lib::Matrix<double,2> & A, int j) {
    Numeric_lib::Matrix<double,2> colVec(A.dim1(),1);
    for (int i{0}; i < A.dim1(); ++i) {
        colVec(i,0) = A(i,j);
    }
    return colVec;
}



Numeric_lib::Matrix<double,2> diagonalVector(const Numeric_lib::Matrix<double,2> & vec) {
    Numeric_lib::Matrix<double,2> diagonalMatrix(vec.dim1(),vec.dim1());
    for (int i{0}; i < vec.dim1(); ++i) {
        diagonalMatrix(i,i) = vec(i,0);
    }
    return diagonalMatrix;
}

