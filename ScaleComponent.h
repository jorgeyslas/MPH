// Scale Distribution for the NPH

#ifndef ScaleComponent_h
#define ScaleComponent_h

#include <iostream>
#include <fstream>
#include <cmath>

class ScaleDistribution {
public:
    virtual double si (int i) = 0;
    virtual double density(int i) = 0;
    virtual double firstDerivative(int i) = 0;
    virtual double secondDerivative(int i) = 0;
    virtual double getParameter() = 0;
    virtual void updateParameter(double parameter) = 0;
    virtual int getid() = 0;
    virtual double auxiliarFn() = 0;
    virtual void print() = 0;
    virtual void exportToFile(std::ofstream & outputFile) = 0;
    virtual ~ScaleDistribution(){}
};


class RiemannZeta: public ScaleDistribution {
private:
    int m_id;
    double m_theta;
    double m_epsilon;
    double m_zeta;
    double m_zeta_fd;
    double m_zeta_sd;
public:
    RiemannZeta(double theta, double epsilon);
    virtual double si(int i);
    virtual double density(int i);
    virtual double firstDerivative(int i);
    virtual double secondDerivative(int i);
    virtual double getParameter();
    virtual void updateParameter(double parameter);
    virtual double auxiliarFn();
    virtual int getid();
    virtual void print();
    virtual void exportToFile(std::ofstream & outputFile);
    virtual ~RiemannZeta(){}
};



class DiscretePareto: public ScaleDistribution {
private:
    int m_id;
    double m_theta;
    int m_progressionType;
    double m_initial;
    double m_ratio;
public:
    DiscretePareto(double theta, int type, double initial, double ratio);
    virtual double si (int i);
    virtual double density(int i);
    virtual double firstDerivative(int i);
    virtual double secondDerivative(int i);
    virtual double getParameter();
    virtual void updateParameter(double parameter);
    virtual double auxiliarFn();
    virtual int getid();
    virtual void print();
    virtual void exportToFile(std::ofstream & outputFile);
    virtual ~DiscretePareto(){}
};


class DiscreteWeibull: public ScaleDistribution {
private:
    int m_id;
    double m_theta;
    int m_progressionType;
    double m_initial;
    double m_ratio;
public:
    DiscreteWeibull(double theta, int type, double initial, double ratio);
    virtual double si (int i);
    virtual double density(int i);
    virtual double firstDerivative(int i);
    virtual double secondDerivative(int i);
    virtual double getParameter();
    virtual void updateParameter(double parameter);
    virtual double auxiliarFn();
    virtual int getid();
    virtual void print();
    virtual void exportToFile(std::ofstream & outputFile);
    virtual ~DiscreteWeibull(){}
};


class DiscreteLogNormal: public ScaleDistribution {
private:
    int m_id;
    double m_theta;
    int m_progressionType;
    double m_initial;
    double m_ratio;
public:
    DiscreteLogNormal(double theta, int type, double initial, double ratio);
    virtual double si (int i);
    virtual double density(int i);
    virtual double firstDerivative(int i);
    virtual double secondDerivative(int i);
    virtual double getParameter();
    virtual void updateParameter(double parameter);
    virtual double auxiliarFn();
    virtual int getid();
    virtual void print();
    virtual void exportToFile(std::ofstream & outputFile);
    virtual ~DiscreteLogNormal(){}
};



#endif /* ScaleComponent_h */
