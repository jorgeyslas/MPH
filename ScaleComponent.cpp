//

#include "ScaleComponent.h"

RiemannZeta::RiemannZeta(double theta, double epsilon) : m_epsilon(epsilon) {
    m_id = 1;
    
    updateParameter(theta);
    
};

double RiemannZeta::si (int i) {
    return i;
}

double RiemannZeta::density(int i) {
    return pow(i, -m_theta) / m_zeta;
}

double RiemannZeta::firstDerivative(int i) {
    return m_zeta_fd;
}

double RiemannZeta::secondDerivative(int i) {
    return m_zeta_sd;
}

double RiemannZeta::getParameter() {
    return m_theta;
}

int RiemannZeta::getid() {
    return m_id;
}

double RiemannZeta::auxiliarFn() {
    return m_zeta;
}

void RiemannZeta::updateParameter(double parameter) {
    
    m_theta = parameter;
    
    int minIterations{50};
    int i{1};
    double sum{1};
    double sum_ant{};
    double sumfd{};
    double sumfd_ant{};
    double sumsd{};
    double sumsd_ant{};
    double relDif{};
        
    do {
        ++i;
        sum_ant = sum;
        sumfd_ant = sumfd;
        sumsd_ant = sumsd;
        sum += 1 / pow(i, m_theta);
        sumfd -= log(i) / pow(i, m_theta); // First derivative
        sumsd += log(i) * log(i) / pow(i, m_theta); // Second derivative
        relDif = fmax(abs((sum - sum_ant) / sum_ant), fmax(abs((sumfd - sumfd_ant) / sumfd_ant), abs((sumsd - sumsd_ant) / sumsd_ant)));
    } while (relDif > m_epsilon || i < minIterations);
       
    m_zeta = sum;
    m_zeta_fd = sumfd;
    m_zeta_sd = sumsd;
}

void RiemannZeta::print() {
    std::cout << "\nRiemann Zeta with parameter " << m_theta << '\n';
}

void RiemannZeta::exportToFile(std::ofstream & outputFile) {
    int initial{1};
    int increment{1};
    
    outputFile << m_id << '\n';
    outputFile << m_theta << '\n';
    outputFile << initial << '\n';
    outputFile << increment;
}



// Pareto

DiscretePareto::DiscretePareto(double theta, int type, double initial, double ratio) : m_theta(theta), m_progressionType(type), m_initial(initial), m_ratio(ratio) {
    
    if (type == 1) { //Arithmetic
        m_id = 3;
    } else { //Geometric
        m_id = 2;
    }
    
}
double DiscretePareto::si (int i) {
    double si{};
    if (m_progressionType == 1) { //Arithmetic
        si = m_initial + (i - 1) * m_ratio;
    }
    else {
        si = m_initial * pow(m_ratio, (i - 1));
    }
    return si;
}

double DiscretePareto::density(int i) {
    return pow(si(i), -m_theta) - pow(si(i + 1), -m_theta);
}

double DiscretePareto::firstDerivative(int i) {
    return -log(si(i)) * pow(si(i), -m_theta) + log(si(i + 1)) * pow(si(i + 1), -m_theta);
}

double DiscretePareto::secondDerivative(int i) {
    return log(si(i)) * log(si(i)) * pow(si(i), -m_theta)  - log(si(i + 1)) * log(si(i + 1)) * pow(si(i + 1), -m_theta);
}

double DiscretePareto::getParameter() {
    return m_theta;
}

void DiscretePareto::updateParameter(double parameter) {
    m_theta = parameter;
}

double DiscretePareto::auxiliarFn() {
    return m_ratio;
}

int DiscretePareto::getid() {
    return m_id;
}

void DiscretePareto::print() {
    if (m_id == 3) {
        std::cout << "\nDiscrete Pareto on arithmetic progresion with parameter " << m_theta << ", initial value " << m_initial << " and increment " << m_ratio << '\n';
    }
    else {
        std::cout << "\nDiscrete Pareto on geometric progresion with parameter " << m_theta << ", initial value " << m_initial << " and ratio " << m_ratio << '\n';
    }
}

void DiscretePareto::exportToFile(std::ofstream & outputFile) {
    outputFile << m_id << '\n';
    outputFile << m_theta << '\n';
    outputFile << m_initial << '\n';
    outputFile << m_ratio;
}



//Weibull

DiscreteWeibull::DiscreteWeibull(double theta, int type, double initial, double ratio) : m_theta(theta), m_progressionType(type), m_initial(initial), m_ratio(ratio) {
    
    if (type == 1) { //Arithmetic
        m_id = 5;
    } else { //Geometric
        m_id = 4;
    }
    
}

double DiscreteWeibull::si (int i) {
    double si{};
    
    if (i == 0) {
        return 0;
    }
    
    if (m_progressionType == 1) { //Arithmetic
        si = m_initial + (i - 1) * m_ratio;
    }
    else {
        si = m_initial * pow(m_ratio, (i - 1));
    }
    return si;
}

double DiscreteWeibull::density(int i) {
    return exp(-pow(si(i - 1), m_theta)) - exp(-pow(si(i), m_theta));
}

double DiscreteWeibull::firstDerivative(int i) {
    if (i == 1) {
        return log(si(i)) * pow(si(i), m_theta) * exp(-pow(si(i), m_theta));
    }
    else {
        return -log(si(i - 1)) * pow(si(i - 1), m_theta) * exp(-pow(si(i - 1), m_theta)) + log(si(i)) * pow(si(i), m_theta) * exp(-pow(si(i), m_theta));
    }
}

double DiscreteWeibull::secondDerivative(int i) {
    if (i == 1) {
        return - log(si(i)) * log(si(i)) * pow(si(i), m_theta) * exp(-pow(si(i), m_theta)) * (pow(si(i), m_theta) - 1);
    }
    else {
        return  log(si(i - 1)) * log(si(i - 1)) * pow(si(i - 1), m_theta) * exp(-pow(si(i - 1), m_theta)) * (pow(si(i - 1), m_theta) - 1) - log(si(i)) * log(si(i)) * pow(si(i), m_theta) * exp(-pow(si(i), m_theta)) * (pow(si(i), m_theta) - 1);
    }
}

double DiscreteWeibull::getParameter() {
    return m_theta;
}

void DiscreteWeibull::updateParameter(double parameter) {
    m_theta=parameter;
}

double DiscreteWeibull::auxiliarFn() {
    return 0;
}

int DiscreteWeibull::getid() {
    return m_id;
}


void DiscreteWeibull::print() {
    if (m_id == 5) {
        std::cout << "\nDiscrete Weibull on arithmetic progresion with parameter " << m_theta << ", initial value " << m_initial << " and increment " << m_ratio << '\n';
    }
    else {
        std::cout << "\nDiscrete Weibull on geometric progresion with parameter " << m_theta << ", initial value " << m_initial << " and ratio " << m_ratio << '\n';
    }
}

void DiscreteWeibull::exportToFile(std::ofstream & outputFile) {
    outputFile << m_id << '\n';
    outputFile << m_theta << '\n';
    outputFile << m_initial << '\n';
    outputFile << m_ratio;
}




//Lognormal

DiscreteLogNormal::DiscreteLogNormal(double theta, int type, double initial, double ratio) : m_theta(theta), m_progressionType(type), m_initial(initial), m_ratio(ratio) {
    
    if (type == 1) { //Arithmetic
        m_id = 7;
    } else { //Geometric
        m_id = 6;
    }
    
}

double DiscreteLogNormal::si (int i) {
    double si{};
    
    if (i == 0) {
        return 0;
    }
    
    if (m_progressionType == 1) { //Arithmetic
        si = m_initial + (i - 1) * m_ratio;
    }
    else {
        si = m_initial * pow(m_ratio, (i - 1));
    }
    return si;
}

double DiscreteLogNormal::density(int i) {
    if (i == 1) {
        return 1 - exp(-pow(log(si(i)), m_theta));
    }
    else {
        return exp(-pow(log(si(i - 1)), m_theta)) - exp(-pow(log(si(i)), m_theta));
    }
}

double DiscreteLogNormal::firstDerivative(int i) {
    if (i == 1) {
        return 0;
    }
    else {
        return 0;
    }
}

double DiscreteLogNormal::secondDerivative(int i) {
    if (i == 1) {
        return 0;
    }
    else {
        return  0;
    }
}

double DiscreteLogNormal::getParameter() {
    return m_theta;
}

void DiscreteLogNormal::updateParameter(double parameter) {
    m_theta=parameter;
}

double DiscreteLogNormal::auxiliarFn() {
    return 0;
}

int DiscreteLogNormal::getid() {
    return m_id;
}


void DiscreteLogNormal::print() {
    if (m_id == 7) {
        std::cout << "\nDiscrete Logrnomal on arithmetic progresion with parameter " << m_theta << ", initial value " << m_initial << " and increment " << m_ratio << '\n';
    }
    else {
        std::cout << "\nDiscrete Logrnomal on geometric progresion with parameter " << m_theta << ", initial value " << m_initial << " and ratio " << m_ratio << '\n';
    }
}

void DiscreteLogNormal::exportToFile(std::ofstream & outputFile) {
    outputFile << m_id << '\n';
    outputFile << m_theta << '\n';
    outputFile << m_initial << '\n';
    outputFile << m_ratio;
}
