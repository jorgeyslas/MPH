// Interface with the user

#include "Interface.h"



// Initial
int TypeOfDistribution() {
    int option{};
    do {
        std::cout << "\nType of distribution:\n";
        std::cout << "      1. PH - Phase-type\n";
        std::cout << "      2. MPH - Multivariate phase-type\n";
        std::cout << "      3. IPH - Inhomogeneous phase-type\n";
        std::cout << "      4. IMPH - Inhomogeneous MPH\n";
        //std::cout << "      5. NPH - Inifite dimensional phase-type \n";
        //std::cout << "      6. MNPH - Multivariate inifite dimensional phase-type\n";
        std::cout << "Select(1-4): ";
        std::cin >> option;

        if (option < 1 || option > 4) {
            std::cout << "Please enter a valid option\n";
        }
        
    } while (option < 1 || option > 4);
    return option;
}


int TypeOfTransformationMult() {
    int option{};
    do {
        std::cout << "\nType of transformation:\n";
        std::cout << "      1. Matrix-Pareto - (e^X - 1)\n";
        std::cout << "      2. Parameter dependent Matrix-Pareto - beta(e^X - 1)\n";
        std::cout << "      3. Matrix-Weibull - X^(1/beta)\n";
        std::cout << "Select(1-3): ";
        std::cin >> option;

        if (option < 1 || option > 3) {
            std::cout << "Please enter a valid option\n";
        }
        
    } while (option < 1 || option > 3);
    return option;
}

int TypeOfTransformation() {
    int option{};
    do {
        std::cout << "\nType of transformation:\n";
        std::cout << "      1. Matrix-Pareto - (e^X - 1)\n";
        std::cout << "      2. Parameter dependent Matrix-Pareto - beta(e^X - 1)\n";
        std::cout << "      3. Matrix-Weibull - X^(1/beta)\n";
        std::cout << "      4. Matrix-Gompertz - log(beta * X + 1) / beta \n";
        std::cout << "      5. Exponential-PH - (mu - sigma * log(X))\n";
        std::cout << "      6. Shifted-power-PH - (mu + (sigma / xi) * (X^(-xi) - 1))\n";
        std::cout << "Select(1-5): ";
        std::cin >> option;

        if (option < 1 || option > 6) {
            std::cout << "Please enter a valid option\n";
        }
        
    } while (option < 1 || option > 6);
    return option;
}


void SimulationOption() {
    int option{TypeOfDistribution()};
    
    if (option == 1) {
        SimulationForPH();
    }
    else if (option == 2) {
        SimulationForMPH();
    }
    else if (option == 3) {
        int transformation{TypeOfTransformation()};
        SimulationForPH(transformation);
    }
    else if (option == 4) {
        int transformation{TypeOfTransformationMult()};
        SimulationForMPH(transformation);
    }
}


int TypeMultDst() {
    int type{};
    do {
        std::cout << "\nType of multivariate distribution:\n";
        std::cout << "      1. General\n";
        std::cout << "      2. Bivariate with explicit density\n";
        std::cout << "Select(1-2): ";
        std::cin >> type;
               
        if (type != 1 && type != 2) {
            std::cout << "Please enter a valid option\n";
        }
               
    } while (type != 1 && type != 2);
    return type;
}


void FittingOption() {
    int option{TypeOfDistribution()};
    
    if (option == 1) {
        EMforPH();
    }
    else if (option == 2) {
        int MPHType{TypeMultDst()};
        
        if (MPHType == 1) {
            EMforMPH();
        }
        else if (MPHType == 2) {
            EMforBivariatePH();
        }
    }
    else if (option == 3) {
        int transformation{TypeOfTransformation()};
        EMforPH(transformation);
    }
    else if (option == 4) {
        int transformation{TypeOfTransformationMult()};
        
        if (transformation == 1) {
            int MultType{TypeMultDst()};
            
            if (MultType == 1) {
                EMforMPH(transformation);
            }
            else if (MultType == 2) {
                EMforBivariatePH(transformation);
            }
        }
        else {
            EMforBivariatePH(transformation);
        }
    }
}

void option() {
    int option{};
    do {
        std::cout << "\nWhat do you want to do?:\n";
        std::cout << "      1. Fit distribution\n";
        std::cout << "      2. Simulate distribution\n";
        std::cout << "Select(1-2): ";
        std::cin >> option;
        
        if (option == 1) {
            FittingOption();
        }
        else if (option == 2) {
            SimulationOption();
        }
        else {
            std::cout << "Please enter a valid option\n";
        }
        
    } while (option != 1 && option != 2);
}


/* ************************************
    General Interface functions
************************************ */
int askNumberOfSimulations() {
    int numSimulations{};
    do {
        std::cout << "\nNumber of simulations: ";
        std::cin >> numSimulations; //Seed for simulation
        if (numSimulations <= 0) {
            std::cout << "Please enter a positive integer\n";
        }
    } while (numSimulations <= 0);
    return numSimulations;
}


int askIfCensored(){
    int censorIndicator{};
    do {
        std::cout << "\nThe sample contains:\n";
        std::cout << "   1. No censored observations \n";
        std::cout << "   2. Some right-censored observations\n";
        std::cout << "Select 1 or 2: ";
        std::cin >> censorIndicator;
        if (censorIndicator != 1 && censorIndicator != 2) {
            std::cout<<"Please enter a valid option\n";
        }
    } while (censorIndicator != 1 && censorIndicator != 2);
    return censorIndicator;
}

int askIfWeighted() {
    int weightOption{};
    do {
        std::cout << "\nType of sample:\n";
        std::cout << "      1. Unweighted\n";
        std::cout << "      2. Weighted\n";
        std::cout << "Select (1-2): ";
        std::cin >> weightOption;
        if (weightOption != 1 && weightOption != 2) {
            std::cout << "Please enter a valid option\n";
        }
    } while (weightOption != 1 && weightOption != 2);
    return weightOption;
}


int askNumberOfPhases(int indicator = 0) {
    int p{};
    do {
        if (indicator == 0) {
            std::cout << "Number of phases of the distribution to be fitted, (p): ";
        }
        else {
            std::cout << "Number of phases of the distribution to be fitted, (p" << indicator << "): ";
        }
        std::cin >> p;
        if (p <= 0) {
            std::cout << "Please enter a positive integer\n";
        }
    } while (p <= 0);
    return p;
}


int askNumberOfEMSteps() {
    int stepsEM{};
    do {
        std::cout << "\nNumber of steps for the EM: ";
        std::cin >> stepsEM;
        if (stepsEM <= 0) {
            std::cout << "\nThe number of steps has to be non-negative.";
        }
    } while (stepsEM <= 0);
    return stepsEM;
}

int askNumberOfEMStepsFirst() {
    int stepsEM{};
    do {
        std::cout << "\nNumber of steps for the first EM: ";
        std::cin >> stepsEM;
        if (stepsEM <= 0) {
            std::cout << "\nThe number of steps has to be non-negative.";
        }
    } while (stepsEM <= 0);
    return stepsEM;
}

int askNumberOfEMStepsSecond() {
    int stepsEM{};
    do {
        std::cout << "\nNumber of steps for the second EM: ";
        std::cin >> stepsEM;
        if (stepsEM <= 0) {
            std::cout << "\nThe number of steps has to be non-negative.";
        }
    } while (stepsEM <= 0);
    return stepsEM;
}


int askForMatrixExponentialMethod() {
    int emMethod{};
    do {
        std::cout << "\nMethod for matrix exponential:\n";
        std::cout << "      1. Pade approximation\n";
        std::cout << "      2. Runge-Kutta\n";
        std::cout << "      3. Uniformisation\n";
        std::cout << "Select(1-3): ";
        std::cin >> emMethod;
        if (emMethod < 1 || emMethod > 3) {
            std::cout << "Please enter a valid option\n";
        }
    } while (emMethod < 1 || emMethod > 3);
    return emMethod;
}


int askIfMoreIterations() {
    int stepsIndicator{};
    do {
        std::cout << "\nGenerate more iterations:\n";
        std::cout << "      1. Yes\n";
        std::cout << "      2. No\n";
        std::cout << "Select(1-2): ";
        std::cin >> stepsIndicator;
        if (stepsIndicator != 1 && stepsIndicator != 2) {
            std::cout << "Please enter a valid option\n";
        }
    } while (stepsIndicator != 1 && stepsIndicator != 2);
    return stepsIndicator;
}


int askIfSaveFit(){
    int saveInFile{};
    do {
        std::cout << "\nSave parameters in file:\n";
        std::cout << "      1. Yes\n";
        std::cout << "      2. No\n";
        std::cout << "Select(1-2): ";
        std::cin >> saveInFile;
        if (saveInFile != 1 && saveInFile != 2) {
            std::cout << "Please enter a valid option\n";
        }
    } while (saveInFile != 1 && saveInFile != 2);
    return  saveInFile;
}


int askSimulateFromFit() {
    int stepsIndicator{};
    do {
        std::cout << "\nSimulate from fitted distribution:\n";
        std::cout << "      1. Yes\n";
        std::cout << "      2. No\n";
        std::cout << "Select(1-2): ";
        std::cin >> stepsIndicator;
        if (stepsIndicator != 1 && stepsIndicator != 2) {
            std::cout << "Please enter a valid option\n";
        }
    } while (stepsIndicator != 1 && stepsIndicator != 2);
    return stepsIndicator;
}


int askMaxMethod() {
    int maxMethod{};
    do {
        std::cout << "\nMethod for maximization:\n";
        std::cout << "      1. Gradient ascent\n";
        std::cout << "      2. Newton–Raphson\n";
        std::cout << "Select(1-2): ";
        std::cin >> maxMethod;
        if (maxMethod != 1 && maxMethod != 2) {
            std::cout << "Please enter a valid option\n";
        }
    } while (maxMethod != 1 && maxMethod != 2);
    return  maxMethod;
}


double askLearningRate() {
    double lambda{};
    do {
        std::cout << "\nLearning rate of gradient ascent: ";
        std::cin >> lambda;
        if (lambda <= 0) {
            std::cout << "\nThe rate has to be positive.";
        }
    } while (lambda <= 0);
    return lambda;
}

/* *********
    PH
******** */

void SimulationForPH(int transformation) { // transformation = 1 applies e^x - 1 to the simulated path of the PH to generate a Matrix-Pareto
    // Input and output files
    const std::string phasesFile{"phases.txt"}; // Input file - parameters of the PH
    const std::string simulationFile{"simSample.txt"}; // Output file - simulations
    
    // Determinates the dimension of the distribution based on the input file
    std::ifstream inPhases{phasesFile};
    long p{};
    p = dimensionsOfFile(inPhases).m_p; // Subtract the implicit size of p
    std::cout << p << " phases in the file\n";
    
    // Parameters of the distribution
    PhaseType PH(p);
    PH.readParametersFromFile(inPhases);
    std::cout << "\nParameters of the distribution\n";
    PH.print(); // Prints the values
    
    inPhases.close();
    
    std::vector<double> beta{};
    double aux{};
    if (transformation == 2 || transformation == 3 || transformation == 4) {
        std::cout << "\nBeta: ";
        std::cin >> aux;
        beta.push_back(aux);
    }
    else if (transformation == 5) {
        std::cout << "\nMu: ";
        std::cin >> aux;
        beta.push_back(aux);
        std::cout << "\nSigma: ";
        std::cin >> aux;
        beta.push_back(aux);
    }
    else if (transformation == 6) {
        std::cout << "\nMu: ";
        std::cin >> aux;
        beta.push_back(aux);
        std::cout << "\nSigma: ";
        std::cin >> aux;
        beta.push_back(aux);
        std::cout << "\nXi: ";
        std::cin >> aux;
        beta.push_back(aux);
    }
    
    // Seed for simulation
    int seed{};
    std::cout << "\nSeed: ";
    std::cin >> seed;

    MTRNG randomSimulation(seed);
    
    int numSimulations{askNumberOfSimulations()};
    
    std::ofstream outFile{simulationFile};
    
    simulatePH(PH, numSimulations, randomSimulation, outFile, transformation, beta);
    
    outFile.close();
}

void EMforPH(int transformation) { // transformation 
    
    // Names of the I/O files - Change if needed
    //  Input files
    const std::string unweightedSampleFile{"sample.txt"}; // Unweighted uncensored Sample
    const std::string weightedSampleFile{"weighted.txt"}; // Weighted uncensored Sample
    const std::string unweightedCensoredFile{"censored.txt"}; // Unweighted uncensored Sample
    const std::string weightedCensoredFile{"weightedCensored.txt"}; // Weighted uncensored Sample
    const std::string structureFile{"legal.txt"}; // Structure of the PH
    const std::string initialPhasesFile{"phases.txt"}; //Initial values of the PH
    //  Output files
    const std::string fittedPhasesFile{"phasesFitted.txt"}; // Fitted values of the PH
    const std::string fittedBeta{"Beta.txt"}; // Fitted value of the parameter-dependent distribution
    
    // Constant that determinates when to print the Loglikelihood
    const int printLikehood{10};
    
    // Epsilon for the Uniformization method
    const double epsilon{1e-5};
    
    //Stopping criteria for gradient ascent
    //double epsilonGA{0.1};
    double epsilonGA{0.001};
    
    // Vector that will contain the uncensored sample
    std::vector<Sample> observations;
    // Vector that will contain the right-censored sample
    std::vector<Sample> rightCensored;
     
    
    int sampleType{askForSampleType()};
    
//    // Type of sample - (1)From file (2) generated from a continuous distribution - available only for a PH fit
//    int sampleType{};
//
//    if (transformation == 0) {
//        sampleType = askForSampleType();
//    } else{
//        sampleType = 1;
//    }
    
    
    if (sampleType == 1) { // File option
        // Ask if the data has censored values: (1) Uncensores (2) censored
        int censorIndicator{askIfCensored()};
        
        // Choice between file with weights (2) and without weights (1)
        int weightOption{askIfWeighted()};
       
        // Assign the corresponding file name
        std::string sampleFile;
        if (weightOption == 1) {
            sampleFile = (censorIndicator == 1) ? unweightedSampleFile : unweightedCensoredFile;
        }
        else {
            sampleFile = (censorIndicator == 1) ? weightedSampleFile : weightedCensoredFile;
        }
        std::ifstream infil{sampleFile};
        
        // Read data from file
        readDataForPH(infil, observations, rightCensored, weightOption, censorIndicator);
        
        // Disply basic information and sort data
        if (censorIndicator == 2) { // Right censored
            infoDataPH(observations, rightCensored);
            sortObservations(rightCensored);
        }
        else { // Uncensored
            sampleStatisticsPH(observations);   // Print basic information of the data
        }
        sortObservations(observations); // Sort data based on obs
        
        infil.close();
        
    }
    else if (sampleType == 2) { // Continuous distribution
        askContinousDist(observations); // Gives list of options and puts sample in vector
        sampleStatisticsPH(observations);  // Prints basic information - not sorted since it is created in this way
    }
    
     
    int structureType{askForStructure()}; // Type of initial structure for the PH (3 options): (1) General and random (2) From file and random (3) Initial from file
     
    long p{}; // Number of phases
    int seed{1};
     
    // Temporal PH object
    PhaseType *tempPH;
    
    // Determintes the seed (if needed) and value of p - Needed to initializates the PhaseType object
    if (structureType == 1) { // General and random
        // Get the number of phases
        p = askNumberOfPhases();
        
        std::cout << "\nSeed: ";
        std::cin >> seed; // Seed for random initialitation
        
        // Creates the temporal object
        tempPH = new PhaseType(p);
        
        // Puts initial values in temporal object
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> piLegal(1,p);
        piLegal += 1;
        Numeric_lib::Matrix<double,2> TLegal(p,p);
        TLegal += 1;
          
        double scale{1.0};
        scale = observations[(observations.size() - 1) / 2].obs * 15; // change if needed , can help on performance
          
        tempPH->randomPhase(piLegal, TLegal, scale, myRandom);
    }
    else if (structureType == 2) { // Structure from file and random
        std::ifstream inLegalFile{structureFile};
        p = dimensionsOfFile(inLegalFile).m_p; // Subtract the implicit size of p
        
        std::cout << p << " phases in the structure of the file";
         
        std::cout << "\nSeed: ";
        std::cin >> seed; // Seed for random initialitation
        
        // Creates the temporal object
        tempPH = new PhaseType(p);
          
        // Put initial values in temporal object
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> piLegal(1,p);
        Numeric_lib::Matrix<double,2> TLegal(p,p);
        
        readStructureForPH(inLegalFile, piLegal, TLegal);
          
        double scale{1.0};
        scale = observations[(observations.size() - 1) / 2].obs;
          
        tempPH->randomPhase(piLegal, TLegal, scale, myRandom);
        
        inLegalFile.close();
    }
    else { // Initial phases from file
        std::ifstream inPhases{initialPhasesFile};
        p = dimensionsOfFile(inPhases).m_p; // Subtract the implicit size of p
         
        std::cout << p << " phases in the structure of the file";
        
        // Creates temporal object
        tempPH = new PhaseType(p);
        
        // Put initial values in temporal object
        tempPH->readParametersFromFile(inPhases);
        
        inPhases.close();
    }
    
    // Initializates the phase-type object
    PhaseType PH{*tempPH};
    
    // Deletes temporal object
    delete tempPH;
    
     
    std::cout << "\nIntilian parameters of the distribution\n";
    PH.print(); // Prints initial values of the PH
    
    
    std::vector<double> beta{};
    double auxBeta{};
    double lambda{};
    
    if (transformation == 2 || transformation == 3 || transformation == 4) {
        std::cout << "\nInitial beta: ";
        std::cin >> auxBeta;
        beta.push_back(auxBeta);
        
        lambda = askLearningRate();
    }
    else if (transformation == 5) {
        std::cout << "\nInitial mu: ";
        std::cin >> auxBeta;
        beta.push_back(auxBeta);
        std::cout << "\nInitial sigma: ";
        std::cin >> auxBeta;
        beta.push_back(auxBeta);
        
        lambda = askLearningRate();
    }
    else if (transformation == 6) {
        std::cout << "\nInitial mu: ";
        std::cin >> auxBeta;
        beta.push_back(auxBeta);
        std::cout << "\nInitial sigma: ";
        std::cin >> auxBeta;
        beta.push_back(auxBeta);
        std::cout << "\nInitial xi: ";
        std::cin >> auxBeta;
        beta.push_back(auxBeta);
        
        std::cout << "\nmu - sigma / xi: " << beta[0] - beta[1] / beta[2];
        std::cout << "\nMinimum observation: " << observations[0].obs << '\n';
        
        lambda = askLearningRate();
    }
    
    
    // Selection of the method for the matrix exponential
    int emMethod{askForMatrixExponentialMethod()};
    
    
    // Indicates if one wants more iterations of the EM
    int stepsIndicator{1};
     
    // EM algorithm
    while (stepsIndicator == 1) {
        int stepsEM{askNumberOfEMSteps()};
            
        Timer myTime; // Object to measure running time
        
        if (emMethod == 1) { // Matlab algorithm
            EMIterate(stepsEM, observations, rightCensored, PH, printLikehood, transformation, beta, lambda, epsilonGA);
        }
        else if (emMethod == 2) { // Runge-Kutta
            EMIterate_RK(stepsEM, observations, rightCensored, PH, printLikehood, transformation, beta, lambda, epsilonGA);

        }
        else if (emMethod == 3) { // Uniformisation
            EMIterate_Uni(stepsEM, epsilon, observations, rightCensored, PH, printLikehood, transformation, beta, lambda, epsilonGA);
        }
        
        // Prints the parameters of the fitted distribution and some basic information to compare with the sample
        std::cout << "\nParameters of the fitted distribution\n";
        PH.print();
        
        if (transformation == 0) {
            std::cout << "\nExpected value of fitted distribution: " << PH.mean() << '\n';
            std::cout << "\nStandard deviation of fitted distribution: " << PH.sd() << '\n';
        }
        else if (transformation == 2 || transformation == 3 || transformation == 4) {
            std::cout << "\nBeta: " << beta[0] <<'\n';
        }
        else if (transformation == 5) {
            std::cout << "\nMu: " << beta[0] << ". Sigma: " << beta[1] <<'\n';
        }
        else if (transformation == 6) {
            std::cout << "\nMu: " << beta[0] << ". Sigma: " << beta[1] << ". Xi: " << beta[2] <<'\n';
        }
        
        
        // Running time
        std::cout << "\nTime elapsed: " << myTime.elapsed() << " seconds\n";
        
        // Ask if one wants more Iterations
        stepsIndicator = askIfMoreIterations();
         
    }
     
    // Option to save parameters of the fitted distribution
    int saveInFile{askIfSaveFit()};
    
    if (saveInFile == 1) {
        std::ofstream outFile{fittedPhasesFile};
        PH.exportToFile(outFile);
        
        outFile.close();
        
        if (transformation > 1) {
            std::ofstream outFileBeta{fittedBeta};
            outFileBeta << beta[0];
            
            if (transformation == 5) {
                outFileBeta << '\n' << beta[1];
            }
            else if (transformation == 6) {
                outFileBeta << '\n' << beta[1] << '\n' << beta[2];
            }
        
            outFileBeta.close();
        }
        
    }
}



/* *********
    MPH
******** */


void SimulationForMPH(int transformation) { // transformation = 1 applies e^x - 1 to the simulated sample of the MPH to generate a multivariate Matrix-Pareto
    
    // Names of Input and output files
    const std::string phasesFile{"MPHpar.txt"}; // Input file - parameters of the MPH
    const std::string simulationFile{"simSampleMPH.txt"}; // Output file - simulations
    
    // Determinates the number of phases and dimension of the distribution base on the input file
    std::ifstream inPhases{phasesFile};
    DimensionsMPH dimensions{dimensionsOfFile(inPhases)};
    long p{dimensions.m_p};
    long d{dimensions.m_dim};
    std::cout << p << " phases in the file and dimension"  << d << '\n';
    
    // Parameters of the distribution
    MPH mph(p,d);
    mph.readParametersFromFile(inPhases);
    std::cout << "\nParameters of the distribution\n";
    mph.print(); // Prints parameters of the MPH
    
    inPhases.close();
    
    std::vector<double> beta{};
    if (transformation == 2 || transformation == 3 || transformation == 4) {
        double aux{};
        for (int i{0}; i < d; ++i) {
            std::cout << "\nBeta"<< i + 1 <<": ";
            std::cin >> aux;
            beta.push_back(aux);
        }
    }
    else if (transformation == 5) {
        double aux{};
        for (int i{0}; i < d; ++i) {
            std::cout << "\nMu"<< i + 1 <<": ";
            std::cin >> aux;
            beta.push_back(aux);
            std::cout << "\nSigma"<< i + 1 <<": ";
            std::cin >> aux;
            beta.push_back(aux);
        }
    }
    else if (transformation == 6) {
        double aux{};
        for (int i{0}; i < d; ++i) {
            std::cout << "\nMu"<< i + 1 <<": ";
            std::cin >> aux;
            beta.push_back(aux);
            std::cout << "\nSigma"<< i + 1 <<": ";
            std::cin >> aux;
            beta.push_back(aux);
            std::cout << "\nXi"<< i + 1 <<": ";
            std::cin >> aux;
            beta.push_back(aux);
        }
    }
    
    
    // Seed for simulation
    int seed{};
    std::cout << "\nSeed: ";
    std::cin >> seed;
    
    MTRNG randomSimulation(seed);
    
    int numSimulations{askNumberOfSimulations()};
    
    std::ofstream outFile{simulationFile};
    
    simulateMPH(mph, numSimulations, randomSimulation, outFile, transformation, beta);
    
    outFile.close();
}


void tailDependence() {
    // Relative error for quantile calculation
    double epsilon{0.01};
    
    // Names of Input and output files
    const std::string phasesFile{"MPHFitted_6p_s7_5000_f4_MV2.txt"}; //Input file - parameters of the MPH
    
    // Determinates the number of phases and dimension of the distribution base on the input file
    std::ifstream inPhases{phasesFile};
    DimensionsMPH dimensions{dimensionsOfFile(inPhases)};
    long p{dimensions.m_p};
    long d{dimensions.m_dim};
    std::cout << p << " phases in the file and dimension" << d << '\n';
    
    //Parameters of the distribution
    MPH mph(p,d);
    mph.readParametersFromFile(inPhases);
    std::cout << "\nParameters of the distribution\n";
    mph.print(); //Prints parameters of the MPH
    
    inPhases.close();
    
    
    double u{0.99};
    std::cout << "\nu for upper tail dependence: ";
    std::cin >> u;
    
    double ul{0.01};
    std::cout << "\nu for lower tail dependence: ";
    std::cin >> ul;
    
    
    // Seed for simulation
    int seed{};
    std::cout << "\nSeed: ";
    std::cin >> seed;
    
    MTRNG randomSimulation(seed);
    
    int numSimulations{askNumberOfSimulations()};
    
    double lambdaU{};
    double lambdaL{};
    
    lambdaU = upperTailDependence(mph, numSimulations, u, epsilon, randomSimulation);
    lambdaL = lowerTailDependence(mph, numSimulations, 1 - u, epsilon, randomSimulation);
    std::cout << "\nLambda_U: " << lambdaU << '\n' << "Lambda_L: " << lambdaL;
    
}


void tailDependence3d() {
    // Relative error for quantile calculation
    double epsilon{0.01};
    
    // Names of Input and output files
    const std::string phasesFile{"MPHFitted_6p_s7_5000_f4_MV2.txt"}; //Input file - parameters of the MPH
    
    // Determinates the number of phases and dimension of the distribution base on the input file
    std::ifstream inPhases{phasesFile};
    DimensionsMPH dimensions{dimensionsOfFile(inPhases)};
    long p{dimensions.m_p};
    long d{dimensions.m_dim};
    std::cout << p << " phases in the file and dimension" << d << '\n';
    
    //Parameters of the distribution
    MPH mph(p,d);
    mph.readParametersFromFile(inPhases);
    std::cout << "\nParameters of the distribution\n";
    mph.print(); //Prints parameters of the MPH
    
    inPhases.close();
    
    
    double u{0.99};
    std::cout << "\nu for upper tail dependence: ";
    std::cin >> u;
    
    double ul{0.01};
    std::cout << "\nu for lower tail dependence: ";
    std::cin >> ul;
    
    
    // Seed for simulation
    int seed{};
    std::cout << "\nSeed: ";
    std::cin >> seed;
    
    MTRNG randomSimulation(seed);
    
    int numSimulations{askNumberOfSimulations()};
    
    std::vector<double> lambdaU{};
    
    lambdaU = upperTailDependence3d(mph, numSimulations, u, epsilon, randomSimulation);
    
    std::cout << "\nLambda_U12: " << lambdaU[0] << '\n' << "Lambda_U13: " << lambdaU[1] << "Lambda_U23: " << lambdaU[2];
    
}




void EMforMPH(int transformation) {
    // Names of I/O files - Change if needed
    //  Input files
    const std::string unweightedSampleFile{"sampleMPH.txt"}; // Unweighted uncensored Sample
    const std::string unweightedCensoredFile{"censoredMPH.txt"}; // Unweighted uncensored Sample
    const std::string structureFile{"legalMPH.txt"}; // Structure of the PH
    const std::string initialPhasesFile{"initialMPH.txt"}; // Initial values of the PH
    //  Output files
    const std::string fittedPhasesFile{"MPHFitted.txt"}; // Fitted values
    const std::string simulationFile{"simMPHFitted.txt"}; // Simulation file
    
    // Constant that determinates when to print the Loglikelihood
    const int printLikehood{10};
    
    // Epsilon for the Uniformization method
    const double epsilon{0.0001};
    
    // Vector that will contain the uncensored sample
    std::vector<std::vector<Sample>> observations;
    std::vector<Sample> sumedObervations;
    
    // Vector that will contain the right-censored sample
    std::vector<std::vector<Sample>> rightCensored;
    std::vector<Sample> sumedRightCensored;
     
     // Type of sample - (1)From file (2) generated from a continuous distribution - available only for a MPH fit
    int sampleType{};

    if (transformation == 0) {
        sampleType = askForSampleType();
    } else{
        sampleType = 1;
    }
    
    long dim{0};
    if (sampleType == 1) { // File option
        // Ask if the data has censored values
        int censorIndicator{askIfCensored()};
        
        std::string sampleFile;
        
        sampleFile = (censorIndicator == 1) ? unweightedSampleFile : unweightedCensoredFile;

        std::ifstream infil{sampleFile};

        dim = readDataForMPH(infil, observations, sumedObervations, rightCensored, sumedRightCensored, censorIndicator, transformation);   // Read data from file
        

        if (censorIndicator == 2) { // Censored
            std::cout << "\nSummary of the sample\n";
            infoDataMPH(observations, rightCensored);
            
            std::cout << "\nSummary of the sumed sample\n";
            infoDataPH(sumedObervations, sumedRightCensored);
            
            // Sort data for RK
            for (int j{0}; j < dim; ++j) {
                sortObservations(rightCensored[j]);
            }
            sortObservations(sumedRightCensored);
        }
        else { // Uncensored
            std::cout << "\nSummary of the sample\n";
            sampleStatisticsMPH(observations);   // Print basic info of data
            
            std::cout << "\nSummary of the sumed sample\n";
            sampleStatisticsPH(sumedObervations);
        }
        // Sort data based on obs
        for (int j{0}; j < dim; ++j) {
            sortObservations(observations[j]);
        }
        sortObservations(sumedObervations);
        
        infil.close();
    }
    else if (sampleType == 2) { // Continuous distribution
        // One needs to initializate the samples
        dim = 2; // Only includes bivariate distributions
        std::vector<Sample> auxVector;
        for (int i{0}; i < dim; ++i) {
            observations.push_back(auxVector);
            rightCensored.push_back(auxVector);
        }
        
        askMultivariateDist(observations, sumedObervations); // Gives list of options and put sample in vector
        
        std::cout << "\nSummary of marginals\n";
        sampleStatisticsMultDist(observations);   // Print basic info of data
        
        std::cout << "\nSummary of the sumed marginals\n";
        sampleStatisticsPH(sumedObervations);
    }
    
    // Type of initial structure for the MPH
    int structureType{askForStructure()};

    long p{};
    int seed{1};
    
    MPH *tempMPH;
    
    if (structureType == 1) { // General and random
        p = askNumberOfPhases();
         
        std::cout << "\nSeed: ";
        std::cin >> seed; // Seed for random initialitation
        
        tempMPH = new MPH(p,dim);
        
        // Put initial values in temporal object
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> piLegal(1,p);
        piLegal += 1;
        Numeric_lib::Matrix<double,2> TLegal(p,p);
        TLegal += 1;
        Numeric_lib::Matrix<double,2> RLegal(p,p);
        RLegal += 1;
          
        double scale{1.0};
        scale = sumedObervations[(sumedObervations.size() - 1) / 2].obs;
          
        tempMPH->randomPhase(piLegal, TLegal, RLegal, scale, myRandom);
    }
    else if (structureType == 2) { // Structure from file and random
        std::ifstream inLegalFile{structureFile};
        p = dimensionsOfFile(inLegalFile).m_p; // Subtract the implicit size of p
        std::cout << p << " phases in the structure of the file";
         
        std::cout << "\nSeed: ";
        std::cin >> seed; // Seed for random initialitation
        
        tempMPH = new MPH(p,dim);
        
        // Put initial values in the temporal object
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> piLegal(1,p);
        Numeric_lib::Matrix<double,2> TLegal(p,p);
        Numeric_lib::Matrix<double,2> RLegal(p,dim);
          
        readStructureForMPH(inLegalFile, piLegal, TLegal, RLegal);
        
        double scale{1.0};
        scale = sumedObervations[(sumedObervations.size() - 1) / 2].obs;
          
        tempMPH->randomPhase(piLegal, TLegal,RLegal, scale, myRandom);
        
        inLegalFile.close();
    }
    else { // Initial phases from file
        std::ifstream inPhases{initialPhasesFile};
        p = dimensionsOfFile(inPhases).m_p; // Subtract the implicit size of p
         
        std::cout << p << " phases in the structure of the file";
        
        tempMPH = new MPH(p,dim);
          
        tempMPH->readParametersFromFile(inPhases);
        
        inPhases.close();
    }
    
    // Initializates the MPH object
    MPH mph{*tempMPH};
    
    delete tempMPH;
    
     
    std::cout << "\nIntilian parameters of the distribution\n";
    mph.print(); // Prints initial values of the PH
    
    // Selection of the method for the matrix exponential
    std::cout << "\n(RK not recomended if the reward matrix has zeros)";
    int emMethod{askForMatrixExponentialMethod()};
     
    // Indicates if one wants more iterations of the EM
    int stepsIndicator{1};
    
    // EM algorithm
    while (stepsIndicator == 1) {
        int stepsFirstEM{askNumberOfEMStepsFirst()};
        int stepsSecondEM{askNumberOfEMStepsSecond()};
            
        Timer myTime; // Object to measure running time
        
        std::vector<double> w(dim,1);
        std::vector<int> newStates;
        
        PhaseType PHsumOfMarginals(mph.linearCombination(w, newStates));
        
        
        if (emMethod == 1) { // Matlab algorithm
            EMIterateMPH(stepsFirstEM, stepsSecondEM, observations, sumedObervations, rightCensored, sumedRightCensored, mph, PHsumOfMarginals, printLikehood);
        }
        else if (emMethod == 2) { // Runge-Kutta
            EMIterateMPH_RK(stepsFirstEM, stepsSecondEM, observations, sumedObervations, rightCensored, sumedRightCensored, mph, PHsumOfMarginals, printLikehood);

        } else if (emMethod==3) { // Uniformization
            EMIterateMPH_Uni(stepsFirstEM, stepsSecondEM, epsilon ,observations, sumedObervations, rightCensored, sumedRightCensored, mph, PHsumOfMarginals, printLikehood);
        }

        // Prints the parameters of the fitted distribution and some basic information to compare with the sample
        std::cout << "Parameters of the fitted distribution\n";
        mph.print();
        std::cout << "\nExpected value of fitted distribution: " << '\n';
        printVector(mph.mean());
        std::cout << "\nStandard deviation of fitted distribution: " << '\n';
        printVector(mph.sd());
        std::cout << "\nCorrelations: " << '\n';
        printMatrix(mph.correlationMatrix());
            
        std::cout << "\nParameters of the sum of the fitted distribution\n";
        std::cout << "Expected value of sum: " << PHsumOfMarginals.mean() << '\n';
        std::cout << "Standard deviation of sum: " << PHsumOfMarginals.sd() << '\n';
        
        
        // Running time
        std::cout << "\nTime elapsed: " << myTime.elapsed() << " seconds\n";
        
        // Ask if one wants more Iterations
        stepsIndicator = askIfMoreIterations();
    }
     
    // Option to save parameters of the fitted distribution
    int saveInFile{askIfSaveFit()};
    
    if (saveInFile == 1) {
        std::ofstream outFile{fittedPhasesFile};
        mph.exportToFile(outFile);
        
        outFile.close();
    }
    
    
    int simulateFromFit{askSimulateFromFit()};
    
    if (simulateFromFit == 1) {
        // Seed for simulation
        int seed{};
        std::cout << "\nSeed: ";
        std::cin >> seed;
        
        MTRNG randomSimulation(seed);
        
        int numSimulations{askNumberOfSimulations()};
        
        std::ofstream outFile{simulationFile};
        
        // Temp
        std::vector<double> beta{};
        
        simulateMPH(mph, numSimulations, randomSimulation, outFile, transformation, beta);
        
        outFile.close();
    }
    
}



void EMforBivariatePH(int transformation) {
    // Names of I/O files - Change if needed
    //  Input files
    const std::string unweightedSampleFile{"sampleBiv.txt"}; // Unweighted uncensored Sample
    const std::string structureFile{"legalBiv.txt"}; // Structure of the PH
    const std::string initialPhasesFile{"initialBiv.txt"}; // Initial values of the PH
    //  Output files
    const std::string fittedPhasesFile{"BivFitted.txt"}; // Fitted values
    const std::string simulationFile{"simBivFitted.txt"}; // Simulation file
    const std::string fittedBeta{"BetaBiv.txt"}; // Fitted value of the parameter-dependent distribution
    
    // Constant that determinates when to print the Loglikelihood
    const int printLikehood{10};
    
    // Epsilon for stoping criteria of gradient ascent
    const double epsilon{0.1};
    
    // Vector that will contain the uncensored sample - An algotithm with censored data is not available yet
    std::vector<BivariateSample> observations;
    
    // Type of sample - (1)From file (2) generated from a continuous distribution - available only for a MPH fit
    int sampleType{};

    if (transformation == 0) {
        sampleType = askForSampleType();
    } else{
        sampleType = 1;
    }

    if (sampleType == 1) { // File option
        std::string sampleFile{unweightedSampleFile};
            
        std::ifstream infil{sampleFile};

        readDataForBivPH(infil, observations, transformation);  //Read data from file

        std::cout << "\nSummary of the sample\n";
        sampleStatisticsBivPH(observations);   // Print basic info of data
        
        infil.close();
    }
    else if (sampleType == 2) { // Continuous distribution
        askBivariateDist(observations); // Gives list of options and put sample in vector
        
        double sumofWeights{0.0};
        for (int k{0}; k < observations.size(); ++k) {
            if (observations[k].weight > 1) {
                std::cout << k << " " << observations[k].x1 << " " << observations[k].x2;
            }
            sumofWeights += observations[k].weight;
        }
        std::cout << "\nSample size: " << observations.size() << '\n';
        std::cout << "Sum of Weights: " << sumofWeights << '\n';
    }

    //Type of initial structure for the bivariate PH
    int structureType{askForStructure()};

    long p1{};
    long p2{};
    int seed{1};
     
    BivariatePH *tempBivph;
    
    // Initializates the Bivariate PH object
    if (structureType == 1) { // General and random
        p1 = askNumberOfPhases(1);
        p2 = askNumberOfPhases(2);
        
        std::cout << "\nSeed: ";
        std::cin >> seed; // Seed for random initialitation
        
        // Allocate memory for the temporal object
        tempBivph = new BivariatePH(p1,p2);
        
        // Initiali values of the temporal object
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> alphaLegal(1,p1);
        alphaLegal += 1;
        Numeric_lib::Matrix<double,2> T11Legal(p1,p1);
        T11Legal += 1;
        Numeric_lib::Matrix<double,2> T12Legal(p1,p2);
        T12Legal += 1;
        Numeric_lib::Matrix<double,2> T22Legal(p2,p2);
        T22Legal += 1;
          
        double scale{1.0};
        
        tempBivph->randomPhase(alphaLegal, T11Legal, T12Legal,T22Legal ,scale, myRandom);
        
    }
    else if (structureType == 2) { // Structure from file and random
        
        std::ifstream inLegalFile{structureFile};
        DimensionsMPH auxiliarDim{dimensionsOfBivFile(inLegalFile)};
        
        p1 = auxiliarDim.m_p;
        p2 = auxiliarDim.m_dim;
        
        std::cout << p1 + p2 << " phases in the structure of the file, separated in p1 = " << p1 << " and p2 = " << p2 << '\n';
          
        std::cout << "\nSeed: ";
        std::cin >> seed; // Seed for random initialitation
        
        
        tempBivph = new BivariatePH(p1,p2);
          
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> alphaLegal(1,p1);
        Numeric_lib::Matrix<double,2> T11Legal(p1,p1);
        Numeric_lib::Matrix<double,2> T12Legal(p1,p2);
        Numeric_lib::Matrix<double,2> T22Legal(p2,p2);
          
        readStructureForBivPH(inLegalFile, alphaLegal, T11Legal, T12Legal, T22Legal);
        
        double scale{1.0};
        tempBivph->randomPhase(alphaLegal, T11Legal, T12Legal, T22Legal ,scale, myRandom);
        
        inLegalFile.close();
        
    }
    else { // Initial phases from file
        std::ifstream inPhases{initialPhasesFile};
        DimensionsMPH auxiliarDim{dimensionsOfBivFile(inPhases)};
                
        p1 = auxiliarDim.m_p;
        p2 = auxiliarDim.m_dim;

        std::cout << p1 + p2 << " phases in the structure of the file, separated in p1 = " << p1 << " and p2 = " << p2 << '\n';
        
        tempBivph = new BivariatePH(p1,p2);
          
        tempBivph->readParametersFromFile(inPhases);
        
        inPhases.close();
        
    }
     
    // Initializates the Bivariate PH object
    BivariatePH bivph{*tempBivph};
    
    delete tempBivph;
    
     
    std::cout << "\nIntilian parameters of the distribution\n";
    bivph.print(); // Prints initial values of the PH
    
    double beta1{};
    double beta2{};
    double lambda{};
    
    if (transformation == 2 || transformation == 3) {
        std::cout << "\nInitial beta1: ";
        std::cin >> beta1;
        std::cout << "\nInitial beta2: ";
        std::cin >> beta2;
        
        lambda = askLearningRate();
    }
    
    //Selection of the method for the matrix exponential
//    std::cout<<"\n(RK not recomended if the reward matrix has zeros)";
//    int emMethod{ askForMatrixExponentialMethod() };
//
    // Indicates if one wants more iterations of the EM
    int stepsIndicator{1};
    
    
    // EM algorithm
    while (stepsIndicator == 1) {
        std::cout << "\nMatrix exponential is computed using a Pade approximation\n";
        int stepsEM{askNumberOfEMSteps()};
        
        Timer myTime; // Object to measure running time

        EMIterateBivPH(stepsEM, observations, bivph, printLikehood, transformation, beta1, beta2, lambda, epsilon);

        // Prints the parameters of the fitted distribution and some basic information to compare with the sample
        std::cout << "Parameters of the fitted distribution\n";
        bivph.print();
        if (transformation == 0) {
            std::cout << "\nExpected value of fitted distribution: " << '\n';
            printVector(bivph.mean());
            std::cout << "\nStandard deviation of fitted distribution: " << '\n';
            printVector(bivph.sd());
            std::cout << "\nCorrelations: " << '\n' << bivph.correlation() << '\n';
        }
        else if (transformation > 1) {
            std::cout << "\nBeta1: " << beta1 << '\t' << "\nBeta2: " << beta2 << '\n';
        }
        
        
        // Running time
        std::cout << "\nTime elapsed: " << myTime.elapsed() << " seconds\n";
        
        // Ask if one wants more Iterations
        stepsIndicator = askIfMoreIterations();
         
    }
    
    // Option to save parameters of the fitted distribution
    int saveInFile{askIfSaveFit()};
    
    if (saveInFile == 1) {
        std::ofstream outFile{fittedPhasesFile};
        bivph.exportToFile(outFile);
        
        outFile.close();
        
        if (transformation > 1) {
            std::ofstream outFileBeta{fittedBeta};
            outFileBeta << beta1 << '\n' << beta2;
        
            outFileBeta.close();
        }
    }
    
    // Option to simulate from fitted distribution
    int simulateForFit{askSimulateFromFit()};
    
    if (simulateForFit == 1) {
        MPH mph(bivph.getPi(), bivph.getT(), bivph.getR());
    
        // Seed for simulation
        int seed{};
        std::cout << "\nSeed: ";
        std::cin >> seed;
        
        MTRNG randomSimulation(seed);
        
        int numSimulations{askNumberOfSimulations()};
        
        std::ofstream outFile{simulationFile};
        
        //Temp
        std::vector<double> beta{};
        beta.push_back(beta1);
        beta.push_back(beta2);
        
        simulateMPH(mph, numSimulations, randomSimulation, outFile, transformation, beta);
        
        outFile.close();
    }
    
}



/* *********
    NPH
******** */


int readOptionScale() {
    int readScale{};
    do {
        std::cout << "\nScale function:\n";
        std::cout << "   1. Read from file \n";
        std::cout << "   2. Select from list\n";
        std::cout << "Select 1 or 2: ";
        std::cin >> readScale;
        if (readScale != 1 && readScale != 2) {
            std::cout<<"Please enter a valid option\n";
        }
    } while (readScale != 1 && readScale != 2);
    
    return readScale;
}


// List of scale distribution available - Initializates the object and retuns a pointer to it
ScaleDistribution * askTypeOfScale() {
    double epsilon{1e-8}; // For Zeta Riemann
    double initial{1}; // Initial value of the scale distribution
    
    int typeScaleDist{};
    do {
        std::cout << "\nScale distribution:\n";
        std::cout << "      1. Riemann Zeta\n";
        std::cout << "      2. Discretizate Pareto over geometric progression\n";
        std::cout << "      3. Discretizate Pareto over arithmetic progression\n";
        std::cout << "      4. Discretizate Weibull over geometric progression\n";
        std::cout << "      5. Discretizate Weibull over arithmetic progression\n";
        std::cout << "Select(1-5): ";
        std::cin >> typeScaleDist;
        if (typeScaleDist < 1 || typeScaleDist > 5) {
            std::cout << "Please enter a valid option\n";
        }
    } while (typeScaleDist < 1 || typeScaleDist > 5);
    
    
    ScaleDistribution *scaleComponent_ptr;
       
    switch (typeScaleDist) {
        case 1: {
            std::cout << "\nRiemann Zeta\n";
            std::cout << "Theta: ";
            double theta{};
            std::cin >> theta;
            scaleComponent_ptr = new RiemannZeta(theta, epsilon);
            break;
        }
        case 2: {
            std::cout << "\nDiscretizate Pareto over geometric progression\n";
            std::cout << "Theta: ";
            double theta{};
            std::cin >> theta;
            std::cout << "Ratio (larger than 1): ";
            double ratio{};
            std::cin >> ratio;
            scaleComponent_ptr = new DiscretePareto(theta, 2, initial, ratio);
            break;
        }
        case 3: {
            std::cout << "\nDiscretizate Pareto over arithmetic progression\n";
            std::cout << "Theta: ";
            double theta{};
            std::cin >> theta;
            std::cout << "Increment (positive): ";
            double ratio{};
            std::cin >> ratio;
            scaleComponent_ptr = new DiscretePareto(theta, 1, initial, ratio);
            break;
        }
        case 4: {
            std::cout << "\nDiscretizate Weibull over geometric progression\n";
            std::cout << "Tau: ";
            double tau{};
            std::cin >> tau;
            std::cout << "Ratio (larger than 1): ";
            double ratio{};
            std::cin >> ratio;
            scaleComponent_ptr = new DiscreteWeibull(tau, 2, initial, ratio);
            break;
        }
        case 5: {
            std::cout << "\nDiscretizate Weibull over arithmetic progression\n";
            std::cout << "Tau: ";
            double tau{};
            std::cin >> tau;
            std::cout << "Increment (positive): ";
            double ratio{};
            std::cin >> ratio;
            scaleComponent_ptr = new DiscreteWeibull(tau, 1, initial, ratio);
            break;
        }
        default:
            scaleComponent_ptr = nullptr;
            break;
    }
    
    return scaleComponent_ptr;
}

// Read scale distribution from file, initializates the object and retuns a pointer to it
ScaleDistribution* readScaleParameters(std::ifstream & inputFile) {
    double epsilon{1e-8}; // For Zeta Riemann
    int type{};
    double parameter{};
    double initial{};
    double increment{};
    
    inputFile >> type;
    inputFile >> parameter;
    inputFile >> initial;
    inputFile >> increment;
    
    ScaleDistribution *scaleComponent_ptr;
    
    switch (type) {
        case 1:
            std::cout << "Riemann Zeta with parameter " << parameter << '\n';
            scaleComponent_ptr = new RiemannZeta(parameter, epsilon);
            break;
        case 2:
                std::cout << "Discrete Pareto on geometric progresion with parameter " << parameter << ", initial value " << initial << " and ratio " << increment << '\n';
                scaleComponent_ptr = new DiscretePareto(parameter, 2, initial, increment);
                break;
        case 3:
            std::cout << "Discrete Pareto on arithmetic progresion with parameter " << parameter << ", initial value " << initial << " and increment " << increment << '\n';
            scaleComponent_ptr = new DiscretePareto(parameter, 1, initial, increment);
            break;
        case 4:
            std::cout << "Discrete Weibull on geometric progresion with parameter " << parameter << ", initial value " << initial << " and ratio " << increment << '\n';
            scaleComponent_ptr = new DiscreteWeibull(parameter, 2, initial, increment);
            break;
        case 5:
            std::cout << "Discrete Weibull on arithmetic progresion with parameter " << parameter << ", initial value " << initial << " and increment " << increment << '\n';
            scaleComponent_ptr = new DiscreteWeibull(parameter, 1, initial, increment);
            break;
        default:
            scaleComponent_ptr = nullptr;
            break;
    }
    
    return scaleComponent_ptr;
}


void SimulationForNPH() {
    // Input and output files
    //  Input files
    const std::string phasesFile{"phases.txt"}; // Parameters of the PH
    const std::string scaleFile{"scaleComponent.txt"}; // Parameters of the scale component
    //  Output files
    const std::string simulationFile{"simSampleNPH.txt"}; // Simulations
    
    // Determinates the dimension of the distribution based on the input file
    std::ifstream inPhases{phasesFile};
    long p{};
    p = dimensionsOfFile(inPhases).m_p; //Subtract the implicit size of p
    std::cout << p << " phases in the file\n";
    
    // Parameters of the distribution
    PhaseType PH(p);
    PH.readParametersFromFile(inPhases);
    std::cout << "\nParameters of the PH component\n";
    PH.print(); // Prints the values
    
    inPhases.close();
    
    int readScale{readOptionScale()};
    
    ScaleDistribution *scaleComponent_ptr;
    
    if (readScale == 1) {
        std::ifstream inScale{scaleFile};
        
        scaleComponent_ptr = readScaleParameters(inScale);
        
        inScale.close();
    }
    else {
        scaleComponent_ptr = askTypeOfScale();
    }
    
    // Seed for simulation
    int seed{};
    std::cout << "\nSeed: ";
    std::cin >> seed;

    MTRNG randomSimulation(seed);
    
    int numSimulations{askNumberOfSimulations()};
    
    std::ofstream outFile{simulationFile};
    
    simulateNPH(PH, scaleComponent_ptr, numSimulations, randomSimulation, outFile);
    
    outFile.close();
    
    delete scaleComponent_ptr;
}




int FixOrEstimate() {
    int readScale{};
    do {
        std::cout << "\nParameter of scale distribution:\n";
        std::cout << "   1. Estimate\n";
        std::cout << "   2. Fix\n";
        std::cout << "Select 1 or 2: ";
        std::cin >> readScale;
        if (readScale != 1 && readScale != 2) {
            std::cout<<"Please enter a valid option\n";
        }
    } while (readScale != 1 && readScale != 2);
    
    return readScale;
}



void EMforNPH() {
    //Names of the I/O files - Change if needed
    //  Input files
    const std::string unweightedSampleFile{"sample.txt"}; // Unweighted uncensored Sample
    const std::string weightedSampleFile{"weighted.txt"}; // Weighted uncensored Sample
    const std::string unweightedCensoredFile{"censored.txt"}; // Unweighted uncensored Sample
    const std::string weightedCensoredFile{"weightedCensored.txt"}; // Weighted uncensored Sample
    const std::string structureFile{"legal.txt"}; // Structure of the PH
    const std::string initialPhasesFile{"phases.txt"}; // Initial values of the PH
    const std::string initialScaleFile{"scale.txt"}; // Initial values of the scale component
    
    //  Output files
    const std::string fittedPhasesFile{"phasesFitted.txt"}; // Fitted values of the PH component
    const std::string fittedScaleFile{"scaleFitted.txt"}; // Fitted Scale distribution
    
    // Constant that determinates when to print the Loglikelihood
    const int printLikehood{10};
    
    // Epsilon for the Uniformization method
    //const double epsilon{0.0001};
    
    // Vector that will contain the uncensored sample
    std::vector<Sample> observations;
    // Vector that will contain the right-censored sample
    std::vector<Sample> rightCensored;
     
    int sampleType{};
    sampleType = askForSampleType(); // Type of sample - From file or generated from a continuous distribution
     
    if (sampleType == 1) { // File option
        // Ask if the data has censored values
        int censorIndicator{askIfCensored()};
        
        // Choice between file with weights or without weights
        int weightOption{askIfWeighted()};
       
        // Assigns the corresponding file name
        std::string sampleFile;
        if (weightOption == 1) {
            sampleFile = (censorIndicator == 1) ? unweightedSampleFile : unweightedCensoredFile;
        }
        else {
            sampleFile = (censorIndicator == 1) ? weightedSampleFile : weightedCensoredFile;
        }
        std::ifstream infil{sampleFile};
        
        // Read data from file - 0 at the end to not transform data
        readDataForPH(infil, observations, rightCensored, weightOption, censorIndicator);
        // Disply basic information and sort data
        if (censorIndicator == 2) { // Right censored
            infoDataPH(observations, rightCensored);
            sortObservations(rightCensored);
        }
        else { // Uncensored
            sampleStatisticsPH(observations);   // Print basic info of data
        }
        sortObservations(observations); // Sort data base on obs
        
        infil.close();
        
    }
    else if (sampleType == 2) { // Continuous distribution
        askContinousDist(observations); // Gives list of options and put sample in vector
        sampleStatisticsPH(observations);  // Prints basic info - not sorted since it is created in that way
    }
     
     
    int structureType{askForStructure()}; // Type of initial structure for the PH (3 options)
     
    long p{}; // Number of phases
    int seed{1};
     
    // Temporal PH object
    PhaseType *tempPH;
    
    // Determintes the seed (if needed) and value of p - Needed to initializates the PhaseType object
    if (structureType == 1) { // General and random
        // Get the number of phases
        p = askNumberOfPhases();
        
        std::cout << "\nSeed: ";
        std::cin >> seed; // Seed for random initialitation
        
        // Creates the temporal object
        tempPH = new PhaseType(p);
        
        // Put initial values in temporal object
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> piLegal(1,p);
        piLegal += 1;
        Numeric_lib::Matrix<double,2> TLegal(p,p);
        TLegal += 1;
          
        double scale{1.0};
        scale = observations[(observations.size() - 1) / 2].obs;
          
        tempPH->randomPhase(piLegal, TLegal, scale, myRandom);
    }
    else if (structureType == 2) { // Structure from file and random
        std::ifstream inLegalFile{structureFile};
        p = dimensionsOfFile(inLegalFile).m_p; // Subtract the implicit size of p
        
        std::cout << p << " phases in the structure of the file";
         
        std::cout << "\nSeed: ";
        std::cin >> seed; //Seed for random initialitation
        
        //Creates the temporal object
        tempPH = new PhaseType(p);
          
        //Put initial values in temporal object
        MTRNG myRandom(seed);
        Numeric_lib::Matrix<double,2> piLegal(1,p);
        Numeric_lib::Matrix<double,2> TLegal(p,p);
        
        readStructureForPH(inLegalFile, piLegal, TLegal);
          
        double scale{1.0};
        scale = observations[(observations.size() - 1) / 2].obs;
          
        tempPH->randomPhase(piLegal, TLegal, scale, myRandom);
        
        inLegalFile.close();
    }
    else { // Initial phases from file
        std::ifstream inPhases{initialPhasesFile};
        p = dimensionsOfFile(inPhases).m_p; // Subtract the implicit size of p
         
        std::cout << p << " phases in the structure of the file";
        
        // Creates temporal object
        tempPH = new PhaseType(p);
        
        // Put initial values in temporal object
        tempPH->readParametersFromFile(inPhases);
        
        inPhases.close();
    }
    
    // Initializates the phase-type object
    PhaseType PH{*tempPH};
    
    // Deletes temporal object
    delete tempPH;
    
     
    std::cout << "\nIntilian parameters of the distribution\n";
    PH.print(); // Prints initial values of the PH component
    
    
    // Inital values of the scale component
    int readScale{readOptionScale()};
    
    // Pointer to the scale component
    ScaleDistribution *scaleComponent_ptr;
    
    if (readScale == 1) { // Read initial values from file
        std::ifstream inScale{initialScaleFile};
        
        scaleComponent_ptr = readScaleParameters(inScale);
        
        inScale.close();
    }
    else { // Read initial values from list
        scaleComponent_ptr = askTypeOfScale();
    }
    
    // Ask if estimation of the parameter of the scale component has to be done
    int fixParameter{FixOrEstimate()};
    
    //Selection method for matrix exponential
    int emMethod{askForMatrixExponentialMethod()};
    
    // Selection of maximization method
    int maxMethod{};
    
    // Learning rate of gradient ascent
    double lambda{};
    
    if (fixParameter == 1 && scaleComponent_ptr->getid() > 2 ) {
        maxMethod = askMaxMethod();
        
        if (maxMethod == 1) { //Gradient ascent
            lambda = askLearningRate();
        }
    }
     
    //Indicates if one wants more iterations of the EM
    int stepsIndicator{1};
    
    // EM algorithm
    while (stepsIndicator == 1) {
        int stepsEM{askNumberOfEMSteps()};
            
        Timer myTime; // Object to measure running time
        
        if (emMethod == 1) { // Matlab algorithm
            EMIterateNPH(stepsEM, observations, rightCensored, PH, scaleComponent_ptr, fixParameter, printLikehood, maxMethod, lambda);
        }
        else if (emMethod == 2) { // Runge-Kutta
            //EMIterate_RK(stepsEM, observations, rightCensored, PH, printLikehood);

        }
        else if (emMethod == 3) { // Uniformisation
            //EMIterate_Uni(stepsEM, epsilon, observations, rightCensored, PH, printLikehood);
        }
        
        // Prints the parameters of the fitted distribution and some basic information to compare with the sample
        std::cout << "\nParameters of the fitted distribution:";
        std::cout << "\nPhase type component\n";
        PH.print();
        
        std::cout << "\nScale component";
        scaleComponent_ptr->print();
        
        // Running time
        std::cout << "\nTime elapsed: " << myTime.elapsed() << " seconds\n";
        
        // Ask if one wants more Iterations
        stepsIndicator = askIfMoreIterations();
         
    }
     
    // Option to save parameters of the fitted distribution
    int saveInFile{askIfSaveFit()};
    
    if (saveInFile == 1) {
        std::ofstream outFile{fittedPhasesFile};
        PH.exportToFile(outFile);
    
        std::ofstream outFileScale{fittedScaleFile};
        scaleComponent_ptr->exportToFile(outFileScale);
        
        outFile.close();
        outFileScale.close();
    }

    delete scaleComponent_ptr;
}
