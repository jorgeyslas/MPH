// Interface with the user

#ifndef Interface_h
#define Interface_h

#include "Matrix.h"
#include "MatrixOperations.h"
#include "Timer.h"
#include "PhaseType.h"
#include "RandomNumbers.h"
#include "Simulation.h"
#include "EM.h"
#include "Distributions.h"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>

//General interface
void option();

//Univariate PH
void SimulationForPH(int transformation = 0);
void EMforPH(int transformation = 0);

//MPH
void SimulationForMPH(int transformation = 0);
void EMforMPH(int transformation = 0);

void EMforBivariatePH(int transformation = 0);

void tailDependence();
void tailDependence3d();

//NPH
void SimulationForNPH();
void EMforNPH();

#endif /* Interface_h */
