/*
 * model050-theoretical.h
 */

#include <stdlib.h> // has abs function
#include <math.h>
#include <iostream>
#include <fstream> // file in & out
//#include <string.h>
#include <stdio.h>
#include <dirent.h> // For reading directories
#include <time.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
//using namespace std;

#define PI 3.1415926535897932384626433832795

typedef boost::numeric::ublas::vector<double> dvec;

const int inix = -9999, nSpecies_max = 100;
const double dnix = -9999.0, bmn = -0.001, bsd = 0.05, sigma = 1.0;

double calcR2 (double rho, int n);
