/*
 * trophic_levels.h
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
#include "boost/multi_array.hpp"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
//using namespace std;

#define PI 3.1415926535897932384626433832795

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;
typedef boost::numeric::ublas::vector<int> ivec;
typedef boost::numeric::ublas::vector<double> dvec;
typedef boost::numeric::ublas::matrix<double> dmat;
typedef boost::numeric::ublas::matrix<double> darr;

const int inix = -9999, nSpecies = 20;
const double dnix = -9999.0, bmn = -0.001, bsd = 0.05;

struct SpeciesData{
	dmat compMat;
};

int makeCommunity (int nSpecies, SpeciesData & speciesData, 
        base_generator_type & generator);
