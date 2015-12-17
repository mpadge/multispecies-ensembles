/***************************************************************************
 *  Project:    multispecies-ensembles
 *  File:       utils.h
 *  Language:   C++
 *
 *  multispecies-ensembles is free software: you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  multispecies-ensembles is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 *  Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  multispecies-ensembles.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright   Mark Padgham December 2015
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    Simulates ensembles of multiple species responding to 
 *                  partially shared stochastic environmental variation.
 *
 *  Project Structure:  
 *      1. model        
 *              The full simulation model        
 *      2. model050     
 *              Estimates the limits of degree of sharing (rho^2) above which
 *              correlations with aggregate abundance exceed R^2=50%.
 *      3. model050-theoretical
 *              As for model050, but using analytic expresssions for species
 *              interactions only, neglecting other effects.
 *      4. model-envcor
 *              Examines additional effects of (i) correlations between
 *              otherwise independent parts of each species' environmental
 *              variation, and (ii) correlations between interaction strengths
 *              and degrees of environmental sharing.
 *      5. trophic-levels
 *              Separate routine to estimate number of equivalent trophic levels
 *              from random community matrices.
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options
 ***************************************************************************/


#ifndef UTILS_H
#define UTILS_H

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

#include <boost/program_options.hpp>

#define PI 3.1415926535897932384626433832795

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;
typedef boost::numeric::ublas::vector<int> ivec;
typedef boost::numeric::ublas::vector<double> dvec;
typedef boost::numeric::ublas::matrix<double> dmat;

const double DOUBLE_MAX = std::numeric_limits<double>::max (),
      DOUBLE_MIN = -DOUBLE_MAX;

struct RegrResults {
    double r2, slope, intercept, SS, tval;	};

struct myTime{
    int hh, mm;
    float ss;	};

RegrResults regression (dvec x, dvec y);
void timeout (double tseconds);

#endif
