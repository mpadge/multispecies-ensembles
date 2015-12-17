/***************************************************************************
 *  Project:    multispecies-ensembles
 *  File:       pop-fns.h
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


#include "utils.h"

#ifndef POPFNS_H
#define POPFNS_H

const int len_t = 1000;
/*
 * NOTE that setting len_t = 100 produces notably noisier results, so it
 * needs to be greater, hence the value of 1,000. 
 * 
 * Also NOTE that the sigma values are actually SDs, not variances!
 */
const double bmn = -0.001, bsd = 0.05, sigma_a = 0.5, sigma_rho = 0.5;

struct SpeciesPars{
    double ar [4], rho;
};
/* SpeciesPars is a vector of (nSpecies) with:
 * 	a		= AR coefficients around means of [1,0,0], with
 * 				SDs of sigma_a * [1, 1/2. 1/4]
 * 	rho	    = degree of environment sharing (mean = rho)
 */
typedef boost::numeric::ublas::vector<SpeciesPars> species_vec;
struct SpeciesData{
    species_vec spvec;
    dmat compMat;
};

void makeCommunity (int nSpecies, double sync, SpeciesData & speciesData, 
        base_generator_type & generator, int adj_rho);
void makeNoise (int nSpecies, dmat & noise, base_generator_type & generator);
void runPop (int nSpecies, dmat noise, SpeciesData speciesData, dvec & p_vec);

#endif
