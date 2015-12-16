/*
 * pop-fns.h
 */

#include "utils.h"

#ifndef POPFNS_H
#define POPFNS_H

const int len_t = 1000;
/*
 * NOTE that setting len_t = 100 produces notably noisier results, so it
 * needs to be greater, hence the value of 1,000, even through computation
 * then takes some time.
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
