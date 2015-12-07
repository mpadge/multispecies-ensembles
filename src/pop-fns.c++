/*
 * pop-fns.c++
 *
 * Contains the core functions used to generate the population simulations.
 * These depend on the following parameters:
 * 	1.	bmn = -0.001; mean of competition matrix; represents predation
 * 	2.	bsd = 0.05; represents strength of competition
 * 	3.	sigma_a = 0.1; SD of variation of autocorrelation coefficients
 * 	4.	sigma_rho = 0.1; SD of variation of degree of environmental sharing.
 */

#include "pop-fns.h"

void makeCommunity (int nSpecies, double sync, SpeciesData & speciesData, 
        base_generator_type & generator, int adj_rho)
{
	double tempd;

	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);

	double q = 0.7;
	for (int i=0; i<nSpecies; i++) 
    {
		speciesData.spvec (i).ar [0] = 1.0 + sigma_a * rnorm ();
		speciesData.spvec (i).ar [1] = sigma_a * rnorm () * q;
		speciesData.spvec (i).ar [2] = sigma_a * rnorm () * q * q;
		speciesData.spvec (i).ar [3] = sigma_a * rnorm () * q * q * q;
		tempd = 0.0;
		for (int j=0; j<4; j++) 
            tempd += fabs (speciesData.spvec (i).ar [j]);
		for (int j=0; j<4; j++) 
			speciesData.spvec (i).ar [j] = speciesData.spvec (i).ar [j] / tempd;

        /* NOTE the importance of the following condition. The upper limit has
         * to be the opposite of the lower limit, so the distribution stays at
         * all times symmetrical. Without this, the distribution becomes biased
         * with increasing sigma_rho, producing greater *mean* values, and
         * skewing everything.
         *
         * NOTE also the important implication that values of rho are normally
         * distributed, and *NOT* values of rho^2.
         */
		speciesData.spvec (i).rho = sync + sigma_rho * rnorm ();
		while (speciesData.spvec (i).rho <= 0.0 || 
                speciesData.spvec (i).rho >= (2.0 * sync))
			speciesData.spvec (i).rho = sync + sigma_rho * rnorm ();
	}

	for (int i=0; i<nSpecies; i++) 
		speciesData.compMat (i, i) = 0.0;
	for (int i=0; i<(nSpecies - 1); i++) 
		for (int j=(i + 1); j<nSpecies; j++) {
			tempd = rnorm ();
			speciesData.compMat (i, j) = bmn + bsd * tempd;
			speciesData.compMat (j, i) = bmn - bsd * tempd;
        }
    // NOTE that with bmn = -0.001 and bsd = 0.05, values of compMat still
    // generally reach extrema of around +/- 0.2.

	if (adj_rho > 0) {
		// Then sort the rho_i values according to compVec strengths. Start by
        // calculating the mean interaction strength of species:
		boost::numeric::ublas::vector<double> compVec (nSpecies);
		for (int i=0; i<nSpecies; i++) 
        {
			compVec (i) = 0.0;
			for (int j=0; j<nSpecies; j++)
				if (j != i) 
                    compVec (i) += fabs (speciesData.compMat (i, j));
			compVec (i) = compVec (i) / ((double) nSpecies - 1.0);
		}

		//static int index = 0;
		int index = 0;
		std::vector <std::pair<double, int> > cmatVector;
		for (int i=0; i<nSpecies; i++)
			cmatVector.push_back (std::pair<double, int>(compVec (i), index++));
		std::sort (cmatVector.begin(), cmatVector.end());
		std::vector <std::pair <double, int> >::const_iterator itr;
		std::vector <double> rhoVector;
		for (itr = cmatVector.begin(); itr != cmatVector.end(); ++itr)
			rhoVector.push_back (speciesData.spvec ((*itr).second).rho);
		for (int i=0; i<nSpecies; i++)
			speciesData.spvec (i).rho = rhoVector [i];
		cmatVector.resize (0);
		rhoVector.resize (0);
	}
} // end SpeciesData makeCommunity


void makeNoise (int nSpecies, dmat & noise, base_generator_type & generator)
{
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::uniform_real<> > runif(generator, uni_dist);
	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);

	for (int t=0; t<len_t; t++)
		for (int i=0; i<=nSpecies; i++)
			noise (t, i) = rnorm();
} // end dmat makeNoise

void runPop (int nSpecies, dmat noise, SpeciesData speciesData, dvec & pop_t)
{
	double tempd;

	boost::numeric::ublas::vector<double> spvec0 (nSpecies), spvec1 (nSpecies), 
        spvec2 (nSpecies), spvec3 (nSpecies), tempvec (nSpecies);

	// Add shared signal to species noise
	for (int t=0; t<len_t; t++) 
    {
		pop_t (t) = 0.0;
		for (int i=0; i<nSpecies; i++) 
        {
            tempd = speciesData.spvec (i).rho;
			noise (t, i) = tempd * noise (t, nSpecies) +
				sqrt (1.0 - tempd * tempd) * noise (t, i);
		}
	}

	// Run in first three AR steps
	for (int i=0; i<nSpecies; i++) 
    {
		spvec0 (i) = noise (0, i);
		for (int j=0; j<nSpecies; j++)
			spvec0 (i) += noise (0, j) * speciesData.compMat (j, i);
		pop_t (0) += spvec0 (i);
		spvec1 (i) = speciesData.spvec (i).ar [0] * noise (1, i) + 
            speciesData.spvec (i).ar [1] * spvec0 (i);
	}
	tempvec = spvec1;
	for (int i=0; i<nSpecies; i++)
		for (int j=0; j<nSpecies; j++)
			tempvec (i) += spvec1 (j) * speciesData.compMat (j, i);
	
	spvec1 = tempvec;
	for (int i=0; i<nSpecies; i++) 
    {
		pop_t (1) += spvec1 (i);
		spvec2 (i) = speciesData.spvec (i).ar [0] * noise (2, i) + 
            speciesData.spvec (i).ar [1] * spvec1 (i) +
			speciesData.spvec (i).ar [2] * spvec0 (i);
	}
	tempvec = spvec2;
	for (int i=0; i<nSpecies; i++)
		for (int j=0; j<nSpecies; j++)
			tempvec (i) += spvec2 (j) * speciesData.compMat (j, i);
	spvec2 = tempvec;
	for (int i=0; i<nSpecies; i++) 
        pop_t (2) += spvec2 (i);

	// Then the full run
	for (int t=3; t<len_t; t++) 
    {
		for (int i=0; i<nSpecies; i++)
			spvec3 (i) = speciesData.spvec (i).ar [0] * noise (t, i) + 
                speciesData.spvec (i).ar [1] * spvec2 (i) +
				speciesData.spvec (i).ar [2] * spvec1 (i) + 
                speciesData.spvec (i).ar [3] * spvec0 (i);
		tempvec = spvec3;
		for (int i=0; i<nSpecies; i++)
			for (int j=0; j<nSpecies; j++)
				tempvec (i) += spvec3 (j) * speciesData.compMat (j, i);
		for (int i=0; i<nSpecies; i++) 
            pop_t (t) += tempvec (i);
		spvec0 = spvec1;
		spvec1 = spvec2;
		spvec2 = tempvec;
	} // end for i over len_t
} // end dev runPop
