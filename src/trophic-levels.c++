/*
 * trophic_levels.c++
 *
 * Constructs a community matrix, and then counts equivalent trophic levels,
 * defining an inter-level connection as one exceeding a given threshold, expressed
 * as a multiple of the SD of the community matrix. Numbers of levels depend on
 * this threshold, but for reasonably large values (> 2 * SD), numbers of levels
 * range between 2 and 6. Proportions of actual networks generated for various
 * numbers of trophic levels are given below for thresholds defined by several 
 * multiples of SD.
 */

#include "trophic-levels.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[])
{
	int nLevels;
	double tempd, sumr, sumr2;
	SpeciesData speciesData;
	time_t seed;
	base_generator_type generator(42u);

	//std::cout.setf (std::ios::fixed, std::ios::floatfield);   // floatfield set to fixed
	//std::cout.precision(2);
	int nTrials = atoi (argv [1]);
	if (nTrials < 1) 
        nTrials = 1000;
	std::cout << "nTrials = " << nTrials << std::endl;

	time (&seed);
	generator.seed(static_cast<unsigned int>(seed));
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::uniform_real<> > runif(generator, uni_dist);
	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);
	// Burn both generators in
	for (int i=0; i<20; i++) {
		tempd = runif();
		tempd = rnorm();
	}

	// Set up data to fixed size of nSpecies
	boost::numeric::ublas::matrix<double> tempdmat (nSpecies, nSpecies);
	speciesData.compMat = tempdmat;
	boost::numeric::ublas::vector<int> troph_levels (10);
	for (int i=0; i<10; i++) 
        troph_levels (i) = 0;

	for (int i=0; i<nTrials; i++) 
    {
		nLevels = makeCommunity (nSpecies, speciesData, generator);
		troph_levels (nLevels - 1)++;
	}
	nLevels = 0;
	for (int i=0; i<10; i++) 
        nLevels += troph_levels (i);
	std::cout << "Proportion of networks with numbers of trophic levels:" << 
        std::endl;
	for (int i=0; i<10; i++)
		if (troph_levels (i) > 0) 
			std::cout << "(" << i + 1 << ", " << 
                (double) troph_levels (i) / (double) nLevels << ")" << std::endl;

	return 0;
}; // end main

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         SUBFUNCTIONS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


int makeCommunity (int nSpecies, SpeciesData & speciesData, 
        base_generator_type & generator)
{
	int nLevels;
	double tempd;
	bool search_done;

	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);

	for (int i=0; i<nSpecies; i++)
		speciesData.compMat (i, i) = 0.0;
	for (int i=0; i<(nSpecies - 1); i++)
		for (int j=(i + 1); j<nSpecies; j++) 
        {
			tempd = rnorm ();
			speciesData.compMat (i, j) = bmn + bsd * tempd;
			speciesData.compMat (j, i) = bmn - bsd * tempd;
        }

	// Then count number of trophic levels, defining an inter-level connection
	// as having strength >
	double clim = 2.0 * bsd;
	/* The following results are from 1,000 trials, with 2 columns representing
	 * 50 (left) and 100 (right) species
		Results for 1.5 * bsd:
			1,
			2,	0.149,	0.413
			3,	0.811,	0.586
			4,	0.038	0.001
			5,	0.002
		Results for 2 * bsd:
			1,
			2,	0.142,	0.008
			3,	0.695,	0.751
			4,	0.152,	0.232
			5,	0.02,	0.008
			6,	0.001	0.001
		Results for 2.5 * bsd:
			1,	0.206	0.001
			2,	0.639,	0.215
			3,	0.137,	0.669
			4,	0.014,	0.103
			5,	0.003,	0.011
			6,	0.001,	0.001
		Results for 3 * bsd:
			1,			0.479,	
			2,			0.476,	
			3,			0.043,	
			4,			0.002,	
	*/

	boost::numeric::ublas::vector <int> troph_levels (nSpecies);
	for (int i=0; i<nSpecies; i++) 
    {
		troph_levels (i) = INT_MIN;
		for (int j=0; j<nSpecies; j++) 
        {
			if (j != i) {
				if (speciesData.compMat (i, j) > (bmn + clim))
					troph_levels (i) = 1;
				else if (speciesData.compMat (i, j) < (bmn - clim))
					troph_levels (i) = 0;
		}	} // end if j != i & for j
	} // end for i
	// This is useful for checking the effect of different threshold limits. For
	// the whole thing to work, there need to be "*" values.
	/*
	for (int i=0; i<nSpecies; i++) {
		if (troph_levels (i) == INT_MIN) { std::cout << "*";	}
		else { std::cout << troph_levels (i);	}
		std::cout.flush();
	}
	std::cout << std::endl;
	*/

	search_done = false;
	boost::numeric::ublas::vector <bool> bump_level (nSpecies);
	while (!search_done) {
		for (int i=0; i<nSpecies; i++) 
            bump_level (i) = false;
		nLevels = 0;
		for (int i=0; i<nSpecies; i++) 
        {
			if (troph_levels (i) == 0) 
            {
				for (int j=0; j<nSpecies; j++) 
                {
					if (j != i && 
                            speciesData.compMat (i, j) < (bmn - clim) && 
                            troph_levels (j) == 1) 
                    {
						for (int k=0; k<nSpecies; k++) 
                        {
							if (k != i && k != j && 
                                    speciesData.compMat (i, k) > (bmn + clim)) 
                            {
								bump_level (i) = true;
								nLevels++;
						}	} // end if & for k
				}	} // end if & for j
		}	} // end if & for i
		if (nLevels > 0) 
        {
			for (int i=0; i<nSpecies; i++) 
            {
				if (bump_level (i))
					troph_levels (i) = 1;
				else if (troph_levels (i) > 0) 
                    troph_levels (i) ++;
			}
		}
		else 
            search_done = true;
	} // end while (!search_done)
	nLevels = 0;
	for (int i=0; i<nSpecies; i++)
		if (troph_levels (i) > nLevels) 
            nLevels = troph_levels (i);

	return nLevels;
} // end SpeciesData makeCommunity
