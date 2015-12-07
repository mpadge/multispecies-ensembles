/*
 * model-envcor.cc
 *
 * Adapted to examine the 2 effects of:
 * 	(i)     Correlating the otherwise independent parts of each species'
 * 	        environmental variation (the x_i) according to strengths of
 * 	        interaction; and 
 * 	(ii)    Increasing values of alpha_i in proportion to average strength of
 * 	        interaction.  The first effect requires a scaling coefficient,
 * 	        called nscale, while the second is examined simply by sorting the
 * 	        allocated values of alpha_i.
 * 
 * Unsuprisingly, none of these really make much difference, as illustrated by
 * the following output, produced with fixed seed of 111, averaged over 100
 * trials, with 50 species and sync = (0.1, 0.2, 0.3, 0.5) for the four R2
 * values:
 * 
 * 	alpha adjusted		nscale	R2
 * 			sync values	->	0.1			0.2			0.3			0.5
 * 	--------------		------	-----------------------------------------------------
 * 		no			  0.0	0.28		0.66		0.83		0.93
 * 		yes			  0.0	0.28		0.65		0.83		0.93
 * 		no			  0.2	0.29		0.67		0.83		0.94
 * 		yes			  0.2	0.30		0.67		0.84		0.95
 * 	--------------		------	-----------------------------------------------------
 * The results for sync = 0.2 are the only ones that ever so slightly differ
 * from the general pattern of both effects increasing correlations. The
 * respective percentage increases w.r.t. the first value with neither effect
 * are:
 * 	alpha adjusted		nscale	R2
 * 			sync values	->	0.1			0.2			0.3			0.5
 * 	--------------		------	-----------------------------------------------------
 * 		yes			  0.0	2.36			-0.12		0.19			0.08
 * 		no			  0.2	3.98			1.91			0.95			0.32
 * 		yes			  0.2	6.35			1.79			1.13			0.40
 * 	--------------		------	-----------------------------------------------------
 * and the ratios of effect sizes of non-independent x_i (nscale=0.2) to
 * adjusting alpha_i values are then (1.69, 15.91, 5, 4) for sync = (0.1, 0.2,
 * 0.3, 0.5).
 */

#include "utils.h"
#include "pop-fns.h"
#include "model-envcor.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[])
{
	int nTrials, nSpecies, adj_rho;
    double sync, nscale, tempd [3], r2;
	SpeciesData speciesData;
	RegrResults regr;
	std::ifstream in_file;
	std::string linetxt, fname;
	std::stringstream ss_code;
	clock_t timer[2];
	time_t seed;
	base_generator_type generator(42u);

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("nTrials,n", boost::program_options::value <int> 
                (&nTrials)->default_value (10), "nTrials")
            ("sync,s", boost::program_options::value <double>
                (&sync)->default_value (0.5), "sync")
            ("nSpecies,p", boost::program_options::value <int> 
                (&nSpecies)->default_value (100), "nSpecies")
            ("nscale,l", boost::program_options::value <double>
                (&nscale)->default_value (0.2), "nscale")
            ("adj_rho,a", boost::program_options::value <int> 
                (&adj_rho)->default_value (0), "adj_rho")
            ;

        boost::program_options::options_description cmdline_options;
        cmdline_options.add(generic).add(config);

        boost::program_options::options_description visible("Allowed options");
        visible.add(generic).add(config);

        boost::program_options::variables_map vm;
        store(boost::program_options::command_line_parser(argc, argv).
                options(cmdline_options).run(), vm);

        notify(vm);

        if (vm.count("help")) {
            std::cout << visible << std::endl;
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "model_dependent_xi, version 1.0" << std::endl;
            return 0;
        }

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }    

	std::cout << "nTrials = " << nTrials << "; sync = " << sync <<
		"; nSpecies = " << nSpecies << "; nscale = " << nscale;
	if (adj_rho == 0) 
        std::cout << "; alpha will NOT be adjusted." << std::endl;
	else 
        std::cout << "; alpha WILL be adjusted." << std::endl;

	time (&seed);
	//generator.seed(static_cast<unsigned int>(seed));
	generator.seed(static_cast<unsigned int>(111));
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::uniform_real<> > runif(generator, uni_dist);
	boost::normal_distribution<> norm_dist (0.0, 1.0);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm(generator, norm_dist);
	// Burn both generators in
	for (int i=0; i<20; i++) 
    {
		tempd [0] = runif();
		tempd [0] = rnorm();
	}

	// Set up data to fixed size of nSpecies_max
	boost::numeric::ublas::vector<SpeciesPars> tempvec (nSpecies);
	boost::numeric::ublas::matrix<double> tempdmat (nSpecies, nSpecies);
	speciesData.spvec = tempvec;
	speciesData.compMat = tempdmat;
	boost::numeric::ublas::matrix<double> noise (len_t, nSpecies + 1);
	boost::numeric::ublas::vector<double> noise_global (len_t);
	boost::numeric::ublas::vector<double> pop_t (len_t);

	r2 = 0.0;
	for (int i=0; i<nTrials; i++) 
    {
		makeCommunity (nSpecies, sync, speciesData, generator, adj_rho);
		makeNoise (nSpecies, noise, generator);
		// Then adjust the noise according to strengths of species' interactions.
		// The adjustment is scaled by the above parameter, nscale.
		for (int j=0; j<(nSpecies - 1); j++) 
        {
			for (int k=(j+1); k<nSpecies; k++) 
            {
				tempd [0] = fabs (speciesData.compMat (j, k)) + 
                    fabs (speciesData.compMat (k, j));
				tempd [0] = (double) nscale * tempd [0] / (200.0 * bsd);
				for (int t=0; t<len_t; t++) 
                {
					tempd [1] = (1.0 - tempd [0]) * noise (t, j) + tempd [0] * noise (t, k);
					tempd [2] = (1.0 - tempd [0]) * noise (t, k) + tempd [0] * noise (t, j);
					noise (t, j) = tempd [1];
					noise (t, k) = tempd [2];
				} // end for t
            }
        } // end for j
		for (int j=0; j<len_t; j++)
			noise_global (j) = noise (j, nSpecies);
		runPop (nSpecies, noise, speciesData, pop_t);
		regr = regression (noise_global, pop_t);
		r2 += regr.r2;
		std::cout << "."; std::cout.flush();
	}
	std::cout << std::endl << "R2 = " << r2 / (double) nTrials << std::endl;

	return 0;
}; // end main
