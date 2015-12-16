/*
 * model050.c++
 * 
 * Estimates limits of rho^2 above which correlations exceeded R^2=50%
 * 
 * Dumps the file "rholims_for_r050.txt" containing ***R^2*** (not R!) values
 */

#include "model050.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char *argv[])
{
    int nTrials, count = 0;
    double sync, rho [2], r2, tempd;
    std::ofstream out_file;
    SpeciesData speciesData;
    RegrResults regr;
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
    std::cout << "nTrials = " << nTrials << std::endl;

    time (&seed);
    generator.seed(static_cast<unsigned int>(seed));
    //generator.seed(static_cast<unsigned int>(2));
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

    boost::numeric::ublas::matrix<double> noise (len_t, nSpecies_max + 1);
    boost::numeric::ublas::vector<double> noise_global (len_t);
    boost::numeric::ublas::vector<double> pop_t (len_t);
    boost::numeric::ublas::vector<SpeciesPars> tempvec (nSpecies_max);
    boost::numeric::ublas::matrix<double> tempdmat (nSpecies_max, nSpecies_max);
    speciesData.spvec = tempvec;
    speciesData.compMat = tempdmat;

    std::stringstream ss;
    ss.str ("");
    ss<< nTrials;
    std::string fname = "rholims_for_r050_ntrials" + ss.str () + ".txt";
    out_file.open (fname);
    out_file << "n,\trho" << std::endl;
    for (int nSpecies = 2; nSpecies <= nSpecies_max; nSpecies++) {
        // This searche uses the list of 3 parameters sorted as [low, high, middle]
        rho [0] = 0.01; rho [1] = 0.99;

        while (fabs (rho [1] - rho [0]) > 0.001) {
            sync = (rho [0] + rho [1]) / 2.0;

            r2 = 0.0;
            for (int i=0; i<nTrials; i++) {
                makeCommunity (nSpecies, sync, speciesData, generator, 0);
                makeNoise (nSpecies, noise, generator);
                for (int j=0; j<len_t; j++) {
                    noise_global (j) = noise (j, nSpecies);
                }
                runPop (nSpecies, noise, speciesData, pop_t);
                regr = regression (noise_global, pop_t);
                r2 += regr.r2;
            }
            r2 = r2 / (double) nTrials;
            if (r2 > 0.5) 
                rho [1] = sync;
            else 
                rho [0] = sync;
        } // end while
        sync = sync * sync;
        out_file << nSpecies << ",\t" << sync << std::endl;

        tempd = 100.0 * ((double) nSpecies - 1.0) / 
            ((double) nSpecies_max - 1.0);
        std::cout << "\r" << nSpecies << " / " << nSpecies_max <<
            " = " << tempd << "%";
        std::cout.flush ();
    }
    out_file.close();
    std::cout << " done." << std::endl;

    return 0;
}; // end main
