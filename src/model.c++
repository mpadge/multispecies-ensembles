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


#include "model.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************
 *
 * NOTE that the whole program is designed to be compared to the theoretical
 * results which don't presume any kind of re-scaling, therefore all random
 * values, including noise and AC coefficients, have to just be taken as
 * generated, and may not be re-scaled in any way at all. This does make
 * competition rather odd, because a predator with negative abundance will
 * actually benefit prey. But that's also precisely the way the theory is
 * written, so has to be accepted.
 */
int main(int argc, char *argv[])
{
    int nTrials;
	double sync, sumr, sumr2, tempd;
	SpeciesData speciesData;
	RegrResults regr;
	std::ifstream in_file;
	std::string linetxt, fname;
	std::stringstream ss_code;
	std::ofstream out_file;
	clock_t timer[2];
	time_t seed;
	base_generator_type generator(42u);

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "version 1.0.1")
            ("help", "this is not helpful")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("nTrials,n", boost::program_options::value <int> 
                (&nTrials)->default_value (10), "nTrials")
            ("sync,s", boost::program_options::value <double>
                (&sync)->default_value (0.5), "sync")
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

    std::cout << "nTrials = " << nTrials << "; sync = " << sync << std::endl;

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

    // Set up data to fixed size of nSpecies_max
    boost::numeric::ublas::vector<SpeciesPars> tempvec (nSpecies_max);
    boost::numeric::ublas::matrix<double> tempdmat (nSpecies_max, nSpecies_max);
    speciesData.spvec = tempvec;
    speciesData.compMat = tempdmat;
    dvec noise_global (len_t);
    dvec pop_t (len_t);
    boost::numeric::ublas::matrix<double> noise (len_t, nSpecies_max + 1);
    boost::numeric::ublas::matrix<double> r2 (nSpecies_max, nTrials);
    // Averaging correlations is kinda wrong to being with, and there's no
    // sure way of saying whether one should average R or R2. At the moment,
    // it's R2 values, but individuals ones are stored in a matrix to allow easy
    // conversion to the alternative if desired.

    std::stringstream ss_ntrials, ss_sync;
    ss_ntrials.str (""); ss_ntrials << nTrials;
    ss_sync.str (""); ss_sync << floor (100.0 * sync);
    fname = "results_" + ss_ntrials.str() + "trials_sync";
    if (sync < 1.0) 
        fname += "0";
    if (sync < 0.1) 
        fname += "0";
    fname += ss_sync.str() + ".txt";

    out_file.open (fname.c_str(), std::ofstream::out);
    out_file << "nspecies,\tr2mn,\tr2sd" << std::endl;
    for (int nSpecies = 2; nSpecies <= nSpecies_max; nSpecies++) {
        sumr = 0.0;
        sumr2 = 0.0;
        for (int i=0; i<nTrials; i++) {
            makeCommunity (nSpecies, sync, speciesData, generator, 0);
            makeNoise (nSpecies, noise, generator);
            for (int j=0; j<len_t; j++) 
                noise_global (j) = noise (j, nSpecies);
            runPop (nSpecies, noise, speciesData, pop_t);
            regr = regression (noise_global, pop_t);
            r2 (nSpecies - 1, i) = regr.r2;
            sumr += r2 (nSpecies - 1, i);
            sumr2 += r2 (nSpecies - 1, i) * r2 (nSpecies - 1, i);
        }
        sumr = sumr / (double) nTrials;
        sumr2 = sumr2 / (double) nTrials - sumr * sumr;
        sumr2 = sqrt (sumr2);
        out_file << nSpecies << ",\t" << sumr << ",\t" << sumr2 << std::endl;

        tempd = 100.0 * (nSpecies - 1) / nSpecies_max;
        std::cout << "\r" << nSpecies << " / " <<
            nSpecies_max << " = " << tempd << "%";
        std::cout.flush ();
    }
    std::cout << " done." << std::endl;
    out_file.close();

    return 0;
}; // end main
