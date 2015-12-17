/***************************************************************************
 *  Project:    multispecies-ensembles
 *  File:       model-envcor.c++
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
 *
 *  Note that effect (ii) is implemented directly in pop-fns, while this code
 *  implemenets (i). This requires a scaling coefficient, called nscale.
 *  Although it works, it only generates positive increases for very small
 *  values of nscale (<0.001), with resultant R2 values strongly dependent on
 *  the chosen value. 
 *
 ***************************************************************************/

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
             (&nTrials)->default_value (1), "nTrials")
            ("sync,s", boost::program_options::value <double>
             (&sync)->default_value (0.1), "sync")
            ("nSpecies,p", boost::program_options::value <int> 
             (&nSpecies)->default_value (10), "nSpecies")
            ("nscale,l", boost::program_options::value <double>
             (&nscale)->default_value (0.0), "nscale")
            ("adj_rho,a", boost::program_options::value <int> 
             (&adj_rho)->default_value (0), 
             "adj_rho (non-zero values flag adjustment)")
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
        std::cout << "; rho will NOT be adjusted." << std::endl;
    else 
        std::cout << "; rho WILL be adjusted." << std::endl;

    time (&seed);
    generator.seed(static_cast<unsigned int>(seed));
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
    boost::numeric::ublas::matrix<double> noise_adj (len_t, nSpecies);
    boost::numeric::ublas::vector<double> noise_global (len_t);
    boost::numeric::ublas::vector<double> pop_t (len_t);

    r2 = 0.0;
    for (int i=0; i<nTrials; i++) 
    {
        makeCommunity (nSpecies, sync, speciesData, generator, adj_rho);
        makeNoise (nSpecies, noise, generator);
        // Then adjust the noise according to strengths of species' interactions:
        for (int j=0; j<nSpecies; j++) 
            for (int k=0; k<nSpecies; k++) 
                if (k != j)
                {
                    tempd [0] = fabs (speciesData.compMat (j, k)) + 
                        fabs (speciesData.compMat (k, j));
                    tempd [0] = nscale * tempd [0] / (2.0 * bsd);
                    for (int t=0; t<len_t; t++)
                        noise_adj (t, j) = tempd [0] * noise (t, k) +
                            sqrt (1.0 - tempd [0] * tempd [0]) * noise (t, j);
                } // end if k != j
        for (int j=0; j<nSpecies; j++) 
            for (int t=0; t<len_t; t++)
                noise (t, j) = noise_adj (t, j);

        for (int j=0; j<len_t; j++)
            noise_global (j) = noise (j, nSpecies);
        runPop (nSpecies, noise, speciesData, pop_t);
        r2 += regression (noise_global, pop_t).r2;
    }
    r2 = r2 / (double) nTrials;
    std::cout << "R2 = " << r2 << std::endl;

    return 0;
}; // end main
