/***************************************************************************
 *  Project:    multispecies-ensembles
 *  File:       model050-theoretical.c++
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


#include "model050-theoretical.h"

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


// Dumps the file "rholims_for_r050_theoretical.txt"

int main()
{
    double rho [3], r2;
    std::ofstream out_file;

    int count = 0;

    out_file.open ("rholims_for_r050_theoretical.txt");
    out_file << "n,\trho" << std::endl;
    for (int nSpecies = 2; nSpecies <= nSpecies_max; nSpecies++) {
        // This searche uses the list of 3 parameters sorted as [low, high, middle]
        rho [0] = 0.01; rho [1] = 0.99;

        while (fabs (rho [1] - rho [0]) > 0.001) {
            rho [2] = (rho [0] + rho [1]) / 2.0;
            r2 = calcR2 (rho [2], nSpecies);
            if (r2 > 0.5) 
                rho [1] = rho [2];
            else 
                rho [0] = rho [2];
        } // end while
        out_file << nSpecies << ",\t" << rho [2] << std::endl;
    }
    out_file.close();

    return 0;
}; // end main

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         SUBFUNCTIONS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


double calcR2 (double rho, int n)
{
    double nd = (double) n;

    double A = sigma * sigma + 
        2.0 * (nd - 1.0) * rho * rho * bmn * sigma * sigma / 
        (1.0 - (nd - 1.0) * bmn) +
        2.0 * (nd - 1.0) * (1.0 - rho * rho) * sigma * sigma * bmn * bmn / 
        (1.0 - (nd - 2.0) * bmn - (nd - 1.0) * bmn * bmn);
    double B = (nd - 1.0) * (nd - 2.0) * bmn * bmn;
    double C = 1.0 - (nd - 1.0) * (bmn * bmn + bsd * bsd);

    double D = rho * rho * sigma * sigma +
        2.0 * (nd - 1.0) * rho * rho * bmn * sigma * sigma / 
        (1.0 - (nd - 1.0) * bmn) +
        2.0 * (nd - 1.0) * (1.0 - rho * rho) * bmn * bmn * sigma * sigma /
        (1.0 - (nd - 2.0) * bmn - (nd - 1.0) * bmn * bmn);
    double E = (nd - 2.0) * bmn * bmn;
    double F = 1.0 - (nd - 1.0) * (nd - 2.0) * bmn * bmn - bmn * bmn + bsd * bsd;

    double ni2 = (A * F + B * D) / (C * F - B * E);
    double ninj = (D + E * ni2) / F;

    double s0n = nd * rho * sigma * sigma / (1.0 - (nd - 1.0) * bmn);
    double sum_ni2 = nd * nd * rho * rho * sigma * sigma + 
        nd * (1.0 - rho * rho) * sigma * sigma + 
        2.0 * nd * nd * (nd - 1.0) * rho * rho * sigma * sigma * bmn / 
        (1.0 - (nd - 1.0) * bmn) + 
        2.0 * nd * (nd - 1.0) * (1.0 - rho * rho) * sigma * sigma * bmn * 
        (1.0 + bmn) / (1.0 - (nd - 2.0) * bmn - (nd - 1.0) * bmn * bmn) + 
        nd * (nd - 1.0) * ((nd - 1.0) * bmn * bmn + bsd * bsd) * ni2 + 
        nd * (nd - 1.0) * ((nd - 1.0) * bmn * bmn - bsd * bsd) * ninj;

    double R2 = s0n * s0n / (sigma * sigma * sum_ni2);

    return R2;
}
