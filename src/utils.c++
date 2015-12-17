/***************************************************************************
 *  Project:    multispecies-ensembles
 *  File:       utils.c++
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

RegrResults regression (dvec x, dvec y)
{
    double sx, sx2, sy, sy2, sxy, t1, t2, xmn, ymn;
    RegrResults regr_results;

    sx = 0.0; sx2 = 0.0;
    sy = 0.0; sy2 = 0.0; sxy = 0.0;
    int count = 0, n = x.size();
    for (int i=0; i<n; i++) 
    {
        if (x(i) > DOUBLE_MIN && y(i) > DOUBLE_MIN) 
        {
            count++;
            sx += x(i);
            sx2 += x(i) * x(i);
            sy += y(i);
            sy2 += y(i) * y(i);
            sxy += x(i) * y(i);
        }
    }
    xmn = sx / (double) count;
    ymn = sy / (double) count;
    if (count > 0) 
    {
        t1 = (sxy - sx * sy / (double) count);
        t2 = (sx2 - sx * sx / (double) count) * (sy2 - sy * sy / (double) count);
        regr_results.r2 = t1 / sqrt(t2); // the R-value
        regr_results.slope = t1 / (sx2 - sx * sx / (double) count); // Slope
        regr_results.intercept = sy / (double) count - 
            regr_results.slope * sx / (double) count; // Intercept

        regr_results.r2 = regr_results.r2 * regr_results.r2;
        if (regr_results.slope < 0) 
            regr_results.r2 = -regr_results.r2;

        // Then calculate SS and tval
        sy2 = 0.0; sx2 = 0.0; count = 0;
        regr_results.SS = 0.0;
        for (int i=0; i<n; i++) 
        {
            if (x(i) > DOUBLE_MIN && y(i) > DOUBLE_MIN) 
            {
                count++;
                t1 = regr_results.slope * x(i) + regr_results.intercept;
                regr_results.SS += (y(i) - t1) * (y(i) - t1);
                sx2 += (x(i) - xmn) * (x(i) - xmn);
                sy2 += (y(i) - ymn) * (y(i) - ymn);	
            }
        }
        if (count > 0) // tval calculation
        { 
            regr_results.SS = regr_results.SS / (double) count;
            regr_results.tval = sqrt(sy2 / ((double) count - 2)) / sqrt(sx2);
            regr_results.tval = regr_results.slope / regr_results.tval;
        } else {
            regr_results.SS = DOUBLE_MIN; 
            regr_results.tval = DOUBLE_MIN;	
        }
    } else { // if count == 0
        regr_results.r2 = DOUBLE_MIN;
        regr_results.slope = DOUBLE_MIN;
        regr_results.intercept = DOUBLE_MIN;
        regr_results.SS = DOUBLE_MIN;
        regr_results.tval = DOUBLE_MIN;
    }

    return regr_results;
} // end function regression


void timeout(double tseconds)
{
    int hh = floor(tseconds / 3600.0);
    if (hh == 0) 
        std::cout<<"00:";
    else if (hh < 10) 
        std::cout<<"0"<<hh<<":";
    else 
        std::cout<<hh<<":";
    double trem = tseconds - (double) hh * 3600.0;
    int mm = floor(trem / 60.0);
    if (mm == 0) 
        std::cout<<"00:";
    else if (mm < 10) 
        std::cout<<"0"<<mm<<":";
    else 
        std::cout<<mm<<":";
    double ss = trem - (double) mm * 60.0;
    if (ss == 0.0) 
        std::cout<<"00:";
    else if (ss < 10) 
        std::cout<<"0"<<ss;
    else 
        std::cout<<ss;
} // end function timeout
