/*
 * utils.c++
 *
 */

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
		if (x(i) > dnix && y(i) > dnix) 
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
			if (x(i) > dnix && y(i) > dnix) 
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
			regr_results.SS = dnix; regr_results.tval = dnix;	
        }
    } else { // if count == 0
		regr_results.r2 = dnix;
		regr_results.slope = dnix;
		regr_results.intercept = dnix;
		regr_results.SS = dnix;
		regr_results.tval = dnix;
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
