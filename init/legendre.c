/* Wrapper routine to call GSL library function to compute Associated Legendre */
/* polynomial function $P_l^m(x)$ to be used in calculating corresponding      */
/* spherical harmonic $Y_l^m(x,y)$                                             */

# include <stdio.h>
# include <gsl/gsl_sf_legendre.h>
# include "legendre.h"

void getsphplm(const int *l, const int *m, const double *x, double *Plm)
{

    int l_int, m_int;
    double x_dbl, Plm_dbl;

    l_int = *l;
    m_int = *m;
    x_dbl = *x;   

    Plm_dbl = *Plm;

    Plm_dbl = gsl_sf_legendre_sphPlm(l_int, m_int, x_dbl);

//    printf("%d %d %g \n",l_int,m_int,Plm_dbl);

    *Plm = Plm_dbl;

    return;
}
