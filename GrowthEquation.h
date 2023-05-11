#ifndef GROWTH_EQUATION_H
#define GROWTH_EQUATION_H

#include <vector>

/*
#include "/home/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/Common.h"
#include "/home/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/MyCosmology.h"
*/




/*

#include "/Users/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/Common.h"
#include "/Users/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/MyCosmology.h"
*/



#include "/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/Common.h"
#include "/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/MyCosmology.h"


//#include "MySpline.cpp"


class GrowthEquation {
public:
    GrowthEquation();
  void  InitialzeFloop(void);
    //GrowthEquation(const Cosmology& C, double a_i, double Dp_i);

    double Evaluate(double z) const;
    double operator()(double z) const { return Evaluate(z); }

    /* Logarithmic growth rate: G = dlogF/dloga */
    double CurlyF(double z);
    double CurlyG(double z) ;

private:
   // Cosmology C;

    //double F;           // F(a)
   // double dFda;        // dF/da

    /* X = (D, dD/da) */
    std::vector<double> dXda(double a, const std::vector<double>& X);
};










#endif
