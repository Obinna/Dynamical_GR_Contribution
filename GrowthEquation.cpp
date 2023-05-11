#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <vector>
using std::vector;

#include <boost/bind.hpp>

#include "GrowthEquation.h"

//#include "/home/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/ODE.h"

//#include "/Users/obinna/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/ODE.h"

 #include "/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/ODE.h"
//#include "array.h"

//#include "MySpline.cpp"


GrowthEquation::GrowthEquation() 
{
}

double Omegag(double a)
{
    Cosmology C;
   double z = 1.0/a -1.0;
   
    double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
    return tmpz;
}
double CHg(double a)
{
    Cosmology C;
    double z = 1.0/a -1.0;
    double c_light = C.H_0*3000.0/C.h;
    // double c_light =  C.H_0*3000.0;
    return C.H(z)/(c_light*(1.0+z));
}


double dHda(double a)
{
        double delta = 1.0e-4;
	double x1 = a -delta;
	double x2 = a + delta ;
	double y1 = CHg(x1);
	double y2 = CHg(x2);
	double deri = (y2-y1)/(x2-x1);
        return deri;
}

vector<double> GrowthEquation::dXda(double a, const vector<double>& X) 
{   
    Cosmology C;
    double z = 1.0/a -1.0;
    double DD = C.D_growth_Tabulated(z);
    vector<double> dX(2);
    dX[0] = X[1];
    dX[1] = -(2/a + dHda(a)/CHg(a)) * X[1] + 1.5*Omegag(1.0)*pow(CHg(1.0),2)/(pow(a,3)*pow(CHg(a),2)) * (X[0] + DD*DD);

//cout<<"\t"<<"\t"<<"dX[0]"<<"\t"<<dX[0]<<"\t"<<"dX"<<"\t"<<dX[1]<<endl;
    return dX;
}
const int NM = 500;
vector<double> A(NM+1);
vector<vector<double> > X;

void GrowthEquation::InitialzeFloop(void) {
    //C = C_;
    
    const  double a_i = 1.0/NM;
    const double a_f = 1.0 + 1.0/NM;;

    /* Solve for CurlyF(a) between a_i and a_f, assuming F(a) ~ a for small a */
    
    for(int i = 0; i <= NM; i++)
        A[i] = a_i + i*(a_f - a_i)/NM;
    vector< double> x0(2);
    x0[0] = a_i;
    x0[1] = 1;
     X = RungeKutta4(a_i, a_f, x0, bind(&GrowthEquation::dXda, this, _1, _2), NM);

    /* Interpolate to a = 0 */
    A.insert(A.begin(), 0.);
    X[0].insert(X[0].begin(), 0.);
    X[1].insert(X[1].begin(), 0.);

    /* Normalize to D=1 at a=1 */
     //F = CubicSpline(A, X[0]);
     Spline<double, double> CubicSpline_f(A,X[0]);
    
    double F0 =  CubicSpline_f.interpolate(1);///F(1);
    for(int i = 0; i <= NM+1; i++) 
    {
        X[0][i] /= F0;
        X[1][i] /= F0;
   /// cout<<"\t i"<<"\t"<<i<<"\t"<<"F"<<"\t"<<X[0][i]<<"\t"<<"G"<<"\t"<<X[1][i]<<endl;
    }

    //F = CubicSpline(A, X[0]);
    //dFda = CubicSpline(A, X[1]);
}


double SplineCurlyF(double a)
{
	Spline<double, double> CubicSpline_curlyF(A,X[0]);
	return CubicSpline_curlyF.interpolate(a);
   	
}

double DcurlyFda(double a)
{
	Spline<double, double> CubicSpline_deri(A,X[1]);
	return CubicSpline_deri.interpolate(a);
   	
}

/*
GrowthEquation::GrowthEquation(const Cosmology& C_, double a_i, double dDda_i) {
    C = C_;
    const int N = 10000;
    const double a_f = 1.01;

    vector<double> A(N+1);
    for(int i = 0; i <= N; i++)
        A[i] = a_i + i*(a_f - a_i)/N;
    vector<double> x0(2);
    x0[0] = a_i;        // D_i = a_i
    x0[1] = dDda_i;
    vector<vector<double> > X = RungeKutta4(a_i, a_f, x0, bind(&GrowthEquation::dXda, this, _1, _2), N);

    
    D = CubicSpline(A, X[0]);
    double D0 = D(1);
    for(int i = 0; i <= N; i++) {
        X[0][i] /= D0;
        X[1][i] /= D0;
    }

    D = CubicSpline(A, X[0]);
    dDda = CubicSpline(A, X[1]);
}
*//* Normalize to D=1 at a=1 */


double GrowthEquation::CurlyF(double z) 
{
    double a = 1/(1+z);
    return SplineCurlyF(a)/SplineCurlyF(1.0);
}

double GrowthEquation::CurlyG(double z) 
 {
    double a = 1/(1+z);
    return a/SplineCurlyF(a) * DcurlyFda(a);
}


