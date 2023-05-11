#include "DynamicalGRBispectrum.h"

//#include "GrowthEquation.h"
//#include "MyBessel.cpp"

#include <boost/bind.hpp>
#include <stdio.h>
#include <stdlib.h>



#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>



DynamicalGR::DynamicalGR()
{
    //SetPowerSpectrum();
}

DynamicalGR::~DynamicalGR()
{
    //SetPowerSpectrum();
}




double HH(double z)
{
    Cosmology C;
    double c_light = C.H_0*3000.0;
    //const  double c_light_km_s =  299792;
    // double c_light =  C.H_0*3000.0;
    return C.H(z)/(c_light*(1.0+z));
}

double Hpr_in(double z)
{
  //Cosmology C;
	double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = HH(x1);
	double y2 = HH(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*HH(z)*deri;
	return tmp;
}


double Omegam_in(double z)//Omegam
{
    Cosmology C;
    double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
    return tmpz;
}

double fg(double z)
{
  //Cosmology C;
  //double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
   //return pow(tmpz,0.55);
  return pow(Omegam_in(z),0.55);
}
////power spectrum from camb////////////////////////////////

////Power spectrum from Class


const int N = 587; ///You have to make sure this number corresponds to the number of rows in class//cat explanatory13_pk.dat | wr -l
std::vector <double> k_array(N); // define vectors to k-values
std::vector <double> pk_array(N);  // define vectors to hold pofk


void DynamicalGR::ReadWrite(void)
{
  ifstream inFile; //class that reads the file
    string line; ///to take care of the header comments in .dat
  inFile.open("explanatory02_pk.dat");///load class file
 
  if(!inFile)
    {
      cerr << "could not open file" << endl; ///check whether the file is readable
    }
 while (getline(inFile, line))
    {
      if (line[0] == '#') continue;///Tell compiler to ignore comments in class's .dat file

        for(int i = 1; i < N; i++)
           {
  inFile >> k_array[i] >> pk_array[i]; /// Save the content of the file in k_array and pk_array

 //cout<<"\t"<<"k="<<"\t"<<k_array[i]<<"\t"<<""<<"pk="<<"\t"<<pk_array[i]<<endl;//check whether it is doing the right thing.
          }
    }
    inFile.close();
}

double MySplinePofk(double k)
{
  Spline<double, double> CubicSpline_num(k_array,pk_array);
  double res = CubicSpline_num.interpolate(k);
  return res;
}


////////////////




/////////////////////Simple Bias formular////////////////////


double bias10(double z)
{ 
  return sqrt(1.0+z);
}

double bias10pr(double z)
{
  double delta = 1.0e-4;
  double x1 = z - delta;
  double x2 = z + delta ;
  double y1 = bias10(x1);
  double y2 = bias10(x2);
  double deri = (y2-y1)/(x2-x1);
  double tmp = deri;
  return tmp;
}

double bias20(double z)
{
  
      double b0 = 1.0;
      double tmp = -0.25*b0 - 0.13*sqrt(1.0+z)*sin(0.8*z);
      return tmp;
        //return -bias1(z)/2.0;
}

double TidalBias(double z)
 {  
    double tmp = -4.0*(bias10(z)-1.0)/7.0;
    return tmp;
 }

 
 double bias01(double z, double fnl)
{ 
  double delta_c = 1.686;
  return 2.0*fnl*delta_c*(bias10(z) - 1.0);
}
double bias02(double z, double fnl)
{
    double delta_c = 1.686;
    double tmp0 = (bias10(z) - 1.0);
    double tmp  = 4.0*pow(fnl,2)*delta_c*(delta_c*bias20(z) - 2.0*((4.0/21.)*delta_c + 1.0)*tmp0);
    return tmp;
}

double hatbias11(double z, double fnl)
{
  double delta_c = 1.686;
    double tmp0 = (bias10(z) - 1.0);
    double tmp  = 4.0*pow(fnl,2)*delta_c*(delta_c*bias20(z) + ((11.0/21.)*delta_c - 1.0)*tmp0);
    return tmp;
}

double anl(double fnl)
{
 double tmp = 1.0 + 3.0*fnl/5.0;
 return tmp;
}

double bias11T(double z, double fnl)
{
  double tmp =  2.0* hatbias11(z,fnl) + 20.0*(2.0 -  anl(fnl))/3.0;
  return tmp;
}
double bias11N(double z, double fnl)
{
  double tmp =  2.0* hatbias11(z,fnl) + 20.0*(1.0 -  anl(fnl))/3.0;
  return tmp;
}

double biasvsqN(double z, double fnl)
{
  double tmp =  2.0*bias01(z, fnl);
  return tmp;
}

double biasvsqT(double z, double fnl)
{
  double tmp =  2.0* bias01(z, fnl) + 20.0*(3.0/4.0 -  anl(fnl))/3.0;
  return tmp;
}


double Alpha( double k,double z)
{
     
      Cosmology C;
    // double g0_overg_infty = 3.0/4.0;
    // double c_light =  C.H_0*3000.0;
    // double c_light =  C.H_0*3000.0;
    //const  double c_light_km_s =  299792;
     //double numerator = 2.0*k*k*c_light*c_light*C.D_growth_Tabulated(z)*C.Tk_tab(k);
     //double denominator = 3.0*C.Omega_m*pow(C.H_0,2) ;

     //double tmp= -g0_overg_infty*numerator/denominator;

//cout <<"\t"<<"Alpha="<<tmp<<endl;
     //return tmp;



     double fac = pow((1.0 + (2./3.)*(fg(z)/Omegam_in(z))),-1.0);
     return -(10.0/9.0)*(k*k*C.Tk_tab(k))/(Omegam_in(z)*HH(z)*HH(z))*fac;
}


  double AlphaAlpha2(double z, double  k1, double k2)
{
  return Alpha(k1,z)*Alpha(k2,z);

}

  double b_one(double z, double fnl, double k)
{
    
    Cosmology C;
    //double delta_c = 1.686;
    double HIb10 = bias10(z);
    double DeltaB = bias01(z,fnl)/Alpha(k, z);
  // double tmp  = HIb10+ (HIb01/alpha_z(k, z));
    double tmp = HIb10 + DeltaB;
   return  tmp;
}


double b_twoT(double z,double fnl ,double k1,double k2)
{   
   
     Cosmology C;
    // double delta_c =  1.686;
     //double HIb10 = bias10(z);  
     double HIb20 = bias20(z); 
     double b11 = bias11T(z,fnl);
     double b02 = bias02(z,fnl);

     double DeltaSqB_Ist1 = b11/Alpha(k1, z);
    double DeltaSqB_Ist2 = b11/Alpha(k2, z);
     double DeltaSqB_2nd = b02/AlphaAlpha2(z,  k1, k2);

     double tmp =  HIb20 + (DeltaSqB_Ist1+DeltaSqB_Ist2) + DeltaSqB_2nd;

     return tmp;

} 
double b_twoN(double z,double fnl ,double k1,double k2)
{   
   
     Cosmology C;
    // double delta_c =  1.686;
     //double HIb10 = bias10(z);  
     double HIb20 = bias20(z); 
     double b11 = bias11N(z,fnl);
     double b02 = bias02(z,fnl);

     double DeltaSqB_Ist1 = b11/Alpha(k1, z);
    double DeltaSqB_Ist2 = b11/Alpha(k2, z);
     double DeltaSqB_2nd = b02/AlphaAlpha2(z,  k1, k2);

     double tmp =  HIb20 + (DeltaSqB_Ist1+DeltaSqB_Ist2) + DeltaSqB_2nd;

     return tmp;

} 



///////////////////////////////////Linear galaxy power spectrum//////////////////////////////////

double DynamicalGR::LinearPowerSpectrum(int x, double z,double fnl, double k)
{
     return PowerSpectrum(x,z,fnl, k);
}


double PowerSpectrum(int x, double z, double fnl, double k)
{
      Cosmology C;
     double Dz = C.D_growth_Tabulated(z);
     //double Pk=  pow(Dz,2)*C.P_k_EH_Tabulated(k); //
     double Pk=  pow(Dz,2)*MySplinePofk(k);
   if(x == 1)
       {
         return pow(b_one(z,fnl, k),2)*Pk;
       }
    else if(x == 2)
      {
        double tmp = b_one(z,fnl, k)+  3.0*fg(z)*pow(HH(z),2)/(k*k);
         return  tmp*tmp*Pk;
      }
    else 
     return Pk;
}

  double GalaxyP1P2(int x, double z, double fnl, double  k1, double k2)
  {   
    double tmp =  PowerSpectrum(x, z, fnl, k1)*PowerSpectrum(x, z, fnl , k2);
   return tmp;
   }
  


///////////Tidial bias kernel
double S_kernel(double k1, double k2, double mu12)
{ 
  double tmp = pow(mu12,2) - 1.0/3.0;
   return tmp;
}
double NewtonianF2(double k1, double k2, double mu12)
{
   double first_term = 5.0/7.0;
   double second_term = (k1/k2 + k2/k1)*mu12/2;
   double third_term = (2.0/7.0)*pow(mu12,2);

   double tmp =  first_term + second_term + third_term;
   return 2.0*tmp;
}

double GRF2(double z, double fnl, double k1, double k2, double mu12)
{
  double anl = 3.0*fnl/5.0 + 1.0;
  double coef1 = -3.0*(1 + 2.0*fg(z)/(3.0*Omegam_in(z)) )*Omegam_in(z)*pow(HH(z),2)/2.0;
  double tmp1 = 4.0*(3.0/4.0- anl)*mu12*k1*k2;
  double tmp2 = 2.0*(2.0-anl)*(k1*k1 + k2*k2);
  double tmp = coef1*(tmp1 + tmp2)/pow(k1*k2,2);
  return tmp;
}
double GRF2Newtonian(double z, double fnl, double k1, double k2, double mu12)
{
  double anl = 3.0*fnl/5.0 + 1.0;
  double coef1 = -3.0*(1 + 2.0*fg(z)/(3.0*Omegam_in(z)) )*Omegam_in(z)*pow(HH(z),2)/2.0;
  double tmp1 = 2.0*(1.0-anl)*(k1*k1 + k2*k2);
  double tmp = coef1*(tmp1)/pow(k1*k2,2);
  return tmp;
}
double NewtonianG2(double k1, double k2, double mu12)
{
   double first_term = 3.0/7.0;
   double second_term = (k1/k2 + k2/k1)*mu12/2.0;
   double third_term = (4.0/7.0)*pow(mu12,2);
   double tmp =  first_term+second_term+third_term;
   return 2.0*tmp;
}

double G2GR(double z,double fnl, double k1, double k2, double k3,double muk)
{
   //double k3 = sqrt(k1*k1 + k2*k2 + 2.0*k1dotk2(k1, k2,  muk));
   double first_term = 3.0;
   double second_term = 2.0*(k1/k2 + k2/k1)*muk;
   double third_term = pow(muk,2);

   double tmp1 =  2.0*(first_term + second_term + third_term)/pow(k3,2);
   double coef =  (k3*k3)/(k1*k1*k2*k2);
   double tmp2 = (-3.0 + (6.0* fnl/5.0)* (1.0 + 2.0*fg(z)/(3.0*Omegam_in(z))))*coef;
   double tmp  = (3.0*Omegam_in(z)*pow(HH(z),2)/2.0)*(tmp1 + tmp2);
   return tmp;
}

double boneprz(double z, double fnl, double k1, double k2)
{
  double delta_c = 1.686;
 double tmp =  bias10pr(z)*( (1. + 2.0*fnl*delta_c/Alpha(k1, z))*k2*k2 + (1. + 2.0*fnl*delta_c/Alpha(k1, z))*k1*k1);
  return tmp;
}

double CurlyY(double z, double fnl, double k1, double k2)
{
 double tmp = (b_one(z,fnl,k2)*k1*k1 +  b_one(z,fnl, k1)*k2*k2)/2.0 + (1.+z)*boneprz(z,fnl,k1,k2)/6.0 + bias10(z)*fg(z)*(k1*k1 + k2*k2)/6.0;
 return tmp;
}


double G2vv(double z, double fnl, double k1, double k2, double k3,double muk)
{
   
   double tmp0 = 2.0*k3*k3/(k1*k1*k2*k2);
  double tmp1 = (fg(z)/(Omegam_in(z)))*(1.0+ 0.5*Omegam_in(z))  + 1.0;
   double tmp  = CurlyY(z,fnl, k1, k2) + (1.5*Omegam_in(z)*pow(HH(z),2))*(tmp1);
   return tmp0*tmp;
}

double CurvyN2(double z, double k1, double k2, double muk)
{
  double tmp = k1*k2*muk/(pow(k1,2)*Alpha(k2, z)) + k1*k2*muk/(pow(k2,2)*Alpha(k1, z));
  return -0.5*tmp;

}

double G2vvmatter(double z, double fnl, double k1, double k2, double k3,double muk)
{
   
   double tmp0 = 2.0*k3*k3/(k1*k1*k2*k2);
  double tmp1 =  1.5*Omegam_in(z)*pow(HH(z),2)*(((fg(z)/Omegam_in(z)))*(1.0+ 0.5*Omegam_in(z)) + 1.0);

   double tmp2  = (3.0 + fg(z))*(k1*k1 + k2*k2)/6.0 ;
   double tmp = tmp0*(tmp1 +tmp2);
   return tmp;
}

double SecondOrderMatterKernel(int x, double z, double fnl, double k1, double k2, double k3,double mu12)
{    
    if(x == 1)
       {
         return NewtonianF2(k1,k2,mu12);
       }
       else if(x == 2)
       {
        return GRF2(z, fnl, k1, k2,mu12);
       }
       else if(x == 3)
       {
        return GRF2Newtonian(z,fnl,k1,k2,mu12);
       }
        else if(x == 4)
       {
        return NewtonianF2(k1,k2,mu12) + GRF2(z, fnl, k1, k2,mu12);
       }
       else if(x == 5)
       {
        return NewtonianF2(k1,k2,mu12) + GRF2Newtonian(z, fnl, k1, k2,mu12);
       }
       else if(x ==6)
      {
        double tmp1 = NewtonianG2(k1, k2,mu12) + G2GR(z, fnl,k1,k2, k3,mu12) + G2vvmatter(z, fnl, k1,  k2,  k3, mu12);
        double tmp2 =  3.0*fg(z)*pow(HH(z),2)*tmp1/(k3*k3); 
        return tmp2;
      }
       else 
      {
        double tmp1 = NewtonianG2(k1, k2,mu12) + G2GR(z,fnl, k1,k2, k3,mu12) + G2vvmatter(z,fnl, k1,  k2,  k3, mu12);
        double tmp3 =  3.0*fg(z)*pow(HH(z),2)*tmp1/(k3*k3);  
        double tmp = NewtonianF2(k1,k2,mu12) + GRF2(z, fnl, k1, k2,mu12) + tmp3;//
        return tmp;
      }
}

double LinearOrderMatterKernel(int x, double z, double fnl, double k1)
{    
    if(x == 1)
      {
        return 1.0;
      }
      else if(x == 2)
      {
       return 1.0+ 3.0*fg(z)*pow(HH(z),2)/(k1*k1);
      }
      else
        return 0.0;
}



double SecondOrderKernelRealSpaceTotalMatter(int x, double z, double fnl, double k1, double k2, double k3,double mu12)
{    
    if(x == 1)
       {
         return bias10(z)*NewtonianF2(k1,k2,mu12);
       }
    else if(x ==2)
      {
         return b_twoT(z, fnl,k1,k2);
      }
      else if(x ==3)
      {
         return b_twoN(z, fnl,k1,k2);
      }
    else if(x ==4)
      {
         return TidalBias(z)*S_kernel(k1, k2,mu12);
      }
    else if(x == 5)
      {
         return  biasvsqT(z, fnl)*CurvyN2(z,k1,k2,mu12);
      }
    else if(x ==6)
      {
         return biasvsqN(z, fnl)*CurvyN2(z,k1,k2,mu12);
      }
    else if(x ==7)
      {
         return bias10(z)*NewtonianF2(k1,k2,mu12) + b_twoN(z, fnl,k1, k2) + TidalBias(z)*S_kernel(k1, k2,mu12) + biasvsqN(z, fnl)*CurvyN2(z,k1,k2,mu12);
      }

else 
      {
         return bias10(z)*NewtonianF2(k1,k2,mu12) + b_twoT(z, fnl, k1,k2) + TidalBias(z)*S_kernel(k1, k2,mu12) + biasvsqT(z, fnl)*CurvyN2(z,k1,k2,mu12);//
      }

}


double SecondOrderKernelRealSpacePoisson(int x,double z, double fnl, double k1, double k2, double k3,double mu12)
{    
     if(x == 1)
       {
         return bias10(z)*NewtonianF2(k1,k2,mu12);
       }
    else if(x ==2)
      {
         return b_twoT(z, fnl,k1,k2);
      }
    else if(x ==3)
      {
         return TidalBias(z)*S_kernel(k1, k2,mu12);
      }
    else if(x == 4)
      {
         return  biasvsqT(z, fnl)*CurvyN2(z,k1,k2,mu12);
      }
    else if(x ==5)
      {
        double tmp1 = NewtonianG2(k1, k2,mu12) + G2GR(z, fnl,k1,k2, k3,mu12) + G2vv(z, fnl, k1,  k2,  k3, mu12);
        double tmp2 =  3.0*fg(z)*pow(HH(z),2)*tmp1/(k3*k3); 
        return tmp2;
      }
    else 
      {
        double tmp1 = NewtonianG2(k1, k2,mu12) + G2GR(z,fnl, k1,k2, k3,mu12) ; //+ G2vv(z,fnl, k1,  k2,  k3, mu12);
        double tmp3 =  3.0*fg(z)*pow(HH(z),2)*tmp1/(k3*k3);  
        double tmp = bias10(z)*NewtonianF2(k1,k2,mu12) + b_twoT(z, fnl,k1,k2)+ TidalBias(z)*S_kernel(k1, k2,mu12)+ biasvsqT(z, fnl)*CurvyN2(z,k1,k2,mu12) + tmp3;//
        return tmp;
      }

}


double mu23(double r, double mu12)
{
     ///double r = k2/k1;
     double y = sqrt(r*r+2.0*mu12 *r + 1.0);
     double tmp = -(r + mu12)/(y);
     return tmp;

}

double mu13(double r, double mu12)
{
     ///double r = k2/k1;
     double y = sqrt(r*r+2.0*mu12 *r + 1.0);
     double tmp =- (r* mu12 +1.0)/(y);
     return tmp;

}


double  DynamicalGR::MatterBispectrum(int x,int y,double z, double fnl, double k1, double r, double mu12)
{

  //for y = 2,  > 6

   double Zz = sqrt(r*r + 2.0*mu12 *r + 1.0);
   double k3 = k1*Zz;
   double k2  = r* k1;
  
  double tmpA1 = SecondOrderMatterKernel(x, z,fnl, k1, k2,k3,mu12);
  double tmpA2 = LinearOrderMatterKernel(y, z, fnl,k1);
  double tmpA3 = LinearOrderMatterKernel(y, z, fnl,k2);
  double tmpA =  tmpA1*tmpA2*tmpA3*GalaxyP1P2(4,z,fnl,k1,k2);

  double tmpB1 = SecondOrderMatterKernel(x, z,fnl, k1, k3,k2,mu13(r, mu12));
  double tmpB2 = LinearOrderMatterKernel(y, z, fnl,k1);
  double tmpB3 = LinearOrderMatterKernel(y, z, fnl,k3);
  double tmpB = tmpB1*tmpB2*tmpB3*GalaxyP1P2(4,z,fnl,k1,k3);

  double tmpC1 = SecondOrderMatterKernel(x, z,fnl, k2, k3,k1,mu23(r, mu12));
  double tmpC2 = LinearOrderMatterKernel(y, z, fnl,k2);
  double tmpC3 = LinearOrderMatterKernel(y, z, fnl,k3);
  double tmpC =  tmpC1*tmpC2*tmpC3*GalaxyP1P2(4,z,fnl,k2,k3);
   return (tmpA + tmpB + tmpC);
}


double  DynamicalGR::RealSpaceBispectrumTotalMatter(int x,double z, double fnl, double k1, double r, double mu12)
{

   double y = sqrt(r*r + 2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;
  
  double tmpA1 = SecondOrderKernelRealSpaceTotalMatter(x, z,fnl, k1, k2,k3,mu12);
  double tmpA2 = b_one(z,fnl, k1);
  double tmpA3 = b_one(z,fnl, k2);
  double tmpA =  tmpA1*tmpA2*tmpA3*GalaxyP1P2(4,z,fnl,k1,k2);

  double tmpB1 = SecondOrderKernelRealSpaceTotalMatter(x, z,fnl, k1, k3,k2,mu13(r, mu12));
  double tmpB2 = b_one(z,fnl, k1);
  double tmpB3 = b_one(z,fnl, k3);
  double tmpB = tmpB1*tmpB2*tmpB3*GalaxyP1P2(4,z,fnl,k1,k3);

  double tmpC1 = SecondOrderKernelRealSpaceTotalMatter(x, z,fnl, k2, k3,k1,mu23(r, mu12));
  double tmpC2 = b_one(z,fnl, k2);
  double tmpC3 = b_one(z,fnl, k3);
  double tmpC =  tmpC1*tmpC2*tmpC3*GalaxyP1P2(4,z,fnl,k2,k3);
   return (tmpA + tmpB + tmpC);
}


double  DynamicalGR::RealSpaceBispectrumTotalMatterk1k2k3(int x,double z, double fnl, double k1, double k2, double k3)
  {
     double y = k3/k1;
     double r = k2/k1;
     double mu12 = (y*y -r*r - 1.0)/(2.0*r);
    //double y = sqrt(r*r + 2.0*mu12 *r + 1.0);
   double tmp = RealSpaceBispectrumTotalMatter(x,z, fnl, k1, r,mu12);
    
return tmp;
}

double DynamicalGR::ReducedRealSpaceBispectrumTotalMatter(int x,double z, double fnl,double k1, double r, double mu12)
{
 // double r = k2/k1;
 // double y = sqrt(r*r+2.0*mu12 *r + 1.0);
 // double k2 = r*k1; ///This ratio is the oppositte of what is in BCGS
  // double k3 =  k1*y;

  double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;

double PowerSPec_perm12 = GalaxyP1P2(1,z,fnl,k1,k2);
double PowerSPec_perm32 = GalaxyP1P2(1,z,fnl,k2,k3);
double PowerSPec_perm13 = GalaxyP1P2(1,z,fnl,k3,k1);

double tmp = RealSpaceBispectrumTotalMatter(x, z,fnl,k1, r, mu12)/(PowerSPec_perm12 + PowerSPec_perm32 + PowerSPec_perm13);//See equation 155
return tmp;
}




double  DynamicalGR::RealSpaceBispectrumPoisson(int x,double z, double fnl, double k1, double r, double mu12)
{

     double y = sqrt(r*r + 2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;
  
  double tmpA1 = SecondOrderKernelRealSpacePoisson(x, z,fnl, k1, k2,k3,mu12);
  double tmpA2 = b_one(z,fnl, k1)+  3.0*fg(z)*pow(HH(z),2)/(k1*k1);
  double tmpA3 = b_one(z,fnl, k2)+  3.0*fg(z)*pow(HH(z),2)/(k2*k2);
  double tmpA =  tmpA1*tmpA2*tmpA3*GalaxyP1P2(4,z,fnl,k1,k2);

  double tmpB1 = SecondOrderKernelRealSpacePoisson(x, z,fnl, k1, k3,k2,mu13(r, mu12));
  double tmpB2 = b_one(z,fnl, k1)+  3.0*fg(z)*pow(HH(z),2)/(k1*k1);
  double tmpB3 = b_one(z,fnl, k3)+  3.0*fg(z)*pow(HH(z),2)/(k3*k3);
  double tmpB = tmpB1*tmpB2*tmpB3*GalaxyP1P2(4,z,fnl,k1,k3);

  double tmpC1 = SecondOrderKernelRealSpacePoisson(x, z,fnl, k2, k3,k1,mu23(r, mu12));
  double tmpC2 = b_one(z,fnl, k2)+  3.0*fg(z)*pow(HH(z),2)/(k2*k2);
  double tmpC3 = b_one(z,fnl, k3)+  3.0*fg(z)*pow(HH(z),2)/(k3*k3);
  double tmpC =  tmpC1*tmpC2*tmpC3*GalaxyP1P2(4,z,fnl,k2,k3);


   return (tmpA + tmpB + tmpC);
}



double  DynamicalGR::RealSpaceKernelTram(int y,double z, double fnl, double k1, double k2, double k3)
{

  double mu12 = (k3*k3 -(k1*k1 + k2*k2))/(2.0*k1*k2);
  if(y ==1)
  {
  double tmpA1 = SecondOrderMatterKernel(10, z,fnl, k1, k2,k3,mu12);
  double tmpA2 = 1.0 +  3.0*fg(z)*pow(HH(z),2)/(k1*k1);
  double tmpA3 = 1.0 +  3.0*fg(z)*pow(HH(z),2)/(k2*k2);
  return tmpA1/(2.*tmpA3*tmpA2);
   }
  else if(y ==2)
  {
  double tmpB1 = SecondOrderMatterKernel(4, z,fnl, k1, k2,k3,mu12);
  return tmpB1/2.0;
   }
  else if(y ==3)
  {
  double tmpC1 = SecondOrderMatterKernel(1, z,fnl, k1, k2,k3,mu12);
  return tmpC1/2.0;
  }
  else return 0.0;

}





double DynamicalGR::ReducedRealSpaceBispectrumPoisson(int x,double z, double fnl, double k1, double r, double mu12)
{
 // double r = k2/k1;
 // double y = sqrt(r*r+2.0*mu12 *r + 1.0);
 // double k2 = r*k1; ///This ratio is the oppositte of what is in BCGS
  // double k3 =  k1*y;

  double y = sqrt(r*r+2.0*mu12 *r + 1.0);
   double k3 = k1*y;
   double k2  = r* k1;
//GalaxyP1P2(int x, double z, double  k1, double k2)
double PowerSPec_perm12 = GalaxyP1P2(2,z,fnl,k1,k2);
double PowerSPec_perm32 = GalaxyP1P2(2,z,fnl,k2,k3);
double PowerSPec_perm13 = GalaxyP1P2(2,z,fnl,k3,k1);
double tmp = RealSpaceBispectrumPoisson(x, z, fnl,k1,  r, mu12)/(PowerSPec_perm12+PowerSPec_perm32+PowerSPec_perm13);//See equation 155
return tmp;
}





double DynamicalGR::DiffReducedRealSpaceBispectrum(int y, double z, double fnl,double k1, double r, double mu12)
{
  if(y ==1)
  {
 double tmp1 = ReducedRealSpaceBispectrumTotalMatter(8, z,fnl, k1,  r, mu12) - ReducedRealSpaceBispectrumTotalMatter(5, z,fnl, k1,  r, mu12);
 double tmp = tmp1/ReducedRealSpaceBispectrumTotalMatter(5, z,fnl, k1,  r, mu12);
 return tmp;
  }
  else if(y ==2)
  {
 double tmp2 = ReducedRealSpaceBispectrumPoisson(8, z, fnl, k1,  r, mu12) - ReducedRealSpaceBispectrumTotalMatter(5, z,fnl, k1,  r, mu12);
 double tmp = tmp2/ReducedRealSpaceBispectrumTotalMatter(5, z,fnl, k1,  r, mu12);
 return tmp;
  }
 else
 {
  return 0.0;
 }
}

double DynamicalGR::DiffRealSpaceBispectrum(int y,double z, double fnl,double k1, double r, double mu12)
{
  if(y == 1)
  {
 double tmp1 = RealSpaceBispectrumTotalMatter(8, z,fnl, k1,  r, mu12) - RealSpaceBispectrumTotalMatter(5, z, fnl, k1,  r, mu12);
 double tmp = tmp1/RealSpaceBispectrumTotalMatter(5, z,fnl, k1,  r, mu12);
 return tmp;
  }
   else if(y == 2)
  {
 double tmp1 = RealSpaceBispectrumPoisson(8, z,fnl, k1,  r, mu12) - RealSpaceBispectrumPoisson(5, z,fnl, k1,  r, mu12);
 double tmp = tmp1/RealSpaceBispectrumPoisson(5, z,fnl, k1,  r, mu12);
 return tmp;
  }
  else
  {
  return 0.0;
  }

}









