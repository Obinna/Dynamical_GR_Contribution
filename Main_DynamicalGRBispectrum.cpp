
#include<iostream>
#include <fstream>

#include <time.h>
#include <sys/time.h>

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cmath>
#include <list>



#include <omp.h>
#include <iomanip>
#include <pthread.h>
#include <list>




#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/rk_int.cpp"
#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/routines.cpp"
#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/constants.cpp"
//#include <array>


#include"/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/MyCosmology.cpp"

#include"HoDBiasModel.cpp"
#include"DynamicalGRBispectrum.cpp"



using namespace std;

double returnk2(double k1, double r)
{
  return r*k1;
}
double returnk3(double k1, double r, double mu12)
{
  double k2 =returnk2(k1,r);
  double tmp = sqrt(k1*k1 + k2*k2 + 2.0*k1*k2*mu12);
  return tmp;
}
double returnmu12(double k3, double k1, double r)
{
  double k2 = k1* r;
   double tmp = (k3*k3 - k2* k2 - k1*k1)/(2.0*k1* k2);
  return tmp;

}

double kH(double z)
{
  Cosmology C;
    double c_light = C.H_0*3000.0;
    const  double c_light_km_s =  299792;
    // double c_light =  C.H_0*3000.0;
    return C.H(z)/(c_light*(1.0+z));
}


int main(int argc, char **argv)
{


//const char* output = "PowerSpectrumGRdata1.dat";//Power spectrum

//const char* output = "DynamicalGRdata1.dat";//Poisson Bispectrum

//const char* output = "DynamicalGRdata2.dat";//Bispectrum eq Poisson

//const char* output = "DynamicalGRdata3.dat";//Poisson reduced Poisson 0.01

//const char* output = "DynamicalGRdata4.dat";//Poisson reduced Poisson 0.005

  //const char* output = "DynamicalGRdata5.dat";

  const char* output = "DynamicalGRdata6.dat";

 // const char* output = "DynamicalGRdata7.dat";

  //const char* output = "DynamicalGRdata8.dat";///compare Bernardau

 // const char* output = "Biasdata.dat";

  //const char* output = "nofz.txt";

ofstream myfile;
myfile.setf(ios::scientific);
myfile.open(output);


  Cosmology C;
C.InitializeCosmology();
C.InitializeDz();
C.InitializeSigma_EH();
C.TabulateSigma_EH();
//ComputeBessels();

HoDBiasModel HoD;
//HoD.InitialzieBiasloop();

DynamicalGR DGR;
DGR.ReadWrite();







  double M;
  double z;
  int nM = 1000;
  double za[5] = {10.07,5.72,3.06,1.50,0.0};
  double pa[5];
int N = 500;
vector <double > z_array(N);
vector <double > k_array(N);
vector <double > M_array(N);
vector <double > R_array(N);
vector <double > mu1_array(N);
vector <double > mu2_array(N);
vector <double > fnl_array(N);
vector <long double> Out1(N);
vector <long double> Out2(N);
vector <long double> Out3(N);
vector <long double> Out4(N);
vector <long double> Out5(N);
vector <long double> Out6(N);
vector <long double> Out7(N);
vector <long double> Out8(N);
vector <long double> Out9(N);
vector <long double> Out10(N);
vector <double> theta_array(N);
vector <double> theta_array2(N);
vector <int> ell(N);

std::vector<std::vector< double> > BVtabij(N,std::vector<double>(N));

  const double z1 = 0.5;
  const double z2 = 1.0;



 

  const double rmin = 0.1;
  const double rmax = 10.0;

  const double kvalues1 = 0.05;
  const double kvalues2 = 0.0503;
   const double kvalues3 = 0.005;
    const double kvalues4 = 0.001;

  const double r1 = kvalues2/kvalues1;
  const double r2 = 0.5;

double mumin = -0.9991;
  double mumax = 0.999;
  const double z_min = 0.0002;;
  const double z_max = 2.0;

double thetamin = 0.0001;
  double thetamax = M_PI-0.0001;

const double thetam = 120.0*M_PI/180.0;
const double thetam2 = 177.0*M_PI/180.0;

  const double k_min = abs(kvalues2-kvalues1);
  const double k_max = kvalues2+kvalues1;

  const double fnl_fid = 0.0;



int thread_number81 = 8;
/*

  for(int i=0;i< N; i++)
  {
    k_array[i] = exp(((log(k_max)-log(k_min))*((double) i)/((double) N-1) + log(k_min) ));
    R_array[i] = (rmin+i*(rmax-rmin)/N);
    theta_array[i] = (thetamin+i*(thetamax-thetamin)/N);
  }
       */ 
        



//#pragma omp parallel private(thread_number81)
{
////#pragma omp for schedule(static) nowait
for(int i = 1; i < N;i++)
{
  
k_array[i] = exp( (  (log(k_max)-log(k_min))*((double) i)/((double) N-1) + log(k_min) ));
z_array[i] = z_min + i*(z_max-z_min)/N;//exp( (  (log(z_max)-log(z_min))*((double) i)/((double) N-1) + log(z_min) ));
//theta_array[i] = cos(thetamin+i*(thetamax-thetamin)/N);

theta_array[i] = cos(returnmu12(k_array[i], kvalues1, 1));
theta_array2[i] = cos(returnmu12(k_array[i], kvalues2, 1));

//RealSpaceBispectrum(int x, double z,double k1, double k2, double mu12)
/*
cout<<"\t" << "k values = " << k_array[i]<< " Mpc/h" << endl;
cout <<"\t"<< "mu12 values = " <<theta_array[i]<< " [-1,1]" << endl;
cout <<"\t"<< "acos(x) = " << acos(theta_array[i]) << " radians" << endl;
 cout <<"\t"<< "acos(x) = " << acos(theta_array[i])*180/3.1415 << " degrees" << endl;
*/
/*
  Out1[i] = HoD.Galb1(z_array[i]);
  Out2[i] = HoD.Galb2(z_array[i]);
  Out3[i] = HoD.Galbsq(z_array[i]);
  Out4[i] = HoD.ng(z_array[i]);

  myfile<<z_array[i]<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<" "<<Out4[i]<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"z="<<z_array[i]<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<endl;
*/

/*
Out1[i] = DGR.RealSpaceBispectrumTotalMatter(5,z1,k_array[i],1, cos(thetam2));
Out2[i] = DGR.RealSpaceBispectrumTotalMatter(8,z1,k_array[i],1, cos(thetam2));
Out3[i] = DGR.RealSpaceBispectrumPoisson(8,z1,k_array[i],1, cos(thetam2));
*/
/*
Out1[i] = DGR.RealSpaceBispectrumTotalMatter(5,z1,kvalues1,1, theta_array[i]);
Out2[i] = DGR.RealSpaceBispectrumTotalMatter(8,z1,kvalues1,1, theta_array[i]);
Out3[i] = DGR.RealSpaceBispectrumPoisson(8,z1,kvalues1,1, theta_array[i]);
*/


Out1[i] = DGR.RealSpaceBispectrumTotalMatter(1,z1,fnl_fid,kvalues1,r1,returnmu12(k_array[i], kvalues1, r1));
Out2[i] = DGR.RealSpaceBispectrumTotalMatter(2,z1,fnl_fid,kvalues1,r1, returnmu12(k_array[i], kvalues1, r1));
Out3[i] = DGR.RealSpaceBispectrumTotalMatter(3,z1,fnl_fid,kvalues1,r1, returnmu12(k_array[i], kvalues1, r1));
Out4[i] = DGR.RealSpaceBispectrumTotalMatter(4,z1,fnl_fid,kvalues1,r1, returnmu12(k_array[i], kvalues1, r1));
Out5[i] = DGR.RealSpaceBispectrumTotalMatter(5,z1,fnl_fid,kvalues1,r1, returnmu12(k_array[i], kvalues1, r1));
Out6[i] = DGR.RealSpaceBispectrumTotalMatter(6,z1,fnl_fid,kvalues1,r1, returnmu12(k_array[i], kvalues1, r1));
Out7[i] = DGR.RealSpaceBispectrumTotalMatter(7,z1,fnl_fid,kvalues1,r1, returnmu12(k_array[i], kvalues1, r1));

Out8[i] = DGR.RealSpaceBispectrumTotalMatter(9,z1,fnl_fid,kvalues1,r1, returnmu12(k_array[i], kvalues1, r1));




myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<" "<<abs(Out7[i])<<" "<<abs(Out8[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<", "<<"Out5="<<Out5[i]<<", "<<"Out6="<<Out6[i]<<endl;


/*
 Out1[i] = DGR.RealSpaceBispectrumTotalMatter(1,z1,k_array[i],r1, cos(thetam));
 Out2[i] = DGR.RealSpaceBispectrumTotalMatter(2,z1,k_array[i],r1, cos(thetam));
 Out3[i] = DGR.RealSpaceBispectrumTotalMatter(3,z1,k_array[i],r1, cos(thetam));
 Out4[i] = DGR.RealSpaceBispectrumTotalMatter(4,z1,k_array[i],r1, cos(thetam));
 Out5[i] = DGR.RealSpaceBispectrumTotalMatter(5,z1,k_array[i],r1, cos(thetam));
 Out6[i] = DGR.RealSpaceBispectrumTotalMatter(8,z1,k_array[i],r1, cos(thetam));

myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<", "<<"Out5="<<Out5[i]<<", "<<"Out6="<<Out6[i]<<endl;
*/

 /*

Out1[i] = DGR.RealSpaceBispectrumPoisson(1,z1,k_array[i],r1, cos(thetam2));
Out2[i] = DGR.RealSpaceBispectrumPoisson(2,z1,k_array[i],r1, cos(thetam2));
Out3[i] = DGR.RealSpaceBispectrumPoisson(3,z1,k_array[i],r1, cos(thetam2));
Out4[i] = DGR.RealSpaceBispectrumPoisson(4,z1,k_array[i],r1, cos(thetam2));
Out5[i] = DGR.RealSpaceBispectrumPoisson(5,z1,k_array[i],r1, cos(thetam2));
Out6[i] = DGR.RealSpaceBispectrumPoisson(6,z1,k_array[i],r1, cos(thetam2));
Out7[i] = DGR.RealSpaceBispectrumPoisson(7,z1,k_array[i],r1, cos(thetam2));
Out8[i] = DGR.RealSpaceBispectrumPoisson(8,z1,k_array[i],r1, cos(thetam2));

//Out3[i] = DGR.RealSpaceBispectrumPoisson(1,z2,k_array[i],1, 0.5);
//Out4[i] = DGR.RealSpaceBispectrumPoisson(2,z2,k_array[i],1, 0.5);


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<" "<<abs(Out7[i])<<" "<<abs(Out8[i])<<"\n";

cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<", "<<"Out5="<<Out5[i]<<", "<<"Out6="<<Out6[i]<<endl;

*/
 /*
Out1[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(5,z1,kvalues1,r1, returnmu12(k_array[i], kvalues1, 1));
Out2[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(8,z1,kvalues1,r1, returnmu12(k_array[i], kvalues1, 1));
Out3[i] = DGR.ReducedRealSpaceBispectrumPoisson(8,z1,kvalues1,r1, returnmu12(k_array[i], kvalues1, 1));
//Out4[i] = k_array[i];

myfile<<acos(returnmu12(k_array[i], kvalues1, 1))/M_PI<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"theta/pi="<<"\t"<<acos(returnmu12(k_array[i], kvalues1, 1))/M_PI<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<", "<<"Out5="<<Out5[i]<<", "<<"Out6="<<Out6[i]<<endl;
*/

/*
Out1[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(5,z1,kvalues3,r1, returnmu12(k_array[i], kvalues3, r1));
Out2[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(8,z1,kvalues3,r1, returnmu12(k_array[i], kvalues3, r1));
Out3[i] = DGR.ReducedRealSpaceBispectrumPoisson(8,z1,kvalues3,r1, returnmu12(k_array[i], kvalues3, r1));
//Out4[i] = k_array[i];

myfile<<acos(returnmu12(k_array[i], kvalues1, r1))/M_PI<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"theta/pi="<<"\t"<<acos(returnmu12(k_array[i], kvalues1, 1))/M_PI<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<", "<<"Out5="<<Out5[i]<<", "<<"Out6="<<Out6[i]<<endl;
*/

/*
Out1[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(5,z1,kvalues1,r1, theta_array[i]);
Out2[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(8,z1,kvalues1,r1, theta_array[i]);
Out3[i] = DGR.ReducedRealSpaceBispectrumPoisson(8,z1,kvalues1,r1, theta_array[i]);

//Out4[i] = DGR.ReducedRealSpaceBispectrumPoisson(4,z1,kvalues1,r1, theta_array[i]);

myfile<<acos(theta_array[i])/M_PI<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<"\n";

cout<<"\t" <<"i="<<i<<" "<<"k="<<acos(theta_array[i])/M_PI<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<endl;
*/

/*
Out1[i] = DGR.DiffReducedRealSpaceBispectrum(1,z1,kvalues1,r1, theta_array[i]);
Out2[i] = DGR.DiffReducedRealSpaceBispectrum(2,z1,kvalues1,r1, theta_array[i]);
Out3[i] = DGR.DiffRealSpaceBispectrum(1,z1,kvalues1,r1, theta_array[i]);
Out4[i] = DGR.DiffRealSpaceBispectrum(2,z1,kvalues1,r1, theta_array[i]);

myfile<<acos(theta_array[i])/M_PI<<" "<<Out1[i]<<" "<<Out2[i]<<" "<<Out3[i]<<" "<<Out4[i]<<"\n";

cout<<"\t" <<"i="<<i<<" "<<"k="<<acos(theta_array[i])/M_PI<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<endl;
*/

/*

//ReducedRealSpaceBispectrumTotalMatter(int x,double z,double k1, double r, double mu12)
Out1[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(1,z1,kvalues1,r1, theta_array[i]);

Out2[i] = DGR.ReducedRealSpaceBispectrumPoisson(6,z1,kvalues1,r1, theta_array[i]);
Out3[i] = DGR.ReducedRealSpaceBispectrumPoisson(8,z1,kvalues1,r1, theta_array[i]);

Out4[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(5,z1,kvalues2,r1, theta_array[i]);
Out5[i] = DGR.ReducedRealSpaceBispectrumPoisson(6,z1,kvalues2,r1, theta_array[i]);
Out6[i] = DGR.ReducedRealSpaceBispectrumPoisson(8,z1,kvalues2,r1, theta_array[i]);

Out7[i] = DGR.ReducedRealSpaceBispectrumTotalMatter(5,z1,kvalues3,r1, theta_array[i]);
Out8[i] = DGR.ReducedRealSpaceBispectrumPoisson(6,z1,kvalues3,r1, theta_array[i]);
Out9[i] = DGR.ReducedRealSpaceBispectrumPoisson(8,z1,kvalues3,r1, theta_array[i]);

myfile<<acos(theta_array[i])/M_PI<<" "<<Out1[i]<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<" "<<abs(Out5[i])<<" "<<abs(Out6[i])<<" "<<abs(Out7[i])<<" "<<abs(Out8[i])<<" "<<abs(Out9[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"theta="<<acos(theta_array[i])/M_PI<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<", "<<"Out5="<<Out5[i]<<", "<<"Out6="<<Out6[i]<<", "<<"Out7="<<Out7[i]<<endl;
*/

 //LinearPowerSpectrum(int x, double z, double k)

/*
 Out1[i] = DGR.LinearPowerSpectrum(1,z2, k_array[i]);
 Out2[i] = DGR.LinearPowerSpectrum(2,z2,k_array[i]);
 Out3[i] = DGR.LinearPowerSpectrum(3,z2,k_array[i]);


myfile<<k_array[i] <<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<endl;
*/

//Out1[i] = DGR.DiffReducedRealSpaceBispectrumTotalMatter(z2,0.01,1, cos(theta_array[i]));
//Out2[i] = DGR.DiffReducedRealSpaceBispectrumTotalMatter(z2,0.01,2,cos(theta_array[i]));
//myfile<<theta_array[i]/M_PI<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<"\n";


/*
Out1[i] = DGR.DiffReducedRealSpaceBispectrum(1,z2,kvalues1,1, returnmu12(k_array[i], kvalues1,1 ));
Out2[i] = DGR.DiffReducedRealSpaceBispectrum(2,z2,kvalues1,2,returnmu12(k_array[i],kvalues1,2 ));
Out3[i] = DGR.DiffReducedRealSpaceBispectrum(1,z2,kvalues2,1,returnmu12(k_array[i],kvalues2,1 ));
Out4[i] = DGR.DiffReducedRealSpaceBispectrum(2,z2,kvalues2,2,returnmu12(k_array[i],kvalues2,2 ));


myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<endl;
*/
  /*
//DiffRealSpaceBispectrum(int y,double z,double k1, double r, double mu12)returnmu12(double k3, double k1, double r)
Out1[i] = DGR.DiffRealSpaceBispectrum(1,z2,kvalues1,1,returnmu12(k_array[i],kvalues1,1 ));
Out2[i] = DGR.DiffRealSpaceBispectrum(2, z2,kvalues1,2,returnmu12(k_array[i],kvalues1,2 ));
Out3[i] = DGR.DiffRealSpaceBispectrum(1,z2,kvalues2,1,returnmu12(k_array[i],kvalues2,1 ));
Out4[i] = DGR.DiffRealSpaceBispectrum(2,z2,kvalues2,2,returnmu12(k_array[i],kvalues2,2 ));

myfile<<k_array[i]<<" "<<abs(Out1[i])<<" "<<abs(Out2[i])<<" "<<abs(Out3[i])<<" "<<abs(Out4[i])<<"\n";
cout<<"\t" <<"i="<<i<<" "<<"k="<<k_array[i]<<", "<<"Out1="<<Out1[i]<<", "<<"Out2="<<Out2[i]<<", "<<"Out3="<<Out3[i]<<", "<<"Out4="<<Out4[i]<<endl;

*/

/*
      //#pragma omp parallel for schedule(dynamic,1) collapse(2)
for(int i = 0;i< N;i++)
{
      //karray[i] = exp(((log(kmax)-log(kmin))*((double) i)/((double) Narray-1) + log(kmin) ));
      //phiarray[i] = phimin+ i*(phimax-phimin)/Narray;

   for(int j = 0; j < N; j++)
    {
//BVtab[i][j] = Integrate(bind(&Spectroscopic::GalBispectrumV,this,x,fnl,z,karray[i],karray2[i],mu12,phiarray[j],0.0,_1),-0.998,0.998,EPSREL,EPSABS)/2.0;
//BVtabij[i][j] = integrate1(bind(&Spectroscopic::GalBispectrumV,this,x,fnl,z,karray[i],karray2[i],mu12,phiarray[j],0.0,_1),-0.998,0.998)/2.0;
//BVtabij[i][j] = GalBispectrumV(x,fnl,z,karray[i],x2*karray[i],mu12,phiarray[j],0.0,1);
 //double Spectroscopic::PassbyrefGalV2(int x, void *params, double k1, double phi, double mu1);
      //returnmu12(double k3, double k1, double r)
   Out1[i] = R_array[i];
   Out2[j] = k_array[j]/kvalues1;
   BVtabij[i][j] = DGR.DiffRealSpaceBispectrumTotalMatter(z2,kvalues1,R_array[i],returnmu12(k_array[j],kvalues1,R_array[i] ));

   myfile<<Out1[i]<<" "<<Out2[j]<<" "<<abs(BVtabij[i][j])<<"\n";

cout<<"\t"<<"k1"<<"\t"<<k_array[i]<<"\t"<<"phi="<<"\t"<<R_array[j]<<"\t"<<"bBV"<<"\t"<<"\t"<<BVtabij[i][j]<<endl;
  


   }
}
*/

}
}



myfile.close();






  return 0;
}



