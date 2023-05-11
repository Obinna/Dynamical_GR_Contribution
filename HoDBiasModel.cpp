#include"HoDBiasModel.h"

//#include"MyCosmology.h"

#include <boost/bind.hpp>
 //#include "/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/Quadrature.h"
//#include "/Users/o-admin/Dropbox/UWC_WorkShop/Effective_fnl/Effectivefnl2/Quadrature.cpp"

using namespace std;


const double Mmin = 1.5e13;
const double Mmax = 1.4e13;
HoDBiasModel::HoDBiasModel(void)
{
    //SetPowerSpectrum();
}

HoDBiasModel::~HoDBiasModel()
{
   // SetParameters();
}
/*
const int N = 587; /// 612 You have to make sure this number corresponds to the number of rows in class
std::vector <double> k_array(N); // define vectors to k-values
std::vector <double> pk_array(N);  // define vectors to hold pofk
void ImprovedHIModel::ReadWrite1(void)
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
*/
// Matter Power Spectrum Spline
/*
double MySplinePofk(double k)
{
	Cosmology C;

  //Spline<double, double> CubicSpline_num(k_array,pk_array);
  //double res = CubicSpline_num.interpolate(k);
  return C.P_k_EH_Tabulated(k);
}*/

double WFun(double R, double k)
{
   double x = k*R;
   if(x<1.0e-5)
    {
       return 1.0;
    }
   else if(x>500000.0)
    {
		return 0.0;
    }
   else 
   {
       double top_hat = 3.0*(sin(x)/pow(x,3)-cos(x)/pow(x,2));
       return  top_hat;
    }
}

double DWFunDM(double R, double k)
{
	Cosmology C;
    double x = k*R;
    if(x<1.0e-5)
    {
		return 1.0;
    }
   else if(x>500000.0)
    {
		return 0.0;
    }
   else 
   {
   	double M = 4*M_PI*C.rho_m(0.0)*pow(R,3)/3.0;
   double topbra = (k*k*R -3.0/R);
   double dRdM = R/(3.0*M);
   double top_hatderivtivewrtR = (3.0/pow(x,3))*(topbra*sin(x)+3.0*k*cos(x));

   double top_hatderivtivewrtM  = top_hatderivtivewrtR * dRdM;

   return  top_hatderivtivewrtM;
    }
}



double DSigmaDMintgrd(double R, double k)
{
	 Cosmology C;
   double tmp = k*k* C.P_k_EH_Tabulated(k)*WFun(R,k)*DWFunDM(R,k)/(pow(M_PI,2));
     return tmp;
}


double Sigmaintgrd(double R, double k)
{
    Cosmology C;
    double tmp = k*k*C.P_k_EH_Tabulated(k)* WFun(R,k)* WFun(R,k)/(2*pow(M_PI,2));
   return  tmp;

}

/*
double Mysigma1(double R)
{
      const double EPSREL = 1e-6;
      const double EPSAB = 1e-6;   
      const double QMIN = 1e-4;
      const double QMAX = 1e4;

      double tmp1 = Integrate(bind(Sigmaintgrd,R,_1),QMIN,QMAX,EPSREL,EPSAB);
      //cout<<"\t i"<<"\t"<<"R"<<"\t"<<R <<"\t"<<"Sig_tab"<<"\t"<<"\t"<<tmp1 <<endl;
        return tmp1;
}

double DsigmasqDM(double R)
{
      const double EPSREL = 1e-4;
      const double EPSAB = 1e-4;   
      const double QMIN = 1e-4;
      const double QMAX = 1e4;

      double tmp1 = Integrate(bind(DSigmaDMintgrd,R,_1),QMIN,QMAX,EPSREL,EPSAB);
     // cout<<"\t i"<<"\t"<<"R"<<"\t"<<R <<"\t"<<"Sig_tab"<<"\t"<<"\t"<<tmp1 <<endl;
        return tmp1;
}
*/
///////////////////////Tabulate//////////////////////////////////

/*

 int k_bins2 =100;
 std::vector <double> Rarrayin(k_bins2);
 std::vector <double> Sig_tab1(k_bins2);
 std::vector <double> Sigv_tab1(k_bins2);



void HoDBiasModel::InitialzieSigmaloop(void)
{
	double  EPSREL = 1e-8;
	double Rmin2;
	double Rmax2;
	Rmin2 = log(1e-5);
	Rmax2 = log(100.0);
	const double QMIN = 1e-4;
	const double QMAX = 1e5;	
#pragma omp parallel private(thread_number)
	{
#pragma omp for schedule(static) nowait
		for(int i = 0;i< k_bins2 ;i++)
		{
			Rarrayin[i] = exp(((Rmax2-Rmin2)*((double) i)/((double)  k_bins2-1) + Rmin2 ));
            Sig_tab1[i] =  Mysigma1(Rarrayin[i]);
            Sigv_tab1[i] =  DsigmasqDM(Rarrayin[i]);

            cout<<"\t i"<<"\t"<<"R"<<"\t"<<Rarrayin[i]<<"\t"<<"Sig_tab"<<"\t"<<"\t"<< Sig_tab1[i]<<"\t"<<"DSig_tabDM"<<"\t"<< Sigv_tab1[i]<<endl;
			
			
		}
		
	}
}



double Sigmasq(double z, double R)
{

       Cosmology C;
       double Dz = C.D_growth_Tabulated(z);
        double coef = pow(Dz,2);
	Spline<double, double> CubicSpline_sigma(Rarrayin,Sig_tab1);
	return coef*CubicSpline_sigma.interpolate(R);
   	
}
double DSigmasqDM(double z, double R)
{       Cosmology C;
        double Dz = C.D_growth_Tabulated(z);
        double coef = pow(Dz,2);
	Spline<double, double> CubicSpline_sigmav(Rarrayin,Sigv_tab1);
	return coef*CubicSpline_sigmav.interpolate(R);
   	
}
*/




//////////////////////////////////////////Basic relationship////////////////////////////////////////

double Omegam_inBias(double z)//Omegam
{
    Cosmology C;
    double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
    return tmpz;
}
double delta_cz(double z)
{  Cosmology C;
	//double Dz = C.D_growth_Tabulated(z);
    //double tmp = 3.0*pow((12.0*M_PI),2.0/3.0)*(1.0  + 0.0123* log10(Omegam_inBias(z)));
     double tmp = 1.686*C.D_growth_Tabulated(0.01)/C.D_growth_Tabulated(z);
   return tmp;
    //return 1.686;
}


double RofM(double z, double M)
{
   Cosmology C;
   double tmp = pow((3.0*M )/(4.0*M_PI*C.rho_m(z)),1/3);
   return tmp;
}


////////////////////////////////////////////Basic relationship ends//////////////////////////////////////////////////////////

double nu_Mz(double M, double z)
{

   
  Cosmology C;
  double delta_c = delta_cz(z);// 1.686;
  double Dz      = C.D_growth_Tabulated(z);
  double sig     = C.sigma_EH(M);
  double nu = pow(delta_c/Dz*sig,2);
  //double tmp =  delta_c/(Dz*sig);
  return  nu;
  
}
double dsigmadM(double M)
{
        Cosmology C;
	//dumb derivative       
	double sigma_A = C.sigma_EH(1.01*M);
	double sigma_B = C.sigma_EH(0.99*M);
	return (sigma_A - sigma_B)/(0.02*M);
}


/*! \fn double f_nu_ST(double nu, void *params)
 *  \brief First-Crossing Distribution for Sheth-Tormen*/
double nu_f_nu_ST(double nu)
{
  //equation 3.4 of Toby
  double p  = 0.3;
  double q   = 0.707;
  double A   = 0.322; 
  double nuq = q*nu;
  double tmp = A*(1.0 + pow(nuq,-p))*sqrt(nuq/(2.0*M_PI))*exp(-0.5*nuq);
// cout <<"nufST="<<tmp<<endl;
  return tmp;
}







////Shet Torman
double b_ST_10(double M,double  z)
{
 
   Cosmology C;
  double delta_c = delta_cz(z);// 1.686;
  double nu     = nu_Mz(M,z);
  // double delta_c = 1.686;
  double p  = 0.3;
  double q   = 0.707;
  //double A   = 0.322;
  double nuq = q*nu;
   double dndnu =-(nuq-1)/(2*nu) - p/(nu*(1+pow(nuq,p)));
   // double first_term = -(nuq-1)/delta_c;
   // double second_term = -(2.0*p)/(delta_c*(1+pow(nuq,p)));

  return -(2*nu/delta_c)* dndnu;
}



double b_ST_20(double M,double  z)
{
  //M is in h^-1 Msun
 // double *fp    = (double *) params;
 // double z      = fp[0];
    Cosmology C;
   
  double delta_c = delta_cz(z);// 1.686;
   double nu     = nu_Mz(M,z);
  // double delta_c = 1.686;
   double p  = 0.3;
   double q   = 0.707;
   //double A   = 0.322; 
  // double nuq = q*nu;
   /*
   double ddnddnv =(pow(p,2)+ nuq*p)/(nu*nu*(1.0+pow(nuq,q)))
                  +(pow(nuq,2)-2.0*nuq-1)/(4.0*nu*nu);

   double dndnu =-(nuq-1)/(2.0*nu) - p/(nu*(1.0+pow(nuq,p)));

   double tmp = (4.0*nu*nu*(ddnddnv))/(delta_c*delta_c)
                + (2.0*nu*(dndnu))/(delta_c*delta_c);
                */



     double tmpup1 = 4.0*(pow(p,2) + nu*p*q)- (q*nu-1.0)*(1.0 + pow(q*nu,p))-2.0*p;
     double tmpdown1 = (1.0 + pow(nu*q,p))*(delta_c*delta_c);
     double tmp1 = tmpup1/tmpdown1;
     double tmp2 =  (pow(q*nu,2) - 2.0*nu*q-1.0)/(delta_c*delta_c);
     double tmp = tmp1 + tmp2;

  return tmp;
}




double dlnnudlnM(double M,double z)
{
    Cosmology C;
     // double Dz = C.D_growth_Tabulated(z);
    double delta = 1.0e-4;
    //double M1 = M - delta ;
   // double M2 = M + delta ;
   // double coef = M/nu_Mz(M,z);
    double lnM1 = log(M) - delta ;
    double lnM2 = log(M) + delta ;

    double nu1 = log(nu_Mz(exp(lnM1),z));
    double nu2 = log(nu_Mz(exp(lnM2),z));

   // double nu2 = nu_Mz(1.01*M,z);
	//double  nu1 = nu_Mz(0.99*M,z);
	// return (nu2 - nu1 )/(0.02*M);

    //double nu1 = log(pow(delta_cz(z)/(Dz*C.sigma_EH(exp(lnM1))),2));
   // double nu2 = log(pow(delta_cz(z)/(Dz*C.sigma_EH(exp(lnM2))),2));

    double deri =  (nu2-nu1)/(lnM2-lnM1);
    //double tmp  = deri;
   return deri;
}
/*
double n_m_ST(double M, double z)
{
  //M is in h^-1 Msun
  Cosmology C;
  //double *fp    = (double *) params;
 // double z      = fp[0]; //redshift
 // double rho_m  = fp[1]; //mean background density in h^2 Msun / Mpc^3
  double nu   = nu_Mz(M,z);
  //double dlnnu_dlnM = M*dnudM(M,params)/nu;

  //printf("M %e z %e rho_m %e nu %e dnu_dM %e\n",M,z,rho_m,nu,dnu_dM);

  //number of halos per h^3 Mpc^-3 per dM
 double tmp =  (C.rho_m(z)/(M*M))*nu_f_nu_ST(nu)*dlnnudlnM(M,z);//;
//cout <<"n_mST="<<C.rho_m(z)<<endl;
   return tmp;
}
*/
double NumberCentralGal(double Mmin, double M)
{
	double sigmalogM = 0.45;
	double tmp1 = (log10(M) - log10(Mmin))/sigmalogM;
	double tmp = 0.5*(1.0 + erf(tmp1));
	return tmp;
}
double NumberSatellite(double M)
{
	double Mcut =  1.4e13;
	double MI = 1.3e14;
	double alpha = 1.38;
	double tmp = pow((M-Mcut)/MI,alpha);
	return tmp;
}
 double MeanNumberGal(double Mmin, double M)
 {
 double tmp = NumberCentralGal(Mmin, M)*(1.0 + NumberSatellite(M));
 return tmp;
 }



double n_m_ST(double M, double z)
{
  //M is in h^-1 Msun
  Cosmology C;
  double nu     = nu_Mz(M,z);
  //number of halos per h^3 Mpc^-3 per dM
  double tmp =  (C.rho_m(z)/(M*M))*nu_f_nu_ST(nu)*dlnnudlnM(M,z);//;
//cout <<"n_mST="<<C.rho_m(z)<<endl;
   return tmp;
}

 double Num_b_ST_intgrnd_10(double z, double Mmin, double M)
   {
     double dnudm =   b_ST_10(M,z)*MeanNumberGal(Mmin, M)*n_m_ST(M, z); 
      return dnudm;
   }

   double Num_b_ST_intgrnd_20(double z,double Mmin, double M)
   {
      double dnudm = b_ST_20(M,z)*MeanNumberGal(Mmin, M)* n_m_ST(M, z);
     return dnudm;
    }

   double ngintrd(double z,double Mmin, double M)
   {   
     double MeannST= MeanNumberGal(Mmin, M)* n_m_ST(M, z);
    //cout <<"\tMean="<<"\t"<<MeannST<<endl;
    return MeannST;
    }

double ngint(double z)
{
	const double EPSREL = 1e-10;
    const double EPSABS = 1e-10;
	double tmp  = Integrate(bind(ngintrd,z,Mmin,_1), Mmin,Mmax,EPSREL,EPSABS);
	//cout <<"\ttmp="<<"\t"<<tmp<<endl;
	return tmp;
}

   double b_ST_10_in(double z)
{   
    const double EPSREL = 1e-10; 
   const double EPSABS = 1e-10; 
     double tmp1 = Integrate(bind(Num_b_ST_intgrnd_10,z,Mmin,_1),Mmin,Mmax,EPSREL,EPSABS);
     //double tmp2 = Integrate(bind(ngintrd,z,Mmin,_1),Mmin,Mmax,EPSREL,EPSABS);
     double tmp = tmp1/ngint(z);
    return tmp;
}

double b_ST_20_in(double z)
{ 

     const double EPSREL = 1e-10;
    const double EPSABS = 1e-10;
     double tmp1 = Integrate(bind(Num_b_ST_intgrnd_20,z,Mmin,_1),Mmin,Mmax,EPSREL,EPSABS);
     //double tmp2 = Integrate(bind(ngintrd,z,Mmin,_1), Mmin,Mmax,EPSREL,EPSABS);

     double tmp = tmp1/ngint(z);
    return tmp;
}

double HoDBiasModel::ng(double z)
{
	 return ngint(z);
}




 int k_bins2 =10;
 std::vector <double> bias_array1(k_bins2);
 std::vector <double> bias_array2(k_bins2);
 std::vector <double> z_array(k_bins2);
 

void HoDBiasModel::InitialzieBiasloop(void)
{
	//double  EPSREL = 1e-8;
	double Zmin2;
	double Zmax2;
	Zmin2 = log(1e-5);
	Zmax2 = log(5.0);
	//const double QMIN = 1e-4;
	//const double QMAX = 1e5;	
#pragma omp parallel private(thread_number)
	{
#pragma omp for schedule(static) nowait
		for(int i = 0;i< k_bins2 ;i++)
		{
			z_array[i] = Zmin2 +i *(Zmax2-Zmin2)/k_bins2;//exp(((Zmax2-Zmin2)*((double) i)/((double)  k_bins2-1) + Zmin2 ));
            bias_array1[i] =  b_ST_10_in(z_array[i]);
            bias_array2[i] =  b_ST_20_in(z_array[i]);

            cout<<"\t i"<<"\t"<<"z"<<"\t"<<z_array[i]<<"\t"<<"b1"<<"\t"<<bias_array1[i]<<"\t"<<"b2"<<"\t"<< bias_array2[i]<<endl;
			
			
		}
		
	}
}



double HoDBiasModel::Galb1(double z)
{
	Spline<double, double> CubicSpline_b1(z_array,bias_array1);
	return CubicSpline_b1.interpolate(z)+ 1.0;
   	
}
double HoDBiasModel::Galb2(double z)
{      
	Spline<double, double> CubicSpline_b2(z_array,bias_array2);
	double tmp = 8.0*(Galb1(z)-1)/21.0 +  CubicSpline_b2.interpolate(z);
	return tmp;
   	
}
double HoDBiasModel::Galbsq(double z)
{
	double tmp = -2.0*(Galb1(z) - 1.0)/7.0;
	return tmp;
}





