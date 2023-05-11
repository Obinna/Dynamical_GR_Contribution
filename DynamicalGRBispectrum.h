

#ifndef  _DynamicalGR
#define  _DynamicalGR



class DynamicalGR
{
    public:


   double LinearPowerSpectrum(int x,  double z,double fnl, double k);
   double RealSpaceBispectrumTotalMatter(int x, double z,double fnl,double k1, double r, double mu12);
   double ReducedRealSpaceBispectrumTotalMatter(int x,double z,double fnl,double k1, double r, double mu12); 

   double RealSpaceBispectrumPoisson(int x, double z,double fnl,double k1, double r, double mu12);
   double ReducedRealSpaceBispectrumPoisson(int x,double z,double fnl,double k1, double r, double mu12); 
  double DiffReducedRealSpaceBispectrum(int y,double z,double fnl,double k1, double r, double mu12);
  double DiffRealSpaceBispectrum(int y,double z,double fnl,double k1, double r, double mu12);
double RealSpaceBispectrumTotalMatterk1k2k3(int x,double z, double fnl, double k1, double k2, double k3);


double  RealSpaceKernelTram(int y,double z, double fnl, double k1, double k2, double k3);

double  MatterBispectrum(int x,int y,double z, double fnl, double k1, double r, double mu12);
    DynamicalGR();
 // void SetPowerSpectrum(void);
    void ReadWrite(void);
   ~DynamicalGR();










};

double mu23(double r, double mu12);

double mu13(double r, double mu12);
double PowerSpectrum(int x, double z, double fnl,double k);
double SecondOrderMatterKernel(int x, double z, double fnl, double k1, double k2, double k3,double mu12);
double LinearOrderMatterKernel(int x, double z, double fnl, double k1);



#endif



