

#ifndef  _HoDBiasModel
#define  _HoDBiasModel
//#include<gsl/gsl_spline.h>
//#include"constants.h"


class HoDBiasModel
{
    public:



 double Galb1(double z);
 double Galb2(double z);
 double Galbsq(double z);
  double ng(double z);


  HoDBiasModel();
  void InitialzieBiasloop(void);
  //void ReadWrite1(void);
   ~HoDBiasModel();




        

};
#endif 