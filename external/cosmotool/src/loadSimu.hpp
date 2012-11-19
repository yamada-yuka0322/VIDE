#ifndef __COSMOTOOLBOX_HPP
#define __COSMOTOOLBOX_HPP


namespace CosmoTool
{
  static const int NEED_GADGET_ID = 1;
  static const int NEED_POSITION = 2;
  static const int NEED_VELOCITY = 4;
  static const int NEED_TYPE = 8;

  class SimuData 
  {
  public:
    float BoxSize;
    float time;
    float Hubble;

    float Omega_M;
    float Omega_Lambda;

    long NumPart;
    long TotalNumPart;
    int *Id;
    float *Pos[3];
    float *Vel[3];
    float *uniqueID;
    int *type;
  public:
    SimuData() : Id(0),NumPart(0),type(0), uniqueID(0)  { Pos[0]=Pos[1]=Pos[2]=0; Vel[0]=Vel[1]=Vel[2]=0; uniqueID=0;}
    ~SimuData() 
    {
      for (int j = 0; j < 3; j++)
	{
	  if (Pos[j])
	    delete[] Pos[j];
	  if (Vel[j])
	    delete[] Vel[j];
	}
      if (type)
	delete[] type;
      if (Id)
	delete[] Id;
      if (uniqueID)
  delete[] uniqueID;
    }    
  };

};

#endif
