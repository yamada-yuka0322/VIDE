#ifndef __COSMOTOOLBOX_HPP
#define __COSMOTOOLBOX_HPP

#include <map>
#include <string>

namespace CosmoTool
{
  static const int NEED_GADGET_ID = 1;
  static const int NEED_POSITION = 2;
  static const int NEED_VELOCITY = 4;
  static const int NEED_TYPE = 8;

  struct SimuParticle
  {
    float Pos[3];
    float Vel[3];
    int type;
    int id;

    bool flag_vel, flag_type, flag_id;
  };

  typedef bool (*SimuFilter)(const SimuParticle& p);

  class SimuData 
  {
  public:
    typedef void (*FreeFunction)(void *);
    typedef std::map<std::string, std::pair<void *, FreeFunction> > AttributeMap;

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
    int *type;

    AttributeMap attributes;
    
  public:
    SimuData() : Id(0),NumPart(0),type(0)  { Pos[0]=Pos[1]=Pos[2]=0; Vel[0]=Vel[1]=Vel[2]=0; }
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

      for (AttributeMap::iterator i = attributes.begin();
	   i != attributes.end();
	   ++i)
	{
	  if (i->second.second)
	    i->second.second(i->second.first);
	}
    }    
    
    template<typename T>
    T *as(const std::string& n) 
    {
      AttributeMap::iterator i = attributes.find(n);
      if (i == attributes.end())
	return 0;
      
      return reinterpret_cast<T *>(i->second.first);
    }

    void new_attribute(const std::string& n, void *p, FreeFunction free_func)
    {
      AttributeMap::iterator i = attributes.find(n);
      if (i != attributes.end())
	{
	  if (i->second.second)
	    i->second.second(i->second.first);
	}
      attributes[n] = std::make_pair(p, free_func);
    }
    
  };

};

#endif
