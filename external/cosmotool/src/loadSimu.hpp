/*+
This is CosmoTool (./src/loadSimu.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

guilhem.lavaux@gmail.com

This software is a computer program whose purpose is to provide a toolbox for cosmological
data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
+*/

#ifndef __COSMOTOOLBOX_HPP
#define __COSMOTOOLBOX_HPP

#include <sys/types.h>
#include <map>
#include <string>

namespace CosmoTool
{
  static const int NEED_GADGET_ID = 1;
  static const int NEED_POSITION = 2;
  static const int NEED_VELOCITY = 4;
  static const int NEED_TYPE = 8;
  static const int NEED_MASS = 16;
  static const int NEED_DOUBLE_PRECISION = 32;

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

    bool noAuto;

    float BoxSize;
    float time;
    float Hubble;

    float Omega_M;
    float Omega_Lambda;

    ssize_t NumPart;
    ssize_t TotalNumPart;
    int64_t *Id;
    float *Pos[3];
    float *Vel[3];
    float *Mass;
    int *type;

    AttributeMap attributes;
    
  public:
    SimuData() : Mass(0), Id(0),NumPart(0),type(0),noAuto(false)  { Pos[0]=Pos[1]=Pos[2]=0; Vel[0]=Vel[1]=Vel[2]=0; }
    ~SimuData() 
    {
        if (!noAuto)  {
            for (int j = 0; j < 3; j++) {
                if (Pos[j])
                    delete[] Pos[j];
                if (Vel[j])
                    delete[] Vel[j];
            }
            if (type)
                delete[] type;
            if (Id)
                delete[] Id;
            if (Mass)
                delete[] Mass;
        }
        for (AttributeMap::iterator i = attributes.begin();
             i != attributes.end();
             ++i) {
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
        if (i != attributes.end()) {
            if (i->second.second)
                i->second.second(i->second.first);
        }
        attributes[n] = std::make_pair(p, free_func);
    }
    
  };

};

#endif
