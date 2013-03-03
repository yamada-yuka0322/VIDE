/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/mock/loaders/simulation_loader.hpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 P. M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/



#ifndef _MOCK_SIMULATION_LOADER_HPP
#define _MOCK_SIMULATION_LOADER_HPP

#include <string>
#include <algorithm>
#include <CosmoTool/loadSimu.hpp>

struct SingleParticle
{
  float Pos[3];
  float Vel[3];
  long Index;
  long ID;
};

class SimulationPreprocessor
{
public:
  SimulationPreprocessor() {}
  virtual ~SimulationPreprocessor() {}

  virtual long getEstimatedPostprocessed(long numParticles) { return numParticles; };
  virtual bool accept(const SingleParticle& p) { return true; }
  virtual void reset() {}
};

class SimulationLoader
{
protected:
  bool do_redshift;
  int redshift_axis;

  SimulationLoader() 
  {
    do_redshift = false;
    redshift_axis = 2;
  }

  template<typename T> void reallocArray(T *& a, long newSize, long toCopy)
  {
    T *b = new T[newSize];
    if (a != 0)
     {
       std::copy(a, a+toCopy, b);
       delete[] a;
     }
   a = b;
  } 



  void reallocSimu(CosmoTool::SimuData *s, long newNumPart);
 

  void basicPreprocessing(CosmoTool::SimuData *d, SimulationPreprocessor *preproc);
  void applyTransformations(CosmoTool::SimuData *s);

  void copyParticleToSimu(const SingleParticle& p, CosmoTool::SimuData *s, long index)
  {
    s->Pos[0][index] = p.Pos[0];
    s->Pos[1][index] = p.Pos[1];
    s->Pos[2][index] = p.Pos[2];
    if (s->Vel[0])
      s->Vel[0][index] = p.Vel[0];
    if (s->Vel[1])
      s->Vel[1][index] = p.Vel[1];
    if (s->Vel[2])
      s->Vel[2][index] = p.Vel[2];
    s->Id[index] = p.ID;
  }
public:
  virtual ~SimulationLoader() {}
  
  void doRedshift(bool set = true) { do_redshift = set; }
  void setVelAxis(int axis) { redshift_axis = axis; }

  virtual CosmoTool::SimuData *getHeader() = 0;
  virtual int num_files() = 0;
  virtual CosmoTool::SimuData* loadFile(int id) = 0;

  
  template<typename T>
  void filteredCopy(T *a, bool *accepted, long N)
  {
    long i = 0, j = 0;

    if (a == 0)
      return;

    while (i < N)
      {
	if (!accepted[i])
	  {
	    i++;
	    continue;
	  }

	a[j] = a[i];
	j++;
	i++;
      }
  }
};


template<typename T>
void delete_adaptor(void *ptr)
{
  T *ptr_T = reinterpret_cast<T *>(ptr);

  delete[] ptr_T;
}


// Unit length is the size of one Mpc in the simulation units
SimulationLoader *gadgetLoader(const std::string& snapshot, double Mpc_unitLength, int flags, SimulationPreprocessor *p);
SimulationLoader *flashLoader(const std::string& snapshot, int flags, SimulationPreprocessor *p);
SimulationLoader *multidarkLoader(const std::string& snapshot, SimulationPreprocessor *p);
SimulationLoader *ramsesLoader(const std::string& snapshot, int baseid, bool double_precision, int flags, SimulationPreprocessor *p);
SimulationLoader *sdfLoader(const std::string& snapshot, int flags, int num_splitting, SimulationPreprocessor *p);

/* This is required for some parallel I/O handler (thus MPI beneath it) */
void simulationLoadersInit(int& argc, char **& argv);

#endif
