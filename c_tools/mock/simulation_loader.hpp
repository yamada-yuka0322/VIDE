#ifndef _MOCK_SIMULATION_LOADER_HPP
#define _MOCK_SIMULATION_LOADER_HPP

#include <string>
#include <CosmoTool/loadSimu.hpp>

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
 
  void applyTransformations(CosmoTool::SimuData *s);

public:
  virtual ~SimulationLoader() {}
  
  void doRedshift(bool set = true) { do_redshift = set; }
  void setVelAxis(int axis) { redshift_axis = axis; }

  virtual CosmoTool::SimuData *getHeader() = 0;
  virtual int num_files() = 0;
  virtual CosmoTool::SimuData* loadFile(int id) = 0;

};

template<typename T>
void delete_adaptor(void *ptr)
{
  T *ptr_T = reinterpret_cast<T *>(ptr);

  delete[] ptr_T;
}


// Unit length is the size of one Mpc in the simulation units
SimulationLoader *gadgetLoader(const std::string& snapshot, double Mpc_unitLength, int flags);
SimulationLoader *flashLoader(const std::string& snapshot, int flags);
SimulationLoader *multidarkLoader(const std::string& snapshot, int flags);
SimulationLoader *ramsesLoader(const std::string& snapshot, int baseid, bool double_precision, int flags);


#endif
