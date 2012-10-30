#include <sys/types.h>
#include <regex.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "loadRamses.hpp"
#include "load_data.hpp"
#include "fortran.hpp"
#include <map>

using namespace std;

CosmoTool::GadgetData *CosmoTool::loadRamses(const char *name, bool quiet)
{
  GadgetData *gd = (GadgetData *)malloc(sizeof(GadgetData));
  int id = 1;
  uint32_t totPart = 0;
  
  if (!quiet)
    cout << "Detecting number of files and particles..." << endl;
  while (1)
    {
      ostringstream ss_fname;  
      ss_fname << name << setfill('0') << setw(5) << id;
      string fname = ss_fname.str();

      if (!quiet)
	cout <<  " ... " << fname << endl;

      try
	{
	  UnformattedRead infile(fname);
	  
	  int nCpu, ndim, nPar;
	  
	  infile.beginCheckpoint();
	  nCpu = max(1,infile.readInt32());
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  ndim = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  nPar = infile.readInt32();
	  infile.endCheckpoint();
	  
	  if (!quiet)
	    cout << "    NUMCPU=" << nCpu <<  " NUMDIM=" << ndim << " NumPart=" << nPar << endl;
	  
	  totPart += nPar;
	  id++;
	}
      catch (const NoSuchFileException& e)
	{
	  break;
	}
    }

  assert (totPart <= ((~(size_t)0)/sizeof(ParticleState)));
  size_t memSize = sizeof(ParticleState)*(size_t)totPart;
  if (!quiet)
    cout << " Needing " << memSize / 1024 << " kBytes" << endl;
  gd->particles = (ParticleState *)malloc(memSize);
  assert(gd->particles != 0);
  gd->NumPart = totPart;
  gd->ntot_withmasses = totPart;
  gd->header.npart[0] = totPart;
  for (int k = 1; k < 6; k++)
    {
      gd->header.npart[k] = 0;
      gd->header.mass[k] = 0;
    }
  gd->header.num_files = id;
  gd->header.BoxSize = 200.0 * 1000; // kPc
  gd->header.Omega0 = 0.30; // ????
  gd->header.OmegaLambda = 0.70; // ????
  gd->header.HubbleParam = 0.70; // ????

  if (!quiet)
    cout << " Total number part=" << totPart << endl
	 << "Loading particles ..." << endl;

  uint32_t curPos = 0;
  id = 1;
  while (1)
    {
      ostringstream ss_fname;  
      ss_fname << name << setfill('0') << setw(5) << id;
      string fname = ss_fname.str();

      if (!quiet)
	cout <<  " ... " << id;
      cout.flush();
      
      try {
	UnformattedRead infile(fname);
	
	int nCpu, ndim, nPar;
	
	
	infile.beginCheckpoint();
	nCpu = infile.readInt32();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	ndim = infile.readInt32();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	nPar = infile.readInt32();
	infile.endCheckpoint();
	
	ParticleState *s = &gd->particles[curPos];
	
	for (uint32_t i = nPar; i > 0; i--,s++)
	  s->Mass = 1.0;
	
	s = &gd->particles[curPos];
	infile.beginCheckpoint();
	for (uint32_t i = nPar; i > 0; i--)
	  {
	    s->Pos[0] = infile.readReal32();
	    s->Pos[0] *= gd->header.BoxSize;
	    s++;
	  }
	infile.endCheckpoint();

	s = &gd->particles[curPos];
	infile.beginCheckpoint();
	for (uint32_t i = nPar; i > 0; i--)
	{
	  s->Pos[1] = infile.readReal32();
	  s->Pos[1] *= gd->header.BoxSize;
	  s++;
	}
	infile.endCheckpoint();

	s = &gd->particles[curPos];
	infile.beginCheckpoint();
	for (uint32_t i = nPar; i > 0; i--)
	  {
	    s->Pos[2] = infile.readReal32();
	    s->Pos[2] *= gd->header.BoxSize;
	    s++;
	  }
	infile.endCheckpoint();

 	infile.beginCheckpoint();
	for (uint32_t i = 0; i < nPar; i++)
	  {
	    gd->particles[curPos+i].Vel[0] = infile.readReal32();
	  }
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	for (uint32_t i = 0; i < nPar; i++)
	  {
	    gd->particles[curPos+i].Vel[1] = infile.readReal32();
	  }
 	infile.endCheckpoint();

 	infile.beginCheckpoint();
	for (uint32_t i = 0; i < nPar; i++)
	  {
	    gd->particles[curPos+i].Vel[2] = infile.readReal32();
	  }
	infile.endCheckpoint();

	curPos += nPar;
      }
      catch (const NoSuchFileException& e)
	{
	  break;
	}

 
      id++;
    }

  if (!quiet)
    cout << endl;

  return gd;
}

typedef struct 
{
  double unitLength;
  double aexp;
  double boxSize;
  double unit_t;
  double omega_m;
  double omega_lambda;
} InfoData;

int readInfoFile(const char *basename, int outputId, InfoData& info)
{
  ostringstream ss_fname;  
  ss_fname << basename << "/info_" << setfill('0') << setw(5) << outputId << ".txt";

//  cout << "Opening info file " << ss_fname.str() << endl;
  ifstream infile(ss_fname.str().c_str());
  if (!infile)
    return 0;

  int err;
  regex_t unit_l_rx;

  //  const char *pattern = "^unit_l[ ]*=[ ]*([0-9\\.E+\\-]+)";
  const char *pattern = "^([A-Za-z_]+)[ ]*=[ ]*([0-9\\.E+\\-]+)";

  err = regcomp (&unit_l_rx, pattern, REG_EXTENDED);
  cout << unit_l_rx.re_nsub << endl;
  if (err)
    {
      char errString[255];
      regerror(err, &unit_l_rx, errString, sizeof(errString));
      cout << errString << endl;
      abort();
    }

  map<string,double> infoMap;
  string line;
  while (getline(infile, line))
    {
      regmatch_t allMatch[4];
      if (!regexec(&unit_l_rx, line.c_str(), 4, allMatch, 0))
	{
	  uint32_t start0 = allMatch[1].rm_so, end0 = allMatch[1].rm_eo;
	  uint32_t start1 = allMatch[2].rm_so, end1 = allMatch[2].rm_eo;

	  string keyword = line.substr(start0, end0-start0);
	  istringstream iss(line.substr(start1, end1-start1));
	  double unitLength;

	  iss >> unitLength;

	  infoMap[keyword] = unitLength;
	}
    }
  
  regfree(&unit_l_rx);

  info.unitLength = infoMap["unit_l"];
  info.aexp = infoMap["aexp"];
  info.boxSize = infoMap["boxlen"];
  info.unit_t = infoMap["unit_t"];
  info.omega_m = infoMap["omega_m"];
  info.omega_lambda = infoMap["omega_l"];

  return 1;
}

CosmoTool::SimuData *CosmoTool::loadRamsesSimu(const char *basename, int outputId, int cpuid, bool dp, int flags)
{
  CosmoTool::SimuData *data = new CosmoTool::SimuData();

  int id = 1;
  uint32_t totPart = 0;
  int nCpu = 0;
  InfoData info;

  static const double CM_IN_MPC = 3.08e24;

  if (!readInfoFile(basename, outputId, info))
    return 0;

  double hubble = info.aexp*info.aexp/info.unit_t / (1e5/CM_IN_MPC);
  double L0 = info.boxSize*info.unitLength*hubble/100/CM_IN_MPC/info.aexp;
  double unit_vel = L0*hubble/info.aexp;
  
  while (1)
    {
      ostringstream ss_fname;
      ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << id;
      string fname = ss_fname.str();

      try
        {
          UnformattedRead infile(fname);

          int ndim, nPar;

          infile.beginCheckpoint();
          nCpu = max(1,infile.readInt32());
          infile.endCheckpoint();

          infile.beginCheckpoint();
          ndim = infile.readInt32();
          infile.endCheckpoint();

          infile.beginCheckpoint();
          nPar = infile.readInt32();
          infile.endCheckpoint();

          totPart += nPar;
          id++;
        }
      catch (const NoSuchFileException& e)
        {
          break;
        }
    }


  data->BoxSize = L0;
  data->time = info.aexp;
  data->NumPart = 0;
  data->TotalNumPart = totPart;
  data->Hubble = hubble;
  data->Omega_M = info.omega_m;
  data->Omega_Lambda = info.omega_lambda;

  if (flags == 0)
    return data;

  if (cpuid < 0)
    cpuid = 1;
  else cpuid++;

  uint32_t curPos = 0;
  
  ostringstream ss_fname;
  ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << cpuid;
  
  string fname = ss_fname.str();
    
  try
    {
      UnformattedRead infile(fname);
      
      int nCpu, ndim, nPar;
      
      infile.beginCheckpoint();
      nCpu = infile.readInt32();
      infile.endCheckpoint();
      
      infile.beginCheckpoint();
      ndim = infile.readInt32();
      infile.endCheckpoint();
      
      infile.beginCheckpoint();
      data->NumPart = nPar = infile.readInt32();
      infile.endCheckpoint();

      infile.beginCheckpoint();
      for (int i = 0; i < 4; i++) infile.readInt32();
      infile.endCheckpoint();
      
      infile.beginCheckpoint();
      infile.readInt32();
      infile.endCheckpoint();
      
      infile.beginCheckpoint();
      infile.readReal64();
      infile.endCheckpoint();

      infile.beginCheckpoint();
      infile.readReal64();
      infile.endCheckpoint();
      
      infile.beginCheckpoint();
      infile.readInt32();
      infile.endCheckpoint();
      
      if (flags & NEED_POSITION)
	{
	  data->Pos[0] = new float[nPar];
	  data->Pos[1] = new float[nPar];
	  data->Pos[2] = new float[nPar];
	}
      if (flags & NEED_VELOCITY)
	{
	  data->Vel[0] = new float[nPar];
	  data->Vel[1] = new float[nPar];
	  data->Vel[2] = new float[nPar];
	}
      if (flags & NEED_GADGET_ID)
	{
	  data->Id = new int[nPar];
	}

      for (int k = 0; k < 3; k++)
	{
	  infile.beginCheckpoint();
	  if (flags & NEED_POSITION)
	    {
	      for (uint32_t i = 0; i < nPar; i++)
		{
		  data->Pos[k][i] = dp ? infile.readReal64() : infile.readReal32();
		  data->Pos[k][i] *= data->BoxSize;
		}
	      infile.endCheckpoint();
	    }
	  else
	    infile.endCheckpoint(true);
	}
      
      for (int k = 0; k < 3; k++) {
	infile.beginCheckpoint();
	if (flags & NEED_VELOCITY)
	  {
	    for (uint32_t i = 0; i < nPar; i++)
	      {
		data->Vel[k][i] = dp ? infile.readReal64() : infile.readReal32();
		data->Vel[k][i] *= unit_vel;
	      }
	    infile.endCheckpoint();
	  }
	else
	  infile.endCheckpoint(true);
      }
      
      float minMass = INFINITY;
      infile.beginCheckpoint();
      for (uint32_t i = nPar; i > 0; i--)
	{
	  float dummyF = dp ? infile.readReal64() : infile.readReal32();
	  if (dummyF < minMass) minMass = dummyF;
	}
      infile.endCheckpoint();
      
      infile.beginCheckpoint();
      if (flags & NEED_GADGET_ID)
	{
	  for (uint32_t i = 0; i < nPar; i++)
	    data->Id[i] = infile.readInt32();
	  
	  infile.endCheckpoint();
	}
      else
	infile.endCheckpoint(true);
      
      curPos += nPar;
    }
  catch (const NoSuchFileException& e)  
    {
	cerr << "No such file " << fname << endl;
      delete data;
      return 0;
    }
  
  return data;
}







CosmoTool::PurePositionData *CosmoTool::loadRamsesPosition(const char *basename, int outputId, bool quiet, bool dp)
{
  PurePositionData *gd = (PurePositionData *)malloc(sizeof(PurePositionData));
  int id = 1;
  uint32_t totPart = 0;
  int nCpu = 0;
  InfoData info;
  
  static const double CM_IN_MPC = 3.08e24;

  if (!readInfoFile(basename, outputId, info))
    return 0;

  double hubble = info.aexp*info.aexp/info.unit_t / (1e5/CM_IN_MPC);
  double L0 = info.boxSize*info.unitLength*hubble/100/CM_IN_MPC/info.aexp;
  if (!quiet)
    {
      cout << "L0=" << L0 << " Mpc" << endl;
      cout << "H=" << hubble << " km/s/Mpc" << endl;
    }

  if (!quiet)
    cout << "Detecting number of files and particles..." << endl;
  while (1)
    {
      ostringstream ss_fname;  
      ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << id;
      string fname = ss_fname.str();
      
      if (!quiet)
	cout <<  " ... " << fname << endl;
      

      try
	{
	  UnformattedRead infile(fname);
	  
	  int ndim, nPar;
	  
	  infile.beginCheckpoint();
	  nCpu = max(1,infile.readInt32());
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  ndim = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  nPar = infile.readInt32();
	  infile.endCheckpoint();
	  
	  if (!quiet)
	    cout << "    NUMCPU=" << nCpu <<  " NUMDIM=" << ndim << " NumPart=" << nPar << endl;
	  
	  totPart += nPar;
	  id++;
	}
      catch (const NoSuchFileException& e)
	{
	  break;
	}

    }

  assert (totPart <= ((~(size_t)0)/sizeof(FCoordinates)));
  size_t memSize = sizeof(FCoordinates)*(size_t)(totPart+totPart/nCpu);
  if (!quiet)
    cout << " Needing " << memSize / 1024 << " kBytes" << endl;
  gd->pos = (FCoordinates *)malloc(sizeof(FCoordinates)*totPart);
  assert(gd->pos != 0);
  gd->NumPart = totPart;
  gd->BoxSize = L0*1000; 
  gd->hubble = hubble;

  if (!quiet)
    cout << " Total number part=" << totPart << endl
	 << "Loading particles ..." << endl;
  
  uint32_t curPos = 0;
  id = 1;
  while (1)
    {
      ostringstream ss_fname;  
      ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << id;
     
      string fname = ss_fname.str();
      int *idP;

      if (!quiet)
	(cout <<  " ... " << id).flush();
      
      try
	{
	  UnformattedRead infile(fname);
	  
	  int nCpu, ndim, nPar;
	  
	  infile.beginCheckpoint();
	  nCpu = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  ndim = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  nPar = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
          for (int i = 0; i < 4; i++) infile.readInt32();
	  infile.endCheckpoint();

	  infile.beginCheckpoint();
          infile.readInt32();
	  infile.endCheckpoint();

	  infile.beginCheckpoint();
          infile.readReal64();
	  infile.endCheckpoint();

	  infile.beginCheckpoint();
          infile.readReal64();
	  infile.endCheckpoint();

	  infile.beginCheckpoint();
          infile.readInt32();
	  infile.endCheckpoint();

	  FCoordinates *s = new FCoordinates[nPar];
	  
	  for (int k = 0; k < 3; k++)
	    {
	      infile.beginCheckpoint();
	      for (uint32_t i = 0; i < nPar; i++)
		{
		  s[i][k] = dp ? infile.readReal64() : infile.readReal32();
		  s[i][k] *= gd->BoxSize;
		}
	      infile.endCheckpoint();
	    }

	  // SKIP VELOCITIES
	  for (int k = 0; k < 3; k++) {	
	    infile.beginCheckpoint();
	    for (uint32_t i = nPar; i > 0; i--)
	      {
		(dp ? infile.readReal64() : infile.readReal32());
	      }
	    infile.endCheckpoint();
	  }
	  
	  float minMass = INFINITY;
	  infile.beginCheckpoint();
	  for (uint32_t i = nPar; i > 0; i--)
	    {
	      float dummyF = dp ?  infile.readReal64() : infile.readReal32();
	      if (dummyF < minMass) minMass = dummyF;
	    }
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  for (uint32_t i = 0; i < nPar; i++)
	    {
	      int id = infile.readInt32();
	      id--;
	      assert(id >= 0);
	      assert(id < totPart);
	      memcpy(&gd->pos[id], &s[i], sizeof(FCoordinates));
	    }
	  infile.endCheckpoint();
      
	  delete[] s;
	  
	  curPos += nPar;
	} catch (const NoSuchFileException& e)
	{
	  break;
	}
      id++;
    }

  if (!quiet)
    cout << endl;

  return gd;
}

CosmoTool::PhaseSpaceData *CosmoTool::loadRamsesPhase(const char *basename, int outputId, bool quiet)
{
  PhaseSpaceData *gd = (PhaseSpaceData *)malloc(sizeof(PhaseSpaceData));
  int id = 1;
  uint32_t totPart = 0;
  int nCpu = 0;
  InfoData info;
  
  static const double CM_IN_MPC = 3.08e24;

  if (!readInfoFile(basename, outputId, info))
    return 0;

  double hubble = info.aexp*info.aexp/info.unit_t / (1e5/CM_IN_MPC);
  double L0 = info.boxSize*info.unitLength*hubble/(100*CM_IN_MPC)/info.aexp;
  double unit_vel = 100*L0/info.aexp;
  if (!quiet) {
    cout << "L0=" << L0 << " Mpc" << endl;
    cout << "H=" << hubble << " km/s/Mpc" << endl;
    cout << "unit_vel=" << unit_vel << " km/s" << endl;
  }

  if (!quiet)
    cout << "Detecting number of files and particles..." << endl;
  while (1)
    {
      ostringstream ss_fname;  
      ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << id;
      string fname = ss_fname.str();
      
      if (!quiet)
	cout <<  " ... " << fname << endl;
      

      try
	{
	  UnformattedRead infile(fname);
	  
	  int ndim, nPar;
	  
	  infile.beginCheckpoint();
	  nCpu = max(1,infile.readInt32());
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  ndim = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  nPar = infile.readInt32();
	  infile.endCheckpoint();
	  
	  if (!quiet)
	    cout << "    NUMCPU=" << nCpu <<  " NUMDIM=" << ndim << " NumPart=" << nPar << endl;
	  
	  totPart += nPar;
	  id++;
	}
      catch (const NoSuchFileException& e)
	{
	  break;
	}

    }

  assert (totPart <= ((~(size_t)0)/sizeof(FCoordinates)));
  size_t memSize = sizeof(FCoordinates)*(size_t)(totPart+totPart/nCpu);
  if (!quiet)
    cout << " Needing " << memSize / 1024 << " kBytes" << endl;
  gd->pos = (FCoordinates *)malloc(sizeof(FCoordinates)*totPart);
  assert(gd->pos != 0);
  gd->vel = (FCoordinates *)malloc(sizeof(FCoordinates)*totPart);
  assert(gd->vel != 0);
  gd->NumPart = totPart;
  gd->BoxSize = L0*1000; 
  gd->hubble = hubble;
  
  if (!quiet)
    cout << " Total number part=" << totPart << endl
	 << "Loading particles ..." << endl;
  
  uint32_t curPos = 0;
  id = 1;
  while (1)
    {
      ostringstream ss_fname;  
      ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << id;
     
      string fname = ss_fname.str();
      int *idP;

      if (!quiet)
	(cout <<  " ... " << id).flush();
      
      try
	{
	  UnformattedRead infile(fname);
	  
	  int nCpu, ndim, nPar;
	  
	  infile.beginCheckpoint();
	  nCpu = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  ndim = infile.readInt32();
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  nPar = infile.readInt32();
	  infile.endCheckpoint();
	  

	  FCoordinates *s = new FCoordinates[nPar];
	  FCoordinates *vel = new FCoordinates[nPar];
	  
	  for (int k = 0; k < 3; k++)
	    {
	      infile.beginCheckpoint();
	      for (uint32_t i = 0; i < nPar; i++)
		{
		  s[i][k] = infile.readReal32();
		  s[i][k] *= gd->BoxSize;
		}
	      infile.endCheckpoint();
	    }

	  // SKIP VELOCITIES
	  for (int k = 0; k < 3; k++) {	
	    infile.beginCheckpoint();
	    for (uint32_t i = 0; i < nPar; i++)
	      {
		vel[i][k] = infile.readReal32()*unit_vel;
	      }
	    infile.endCheckpoint();
	  }
	  
	  float minMass = INFINITY;
	  infile.beginCheckpoint();
	  for (uint32_t i = nPar; i > 0; i--)
	    {
	      float dummyF = infile.readReal32();
	      if (dummyF < minMass) minMass = dummyF;
	    }
	  infile.endCheckpoint();
	  
	  infile.beginCheckpoint();
	  for (uint32_t i = 0; i < nPar; i++)
	    {
	      int id = infile.readInt32();
	      id--;
	      assert(id >= 0);
	      assert(id < totPart);
	      memcpy(&gd->pos[id], &s[i], sizeof(FCoordinates));
	      memcpy(&gd->vel[id], &vel[i], sizeof(FCoordinates));
	    }
	  infile.endCheckpoint();
      
	  delete[] vel;
	  delete[] s;
	  
	  curPos += nPar;
	} catch (const NoSuchFileException& e)
	{
	  break;
	}
      id++;
    }

  if (!quiet)
    cout << endl;

  return gd;
}


CosmoTool::PhaseSpaceDataID *CosmoTool::loadRamsesPhase1(const char *basename, int outputId, int cpuid, bool dp, bool quiet)
{
  PhaseSpaceDataID *gd = (PhaseSpaceDataID *)malloc(sizeof(PhaseSpaceDataID));
  int id = 1;
  uint32_t totPart = 0;
  int nCpu = 0;
  InfoData info;
  
  static const double CM_IN_MPC = 3.08e24;

  if (!readInfoFile(basename, outputId, info))
    return 0;

  double hubble = info.aexp*info.aexp/info.unit_t / (1e5/CM_IN_MPC);
  double L0 = info.boxSize*info.unitLength*hubble/100/CM_IN_MPC/info.aexp;
  double unit_vel = 100*L0/info.aexp;
  if (!quiet) {
    cout << "L0=" << L0 << " Mpc" << endl;
    cout << "H=" << hubble << " km/s/Mpc" << endl;
    cout << "unitvel=" << unit_vel << " km/s" << endl;
  }

  if (!quiet)
    cout << "Detecting number of files and particles..." << endl;
  id = cpuid;
  {
    ostringstream ss_fname;  
    ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << id;
    string fname = ss_fname.str();
    
    if (!quiet)
      cout <<  " ... " << fname << endl;
    
    
    try
      {
	UnformattedRead infile(fname);
	
	int ndim, nPar;
	
	infile.beginCheckpoint();
	nCpu = max(1,infile.readInt32());
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	ndim = infile.readInt32();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	nPar = infile.readInt32();
	infile.endCheckpoint();
	
	if (!quiet)
	  cout << "    NUMCPU=" << nCpu <<  " NUMDIM=" << ndim << " NumPart=" << nPar << endl;
	
	totPart += nPar;
	id++;
      }
    catch (const NoSuchFileException& e)
      {
	return 0;
      }
  }
  
  assert (totPart <= ((~(size_t)0)/sizeof(FCoordinates)));
  size_t memSize = sizeof(FCoordinates)*(size_t)(totPart+totPart/nCpu);
  if (!quiet)
    cout << " Needing " << memSize / 1024 << " kBytes" << endl;
  gd->pos = (FCoordinates *)malloc(sizeof(FCoordinates)*totPart);
  assert(gd->pos != 0);
  gd->vel = (FCoordinates *)malloc(sizeof(FCoordinates)*totPart);
  assert(gd->vel != 0);
  gd->ID = (int *)malloc(sizeof(int)*totPart);
  assert(gd->ID != 0);
  gd->NumPart = totPart;
  gd->BoxSize = L0*1000; 
  gd->hubble = hubble;
  
  uint32_t curPos = 0;
  id = cpuid;
  {
    ostringstream ss_fname;  
    ss_fname << basename << "/part_" << setfill('0') << setw(5) << outputId << ".out" << setfill('0') << setw(5) << id;
    
    string fname = ss_fname.str();
    int *idP;
    
    if (!quiet)
      (cout <<  " ... " << id).flush();
    
    try
      {
	UnformattedRead infile(fname);
	
	int nCpu, ndim, nPar;
	
	infile.beginCheckpoint();
	nCpu = infile.readInt32();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	ndim = infile.readInt32();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	nPar = infile.readInt32();
	infile.endCheckpoint();

	infile.beginCheckpoint();
	for (int i = 0; i < 4; i++) infile.readInt32();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	infile.readInt32();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	infile.readReal64();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	infile.readReal64();
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	infile.readInt32();
	infile.endCheckpoint();

	gd->pos = new FCoordinates[nPar];
	gd->vel = new FCoordinates[nPar];
	
	for (int k = 0; k < 3; k++)
	  {
	    infile.beginCheckpoint();
	    for (uint32_t i = 0; i < nPar; i++)
	      {
		gd->pos[i][k] = (dp ? infile.readReal64() : infile.readReal32())*gd->BoxSize;
	      }
	    infile.endCheckpoint();
	  }
	
	// SKIP VELOCITIES
	for (int k = 0; k < 3; k++) {	
	  infile.beginCheckpoint();
	  for (uint32_t i = 0; i < nPar; i++)
	    {
	      gd->vel[i][k] = (dp ? infile.readReal64() : infile.readReal32())*unit_vel;
	    }
	  infile.endCheckpoint();
	}
	
	float minMass = INFINITY;
	infile.beginCheckpoint();
	for (uint32_t i = nPar; i > 0; i--)
	  {
	    float dummyF = (dp ? infile.readReal64() : infile.readReal32());
	    if (dummyF < minMass) minMass = dummyF;
	  }
	infile.endCheckpoint();
	
	infile.beginCheckpoint();
	for (uint32_t i = 0; i < nPar; i++)
	  {
	    gd->ID[i] = infile.readInt32();
	  }
	infile.endCheckpoint();
	
	curPos += nPar;
      } 
    catch (const NoSuchFileException& e) 
      {
	return 0;
      }
  }
  
  if (!quiet)
    cout << endl;

  return gd;
}
