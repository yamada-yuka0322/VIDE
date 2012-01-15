#include <iostream>
#include <fstream>
#include <netcdfcpp.h>
#include <vector>
#include <set>

using namespace std;

int main(int argc, char **argv)
{
  char *n =  argv[1];
  char *n2 = argv[2];
  ifstream f(n);
  ofstream f2(n2);
  int np, nvp;
  int *ids;
  float *vols;
  set<int> mask_id;

  NcFile fp("params.nc");
  NcVar *v_id = fp.get_var("particle_ids");
  long int *e_N = v_id->edges();
  int *pid = new int[e_N[0]];

  v_id->get(pid, e_N[0]);

  for (int i = 0; i < e_N[0]; i++)
    mask_id.insert(pid[i]);
  delete[] pid;
  delete[] e_N;

  f.read((char*)&np, sizeof(int));
  f.read((char*)&nvp,sizeof(int));
  f2.write((char*)&np, sizeof(int));

  ids = new int[nvp];
  vols = new float[nvp];
  f.read((char*)ids, sizeof(int)*nvp);
  f.read((char *)vols, sizeof(float)*nvp);

  vector< vector<int> > adjs;
  vector<int> adj0, new_ids;
  vector<float>  new_vols;
  for (int i = 0; i < nvp; i++)
    {
      int nadj;
      bool  masked_particle = false;


      if (mask_id.find(ids[i]) != mask_id.end())
	{
	  masked_particle = true;
	}

      f.read((char*)&nadj, sizeof(int));
      adj0.resize(nadj);
      
      if (nadj != 0)
	{
	  f.read((char*)&adj0[0], sizeof(float)*nadj);
	  for (int j = 0; j < nadj; j++)
	    if (mask_id.find(adj0[j]) != mask_id.end())
	      masked_particle = true;
	}
      
      if (!masked_particle)
	{
	  new_vols.push_back(vols[i]);
	  new_ids.push_back(ids[i]);
	  adjs.push_back(adj0);
	}
    }

  nvp = new_ids.size();
  f2.write((char*)&nvp, sizeof(int));
  f2.write((char*)&new_ids[0], sizeof(int)*nvp);
  f2.write((char*)&new_vols[0], sizeof(float)*nvp);
  
  for (int i = 0; i < nvp; i++)
    {
      int nadj = adjs[i].size();
      f2.write((char*)&nadj, sizeof(int));
      if (nadj != 0)
	f2.write((char*)&adjs[i], sizeof(int)*nadj);
    }

  return 0;
}
