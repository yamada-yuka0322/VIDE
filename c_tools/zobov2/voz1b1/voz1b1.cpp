#include <cassert>
#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include <CosmoTool/miniargs.hpp>
#include <cmath>
#include "libqhull/qhull_a.h"
#include "voz.h"
#include "voz_io.hpp"

using CosmoTool::MINIARG_STRING;
using CosmoTool::MINIARG_DOUBLE;
using CosmoTool::MINIARG_INT;
using CosmoTool::MINIARG_NULL;
using boost::format;
using namespace std;

#define DL for (d=0;d<3;d++)

bool checkParameters(int *numdiv, int *b)
{
  for (int i = 0; i < 3; i++)
    {
      if (numdiv[i] < 0 || b[i] < 0 || b[i] >= numdiv[i])
        return false;
    }
  return true;
}


struct BoxData
{
  float width[3], totwidth[3];
  float width2[3], totwidth2[3];
  float bf, s[3], g[3], c[3];
  coordT *parts;
  pid_t *orig;
  pid_t nvpall, nvp, nvpbuf;
  bool guard_added;
  double xyz_min[3], xyz_max[3];

  void prepareBox(const PositionData& pdata, int *b_id, int *numdivs);

  void checkParticle(double *xyz, bool& in_main, bool& in_buf);
  void addGuardPoints();
  void findBoundary();
};

void BoxData::checkParticle(double *xyz, bool& in_main, bool& in_buf)
{
  in_main = in_buf = true;
  for (int d = 0; d < 3; d++)
    {
      xyz[d] -= (double)c[d];
      if (xyz[d] > width2[d])
        xyz[d] -= width[d];
      if (xyz[d] < -width2[d])
        xyz[d] += width[d];

      in_buf = in_buf && (fabs(xyz[d]) < totwidth2[d]);
      in_main = in_main && (fabs(xyz[d]) <= width2[d]);
    }
  assert(!in_main || in_buf);
}

void BoxData::prepareBox(const PositionData& pdata, int *b_id, int *numdivs)
{
  double BF = std::numeric_limits<double>::max();
  guard_added = false;

  for (int i = 0; i < 3; i++)
    {
      width[i] = (pdata.xyz_max[i] - pdata.xyz_min[i])/numdivs[i];
      width2[i] = 0.5*width[i];

      totwidth[i] = width[i]+2.*bf;
      totwidth2[i] = width2[i] + bf;

      s[i] = width[i]/(float)NGUARD;
      if ((bf*bf - 2.*s[i]*s[i]) < 0.)
        {
          printf("bf = %f, s = %f.\n",bf,s[i]);
          printf("Not enough guard points for given border.\nIncrease guards to >= %f\n.",
                 sqrt(2.)*width[i]/bf);
          exit(0);
        }
      g[i] = (bf / 2.)*(1. + sqrt(1 - 2.*s[i]*s[i]/(bf*bf)));
      cout << format("s[%d] = %f, bf = %f, g[%d] = %f.") % i%  s[i] % bf % i % g[i] << endl;

      c[i] = b_id[i] * width[i];
    }
  cout.flush();

  cout << format("c: %f,%f,%f") % c[0] % c[1] % c[2] << endl;

  /* Assign temporary array*/
  nvpbuf = 0; /* Number of particles to tesselate, including buffer */
  nvp = 0; /* Without the buffer */

  for (pid_t i = 0; i < pdata.np; i++)
    {
      double xyz[3] = { pdata.xyz[0][i], pdata.xyz[1][i], pdata.xyz[2][i] };
      bool is_it_in_buf, is_it_in_main;

      checkParticle(xyz, is_it_in_main, is_it_in_buf);

      if (is_it_in_buf)
        nvpbuf++;
      if (is_it_in_main)
        nvp++;
      assert(nvp <= nvpbuf);
    }

  nvpbuf += 6*(NGUARD+1)*(NGUARD+1); /* number of guard points */
  nvpall = nvpbuf;

  parts = new coordT[3*nvpall];
  orig = new pid_t[nvpall];

  if (parts == 0)
    {
      cout << "Unable to allocate parts" << endl;
      exit(0);
    }
  if (orig == 0)
    {
      cout << "Unable to allocate orig" << endl;
      exit(0);
    }

  nvp = 0; /* nvp = number of particles without buffer */

  for (int j = 0; j < 3; j++)
    {
      xyz_min[j] = BF;
      xyz_max[j] = -BF;
    }
  for (pid_t i = 0; i < pdata.np; i++)
    {
      bool is_it_in_main, is_it_in_buf;
      double xyz[3] = { pdata.xyz[0][i], pdata.xyz[1][i], pdata.xyz[2][i] };

      checkParticle(xyz, is_it_in_main, is_it_in_buf);

      if (is_it_in_main) {
        assert(nvp < nvpall);
        for (int j = 0; j < 3; j++)
          {
            parts[3*nvp+j] = xyz[j];
            xyz_min[j] = min(xyz_min[j], xyz[j]);
            xyz_max[j] = max(xyz_max[j], xyz[j]);
          }
        orig[nvp] = i;
        nvp++;
      }
  }
  cout << format("nvp = %d") %nvp << endl;
  cout << format("x: %f,%f; y: %f,%f; z:%f,%f") % xyz_min[0] % xyz_max[0] % xyz_min[1] % xyz_max[1] % xyz_min[2] % xyz_max[2] << endl;
  nvpbuf = nvp;
  for (pid_t i = 0; i < pdata.np; i++)
    {
      bool is_it_in_main, is_it_in_buf;
      double xyz[3] = { pdata.xyz[0][i], pdata.xyz[1][i], pdata.xyz[2][i] };

      checkParticle(xyz, is_it_in_main, is_it_in_buf);

      if (is_it_in_buf && !is_it_in_main)
        {
          assert(nvpbuf < nvpall);
          for (int j = 0; j < 3; j++)
            {
              parts[3*nvpbuf+j] = xyz[j];
              xyz_min[j] = min(xyz_min[j], xyz[j]);
              xyz_max[j] = max(xyz_max[j], xyz[j]);
            }
          orig[nvpbuf] = i;

          nvpbuf++;
        }
    }
  cout << format("nvpbuf = %d") % nvpbuf << endl;
  cout << format("x: %f,%f; y: %f,%f; z:%f,%f\n") % xyz_min[0] % xyz_max[0] % xyz_min[1] % xyz_max[1] % xyz_min[2] % xyz_max[2] << endl;

  double predict = 1;
  for (int j = 0; j < 3; j++)
    predict *= (width[j]+2*bf);
  predict *= pdata.np;
  cout << format("There should be ~ %g points; there are %d\n") % predict % nvpbuf << endl;
}

void BoxData::addGuardPoints()
{
  int rot_g[3][3]  = { {0, 0, 1}, {0,1,0}, {1,0,0} };
  int rot_si[3][3] = { {1, 0, 0}, {1,0,0}, {0,1,0} };
  int rot_sj[3][3] = { {0, 1, 0}, {0,0,1}, {0,0,1} };

  pid_t nvpcur = nvpbuf;

  for (int k = 0; k < 3; k++)
    {

      /* Add guard points */
      for (int i=0; i <= NGUARD; i++)
        {
          for (int j=0; j <= NGUARD; j++)
            {
              assert(nvpcur < nvpall);
              /* Bottom */
              for (int a = 0; a < 3; a++)
                parts[3*nvpcur+a] = -width2[a] + (realT)i * s[a] * rot_si[k][a] + (realT)j * s[a] * rot_sj[k][a] - rot_g[k][a] * g[a];

              nvpcur++;

              assert(nvpcur < nvpall);
              /* Top */
              for (int a = 0; a < 3; a++)
                parts[3*nvpcur+a] = (2*rot_g[k][a]-1)*width2[a] + (realT)i * s[a] * rot_si[k][a] + (realT)j * s[a] * rot_sj[k][a] + rot_g[k][a] * g[a];

              nvpcur++;
            }
        }
    }
}

void BoxData::findBoundary()
{
  double BF = std::numeric_limits<double>::max();

  for (int j = 0; j < 3; j++)
    {
      xyz_min[j] = BF;
      xyz_max[j] = -BF;
    }

  for (pid_t i = nvpbuf; i < nvpall; i++) {
    for (int j = 0; j < 3; j++)
      {
        xyz_min[j] = std::min(xyz_min[j], parts[3*i+j]);
        xyz_max[j] = std::max(xyz_max[j], parts[3*i+j]);
      }
  }
}

void saveTesselation(const string& outfile, PositionData& pdata, BoxData& boxdata, PARTADJ *adjs, float *vols)
{
  ofstream out(outfile.c_str());
  if (!out)
    {
      cout << format("Unable to open %s") % outfile << endl;
      exit(0);
    }
  out.write((char *)&pdata.np, sizeof(pid_t));
  out.write((char *)&boxdata.nvp, sizeof(pid_t));
  cout << format("nvp = %d") % boxdata.nvp << endl;

  /* Tell us where the original particles were */
  out.write((char *)boxdata.orig, sizeof(pid_t)*boxdata.nvp);
  /* Volumes*/
  out.write((char *)vols,sizeof(float)*boxdata.nvp);
  /* Adjacencies */
  for (pid_t i = 0; i < boxdata.nvp; i++)
    {
      out.write((char*)&adjs[i].nadj, sizeof(pid_t));
      if (adjs[i].nadj > 0)
        out.write((char *)adjs[i].adj, adjs[i].nadj*sizeof(pid_t));
      else
        (cout << "0").flush();
    }
  out.close();
}

int main(int argc, char *argv[]) {
  PositionData pdata;
  BoxData boxdata;

  coordT deladjs[3*MAXVERVER], points[3*MAXVERVER];
  pointT intpoints[3*MAXVERVER];
  string outfile;
  ofstream out;

  char *suffix, *outDir, *posfile;
  PARTADJ *adjs;
  float *vols;
  double border, boxsize[3];
  int b[3], numdiv[3];
  double totalvol;

  CosmoTool::MiniArgDesc args[] = {
    { "POSITION FILE", &posfile, MINIARG_STRING },
    { "BORDER SIZE", &border, MINIARG_DOUBLE },
    { "BOX_X", &boxsize[0], MINIARG_DOUBLE },
    { "BOX_Y", &boxsize[1], MINIARG_DOUBLE },
    { "BOX_Z", &boxsize[2], MINIARG_DOUBLE },
    { "SUFFIX", &suffix, MINIARG_STRING },
    { "NUM_DIVISION_X", &numdiv[0], MINIARG_INT },
    { "NUM_DIVISION_Y", &numdiv[1], MINIARG_INT },
    { "NUM_DIVISION_Z", &numdiv[2], MINIARG_INT },
    { "B0", &b[0], MINIARG_INT },
    { "B1", &b[1], MINIARG_INT },
    { "B2", &b[2], MINIARG_INT },
    { "OUTPUT DIRECTORY", &outDir, MINIARG_STRING },
    { 0, 0, MINIARG_NULL }
  };

  if (!CosmoTool::parseMiniArgs(argc, argv, args))
    return 1;

  if (!checkParameters(numdiv, b))
    return 2;

  /* Boxsize should be the range in r, yielding a range 0-1 */
  if (!pdata.readFrom(posfile))
    return 3;

  (cout << pdata.np << " particles" << endl).flush();

  pdata.findExtrema();

  (cout << boost::format("np: %d, x: %f,%f; y: %f,%f; z: %f,%f")
    % pdata.np
    % pdata.xyz_min[0] % pdata.xyz_max[0]
    % pdata.xyz_min[1] % pdata.xyz_max[1]
    % pdata.xyz_min[2] % pdata.xyz_max[2]).flush();


  if (border > 0.)
    boxdata.bf = border;
  else
    boxdata.bf = 0.1;

  boxdata.prepareBox(pdata, b, numdiv);
  pdata.destroy();
  boxdata.addGuardPoints();

  adjs = new PARTADJ[boxdata.nvpall];
  if (adjs == 0)
    {
      cout << "Unable to allocate adjs" << endl;
      return 0;
    }

  boxdata.findBoundary();

  cout << format("Added guard points to total %d points (should be %d)")
            % boxdata.nvpall % (boxdata.nvpbuf + 6*(NGUARD+1)*(NGUARD+1)) << endl;
  cout << format("x: %f,%f; y: %f,%f; z:%f,%f") % boxdata.xyz_min[0] % boxdata.xyz_max[0] % boxdata.xyz_min[1] % boxdata.xyz_max[1] % boxdata.xyz_min[2] % boxdata.xyz_max[2] << endl;

  /* Do tesselation*/
  printf("File read.  Tessellating ...\n"); fflush(stdout);
  int exitcode = delaunadj(boxdata.parts, boxdata.nvp, boxdata.nvpbuf, boxdata.nvpall, &adjs);
  if (exitcode != 0)
   {
     printf("Error while tesselating. Stopping here."); fflush(stdout);
     return 4;
   }

  /* Now calculate volumes*/
  printf("Now finding volumes ...\n"); fflush(stdout);
  vols = new float[boxdata.nvp];

  for (pid_t i = 0; i < boxdata.nvp; i++)
    { /* Just the original particles
       * Assign adjacency coordinate array*/
      /* Volumes */
      for (int j = 0; j < adjs[i].nadj; j++)
        {
          for (int d = 0; d < 3; d++)
            {
              deladjs[3*j + d] = boxdata.parts[3*adjs[i].adj[j]+d] - boxdata.parts[3*i+d];

              if (deladjs[3*j+d] < -boxdata.width2[d])
                deladjs[3*j+d] += boxdata.width[d];
              if (deladjs[3*j+d] > boxdata.width2[d])
                deladjs[3*j+d] -= boxdata.width[d];
            }
        }

      exitcode = vorvol(deladjs, points, intpoints, adjs[i].nadj, &(vols[i]));
      vols[i] *= pdata.np/pdata.V0;
      if ((i % 1000) == 0)
        cout << format("%d: %d, %f") % i % adjs[i].nadj % vols[i] << endl;
    }

  /* Get the adjacencies back to their original values */
  for (pid_t i=0; i<boxdata.nvp; i++)
    for (int j=0; j < adjs[i].nadj; j++)
      adjs[i].adj[j] = boxdata.orig[adjs[i].adj[j]];

  totalvol = 0.;
  for (pid_t i=0;i<boxdata.nvp; i++)
    totalvol += vols[i];

  cout << format("Average volume = %g") % (totalvol/boxdata.nvp) << endl;

  /* Now the output!
     First number of particles */
  outfile = str(format("%s/part.%s.%02d.%02d.%02d") % outDir % suffix % b[0] % b[1] % b[2]);

  cout << format("Output to %s") %outfile << endl << endl;
  saveTesselation(outfile, pdata, boxdata, adjs, vols);
  delete[] adjs;

  return(0);
}
