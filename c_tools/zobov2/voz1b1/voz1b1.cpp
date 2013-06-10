#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include <CosmoTool/miniargs.hpp>
#include <cmath>
#include "libqhull/qhull_a.h"
#include "voz.h"
#include "voz_io.hpp"

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

  void prepareBox(const PositionData& pdata);

  void checkParticle(float *xyz, bool& in_main, bool& in_buf);
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
}

void BoxData::prepareBox(const PositionData& pdata, int *b_id)
{
  guard_added = false;

  for (int i = 0; i < 3; i++)
    {
      float s;
      width[i] = (pdata.xyz_max[i] - pdata.xyz_min[i])/(float)numdiv;
      width2[i] = 0.5*width;

      totwidth[i] = width+2.*bf;
      totwidth2[i] = width2 + bf;

      s[i] = width[i]/(float)NGUARD;
      if ((bf*bf - 2.*s[i]*s[i]) < 0.)
        {
          printf("bf = %f, s = %f.\n",bf,s[i]);
          printf("Not enough guard points for given border.\nIncrease guards to >= %f\n.",
                 sqrt(2.)*width/bf);
          exit(0);
        }
      g[i] = (bf / 2.)*(1. + sqrt(1 - 2.*s[i]*s[i]/(bf*bf)));
      cout << format("s[%d] = %f, bf = %f, g[%d] = %f.") % s[i] % i % bf % g[i] % i << endl;

      c[i] = b_id[i] * width[i];
    }
  cout.flush();

  cout << format("c: %f,%f,%f") % c[0] % c[1] % c[2] << endl;

  /* Assign temporary array*/
  nvpbuf = 0; /* Number of particles to tesselate, including buffer */
  nvp = 0; /* Without the buffer */

  for (pid_t i=0; i < np; i++)
    {
      double xyz[3] = { pdata.xyz[0][i], pdata.xyz[1][i], pdata.xyz[2][i] };
      bool is_it_in_buf, is_it_in_main;

      checkParticle(xyz, is_it_in_buf, is_it_in_main);

      if (isitinbuf)
        nvpbuf++;
      if (isitinmain)
        nvp++;
    }

  nvpbuf += 6*(NGUARD+1)*(NGUARD+1); /* number of guard points */

  parts = new coordT[3*nvpbuf];
  orig = new pid_t[nvpuf];

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

  nvp = 0; nvpall = 0; /* nvp = number of particles without buffer */

  for (int j = 0; j < 3; j++)
    {
      xyz_min[j] = -std::numeric_limits<double>::max();
      xyz_max[j] = std::numeric_limits<double>::max();
    }
  for (pid_t i = 0; i < np; i++)
    {
      bool is_it_in_main, is_it_in_buf;
      double xyz[3] = { pdata.xyz[0][i], pdata.xyz[1][i], pdata.xyz[2][i] };

      checkParticle(xyz, is_it_in_main, is_it_in_buf);

      if (is_it_in_main) {
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
  printf("nvp = %d\n",nvp);
  printf("x: %f,%f; y: %f,%f; z:%f,%f\n",xmin,xmax,ymin,ymax,zmin,zmax);
  nvpbuf = nvp;
  for (pid_t i = 0; i < np; i++)
    {
      bool is_it_in_main, is_it_in_buf;
      double xyz[3] = { pdata.xyz[0][i], pdata.xyz[1][i], pdata.xyz[2][i] };

      checkParticle(xyz, is_it_in_main, is_it_in_buf);

      if (isitinbuf && !isitinmain)
        {
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
  nvpall = nvpbuf;

  double predict = 1;
  for (int j = 0; j < 3; j++)
    predict *= (width[j]+2*bf);
  predict *= np;
  cout << format("There should be ~ %g points; there are %d\n") % predict % nvpbuf << endl;
}

void BoxData::addGuardPoints()
{
  int rot_g[3][3]  = { {0, 0, 1}, {0,1,0}, {1,0,0} };
  int rot_si[3][3] = { {1, 0, 0}, {1,0,0}, {0,1,0} };
  int rot_sj[3][3] = { {0, 1, 0}, {0,0,1}, {0,0,1} };

  for (int k = 0; k < 3; k++)
    {

      /* Add guard points */
      for (i=0; i<=NGUARD; i++)
        {
          for (j=0; j<=NGUARD; j++)
            {
              /* Bottom */
              for (int a = 0; a < 3; a++)
                parts[3*nvpall+a] = -width2[a] + (realT)i * s[a] * rot_si[k][a] + (realT)j * s[a] * rot_sj[k][a] - rot_g[k][a] * g[a];

              nvpall++;

              /* Top */
              for (int a = 0; a < 3; a++)
                parts[3*nvpall+a] = (2*rot_g[k][a]-1)*width2[a] + (realT)i * s[a] * rot_si[k][a] + (realT)j * s[a] * rot_sj[k][a] + rot_g[k][a] * g[a];

              nvpall++;
            }
        }
    }
}


int main(int argc, char *argv[]) {
  PositionData pdata;
  BoxData boxdata;

  coordT deladjs[3*MAXVERVER], points[3*MAXVERVER];
  pointT intpoints[3*MAXVERVER];
  string outfile;
  ofstream out;

  char *suffix, *outDir;
  PARTADJ *adjs;
  float *vols;
  pid_t *orig;
  double border, boxsize[3];
  int b[3], numdiv[3];
  double totalvol;

  CosmoTool::MiniArgDesc args[] = {
    { "POSITION FILE", MINIARG_STRING, &posfile },
    { "BORDER SIZE", MINIARG_DOUBLE, &border },
    { "BOX_X", MINIARG_DOUBLE, &boxsize[0] },
    { "BOX_Y", MINIARG_DOUBLE, &boxsize[1] },
    { "BOX_Z", MINIARG_DOUBLE, &boxsize[2] },
    { "SUFFIX", MINIARG_STRING, &suffix },
    { "NUM_DIVISION_X", MINIARG_INT, &numdiv[0] },
    { "NUM_DIVISION_Y", MINIARG_INT, &numdiv[1] },
    { "NUM_DIVISION_Z", MINIARG_INT, &numdiv[2] },
    { "B0", MINIARG_INT, &b[0] },
    { "B1", MINIARG_INT, &b[1] },
    { "B2", MINIARG_INT, &b[2] },
    { "OUTPUT DIRECTORY", MINIARG_STRING, &outDir },
    { 0, MINIARG_NULL, 0 }
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
    % xyz_min[0] % xyz_max[0]
    % xyz_min[1] % xyz_max[1]
    % xyz_min[2] % xyz_max[2]).flush();


  if (border > 0.)
    boxdata.bf = border;
  else
    boxdata.bf = 0.1;

  boxdata.prepareBox(pdata);
  pdata.destroy();
  boxdata.addGuardPoints();

  adjs = new PARTADJ[np];
  if (adjs == 0)
    {
      cout << "Unable to allocate adjs" << endl;
      return 0;
    }

  xmin = BF; xmax = -BF; ymin = BF; ymax = -BF; zmin = BF; zmax = -BF;
  for (i=nvpbuf;i<nvpall;i++) {
    if (parts[3*i] < xmin) xmin = parts[3*i];
    if (parts[3*i] > xmax) xmax = parts[3*i];
    if (parts[3*i+1] < ymin) ymin = parts[3*i+1];
    if (parts[3*i+1] > ymax) ymax = parts[3*i+1];
    if (parts[3*i+2] < zmin) zmin = parts[3*i+2];
    if (parts[3*i+2] > zmax) zmax = parts[3*i+2];
  }

  cout << format("Added guard points to total %d points (should be %d)")
            % boxdata.nvpall % (boxdata.nvpbuf + 6*(NGUARD+1)*(NGUARD+1)) << endl;
  cout << format("x: %f,%f; y: %f,%f; z:%f,%f") % boxdata.xyz_min[0] % boxdata.xyz_max[0] % boxdata.xyz_min[1] % boxdata.xyz_max[1] % boxdata.xyz_min[2] % boxdata.xyz_max[2] << endl;

  /* Do tesselation*/
  printf("File read.  Tessellating ...\n"); fflush(stdout);
  exitcode = delaunadj(boxdata.parts, boxdata.nvp, boxdata.nvpbuf, boxdata.nvpall, &adjs);
  if (exitcode != 0)
   {
     printf("Error while tesselating. Stopping here."); fflush(stdout);
     return 4;
   }

  /* Now calculate volumes*/
  printf("Now finding volumes ...\n"); fflush(stdout);
  vols = new float[nvp];

  for (pid_t i = 0; i < nvp; i++)
    { /* Just the original particles
       * Assign adjacency coordinate array*/
      /* Volumes */
      for (int j = 0; j < adjs[i].nadj; j++)
        {
          for (int d = 0; d < 3; d++)
            {
              deladjs[3*j + d] = parts[3*adjs[i].adj[j]+d] - parts[3*i+d];

              if (deladjs[3*j+d] < -boxdata.width2[d])
                deladjs[3*j+d] += boxdata.width[d];
              if (deladjs[3*j+d] > boxdata.width2[d])
                deladjs[3*j+d] -= boxdata.width[d];
            }
        }

      exitcode = vorvol(deladjs, points, intpoints, adjs[i].nadj, &(vols[i]));
      vols[i] *= np/V0;
      if ((i % 1000) == 0)
        cout << format("%d: %d, %f") % i % adjs[i].nadj % vols[i] << endl;
    }

  /* Get the adjacencies back to their original values */
  for (pid_t i=0; i<nvp; i++)
    for (int j=0; j < adjs[i].nadj; j++)
      adjs[i].adj[j] = orig[adjs[i].adj[j]];

  totalvol = 0.;
  for (pid_t i=0;i<nvp; i++)
    totalvol += vols[i];

  cout << format("Average volume = %g") % (totalvol/nvp) << endl;

  /* Now the output!
     First number of particles */
  outfile = str(format("%s/part.%s.%02d.%02d.%02d") % outDir % suffix % b[0] % b[1] % b[2]);

  cout << format("Output to %s") %outfile << endl << endl;
  out.open(outfile.c_str());
  if (!out)
    {
      cout << format("Unable to open %s") % outfile << endl;
      return 0;
    }
  out.write((char *)&pdata.np, sizeof(pid_t));
  out.write((char *)&boxdata.nvp, sizeof(pid_t));
  cout << format("nvp = %d") % nvp << endl;

  /* Tell us where the original particles were */
  out.write((char *)boxdata.orig, sizeof(pid_t)*boxdata.nvp);
  /* Volumes*/
  out.write((char *)vols,sizeof(float)*nvp);
  /* Adjacencies */
  for (pid_t i = 0; i < nvp; i++)
    {
      out.write((char*)&adjs[i].nadj, sizeof(pid_t));
      if (adjs[i].nadj > 0)
        out.write((char *)adjs[i].adj, adjs[i].nadj*sizeof(pid_t));
      else
        (cout << "0").flush();
    }
  out.close();
  delete[] adjs;

  return(0);
}
