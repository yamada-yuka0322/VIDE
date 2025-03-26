/*+
This is CosmoTool (./src/replicateGenerator.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
#ifndef __REPLICATE_GENERATOR_HPP
#define __REPLICATE_GENERATOR_HPP

#include <algorithm>
#include "algo.hpp"

namespace CosmoTool
{

  template<typename Coord, int N>
  class ReplicateGenerator
  {
  public:
    typedef Coord Coords[N];
    Coords replicate;

    ReplicateGenerator(const Coords& x, Coord shift)
    {
      face = 0; 
      std::fill(replicate, replicate+N, shift);
      numFaces = spower<N,long>(3);
      std::copy(x, x+N, x_base);
      if (!next())
	abort();
    }

    ReplicateGenerator(const Coords& x, Coords& shift)
    {
      face = 0;
      std::copy(shift, shift+N, replicate);
      numFaces = spower<N,long>(3);
      std::copy(x, x+N, x_base);
      if (!next())
        abort();
    }


    bool next()
    {
      if (face == numFaces)
        return false;

      face++;

      bool no_move = true;
      int q_face = face;
      for (int i = 0; i < N; i++)
        {
          int c_face;
          c_face = q_face % 3;
          q_face /= 3;
          x_shifted[i] = x_base[i] + (c_face-1)*replicate[i];
          no_move = no_move && (c_face == 1);
        }
      if (no_move)
        return next();
      return true;
    }

    const Coord *getPosition()
    {
      return x_shifted;
    }

    void getPosition(Coords& x_out)
    {
      std::copy(x_shifted, x_shifted+N, x_out);
    }

  private:
    Coord x_shifted[N], x_base[N];
    long face, numFaces;
  };

};

#endif
