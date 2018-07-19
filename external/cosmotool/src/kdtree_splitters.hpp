/*+
This is CosmoTool (./src/kdtree_splitters.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __KDTREE_SPLITTERS_HPP
#define __KDTREE_SPLITTERS_HPP

#include <algorithm>

namespace CosmoTool
{

  template<int N, typename ValType, typename CType = ComputePrecision> 
  struct KD_homogeneous_cell_splitter
  {
    typedef typename KDDef<N,CType>::KDCoordinates coords; 
    typedef typename KDDef<N,CType>::CoordType ctype; 


    void check_splitting(KDCell<N,ValType,CType> **cells, NodeIntType Ncells, int axis, NodeIntType split_index, ctype midCoord)
    {
      ctype delta = std::numeric_limits<ctype>::max();
      assert(split_index < Ncells);
      assert(axis < N);
      for (NodeIntType i = 0; i < split_index; i++)
	{
	  assert(cells[i]->coord[axis] <= midCoord);
	  delta = min(midCoord-cells[i]->coord[axis], delta);
	}
      for (NodeIntType i = split_index+1; i < Ncells; i++)
	{
	  assert(cells[i]->coord[axis] > midCoord);
	  delta = min(cells[i]->coord[axis]-midCoord, delta);
	}
      assert(delta >= 0);
      assert (std::abs(cells[split_index]->coord[axis]-midCoord) <= delta);
    } 
    
    void operator()(KDCell<N,ValType,CType> **cells, NodeIntType Ncells, NodeIntType& split_index, int axis, coords minBound, coords maxBound)
    {
      if (Ncells == 1)
	{
	  split_index = 0;
	  return;
	}

      ctype midCoord = 0.5*(maxBound[axis]+minBound[axis]);
      NodeIntType below = 0, above = Ncells-1;
      ctype delta_min = std::numeric_limits<ctype>::max();
      NodeIntType idx_min = std::numeric_limits<NodeIntType>::max();

      while (below < above)
	{
          ctype delta = cells[below]->coord[axis]-midCoord;
	  if (delta > 0)
            {
              if (delta < delta_min)
                {
                  delta_min = delta;
                  idx_min = above;
                }
	      std::swap(cells[below], cells[above--]);
            }
	  else
            {
              if (-delta < delta_min)
                {
                  delta_min = -delta;
                  idx_min = below;
                }
	      below++;
            }
	}
      // Last iteration
      {
	ctype delta = cells[below]->coord[axis]-midCoord;	
	if (delta > 0)
	  {
	    if (delta < delta_min)
	      {
		delta_min = delta;
		idx_min = above;
	      }
	  }
	else
	  {
	    if (-delta < delta_min)
	      {
		delta_min = -delta;
		idx_min = above;
	      }
	  }
      }

      if (idx_min != above)
        {
          bool cond1 = cells[idx_min]->coord[axis] > midCoord;
          bool cond2 = cells[above]->coord[axis] > midCoord;
          if ((cond1 && cond2) || (!cond1 && !cond2))
            {
              split_index = above;
              std::swap(cells[above], cells[idx_min]);
            }
          else if (cond2)
            {
              if (above >= 1)
		{
		  split_index = above-1;
		  std::swap(cells[above-1], cells[idx_min]);
		}
	      else
		split_index = 0;
	      assert(split_index >= 0);
            }
          else
            {
	      if (above+1 < Ncells)
		{
		  split_index = above+1;
		  std::swap(cells[above+1], cells[idx_min]);
		}
	      else
		split_index = Ncells-1;

	      assert(split_index < Ncells);
            }
         }
       else split_index = above;


      //      check_splitting(cells, Ncells, axis, split_index, midCoord);
    }
  };


};

#endif
