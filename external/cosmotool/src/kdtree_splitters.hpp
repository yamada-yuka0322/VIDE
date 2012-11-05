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


    void check_splitting(KDCell<N,ValType,CType> **cells, uint32_t Ncells, int axis, uint32_t split_index, ctype midCoord)
    {
      ctype delta = std::numeric_limits<ctype>::max();
      assert(split_index < Ncells);
      assert(axis < N);
      for (uint32_t i = 0; i < split_index; i++)
	{
	  assert(cells[i]->coord[axis] <= midCoord);
	  delta = min(midCoord-cells[i]->coord[axis], delta);
	}
      for (uint32_t i = split_index+1; i < Ncells; i++)
	{
	  assert(cells[i]->coord[axis] > midCoord);
	  delta = min(cells[i]->coord[axis]-midCoord, delta);
	}
      assert(delta >= 0);
      assert (std::abs(cells[split_index]->coord[axis]-midCoord) <= delta);
    } 
    
    void operator()(KDCell<N,ValType,CType> **cells, uint32_t Ncells, uint32_t& split_index, int axis, coords minBound, coords maxBound)
    {
      if (Ncells == 1)
	{
	  split_index = 0;
	  return;
	}

      ctype midCoord = 0.5*(maxBound[axis]+minBound[axis]);
      uint32_t below = 0, above = Ncells-1;
      ctype delta_min = std::numeric_limits<ctype>::max();
      uint32_t idx_min = std::numeric_limits<uint32_t>::max();     

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
