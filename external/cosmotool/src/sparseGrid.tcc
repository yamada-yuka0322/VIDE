#include <cassert>
#include <cstring>

namespace CosmoTool {
   
   template<typename T, int N, typename CType,
	    typename CAccessor>
   SparseGrid<T,N,CType,CAccessor>::SparseGrid(CType L, uint32_t gridRes)
   {
     uint32_t Ntot = 1;
     for (int i = 0; i < N; i++)
       Ntot *= gridRes;

     grid.grid = new GridElement<T>[Ntot];
     grid.boxSize = L;
     grid.gridRes = gridRes;
   }

   template<typename T, int N, typename CType,
	    typename CAccessor>
   SparseGrid<T,N,CType,CAccessor>::~SparseGrid()
   {
     delete[] grid.grid;
   }
    
   template<typename T, int N, typename CType,
	    typename CAccessor>
   void SparseGrid<T,N,CType,CAccessor>::addElementInGrid(T& e)
   {
     uint32_t index[N];
     uint32_t idx = 0;

     for (int i = N-1; i >= 0; i--)
       {
	 CType x = acc(e, i) * grid.gridRes / grid.boxSize;
	 
	 index[i] = (int32_t)std::floor((double)x);
	 assert(index[i] >= 0);
	 assert(index[i] < grid.gridRes);
	 idx = (idx*grid.gridRes) + index[i];
       }
     assert(idx < grid.gridRes*grid.gridRes*grid.gridRes);

     grid.grid[idx].eList.push_back(e);
   }

   template<typename T, int N, typename CType,
	    typename CAccessor>
   SparseIterator<T,N,CType>
      SparseGrid<T,N,CType,CAccessor>::locateElements(coord& c, CType R)
   {
     return SparseIterator<T,N,CType>(grid, c, R);
   }

  template<typename T, int N, typename CType,
	   typename CAccessor>
  void SparseGrid<T,N,CType,CAccessor>::clear()
  {
    uint32_t N3 = grid.gridRes*grid.gridRes*grid.gridRes;

    for (uint32_t i = 0; i < N3; i++)
      grid.grid[i].eList.clear();
  }


  template<typename T, int N, typename CType>
  SparseIterator<T,N,CType>::SparseIterator(GridDesc<T>& d, coord& c, CType R)
  {
    desc = d;
    
    for (uint32_t i = 0; i < N; i++)
      {
	CType x_min = (c[i] - R) * d.gridRes / d.boxSize;
	CType x_max = (c[i] + R) * d.gridRes / d.boxSize;
	
	index_min[i] = (int32_t) std::floor((double)x_min);
	index_max[i] = (int32_t) std::ceil((double)x_max) + 1;
	index_cur[i] = index_min[i];
      }

    eCur = getElement();
    eIterator = eCur->eList.begin();
    if (eIterator == eCur->eList.end())
      operator++();
  }
  
  template<typename T, int N, typename CType>
  SparseIterator<T,N,CType>::~SparseIterator()
  {
  }
  
  template<typename T, int N, typename CType>
  bool SparseIterator<T,N,CType>::last()
  {
    return (index_cur[N-1] == index_max[N-1]);
  }

 
  template<typename T, int N, typename CType>
  GridElement<T> *SparseIterator<T,N,CType>::getElement()
  {
    uint32_t idx = 0;
    for (int i = (N-1); i >= 0; i--)
      {
	int32_t k = (index_cur[i] + desc.gridRes) % desc.gridRes;	
	idx = (idx*desc.gridRes) + k;
      }
    assert(idx < desc.gridRes*desc.gridRes*desc.gridRes);
    
    return &desc.grid[idx];
  }

  template<typename T, int N, typename CType>
  SparseIterator<T,N,CType>& SparseIterator<T,N,CType>::operator++()
  {
    if (last())
      return *this;
    ++eIterator;
    if (eIterator != eCur->eList.end())
      return *this;
      
    do
      {		
	index_cur[0]++;
	for (int i = 0; i < (N-1); i++)
	  if (index_cur[i] == index_max[i])
	    {
	      index_cur[i] = index_min[i];
	      index_cur[i+1]++;
	    }
	
	if (last())
	  return *this;
	
	eCur = getElement();
	eIterator = eCur->eList.begin();
      }
    while (eIterator == eCur->eList.end());

    return *this;
  }
 


};
