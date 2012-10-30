#ifndef __VOID_SPARSE_GRID_HPP
#define __VOID_SPARSE_GRID_HPP

#include <cmath>
#include "config.hpp"
#include <list>

namespace CosmoTool 
{
  
  template <typename T>
  struct bracketAccessor
  {
    typedef typename T::value_type result_type;

    result_type
    operator()(T const& V, size_t const N) const throw ()
    {
      return V[N];
    }
  };

  template<typename T>
  struct GridElement
  {
    std::list<T> eList;
  };

  template<typename T, typename CType = double>
  struct GridDesc
  {
    GridElement<T> *grid;
    uint32_t gridRes;
    CType boxSize;
  };

  template<typename T, int N, typename CType = double>
  class SparseIterator
  {
  protected:
    typedef CType coord[N];

    int32_t index_min[N], index_max[N], index_cur[N];
    GridElement<T> *eCur;
    GridDesc<T> desc;
    typename std::list<T>::iterator eIterator;

    GridElement<T> *getElement();
  public:
    SparseIterator(GridDesc<T>& d, coord& c, CType R);
    SparseIterator();
    ~SparseIterator();
    
    SparseIterator& operator++();

    bool last();

    T& operator*() { return *eIterator; }
    T *operator->() { return eIterator.operator->(); }
  };

  template<typename T, int N, typename CType = double,
	   typename CAccessor = bracketAccessor<T> >
  class SparseGrid
  {
  public:
    typedef CType coord[N];
    typedef SparseIterator<T, N, CType> iterator;
    CAccessor acc;

    SparseGrid(CType L, uint32_t gridRes);
    ~SparseGrid();
    
    void addElementInGrid(T& e);

    iterator locateElements(coord& c, CType R);
    void clear();

  protected:
    GridDesc<T> grid;
  };

};

#include "sparseGrid.tcc"

#endif
