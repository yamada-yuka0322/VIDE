/*+
This is CosmoTool (./src/sparseGrid.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
