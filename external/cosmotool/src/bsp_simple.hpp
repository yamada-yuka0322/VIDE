/*+
This is CosmoTool (./src/bsp_simple.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMOTOOL_SIMPLE_BSP_HPP
#define __COSMOTOOL_SIMPLE_BSP_HPP

#include <stdlib.h>
#include <cmath>
#include "algo.hpp"
#include <queue>
#include <exception>

namespace CosmoTool
{
  
  namespace simple_bsp
  {

    template<typename T, typename PType, int N>
    struct space
    {
      typedef T data_t;
      typedef PType point_t;
      typedef PType coord_t[N];
      static const int dim = N;
    };

    template<typename SType>
    struct Plane
    { 
      typename SType::coord_t n;
      typename SType::point_t d;
    };

    template<typename SType>
    typename SType::point_t dot_product(const typename SType::coord_t& c1,
				      const typename SType::coord_t& c2)
    {
      typename SType::point_t A = 0;

      for(int j = 0; j < SType::dim; j++)
	A += c1[j]*c2[j];
      return A;
    }

    template<typename SType>
    struct Node {
      Plane<SType> plane;
      Node<SType> *minus, *plus;
      typename SType::data_t data;
    };

    template<typename SType>
    void normal2(typename SType::coord_t p[2], typename SType::coord_t& n)
    {
      typename SType::point_t d;
      using CosmoTool::square;

      n[0] = p[1][1]-p[0][1];
      n[1] = -p[1][0]+p[0][0];
      d = std::sqrt(square(n[0])+square(n[1]));
      n[0] /= d;
      n[1] /= d;
    }
    
    template<typename SType>
    void normal3(typename SType::coord_t p[3], typename SType::coord_t& n)
    {
      typename SType::point_t delta0[3] = { p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2] };
      typename SType::point_t delta1[3] = { p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2] };
      typename SType::point_t d;
      using CosmoTool::square;

      n[0] = delta0[1] * delta1[2] - delta0[2] * delta1[1];
      n[1] = delta0[2] * delta1[0] - delta0[0] * delta1[2];
      n[2] = delta0[0] * delta1[1] - delta0[1] * delta1[0];
      d = std::sqrt(square(n[0]) + square(n[1]) + square(n[2]));
      n[0] /= d;
      n[1] /= d;
      n[2] /= d;
    }
    


    template<typename SType>
    struct Facet
    {
      typename SType::coord_t p[SType::dim];
      typename SType::data_t data;

      void center(typename SType::coord_t& c)
      {
	for (int j = 0; j < SType::dim; j++)
	  {
	    c[j] = 0;
	    for (int i = 0; i < SType::dim; i++)
	      {
		c[j] += p[i][j]; 
	      }
	    c[j] /= SType::dim+1;
	  }
      }
      
      void normal(typename SType::coord_t& n)
      {
	if (SType::dim==2)
	  {
	    normal2<SType>(p, n);
	    return;
	  }
	if (SType::dim == 3)
	  {
	    normal3<SType>(p, n);
	    return;
	  }
	abort();
      }

    };


    class InvalidPoint: public std::exception
    {
    };

    template<typename T, typename PType, int N>
    class BSP
    {
    public:
      typedef space<T, PType, N> space_t;
      typedef Plane<space_t> plane_t;
      typedef Node<space_t> node_t;

      node_t *root;
      std::queue<node_t *> allocated;
      
      BSP() throw();
      ~BSP();
      void insert(Facet<space_t>& facet);

      template<typename PIterator>
      void insert(PIterator b, PIterator e, T data)
      {
	Facet<space_t> f;
	int q = 0;

	while (b != e && q < N+1)
	  {
	    for (int j = 0; j < N; j++)
	      f.p[q][j] = (*b)[j];
	    ++b;
	    ++q;
	  }
	if (q != N)
	  throw InvalidPoint();
	
	f.data = data;
	insert(f);
      }

      bool inside(const typename space_t::coord_t& p) const;
    };

  };

};

#include "bsp_simple.tcc"

#endif
