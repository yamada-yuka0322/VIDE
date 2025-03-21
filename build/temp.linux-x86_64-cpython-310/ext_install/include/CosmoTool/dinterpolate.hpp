/*+
This is CosmoTool (./src/dinterpolate.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMO_DINTERPOLATE_HPP
#define __COSMO_DINTERPOLATE_HPP

#include "config.hpp"
#include "mykdtree.hpp"
#include "kdtree_splitters.hpp"
#include <gsl/gsl_eigen.h>

namespace CosmoTool {
 
  template<typename PType, typename IType, int N>
  class DelaunayInterpolate {
  public:
    struct SimplexAccess {
      int32_t *simplex_list;      
    };
 
    typedef KDTree<N, SimplexAccess, PType, KD_homogeneous_cell_splitter<N, SimplexAccess> > QuickTree;
    typedef typename QuickTree::Cell QuickCell;
    typedef PType CoordType[N];
    
    QuickTree *quickAccess;
    QuickCell *cells;
    PType *all_preweight;
    int32_t *point_to_simplex_list_base;    
    IType *values;
    CoordType *positions;    
    uint32_t numPoints;
    uint32_t numSimplex;
    uint32_t *simplex_list;
    gsl_eigen_symmv_workspace *eigen_work;
    bool *disable_simplex;
    
    /**
     * This construct the interpolator. The construction is time consuming so
     * please do it the less possible number of times, especially if you have
     * a large number of points.
     *
     * @param positions list of the positions
     * @param values list of the values taken at each position
     * @param simplex_list list of points for each simplex. The packing
     * is the following:
     * [t(0,1),t(0,2),...,t(0,n+1),t(1,0),t(1,1),...,t(1,n+1),..],
     * with t(i,j) the i-th simplex and j-th point of the simplex. The indexes
     * refer to the previous list of points.
     * @param numPoints the number of points
     */
    DelaunayInterpolate(CoordType *positions, IType *values, uint32_t *simplex_list,
			uint32_t numPoints, uint32_t numSimplex)
	throw (InvalidArgumentException)
    {
      this->positions = positions;
      this->values = values;
      this->simplex_list = simplex_list;
      this->numPoints = numPoints;
      this->numSimplex = numSimplex;
      this->disable_simplex = new bool[numSimplex];
      
      buildPreweight();
      buildQuickAccess();

      eigen_work = gsl_eigen_symmv_alloc(N);

    }
    
    ~DelaunayInterpolate()
    {
      delete[] cells;
      delete quickAccess;
      delete[] point_to_simplex_list_base;      
      delete[] all_preweight;
      delete[] disable_simplex;

      gsl_eigen_symmv_free(eigen_work);
    }

    void buildPreweight()
	throw (InvalidArgumentException);
    void buildQuickAccess();
    void buildHyperplane(const PType *v, CoordType& hyper);

    bool checkPointInSimplex(const CoordType& pos, uint32_t simplex);
    
    uint32_t findSimplex(const CoordType& pos)
      throw (InvalidArgumentException);

    IType computeValue(const CoordType& pos)
      throw (InvalidArgumentException);
  };
};

#include "dinterpolate.tcc"

#endif
