#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

namespace CosmoTool {


  template<typename PType, typename IType, int N>
  void DelaunayInterpolate<PType,IType,N>::buildQuickAccess()
  {
    cells = new QuickCell[numPoints];

    uint32_t point_to_simplex_size = 0;
    uint32_t *numSimplex_by_point = new uint32_t[numPoints];
    uint32_t *index_by_point = new uint32_t[numPoints];

    // First count the number of simplex for each point
    for (uint32_t i = 0; i < numPoints; i++)
      index_by_point[i] = numSimplex_by_point[i] = 0;
    for (uint32_t i = 0; i < (N+1)*numSimplex; i++)
      {
	assert(simplex_list[i] < numPoints);
	if (!disable_simplex[i/(N+1)])
	  numSimplex_by_point[simplex_list[i]]++;
      }

    // Compute the total number and the index for accessing lists.
    for (uint32_t i = 0; i < numPoints; i++)
      {
	index_by_point[i] = point_to_simplex_size;
	point_to_simplex_size += numSimplex_by_point[i]+1;
      }

    // Now compute the real list.
    point_to_simplex_list_base = new int32_t[point_to_simplex_size];
    for (uint32_t i = 0; i < numSimplex; i++)
      {
	for (int j = 0; j <= N; j++)
	  {
	    uint32_t s = (N+1)*i+j;
	    if (disable_simplex[i])
	      continue;

	    uint32_t p = simplex_list[s];
	    assert(index_by_point[p] < point_to_simplex_size);
	    point_to_simplex_list_base[index_by_point[p]] = i;
	    ++index_by_point[p];
	  }
      }

    // Finish the lists
    for (uint32_t i = 0; i < numPoints; i++)
      {
	// check assertion
	assert((i==0 && index_by_point[0]==numSimplex_by_point[0]) 
	       || 
	       ((index_by_point[i]-index_by_point[i-1]) == (numSimplex_by_point[i]+1)));
	assert(index_by_point[i] < point_to_simplex_size);
	point_to_simplex_list_base[index_by_point[i]] = -1;
      }
    
    uint32_t idx = 0;
    for (uint32_t i = 0; i < numPoints; i++)
      {
	cells[i].active = true;
	cells[i].val.simplex_list = &point_to_simplex_list_base[idx];
	// We may have to cast here.
	for (int j = 0; j < N; j++)
	  cells[i].coord[j] = positions[i][j];

	idx += numSimplex_by_point[i]+1;
      }

    // Free the memory allocated for temporary arrays.
    delete[] numSimplex_by_point;
    delete[] index_by_point;

    // Build the kd tree now.
    quickAccess = new QuickTree(cells, numPoints);
  }

  template<typename PType, typename IType, int N>
  void DelaunayInterpolate<PType,IType,N>::buildPreweight()
	throw(InvalidArgumentException)
  {
    double preweight[N*N];
    double preweight_inverse[N*N];
    gsl_permutation *p = gsl_permutation_alloc(N);
    uint32_t numDisabled = 0;

    all_preweight = new PType[N*N*numSimplex];

    for (uint32_t i = 0; i < numSimplex; i++)     
      {
	uint32_t base = i*(N+1);
	uint32_t pref = simplex_list[base];
	// Compute the forward matrix first.
	for (int j = 0; j < N; j++)
	  {
	    PType xref = positions[pref][j];

	    for (int k = 0; k < N; k++)
	      {
		preweight[j*N + k] = positions[simplex_list[k+base+1]][j] - xref;
	      }
	  }

	gsl_matrix_view M = gsl_matrix_view_array(preweight, N, N);
	gsl_matrix_view iM = gsl_matrix_view_array(preweight_inverse, N, N);
	int signum;

	gsl_linalg_LU_decomp(&M.matrix, p, &signum);
        double a = fabs(gsl_linalg_LU_det(&M.matrix, signum));
	if (a < 1e-10)
         {
#ifdef DEBUG
           for (int j = 0; j < N; j++)
            {
              PType xref = positions[pref][j];

              for (int k = 0; k < N; k++)
                {
                  preweight[j*N + k] = positions[simplex_list[k+base+1]][j] - xref;
                }
            }
           std::ofstream f("matrix.txt");
           for (int j = 0; j < N*N; j++)
             f << std::setprecision(12) << preweight[j] << std::endl;
	   throw InvalidArgumentException("Invalid tesselation. One simplex is coplanar.");
#else
           gsl_matrix_set_zero(&iM.matrix);
	   disable_simplex[i] = true;
	   numDisabled++;
#endif
         }
        else {
	  gsl_linalg_LU_invert(&M.matrix, p, &iM.matrix);
	  disable_simplex[i] = false;
        }
 
	for (int j = 0; j < N*N; j++)
	  all_preweight[N*N*i + j] = preweight_inverse[j];
      }

    std::cout << "Number of disabled simplices: " << numDisabled << std::endl;

    gsl_permutation_free(p);
  }

  template<typename PType, typename IType, int N>
  void DelaunayInterpolate<PType,IType,N>::buildHyperplane(const PType *v, CoordType& hyper)
  {
    double M[N][N], eVal[N], eVec[N][N];
    gsl_matrix_view mM, evec;
    gsl_vector_view eval;

    // Construct the symmetric matrix
    for (int k = 0; k < N; k++)
      for (int l = k; l < N; l++)
	{
	  double val = 0;

	  for (int i = 0; i < (N-1); i++)
	    {
	      val += v[i*N+l] * v[i*N+k];
	    }
	  M[l][k] = M[k][l] = val;
	}

    mM = gsl_matrix_view_array(&M[0][0], N, N);
    evec = gsl_matrix_view_array(&eVec[0][0], N, N);
    eval = gsl_vector_view_array(&eVal[0], N);

    // Solve the eigensystem
    gsl_eigen_symmv (&mM.matrix, &eval.vector, &evec.matrix, eigen_work);
    
    double minLambda = INFINITY;
    uint32_t idx = N+1;

    // Look for the smallest eigenvalue
    for (int k = 0; k < N; k++)
      {
	if (minLambda > eVal[k])
	  {
	    minLambda = eVal[k];
	    idx = k;
	  }
      }
    assert(idx != (N+1));

    // Copy the corresponding vector
    for (int k = 0; k < N; k++)
      {
	hyper[k] = eVec[k][idx];
      }
  }

  template<typename PType, typename IType, int N>
  bool DelaunayInterpolate<PType,IType,N>::checkPointInSimplex(const CoordType& pos, uint32_t simplex)
  {
    if (disable_simplex[simplex])
      return false;

    uint32_t *desc_simplex = &simplex_list[simplex*(N+1)];
    CoordType *p[N+1], v[N], hyper;

    for (int k = 0; k <= N; k++)
      p[k] = &positions[desc_simplex[k]];


    for (int i = 0; i <= N; i++)
      {
	// Build vectors 
	for (int k = 1; k <= N; k++)
	  for (int l = 0; l < N; l++)
	    v[k-1][l] = (*p[k])[l] - (*p[0])[l];

	// Build hyperplane.
	buildHyperplane(&v[0][0], hyper);
       
	// Compute the appropriate sign using the last point.
	PType sign = 0;
	for (int k = 0; k < N; k++)
	  sign += hyper[k] * v[N-1][k];

	// Now check the point has the same sign;
	PType pnt_sign = 0;
	for (int k = 0; k < N; k++)
	  pnt_sign += hyper[k] * (pos[k] - (*p[0])[k]);

	if (pnt_sign*sign < 0)
	  return false;


	// Rotate the points.
	for (int k = 1; k <= N; k++)
	  {
	    p[k-1] = p[k];
	  }
	p[N] = &positions[desc_simplex[i]];
      }

    // We checked all possibilities. Return now.
    return true;
  }


  template<typename PType, typename IType, int N>
  uint32_t DelaunayInterpolate<PType,IType,N>::findSimplex(const CoordType& c)
    throw (InvalidArgumentException)
  {
    uint32_t N_ngb = 1;
    QuickCell **cell_Ngb = new QuickCell *[N_ngb];
    typename QuickTree::coords kdc;

    for (int i = 0; i < N; i++)
       kdc[i] = c[i];

    // It may happen that we are unlucky and have to iterate to farther
    // neighbors. It is bound to happen, especially on the boundaries.
    do
      {
	uint32_t i;

	quickAccess->getNearestNeighbours(kdc, N_ngb, cell_Ngb);
	
	for (i = 0; i < N_ngb && cell_Ngb[i] != 0; i++)
	  {
	    int32_t *simplex_list = cell_Ngb[i]->val.simplex_list;
	    uint32_t j = 0;

	    while (simplex_list[j] >= 0)
	      {
		if (checkPointInSimplex(c, simplex_list[j]))
		  {
		    delete[] cell_Ngb;
		    return simplex_list[j];
		  }
	        ++j;
	      }
	  }	
	delete[] cell_Ngb;

	// The point does not belong to any simplex.
	if (i != N_ngb)
	  throw InvalidArgumentException("the given point does not belong to any simplex");

	N_ngb *= 2;
	cell_Ngb = new QuickCell *[N_ngb];
      }
    while (1);

    // Point not reached.
    abort();
    return 0;
  }

  template<typename PType, typename IType, int N>
  IType DelaunayInterpolate<PType,IType,N>::computeValue(const CoordType& c)
    throw (InvalidArgumentException)
  {
    uint32_t simplex = findSimplex(c);
    PType *preweight = &all_preweight[simplex*N*N];
    PType weight[N+1];
    PType p0[N];
    PType sum_weight = 0;

    for (int i = 0; i < N; i++)
      p0[i] = positions[simplex_list[simplex*(N+1) + 0]][i];

    // Now we use the preweight to compute the weight...
    for (int i = 1; i <= N; i++)
      {
	weight[i] = 0;
	for (int j = 0; j < N; j++)
	  weight[i] += preweight[(i-1)*N+j]*(c[j]-p0[j]);

	assert(weight[i] > -1e-7);
	assert(weight[i] < 1+1e-7);
	sum_weight += weight[i];
      }   
    weight[0] = 1-sum_weight;
    assert(weight[0] > -1e-7);
    assert(weight[0] < (1+1e-7));
    
    // We compute the final value by weighing the value at the N+1
    // points by the proper weight.
    IType final = 0;
    for (int i = 0; i <= N; i++)
      final += weight[i] * values[ simplex_list[simplex*(N+1) + i] ];

    return final;
  }
  
};
