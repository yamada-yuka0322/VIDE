/*+
This is CosmoTool (./src/eskow.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __ESKOW_CHOLESKY_HPP
#define __ESKOW_CHOLESKY_HPP

#include <cmath>
#include <vector>
#include "mach.hpp"

/* Implementation of Schnabel & Eskow, 1999, Vol. 9, No. 4, pp. 1135-148, SIAM J. OPTIM. */

template<typename T, typename A>
class CholeskyEskow
{
private:
  static const bool verbose_eskow = true;
  T tau, tau_bar, mu;

  void print_matrix(A& m, int N)
  {
    using std::cout;
    using std::endl;
    using std::setprecision;

    if (verbose_eskow)
      {

	for (int i = 0; i < N; i++)
	  {
	    for (int j = 0; j < N; j++)
	      {
		cout.width(6);
		cout << setprecision(5) << m(i,j) << " ";
	      }
	    cout << endl;
	  }
	cout << endl;
      }
  }

  T max_diag(A& m, int j, int N)
  {
    T maxval = std::abs(m(j,j));

    for (int k = j+1; k < N; k++)
      {
	maxval = std::max(maxval, std::abs(m(k,k)));
      }
    return maxval;
  }

  void minmax_diag(A& m, int j, int N, T& minval, T& maxval, int& i_min, int& i_max)
  {
    i_min = i_max = j;
    minval = maxval = m(j,j);

    for (int k = j+1; k < N; k++)
      {
	maxval = std::max(maxval, m(k,k));
	minval = std::min(minval, m(k,k));
      }
  
    for (int k = j; k < N; k++)
      {
	if (m(k,k) == minval && i_min < 0)
	  i_min = k;
	if (m(k,k) == maxval && i_max < 0)
	  i_max = k;
      }
  }

  void swap_rows(A& m, int N, int i0, int i1)
  {
    for (int r = 0; r < N; r++)
      std::swap(m(r,i0), m(r,i1));
  }

  void swap_cols(A& m, int N, int i0, int i1)
  {
    for (int c = 0; c < N; c++)
      std::swap(m(i0,c), m(i1,c));
  }

  T square(T x)
  {
    return x*x;
  }

  T min_row(A& m, int j, int N)
  {
    T a = 1/m(j,j);
    T v = m(j+1,j+1) - square(m(j+1,j))*a;
    
    for (int i = j+2; i < N; i++)
      {
	v = std::min(v, m(i, i) - square(m(i,j))*a);
      }

    return v;
  }

  int g_max(const std::vector<T>& g, int j, int N)
  {
    T a = g[j];
    int k = j;

    for (int i = j+1; i < N; i++)
      {
	if (a < g[i])
	  {
	    a = g[i];
	    k = i;
	  }
      }
    return k;
  }

public:
  CholeskyEskow()
  {
    tau = std::pow(mach_epsilon<T>(), 1./3);
    tau_bar = std::pow(mach_epsilon<T>(), 2./3);
    mu=0.1;
  }

  void cholesky(A& m, int N, T& norm_E)
  {
    bool phaseone = true;
    T gamma = max_diag(m, 0, N);
    int j;

    norm_E = 0;
    
    for (j = 0; j < N && phaseone; j++)
      {
	T minval, maxval;
	int i_min, i_max;
	
	print_matrix(m, N);
     
	minmax_diag(m, j, N, minval, maxval, i_min, i_max);
	if (maxval < tau_bar*gamma || minval < -mu*maxval)
	  {
	    phaseone = false;
	    break;
	  }
      
	if (i_max != j)	  
	  {
	    std::cout << "Have to swap i=" << i_max << " and j=" << j << std::endl;
	    swap_cols(m, N, i_max, j);
	    swap_rows(m, N, i_max, j);

	  }
      
	if (min_row(m, j, N) < -mu*gamma)
	  {
	    phaseone = false;
	    break;
	  }
      
	T L_jj = std::sqrt(m(j,j));

	m(j,j) = L_jj;
	for (int i = j+1; i < N; i++)
	  {
	    m(i,j) /= L_jj;
	    for (int k = j+1; k <= i; k++)
	      m(i,k) -= m(i,j)*m(k,j);
	  }
      }


    if (!phaseone && j == N-1)
      {
	T A_nn = m(N-1,N-1);
	T delta = -A_nn + std::max(tau*(-A_nn)/(1-tau), tau_bar*gamma);
      
	m(N-1,N-1) = std::sqrt(m(N-1,N-1) + delta);	  	  
      }

    

    if (!phaseone && j < (N-1))
      {
	std::cout << "Phase two ! (j=" << j << ")" << std::endl;

	int k = j-1;
	std::vector<T> g(N);

	for (int i = k+1; i < N; i++)
	  {
	    g[i] = m(i,i);
	    for (int j = k+1; j < i; j++)
	      g[i] -= std::abs(m(i,j));
	    for (int j = i+1; j < N; j++)
	      g[i] -= std::abs(m(j,i));
	  }

	T delta, delta_prev = 0;
	
	for (int j = k+1; j < N-2; j++)
	  {
	    int i = g_max(g, j, N);
	    T norm_j;

	    print_matrix(m, N);

	    if (i != j)
	      {
		swap_cols(m, N, i, j);
		swap_rows(m, N, i, j);		
	      }


	    for (int i = j+1; j < N; j++)
	      {
		norm_j += std::abs(m(i,j));
	      }
	    
	    delta = std::max(delta_prev, std::max((T)0, -m(j,j) + std::max(norm_j,tau_bar*gamma)));
	    if (delta > 0)
	      {
		m(j,j) += delta;
		delta_prev = delta;
	      }
	    
	    if (m(j,j) != norm_j)
	      {
		T temp = 1 - norm_j/m(j,j);
		
		for (int i = j+1; j < N; j++)
		  {
		    g[i] += std::abs(m(i,j))*temp;
		  }
	      }

	    // Now we do the classic cholesky iteration
	    T L_jj = std::sqrt(m(j,j));
	    
	    m(j,j) = L_jj;
	    for (int i = j+1; i < N; i++)
	      {
		m(i,j) /= L_jj;
		for (int k = j+1; k <= i; k++)
		  m(i,k) -= m(i,j)*m(k,j);
	      }
	  }

	// The final 2x2 submatrix is special
	T A00 = m(N-2, N-2), A01 = m(N-2, N-1), A11 = m(N-1,N-1);
	T sq_DELTA = std::sqrt(square(A00-A11) + square(A01));
	T lambda_hi = 0.5*((A00+A11) + sq_DELTA);
	T lambda_lo = 0.5*((A00+A11) - sq_DELTA);
       
	delta = std::max(std::max((T)0, -lambda_lo + std::max(tau*sq_DELTA/(1-tau), tau_bar*gamma)),delta_prev);
	if (delta > 0)
	  {
	    m(N-1,N-1) += delta;
	    m(N,N) += delta;
	    delta_prev = delta;
	  }
	m(N-2,N-2) = A00 = std::sqrt(A00);
	m(N-1,N-2) = (A01 /= A00);
	m(N-1,N-1) = std::sqrt(A11-A01*A01);
	norm_E = delta_prev;
      }
  }

};


#endif
