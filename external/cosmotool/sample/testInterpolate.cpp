/*+
This is CosmoTool (./sample/testInterpolate.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <iostream>
#include "interpolate3d.hpp"

using namespace std;
using namespace CosmoTool;

int main()
{
   VectorField<float,3> *vectors = new VectorField<float,3>[27];

   for (int i = 0; i < 3; i++)
     for (int j = 0; j < 3; j++)
       for (int k = 0; k < 3; k++)
         {
            int idx = i + 3*(j + 3*k);
            vectors[idx].vec[0] = i;
            vectors[idx].vec[1] = j;
	    vectors[idx].vec[2] = k;
         }

   GridSampler<VectorField<float,3> > sampler(vectors, 3, 3, 3, 1);
   Interpolate3D<GridSampler<VectorField<float,3> > > inter(sampler);

   VectorField<float,3> v = inter.get(0.5,0.15,0.5);
   VectorField<float,3> v2 = inter.get(1.5,1.65,1.5);

   cout << v.vec[0] << " " << v.vec[1] << " " << v.vec[2] << endl;
   cout << v2.vec[0] << " " << v2.vec[1] << " " << v2.vec[2] << endl;

   return 0;
}
