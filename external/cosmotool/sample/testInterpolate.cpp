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
