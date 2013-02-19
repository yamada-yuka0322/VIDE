#include <cmath>
#include <iostream>
#include "cosmopower.hpp"

using namespace std;
using namespace CosmoTool;

int main()
{
   CosmoPower cp;
   CosmoPower::CosmoFunction f[] = { CosmoPower::POWER_EFSTATHIOU, CosmoPower::HU_WIGGLES, CosmoPower::HU_BARYON, CosmoPower::POWER_SUGIYAMA };
   int num_F = sizeof(f)/sizeof(f[0]);

   cp.setFunction(f[0]);
   cp.normalize();
   for (int ik = 0; ik < 100; ik++)
    {
      double k = pow(10.0, 4*ik/100.-3);
      
      cout << k << " ";
      for (int q = 0; q < num_F; q++)
        {
          cp.setFunction(f[q]);
          cout << cp.power(k) << " ";
        } 
      cout << endl;
    }
   return 0;
}
