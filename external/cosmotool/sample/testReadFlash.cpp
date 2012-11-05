#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "loadFlash.hpp"

using namespace CosmoTool;
using namespace std;

int main () {
 
  const char* filename = "lss_read_hdf5_chk_0000";
 
  SimuData* data = CosmoTool::loadFlashMulti(filename, 0, 0);

  return 0;
}
