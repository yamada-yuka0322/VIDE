#ifndef __VOZ_IO_HPP
#define __VOZ_IO_HPP

#include <string>

struct PositionData
{
  float *xyz[3];
  pid_t np;
  
  void destroy()
  {
    for (int j = 0; j < 3; j++)
      delete[] xyz[j];
  }

  void findExtrema();

  bool readFrom(const std::string& fname);
};

#endif
