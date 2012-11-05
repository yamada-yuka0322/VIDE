#ifndef _MAK_MINIARGS_HPP
#define _MAK_MINIARGS_HPP

namespace CosmoTool
{
  typedef enum 
    {
      MINIARG_NULL,
      MINIARG_STRING,
      MINIARG_INT,
      MINIARG_DOUBLE,
      MINIARG_FLOAT,
      MINIARG_DOUBLE_3D_VECTOR
    } MiniArgType;
  
  typedef struct
  {
    const char *name;
    void *data;
    MiniArgType argType;
  } MiniArgDesc;

  int parseMiniArgs(int argc, char **argv, MiniArgDesc *desc);
};

#endif
