from .timing import time_block as time_block_orig
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
from contextlib import contextmanager

TIMER_ACTIVE=False

@contextmanager
def time_block_dummy(*args):
  yield

if TIMER_ACTIVE:
  time_block=time_block_orig
else: 
  time_block=time_block_dummy

CIC_PREKERNEL='''
#define NDIM {ndim}
#define CENTERED {centered}
typedef {cicType} BASIC_TYPE;

'''

CIC_KERNEL='''///CL///
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

__kernel void init_pcell(__global int *p_cell, const int value)
{
    int i = get_global_id(0);
    p_cell[i] = value;
}

__kernel void build_indices(__global const BASIC_TYPE *pos,
                            __global int *part_mesh, __global int *part_list, const int N, const BASIC_TYPE delta, const BASIC_TYPE shift_pos)
{
    int i_part = get_global_id(0);
    long shifter = 1;
    long idx = 0;
    int d;

    for (d = 0; d < NDIM; d++) {
        BASIC_TYPE x;
        
        if (CENTERED)
          x = pos[i_part*NDIM + d] - shift_pos;
        else
          x = pos[i_part*NDIM + d];

        int m = (int)floor(x*delta) %% N;
        
        idx += shifter * m;        
        shifter *= N;
    }
    
    // Head of the list
    int initial_elt = atom_xchg(&part_mesh[idx], i_part);
    if (initial_elt == -1) {
        return;
    }
    // Point the next pointer of old_end to i_part
    part_list[i_part] = initial_elt;
}

__kernel void reverse_list(__global int *part_mesh, __global int *part_list)
{
    int mid = get_global_id(0);
   
    int current_part = part_mesh[mid];
    if (current_part >= 0) {
        int next_part = part_list[current_part];
        part_list[current_part] = -1;
        while (next_part != -1) {
          int p = part_list[next_part];
          part_list[next_part] = current_part;
          current_part = next_part;
          next_part = p;
        }
        part_mesh[mid] = current_part;
    }
}

__kernel void dance(__global const BASIC_TYPE *pos,
                    __global BASIC_TYPE *density,
                    __global int *part_mesh, __global int *part_list, const int N, const BASIC_TYPE delta, const BASIC_TYPE shift_pos)
{
    int m[NDIM];
    int shifter = 1;    
    int i;
    int first, i_part;
    int idx = 0;
    
    for (i = 0; i < NDIM; i++) {
        m[i] = get_global_id(i);
        idx += shifter * m[i];
        shifter *= N;
    }

    first = 1;    

//BEGIN LOOPER
%(looperFor)s
//END LOOPER        

            int idx_dance = 0; 
            BASIC_TYPE w = 0;
//LOOPER INDEX
            int r[NDIM] = { %(looperVariables)s };
//END LOOPER

            i_part = part_mesh[idx];
            while (i_part != -1) {
                BASIC_TYPE w0 = 1;
                
                for (int d = 0; d < NDIM; d++) {
                    BASIC_TYPE x;
                    BASIC_TYPE q;
                    BASIC_TYPE dx;

                    if (CENTERED)
                      x = pos[i_part*NDIM + d]*delta - shift_pos;
                    else
                      x = pos[i_part*NDIM + d]*delta;
                    q = floor(x);
                    dx = x - q;
                    
                    w0 *= (r[d] == 1) ? dx : ((BASIC_TYPE)1-dx);
                }

                i_part = part_list[i_part];
                w += w0;
            }
        
            shifter = 1;    
            for (i = 0; i < NDIM; i++) {
                idx_dance += shifter * ((m[i]+r[i])%%N);
                shifter *= N;
            }

            density[idx_dance] += w;
                
            // One dance done. Wait for everybody for the next iteration
            barrier(CLK_GLOBAL_MEM_FENCE);
%(looperForEnd)s
}
'''

class CIC_CL(object):

    def __init__(self, context, ndim=2, ktype=np.float32, centered=False):
        global CIC_PREKERNEL, CIC_KERNEL

        translator = {}
        if ktype == np.float32:
            translator['cicType'] = 'float'
            pragmas = ''
        elif ktype == np.float64:
            translator['cicType'] = 'double'
            pragmas = '#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n'
        else:
            raise ValueError("Invalid ktype")
            
        # 2 dimensions
        translator['ndim'] = ndim
        translator['centered'] = '1' if centered else '0'
        looperVariables = ','.join(['id%d' % d for d in range(ndim)])
        looperFor = '\n'.join(['for (int id{dim}=0; id{dim} < 2; id{dim}++) {{'.format(dim=d) for d in range(ndim)])
        looperForEnd = '}' * ndim
        
        kern = pragmas + CIC_PREKERNEL.format(**translator) + (CIC_KERNEL % {'looperVariables': looperVariables, 'looperFor': looperFor, 'looperForEnd':looperForEnd})
        self.kern_code = kern
        self.ctx = context
        self.queue = cl.CommandQueue(context)#, properties=cl.OUT_OF_ORDER_EXEC_MODE_ENABLE)
        self.ktype = ktype
        self.ndim = ndim
        self.prog = cl.Program(self.ctx, kern).build()
        self.centered = centered

    def run(self, particles, Ng, L):
        assert particles.strides[1] == self.ktype().itemsize  # This is C-ordering
        assert particles.shape[1] == self.ndim

        print("Start again")

        ndim = self.ndim
        part_pos = cl_array.to_device(self.queue, particles)
        part_mesh = cl_array.empty(self.queue, (Ng,)*ndim, np.int32, order='C')
        density = cl_array.zeros(self.queue, (Ng,)*ndim, self.ktype, order='C')
        part_list = cl_array.empty(self.queue, (particles.shape[0],), np.int32, order='C')
        shift_pos = 0.5*L if self.centered else 0
        
        if True:
          delta = Ng/L
    
          with time_block("Init pcell array"):
              e = self.prog.init_pcell(self.queue, (Ng**ndim,), None, part_mesh.data, np.int32(-1))
              e.wait()
        
          with time_block("Init idx array"):
              e=self.prog.init_pcell(self.queue, (particles.shape[0],), None, part_list.data, np.int32(-1))
              e.wait()
        
          with time_block("Build indices"):
              self.prog.build_indices(self.queue, (particles.shape[0],), None,
                      part_pos.data, part_mesh.data, part_list.data, np.int32(Ng), self.ktype(delta), self.ktype(shift_pos))
        
        if True:
          with time_block("Reverse list"):
              lastevt = self.prog.reverse_list(self.queue, (Ng**ndim,), None, part_mesh.data, part_list.data)
        # We require pmax pass, particles are ordered according to part_idx

          with time_block("dance"):
              self.prog.dance(self.queue, (Ng,)*ndim, None, part_pos.data, density.data, part_mesh.data, part_list.data, np.int32(Ng), self.ktype(delta), self.ktype(shift_pos))
    
        self.queue.finish()
        del part_pos
        del part_mesh
        del part_list
        with time_block("download"):
          return density.get()


def cl_CIC_Density(particles, Ngrid, Lbox, context=None, periodic=True, centered=False):
  """
  cl_CIC_Density(particles (Nx3), Ngrid, Lbox, context=None, periodic=True, centered=False)
"""
  if context is None:
    context = cl.create_some_context()

  ktype = particles.dtype
  if ktype != np.float32 and ktype != np.float64:
    raise ValueError("particles may only be float32 or float64")

  if len(particles.shape) != 2 or particles.shape[1] != 3:
    raise ValueError("particles may only be a Nx3 array")

  cic = CIC_CL(context, ndim=3, centered=centered, ktype=ktype)

  return cic.run(particles, Ngrid, Lbox)
