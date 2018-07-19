import numexpr as ne
import numpy as np


def cicParticles(particles, L, N):

    if type(N) not in [int,int]:
      raise TypeError("N must be a numeric type")

    def shifted(i, t):
      a = np.empty(i[0].size, dtype=np.int64)
      return ne.evaluate('(i2+t2)%N + N*((i1+t1)%N + N*((i0+t0)%N) )', local_dict={'i2':i[2], 't2':t[2], 'i1':i[1], 't1':t[1], 'i0':i[0], 't0':t[0], 'N':N}, out=a)


    i =[]
    r = []
    for d in range(3):
        q = ne.evaluate('(p%L)*N/L', local_dict={'p':particles[d], 'L':L, 'N':N })
        o = np.empty(q.size, dtype=np.int64)
        o[:] = np.floor(q)
        i.append(o)
        r.append(ne.evaluate('q-o'))

    D = {'a':r[0],'b':r[1],'c':r[2]}
    N3 = N*N*N

    def accum(density, ss, op):
      d0 = np.bincount(shifted(i, ss), weights=ne.evaluate(op, local_dict=D), minlength=N3)
      ne.evaluate('d + d0', local_dict={'d':density, 'd0':d0}, out=density)

    density = np.empty(N3, dtype=np.float64)

    accum(density, (1,1,1), 'a     *    b  *    c ')
    accum(density, (1,1,0), 'a     *    b  * (1-c)')
    accum(density, (1,0,1), 'a     * (1-b) *    c ')
    accum(density, (1,0,0), 'a     * (1-b) * (1-c)')
    accum(density, (0,1,1), '(1-a) *    b  *    c ')
    accum(density, (0,1,0), '(1-a) *    b  * (1-c)')
    accum(density, (0,0,1), '(1-a) * (1-b) *    c ')
    accum(density, (0,0,0), '(1-a) * (1-b) * (1-c)')

    return density.reshape((N,N,N))
