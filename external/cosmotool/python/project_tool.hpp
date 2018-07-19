
// Only in 3d

template<typename T, typename ProdType>
static T project_tool(T *vertex_value, T *u, T *u0)
{
  T ret0 = 0;
  for (unsigned int i = 0; i < 8; i++)
    {
      unsigned int c[3] = { i & 1, (i>>1)&1, (i>>2)&1 };
      int epsilon[3];
      T ret = 0;
      
      for (int q = 0; q < 3; q++)
        epsilon[q] = 2*c[q]-1;

      for (int q = 0; q < ProdType::numProducts; q++)
        ret += ProdType::product(u, u0, epsilon, q);
      ret *= vertex_value[i];
      ret0 += ret;
    }
    
  return ret0;
}


template<typename T>
static inline T get_u0(const T& u0, int epsilon)
{
  return (1-epsilon)/2 + epsilon*u0;
//  return (epsilon > 0) ? u0 : (1-u0);
}

template<typename T>
struct ProductTerm0
{
  static const int numProducts = 1;
  
  static inline T product(T *u, T *u0, int *epsilon, int q)
  {
    T a = 1;
    
    for (unsigned int r = 0; r < 3; r++)
      a *= get_u0(u0[r], epsilon[r]);
    return a;
  }
};


template<typename T>
struct ProductTerm1
{
  static const int numProducts = 3;
  
  static T product(T *u, T *u0, int *epsilon, int q)
  {
    T a = 1;
    T G[3];
    
    for (unsigned int r = 0; r < 3; r++)
      {
        G[r] = get_u0(u0[r], epsilon[r]);
      }

    T F[3] = { G[1]*G[2], G[0]*G[2], G[0]*G[1] };

    return F[q] * u[q] * epsilon[q];
  }
};

template<typename T>
struct ProductTerm2
{
  static const int numProducts = 3;
  
  static inline T product(T *u, T *u0, int *epsilon, int q)
  {
    T a = 1;
    T G[3];
    
    for (unsigned int r = 0; r < 3; r++)
      {
        G[r] = get_u0(u0[r], epsilon[r]);
      }

    T F[3] = { epsilon[1]*epsilon[2]*u[1]*u[2], epsilon[0]*epsilon[2]*u[0]*u[2], epsilon[0]*epsilon[1]*u[0]*u[1] };    

    return F[q] * G[q];
  }
};



template<typename T>
struct ProductTerm3
{
  static const int numProducts = 1;
  
  static inline T product(T *u, T *u0, int *epsilon, int q)
  {
    return epsilon[0]*epsilon[1]*epsilon[2]*u[0]*u[1]*u[2];
  }
};


template<typename T>
T compute_projection(T *vertex_value, T *u, T *u0, T rho)
{
  T ret;
  
  ret = project_tool<T, ProductTerm0<T> >(vertex_value, u, u0) * rho;
  ret += project_tool<T, ProductTerm1<T> >(vertex_value, u, u0) * rho * rho / 2;
  ret += project_tool<T, ProductTerm2<T> >(vertex_value, u, u0) * rho * rho * rho / 3;
  ret += project_tool<T, ProductTerm3<T> >(vertex_value, u, u0) * rho * rho * rho * rho / 4;
  return ret;
}

template
double compute_projection(double *vertex_value, double *u, double *u0, double rho);

