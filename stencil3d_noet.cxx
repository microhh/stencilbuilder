#include <iostream>
#include <iomanip>
#include <cstdlib>

// Fourth order interpolation function.
inline double interp(const double m2, const double m1, const double p1, const double p2)
{
  return (-1./16)*m2 + (9./16)*m1 + (9./16)*p1 + (-1./16)*p2;
}

// Fourth order gradient function.
inline double grad(const double m2, const double m1, const double p1, const double p2)
{
  return (1./24.)*m2 + (-27./24.)*m1 + (27./24.)*p1 + (-1./24.)*p2;
}

// Test function with a similar structure as the advection operator.
void advection(double * const __restrict__ at, const double * const __restrict__ a,
               const double * const __restrict__ b, const double * const __restrict__ c,
               const int istart, const int iend,
               const int jstart, const int jend,
               const int kstart, const int kend,
               const int icells, const int ijcells)
{
  const int ii1 = 1;
  const int ii2 = 2;
  const int ii3 = 3;
  const int jj1 = 1*icells;
  const int jj2 = 2*icells;
  const int jj3 = 3*icells;
  const int kk1 = 1*ijcells;
  const int kk2 = 2*ijcells;
  const int kk3 = 3*ijcells;

  for (int k=kstart; k<kend; ++k)
    for (int j=jstart; j<jend; ++j)
      #pragma ivdep
      for (int i=istart; i<iend; ++i)
      {
        const int ijk = i + j*jj1 + k*kk1;
        at[ijk] += grad( interp( a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk    ] ) * interp( a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk    ] ),
                         interp( a[ijk-ii2], a[ijk-ii1], a[ijk    ], a[ijk+ii1] ) * interp( a[ijk-ii2], a[ijk-ii1], a[ijk    ], a[ijk+ii1] ),
                         interp( a[ijk-ii1], a[ijk    ], a[ijk+ii1], a[ijk+ii2] ) * interp( a[ijk-ii1], a[ijk    ], a[ijk+ii1], a[ijk+ii2] ),
                         interp( a[ijk    ], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3] ) * interp( a[ijk    ], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3] ))

                 + grad( interp( b[ijk-ii2-jj1], b[ijk-ii1-jj1], b[ijk-jj1], b[ijk+ii1-jj1] ) * interp( a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk    ] ),
                         interp( b[ijk-ii2    ], b[ijk-ii1    ], b[ijk    ], b[ijk+ii1    ] ) * interp( a[ijk-jj2], a[ijk-jj1], a[ijk    ], a[ijk+jj1] ),
                         interp( b[ijk-ii2+jj1], b[ijk-ii1+jj1], b[ijk+jj1], b[ijk+ii1+jj1] ) * interp( a[ijk-jj1], a[ijk    ], a[ijk+jj1], a[ijk+jj2] ),
                         interp( b[ijk-ii2+jj2], b[ijk-ii1+jj2], b[ijk+jj2], b[ijk+ii1+jj2] ) * interp( a[ijk    ], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3] ))

                 + grad( interp( c[ijk-ii2-kk1], c[ijk-ii1-kk1], c[ijk-kk1], c[ijk+ii1-kk1] ) * interp( a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk    ] ),
                         interp( c[ijk-ii2    ], c[ijk-ii1    ], c[ijk    ], c[ijk+ii1    ] ) * interp( a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1] ),
                         interp( c[ijk-ii2+kk1], c[ijk-ii1+kk1], c[ijk+kk1], c[ijk+ii1+kk1] ) * interp( a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2] ),
                         interp( c[ijk-ii2+kk2], c[ijk-ii1+kk2], c[ijk+kk2], c[ijk+ii1+kk2] ) * interp( a[ijk    ], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3] ));
      }
}

// Test function for time integration.
void tendency(double * const __restrict__ at, double * const __restrict__ a,
              const double dt,
              const int istart, const int iend,
              const int jstart, const int jend,
              const int kstart, const int kend,
              const int icells, const int ijcells)
{
  const int ii = 1;
  const int jj = icells;
  const int kk = ijcells;

  for (int k=kstart; k<kend; ++k)
    for (int j=jstart; j<jend; ++j)
      #pragma ivdep
      for (int i=istart; i<iend; ++i)
      {
        const int ijk = i + j*jj + k*kk;
        a[ijk] += dt*at[ijk];
      }

  for (int k=kstart; k<kend; ++k)
    for (int j=jstart; j<jend; ++j)
      #pragma ivdep
      for (int i=istart; i<iend; ++i)
      {
        const int ijk = i + j*jj + k*kk;
        at[ijk] = 0.;
      }
}

int main()
{
  // Test configuration settings.
  const int itot = 256;
  const int jtot = 256;
  const int ktot = 256;
  const int gc   = 4;
  const int iter = 5;

  // Calculate the required variables.
  const int ntot = (itot+2*gc)*(jtot+2*gc)*(ktot+2*gc);
  const int istart = gc;
  const int jstart = gc;
  const int kstart = gc;
  const int iend = itot+gc;
  const int jend = jtot+gc;
  const int kend = ktot+gc;
  const int icells = itot+2*gc;
  const int ijcells = (itot+2*gc)*(jtot+2*gc);

  // Allocate the raw arrays.
  double *a_data  = new double[ntot];
  double *b_data  = new double[ntot];
  double *c_data  = new double[ntot];
  double *at_data = new double[ntot];

  // Initialize the raw arrays.
  for (int n=0; n<ntot; ++n)
  {
    a_data[n] = 0.001 * (std::rand() % 1000) - 0.5;
    b_data[n] = 0.001 * (std::rand() % 1000) - 0.5;
    c_data[n] = 0.001 * (std::rand() % 1000) - 0.5;
    at_data[n] = 0.;
  }

  // Initialize a time step.
  const double dt = 1.e-3;

  // Execute the loop iter times.
  for (int ii=0; ii<iter; ++ii)
  {
    advection(at_data, a_data, b_data, c_data,
              istart, iend,
              jstart, jend,
              kstart, kend,
              icells, ijcells);

    tendency(at_data, a_data,
             dt,
             istart, iend,
             jstart, jend,
             kstart, kend,
             icells, ijcells);
  }

  // Print a value in the middle of the field.
  const int ijk = itot/2 + (jtot/2)*icells + (ktot/2)*ijcells;
  std::cout << std::setprecision(8) << "a = " << a_data[ijk] << std::endl;

  return 0;
}
