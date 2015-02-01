#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "StencilBuilder.h"

using namespace StencilBuilder;

int main()
{
  // Test configuration settings.
  const int itot = 256;
  const int jtot = 256;
  const int ktot = 256;
  const int gc = 4;
  const int iter = 5;

  // Initialize the grid.
  Grid grid(itot, jtot, ktot, gc);

  // Create fields on the grid.
  Field a (grid);
  Field b (grid);
  Field c (grid);
  Field at(grid);

  // Initialize the fields.
  for (int n=0; n<grid.ntot; ++n)
  {
    a[n] = 0.001 * (std::rand() % 1000) - 0.5;
    b[n] = 0.001 * (std::rand() % 1000) - 0.5;
    c[n] = 0.001 * (std::rand() % 1000) - 0.5;

    at[n] = 0.;
  }

  // Initialize the time step.
  const double dt = 1.e-3;

  // Execute the loop iter times.
  for (int n=0; n<iter; ++n)
  {
    // Advection operator.
    at += Gx_h( Ix  (a) * Ix  (a) )
        + Gy  ( Ix_h(b) * Iy_h(a) )
        + Gz  ( Ix_h(c) * Iz_h(a) );

    // Time integration.
    a += dt*at;

    // Tendency reset.
    at = 0.;
  }

  // Print a value in the middle of the field.
  std::cout << std::setprecision(8) << "a = " << a(itot/2, jtot/2, ktot/2) << std::endl;

  return 0;
}
