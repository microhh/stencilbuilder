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
  Field u (grid);
  Field v (grid);
  Field w (grid);
  Field ut(grid);

  // Initialize the fields.
  u.randomize();
  v.randomize();
  w.randomize();

  // Initialize the time step.
  const double dt = 1.e-3;

  // Execute the loop iter times.
  for (int n=0; n<iter; ++n)
  {
    // Advection operator.
    ut += Gx_h( Ix  (u) * Ix  (u) )
        + Gy  ( Ix_h(v) * Iy_h(u) )
        + Gz  ( Ix_h(w) * Iz_h(u) );

    // Time integration.
    u += dt*ut;

    // Tendency reset.
    ut = 0.;
  }

  // Print a value in the middle of the field.
  std::cout << std::setprecision(8) << "u = " << u(itot/2, jtot/2, ktot/2) << std::endl;

  return 0;
}
