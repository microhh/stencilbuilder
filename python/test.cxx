const double ci0 = 1.;
const double ci1 = 1.;
const double ci2 = 1.;
const double ci3 = 1.;

const double cg0 = 1.;
const double cg1 = 1.;
const double cg2 = 1.;
const double cg3 = 1.;

void advecu(double *ut, double *u, double *v, double *w,
            double *dzi4, double *dzhi4, double dxi, double dyi)
{
  const int istart = 0;
  const int iend   = 1;
  const int jstart = 0;
  const int jend   = 1;
  const int kstart = 0;
  const int kend   = 1;
  const int jj = 0;
  const int kk = 0;

  for (int k=kstart; k<kend; ++k)
    for (int j=jstart; j<jend; ++j)
      for (int i=istart; i<iend; ++i)
      {
        const int ijk = i + j*jj + k*kk;
        //$ SBStart ut
        u = Field("u", uloc)
        v = Field("v", vloc)
        w = Field("w", wloc)
        
        dxi = Scalar("dxi")
        dyi = Scalar("dyi")

        dzi4 = Vector("dzi4", zloc)

        ut = gradx( interpx(u) * interpx(u) ) * dxi \
           + grady( interpx(v) * interpy(u) ) * dyi \
           + gradz( interpx(w) * interpz(u) ) * dzi4
        //$ SBEnd
      }
}

void advecv(double *vt, double *u, double *v, double *w,
            double *dzi4, double *dzhi4, double dxi, double dyi)
{
  const int istart = 0;
  const int iend   = 1;
  const int jstart = 0;
  const int jend   = 1;
  const int kstart = 0;
  const int kend   = 1;
  const int jj = 1;
  const int kk = 1;

  for (int k=kstart; k<kend; ++k)
    for (int j=jstart; j<jend; ++j)
      for (int i=istart; i<iend; ++i)
      {
        const int ijk = i + j*jj + k*kk;
        //$ SBStart vt
        u = Field("u", uloc)
        v = Field("v", vloc)
        w = Field("w", wloc)
        
        dxi = Scalar("dxi")
        dyi = Scalar("dyi")

        dzi4 = Vector("dzi4", zloc)

        vt = gradx( interpy(u) * interpx(v) ) * dxi \
           + grady( interpy(v) * interpy(v) ) * dyi \
           + gradz( interpy(w) * interpz(v) ) * dzi4
        //$ SBEnd
      }
}

int main()
{
  double *u, *v, *w;
  double *ut, *vt;
  double *dzi4, *dzhi4;
  double dxi, dyi;

  advecu(ut, u, v, w, dzi4, dzhi4, dxi, dyi);
  advecv(vt, u, v, w, dzi4, dzhi4, dxi, dyi);

  return 0;
}
