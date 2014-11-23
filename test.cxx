void advecu(double *ut, double *u, double *v, double *w,
            double *dzi4, double *dzhi4, double dxi, double dyi)
{
  const int istart = 0;
  const int iend   = 0;
  const int jstart = 0;
  const int jend   = 0;
  const int kstart = 0;
  const int kend   = 0;
  const int jj = 0;
  const int kk = 0;

  for (int k=kstart; k<kend; ++k)
    for (int j=jstart; j<jend; ++j)
      for (int i=istart; i<iend; ++i)
      {
        const int ijk = i + j*jj + k*kk;
        //$ SBStart
        uloc = np.array([1,0,0])
        vloc = np.array([0,1,0])
        wloc = np.array([0,0,1])
        zloc = 0
        
        u = Field("u", uloc)
        v = Field("v", vloc)
        w = Field("w", wloc)
        
        dxi = Scalar("dxi")
        dyi = Scalar("dyi")

        dzi4 = Vector("dzi4", zloc)

        ut = gradx( interpx(u) * interpx(u) ) * dxi \
           + grady( interpx(v) * interpy(u) ) * dyi \
           + gradz( interpx(w) * interpz(u) ) * dzi4

        print("ut[i,j,k] = {0};".format(ut.getString(0,0,0,12)))
        //$ SBEnd
}

int main()
{
  double *u, *v, *w;
  double *dzi4, *dzhi4;
  double dxi, dyi;

  advecu(u, v, w, dzi4, dzhi4, dxi, dyi);

  return 0;
}
