#define restrict __restrict__

namespace StencilBuilder
{
  struct Grid
  {
    Grid(const int itot, const int jtot, const int ktot, const int gc)
      : itot(itot), jtot(jtot), ktot(ktot), gc(gc),
        ntot((itot+2*gc)*(jtot+2*gc)*(ktot+2*gc)),
        istart(gc),
        jstart(gc),
        kstart(gc),
        iend(itot+gc),
        jend(jtot+gc),
        kend(ktot+gc),
        icells(itot+2*gc),
        ijcells((itot+2*gc)*(jtot+2*gc)) {}

    const int itot;
    const int jtot;
    const int ktot;
    const int gc;
    const int ntot;
    const int istart;
    const int jstart;
    const int kstart;
    const int iend;
    const int jend;
    const int kend;
    const int icells;
    const int ijcells;
  };

  // STENCIL OPERATORS.
  // Fourth order interpolation.
  struct Interp
  {
    static double apply_narrow(const double a, const double b) { return ( 9./16.)*(a+b); }
    static double apply_wide  (const double a, const double b) { return (-1./16.)*(a+b); }
  };

  // Fourth order gradient.
  struct Grad
  {
    static double apply_narrow(const double a, const double b) { return (-27./24.)*(a-b); }
    static double apply_wide  (const double a, const double b) { return (  1./24.)*(a-b); }
  };

  // STENCIL NODE CLASS
  // Stencil node in expression tree.
  template<int toCenter, class Inner, class Op>
  struct Stencil
  {
    Stencil(const Inner &inner, const int nn) : inner_(inner), nn_(nn) {}

    const Inner &inner_;
    const int nn_;

    double operator[](const int i) const
    {
      const double wide   = Op::apply_wide  (inner_[i + (-2+toCenter)*nn_], inner_[i + (+1+toCenter)*nn_]);
      const double narrow = Op::apply_narrow(inner_[i + (-1+toCenter)*nn_], inner_[i + (   toCenter)*nn_]);
      return wide + narrow;
    }
  };

  // Stencil generation operator for interpolation.
  template<int toCenter, class Inner>
  Stencil<toCenter, Inner, Interp> interp(const Inner &inner, const int nn)
  {
    return Stencil<toCenter, Inner, Interp>(inner, nn);
  }

  // Stencil generation operator for gradient.
  template<int toCenter, class Inner>
  Stencil<toCenter, Inner, Grad> grad(const Inner &inner, const int nn)
  {
    return Stencil<toCenter, Inner, Grad>(inner, nn);
  }

  // SCALAR OPERATORS
  // Multiplication operator.
  struct Multiply
  {
    static double apply(const double left, const double right) { return left*right; }
  };

  // Addition operator.
  struct Add
  {
    static double apply(const double left, const double right) { return left+right; }
  };

  // OPERATOR NODE CLASS
  // Operator node in expression tree.
  template<class Left, class Op, class Right>
  struct Operator
  {
    Operator(const Left &left, const Right &right) : left_(left), right_(right) {}

    const Left &left_;
    const Right &right_;

    double operator[](const int i) const { return Op::apply(left_[i], right_[i]); }
  };

  // Operator aggregation class, specialization for left scalar multiplication
  template<class Op, class Right>
  struct Operator<double, Op, Right>
  {
    Operator(const double &left, const Right &right) : left_(left), right_(right) {}

    const double &left_;
    const Right &right_;

    double operator[](const int i) const { return Op::apply(left_, right_[i]); }
  };

  // Template classes for the multiplication operator.
  template<class Left, class Right>
  Operator<Left, Multiply, Right> operator*(const Left &left, const Right &right)
  {
    return Operator<Left, Multiply, Right>(left, right);
  }

  // Template classes for the addition operators.
  template<class Left, class Right>
  Operator<Left, Add, Right> operator+(const Left &left, const Right &right)
  {
    return Operator<Left, Add, Right>(left, right);
  }

  // Field class representing the field, whose operations expand compile time.
  class Field
  {
    public:
      Field(const Grid &grid)
        : grid_(grid),
          data_(new double[grid_.ntot]) {}

      ~Field() { delete[] data_; }

      double& operator[](const int i) const { return data_[i]; }

      double& operator()(const int i, const int j, const int k) const
      { return data_[i + j*grid_.icells + k*grid_.ijcells]; }

      // Assignment operator, this operator starts the inline expansion.
      template<class T> void operator= (const T& restrict expression)
      {
        const int jj = grid_.icells;
        const int kk = grid_.ijcells;

        for (int k=grid_.kstart; k<grid_.kend; ++k)
          for (int j=grid_.jstart; j<grid_.jend; ++j)
            #pragma ivdep
            for (int i=grid_.istart; i<grid_.iend; ++i)
            {
              const int ijk = i + j*jj + k*kk;
              data_[ijk] = expression[ijk];
            }
      }

      // Compound assignment operator, this operator starts the inline expansion.
      template<class T> void operator+=(const T& restrict expression)
      {
        const int jj = grid_.icells;
        const int kk = grid_.ijcells;

        for (int k=grid_.kstart; k<grid_.kend; ++k)
          for (int j=grid_.jstart; j<grid_.jend; ++j)
            #pragma ivdep
            for (int i=grid_.istart; i<grid_.iend; ++i)
            {
              const int ijk = i + j*jj + k*kk;
              data_[ijk] += expression[ijk];
            }
      }

    private:
      // Reference to the grid on which the field is created
      const Grid &grid_;
      // Pointer to the raw data.
      double * restrict data_;
  };

  // Specialization for assignment with a constant.
  template<> void Field::operator= (const double& restrict expression)
  {
    const int jj = grid_.icells;
    const int kk = grid_.ijcells;

    for (int k=grid_.kstart; k<grid_.kend; ++k)
      for (int j=grid_.jstart; j<grid_.jend; ++j)
        #pragma ivdep
        for (int i=grid_.istart; i<grid_.iend; ++i)
        {
          const int ijk = i + j*jj + k*kk;
          data_[ijk] = expression;
        }
  }
}
