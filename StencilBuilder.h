#define restrict __restrict__

namespace StencilBuilder
{
  // Unit vectors
  constexpr int iVec[3] = {1, 0, 0};
  constexpr int jVec[3] = {0, 1, 0};
  constexpr int kVec[3] = {0, 0, 1};

  struct Grid
  {
    Grid(const int itot, const int jtot, const int ktot, const int gc) :
      itot(itot), jtot(jtot), ktot(ktot), gc(gc),
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
    static inline double apply(const double a, const double b, const double c, const double d)
    { return (9./16.)*(b+c) - (1./16.)*(a+d); }
  };

  // Fourth order gradient.
  struct Grad
  {
    static inline double apply(const double a, const double b, const double c, const double d)
    { return (27./24.)*(c-b) - (1./24)*(d-a); }
  };

  // STENCIL NODE CLASS
  // Stencil node in expression tree.
  template<int toCenter, class Inner, class Op, const int vec[3]>
  struct Stencil
  {
    Stencil(const Inner& inner) : inner_(inner) {}

    const Inner& inner_;

    // inline double operator[](const int i) const
    // {
    //   return Op::apply(inner_[i + (-2+toCenter)*nn_], inner_[i + (-1+toCenter)*nn_],
    //                    inner_[i + (   toCenter)*nn_], inner_[i + ( 1+toCenter)*nn_]);
    // }

    inline double operator()(const int i, const int j, const int k) const
    {
      return Op::apply(inner_(i + vec[0]*(-2+toCenter), j + vec[1]*(-2+toCenter), k + vec[2]*(-2+toCenter)),
                       inner_(i + vec[0]*(-1+toCenter), j + vec[1]*(-1+toCenter), k + vec[2]*(-1+toCenter)),
                       inner_(i + vec[0]*(   toCenter), j + vec[1]*(   toCenter), k + vec[2]*(   toCenter)),
                       inner_(i + vec[0]*( 1+toCenter), j + vec[1]*( 1+toCenter), k + vec[2]*( 1+toCenter)));
    }
  };

  // Stencil generation operator for interpolation.
  template<int toCenter, class Inner>
  inline Stencil<toCenter, Inner, Interp, iVec> interpx(const Inner& inner)
  {
    return Stencil<toCenter, Inner, Interp, iVec>(inner);
  }

  template<int toCenter, class Inner>
  inline Stencil<toCenter, Inner, Interp, jVec> interpy(const Inner& inner)
  {
    return Stencil<toCenter, Inner, Interp, jVec>(inner);
  }

  template<int toCenter, class Inner>
  inline Stencil<toCenter, Inner, Interp, kVec> interpz(const Inner& inner)
  {
    return Stencil<toCenter, Inner, Interp, kVec>(inner);
  }

  // Stencil generation operator for gradient.
  template<int toCenter, class Inner>
  inline Stencil<toCenter, Inner, Grad, iVec> gradx(const Inner& inner)
  {
    return Stencil<toCenter, Inner, Grad, iVec>(inner);
  }

  template<int toCenter, class Inner>
  inline Stencil<toCenter, Inner, Grad, jVec> grady(const Inner& inner)
  {
    return Stencil<toCenter, Inner, Grad, jVec>(inner);
  }

  template<int toCenter, class Inner>
  inline Stencil<toCenter, Inner, Grad, kVec> gradz(const Inner& inner)
  {
    return Stencil<toCenter, Inner, Grad, kVec>(inner);
  }

  // SCALAR OPERATORS
  // Multiplication operator.
  struct Multiply
  {
    static inline double apply(const double left, const double right) { return left*right; }
  };

  // Addition operator.
  struct Add
  {
    static inline double apply(const double left, const double right) { return left+right; }
  };

  // OPERATOR NODE CLASS
  // Operator node in expression tree.
  template<class Left, class Op, class Right>
  struct Operator
  {
    Operator(const Left& left, const Right& right) : left_(left), right_(right) {}

    const Left& left_;
    const Right& right_;

    // inline double operator[](const int i) const { return Op::apply(left_[i], right_[i]); }
    inline double operator()(const int i, const int j, const int k) const
    { return Op::apply(left_(i, j, k), right_(i, j, k)); }
  };

  // Operator aggregation class, specialization for left scalar multiplication
  template<class Op, class Right>
  struct Operator<double, Op, Right>
  {
    Operator(const double& left, const Right& right) : left_(left), right_(right) {}

    const double& left_;
    const Right& right_;

    // inline double operator[](const int i) const { return Op::apply(left_, right_[i]); }
    inline double operator()(const int i, const int j, const int k) const
    { return Op::apply(left_, right_(i, j, k)); }
  };

  // Template classes for the multiplication operator.
  template<class Left, class Right>
  inline Operator<Left, Multiply, Right> operator*(const Left& left, const Right& right)
  {
    return Operator<Left, Multiply, Right>(left, right);
  }

  // Template classes for the addition operators.
  template<class Left, class Right>
  inline Operator<Left, Add, Right> operator+(const Left& left, const Right& right)
  {
    return Operator<Left, Add, Right>(left, right);
  }

  // Field class representing the field, whose operations expand compile time.
  class Field
  {
    public:
      Field(const Grid& grid) :
        grid_(grid),
        data_(new double[grid_.ntot]) {}

      ~Field() { delete[] data_; }

      inline double& operator[](const int i) const { return data_[i]; }

      inline double& operator()(const int i, const int j, const int k) const
      { return data_[i + j*grid_.icells + k*grid_.ijcells]; }

      // Assignment operator, this operator starts the inline expansion.
      template<class T>
      inline Field& operator= (const T& restrict expression)
      {
        const int jj = grid_.icells;
        const int kk = grid_.ijcells;

        for (int k=grid_.kstart; k<grid_.kend; ++k)
          for (int j=grid_.jstart; j<grid_.jend; ++j)
            #pragma ivdep
            for (int i=grid_.istart; i<grid_.iend; ++i)
            {
              const int ijk = i + j*jj + k*kk;
              data_[ijk] = expression(i, j, k);
            }

        return *this;
      }

      // Overload, NOT specialization, for assignment with a constant.
      inline Field& operator= (const double& restrict expression)
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

        return *this;
      }

      // Compound assignment operator, this operator starts the inline expansion.
      template<class T>
      inline Field& operator+=(const T& restrict expression)
      {
        const int jj = grid_.icells;
        const int kk = grid_.ijcells;

        for (int k=grid_.kstart; k<grid_.kend; ++k)
          for (int j=grid_.jstart; j<grid_.jend; ++j)
            #pragma ivdep
            for (int i=grid_.istart; i<grid_.iend; ++i)
            {
              const int ijk = i + j*jj + k*kk;
              data_[ijk] += expression(i, j, k);
            }

        return *this;
      }

    private:
      // Reference to the grid on which the field is created
      const Grid& grid_;
      // Pointer to the raw data.
      double* restrict data_;
  };
}
