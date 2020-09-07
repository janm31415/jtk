///////////////////////////////////////////////////////////////////////////////
//
// header only matrix library using expression templates
// matrices are stored row major
//
// Author    :  Jan Maes                                            
// Version   :  1.5
// Date      :  04 September 2020
// License   :  MIT License
//
///////////////////////////////////////////////////////////////////////////////

/*
Changelog:
V1.0: 24 September 2019
  - first version
V1.1: 09 January 2020
  - added evaluate_before_assigning flag to automatically avoid aliasing effects
  - added data() method for raw pointer access
V1.2: 31 January 2020
  - added trace()
V1.3: 30 April 2020
  - allowing svd when m < n
  - concatenation of matrices with << operator in horizontal or vertical direction
V1.4: 30 June 2020
  - adding zeros, ones, identity in templatized form, e.g. zeros<float>(m, n).
V1.5: 04 September 2020
  - adding sparse matrices.
*/


#pragma once

#include <vector>
#include <sstream>
#include <array>
#include <algorithm>
#include <stdint.h>
#include <cassert>
#include <ostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <immintrin.h>

#include "timer.h"

namespace jtk
  {

  ///////////////////////////////////////////////////////////////////////////////
  // forward declarations of useful classes and functions
  ///////////////////////////////////////////////////////////////////////////////

  /*
  the basic matrix class. template parameter T should be a numeric type such as float or double.
  template parameter Container should be std::vector<T> for dynamic allocation or std::array<T, size>
  for small sized matrices.
  */
  template <class T, class Container>
  class matrix;

  /*
  the sparse matrix class. template parameter T should be a numeric type such as float or double.
  */
  template <class T>
  class sparse_matrix;

  /*
  given a matrix a, this routine computes its singular value
  decomposition, a = u.w.vt.  the matrix u replaces a on output. the diagonal
  matrix of singular values is output as a vector w. the matrix v (not
  the transpose vt) is output as well.
  */
  template <class T, class Container_mn, class Container_n, class Container_nn>
  bool svd(matrix<T, Container_mn>& a, matrix<T, Container_n>& w, matrix<T, Container_nn>& v);

  /*
  int p = pseudo_inverse(x, a) produces a matrix x of the same dimensions
  as a so that a*x*a = a, x*a*x = x and a*x and x*a
  are hermitian. the computation is based on the svd of a and any
  singular values less than a tolerance are treated as zero.
  the return value p represents the number of non zero singular values.
  */
  template <class T, class Container>
  int pseudo_inverse(matrix<T, Container>& Ainv, matrix<T, Container>& A, T tol);

  /*
  lsd computes the minimum-norm solution to a overdetermined system:
  minimize 2-norm(| b - a*x |)
  using the singular value decomposition (svd) of a.
  the matrix a is destroyed. the output x is stored in vector b.
  */
  template <class T, class Container, class Container2>
  void lsd(matrix<T, Container>& a, matrix<T, Container2>& b, T tol);

  /*
  Given a matrix a, this routine replaces it by the lu decomposition of a rowwise
  permutation of itself. a is input and output.
  permutations is an output vector that records the row permutation effected by the partial
  pivoting; d is output as +-1 depending on whether the number of row interchanges was even
  or odd, respectively. this routine is used in combination with lubksb to solve linear equations
  or invert a matrix.
  */
  template <class T, class Container>
  void ludcmp(matrix<T, Container>& a, std::vector<uint64_t>& permutations, T& d);

  /*
  lubksb solves the set of n linear equations a.x = b. here a is input, not as the matrix
  a but rather as its lu decomposition, determined by the routine ludcmp. permutations is input
  as the permutation vector returned by ludcmp. b is input as the right-hand side vector
  and returns with the solution vector x. a and permutations are not modified by this routine
  and can be left in place for successive calls with different right-hand sides b. This routine takes
  into account the possibility that b will begin with many zero elements, so it is effcient for use
  in matrix inversion.
  */
  template <class T, class Container, class Container2>
  void lubksb(matrix<T, Container2>& b, const matrix<T, Container>& a, const std::vector<uint64_t>& permutations);

  /*
  this method solves a*x = b using lu decomposition.
  a is destroyed, and b is replaced by the solution x.
  */
  template <class T, class Container, class Container2>
  void solve(matrix<T, Container>& a, matrix<T, Container2>& b);

  /*
  This method computes the inverse of a matrix a. it uses the lu decomposition in the background.
  */
  template <class T, class Container>
  void invert(matrix<T, Container>& ainv, const matrix<T, Container>& a);

  /*
  This method computes the determinant of matrix a. it uses the lu decomposition in the background.
  */
  template <class T, class Container>
  T determinant(const matrix<T, Container>& a);

  /*
  computes cholesky factorizations of square hermitian positive definite matrices.
  the algorithm only uses the superdiagonal half of the input matrix.
  the input matrix is replaced by the output of the cholesky algorithm.
  as output we get the upper-triangular factor r of the cholesky factorization
  so that a = (r*)*(r), with a the input matrix.
  */
  template <class T, class Container>
  void cholesky(matrix<T, Container>& a);

  /*
  we solve ax = b where a is upper triangular.
  */
  template <class T, class Container, class Container2>
  void back_substitution(matrix<T, Container2>& x, const matrix<T, Container2>& b, const matrix<T, Container>& a);

  /*
  we solve ax = b where a is lower triangular.
  */
  template <class T, class Container, class Container2>
  void forward_substitution(matrix<T, Container2>& x, const matrix<T, Container2>& b, const matrix<T, Container>& a);

  /*
  This method solves a*x = b using cholesky decomposition.
  a is destroyed, and b is replaced by the solution x.
  */
  template <class T, class Container, class Container2>
  void solve_cholesky(matrix<T, Container>& a, matrix<T, Container2>& b);

  /*
  Constructs the qr decomposition of a rectangular matrix a with or without pivoting, i.e.
  a*p = q*r
  the upper triangular matrix r is returned in a. The matrix qt is the transpose of the orthogonal
  matrix q. if pivot is true, then the permutation matrix p can be constructed by making its j-th column equal to the
  permutations[j]-th column of the identity matrix. this method calls qrfac.
  */
  template <class T, class Container>
  void qrdcmp(matrix<T, Container>& A, matrix<T, Container>& QT, std::vector<uint64_t>& permutations, bool pivot);

  /*
  this method solves a*x = b using qr decomposition with pivoting.
  a is destroyed, and b is replaced by the solution x.
  the system can be overdetermined in which case it is solved in least squares sense.
  this method calls qrfac and qrsolv.
  */
  template <class T, class Container, class Container2>
  void solve_qr(matrix<T, Container>& a, matrix<T, Container2>& b);

  /*
  Returns the l2 norm for a vector and the frobenius norm for a matrix.
  */
  template <class T, class Container >
  T norm(const matrix<T, Container>& a);

  /*
  Returns the trace of the matrix.
  */
  template <class T, class Container >
  T trace(const matrix<T, Container>& a);

  /*
  constructs the qr decomposition of a rectangular mxn matrix a with pivoting, i.e. a*p = q*r.
  on output, a contains the strict upper trapezoidal part of r, while rdiag contains the
  diagonal elements of r. the remaining part of a contains a factored form of q, i.e. the
  non-trivial elements of the u vectors, where the u vectors represent the householder
  transformation of q which, for column k (k = 1,2,...,min(m,n)), is of the form
  i - (1/u(k))*u*ut where u has zeros in the first k-1 positions. with the boolean pivot you
  can switch on/off pivoting. the permutations vector is non-empty when pivot is switched on.
  the permutation matrix p can be constructed by making its j-th column equal to the
  permutations[j]-th column of the identity matrix. acnorm is an output vector of length
  n which contains the norms of the corresponding columns of the input matrix a.
  */
  template <class T, class Container, class Container2>
  void qrfac(matrix<T, Container>& a, bool pivot, std::vector<uint64_t>& permutations,
    matrix<T, Container2>& rdiag, matrix<T, Container2>& acnorm);

  /*
  given an mxn matrix a, and nxn diagonal matrix d, and an mx1 vector b, the system
  a*x = b, d*x = 0
  is solved in the least squares sense.
  the input to this method is obtained from the output of method qrfac.
  r is the nxn matrix that represents the full upper triangle of the qr decomposition of matrix
  a. permutations is the output permutations vector of method qrfac. diag is an input vector
  of length n which must contain the diagonal elements of the matrix d. qtb is an input
  vector of length n which contains the first n elements of the vector (q transpose)*b.
  x is the output vector that represents the least squares solution to the system
  a*x = b, d*x = 0. sdiag is an output vector of length n which contains the diagonal
  elements of the upper triangular matrix s, where s is defined by
  pT*(aT*a + d*d)*p = sT*s
  */
  template <class T, class Container, class Container2>
  void qrsolv(matrix<T, Container>& r, const std::vector<uint64_t>& permutations, const matrix<T, Container2>& diag,
    const matrix<T, Container2>& qtb, matrix<T, Container2>& x, matrix<T, Container2>& sdiag);

  /*
  this method computes a forward-difference approximation to the mxn jacobian matrix associated
  with a with a specified problem of m functions in n variables.
  bool(*f)(const matrix<T, Container>& x, matrix<T, Container> fvec&, void* user_data) is a function pointer that is provided
  by the user. The first argument is a nx1 vector and represents the position at which the functions should be evaluated.
  The second argument is an output mx1 vector containing the values of the m functions evaluated at the input position. The third
  argument can be used to pass user specific data. If the user function returns false, all operations are aborted and
  fdjac will return false (otherwise true).
  x is a vector of length n and fvec is a vector of length m containing the functions evaluated at x. fjac is the output mxn
  matrix representing the jacobian and epsfcn is an input value representing a suitable step size for the finite difference
  computation. user_data is optional and is intented to provide user specific data to function pointer f.
  */
  template <class T, class Container_mn, class Container_n, class Container_m>
  bool fdjac(bool(*f)(const matrix<T, Container_n>&, matrix<T, Container_m>&, void*), matrix<T, Container_n>& x, matrix<T, Container_m>& fvec, matrix<T, Container_mn>& fjac, T epsfcn, void* user_data = nullptr);

  /*
  levenberg-marquardt optimization.
  the method lmdif is intended to minimize the sum of squares of m nonlinear functions
  in n variables. the user must provide a function pointer which calculates the functions.
  the jacobian is then calculated by a forward-difference approximation.
  f is a function pointer that has the declaration
  bool f(const matrix<T, Container>& pos, matrix<T, Container>& fvec, void* user_data) where
  pos is an input vector of length n, and fvec is an output vector of length m representing
  the error functions (sum of squares) of m nonlinear functions evaluated in pos. any extra data needed
  to compute the output vector fvec can be provided with the void pointer user_data.
  If this function returns false, lmdif will stop its computations and return.
  x is a vector of length n. On input x must contain an initial estimate of the solution. On
  output x contains the final estimate of the solution.
  fvec is an output vector of length m containing the functions evaluated at output x.
  ftol is a nonnegative input variable. termination occurs when both the actual and predicated
  relative reductions in the sum of squares are at most ftol. therefore, ftol measures the relative
  error desired in the sum of squares.
  xtol is a nonnegative input variable. termination occurs when the relative error between two consecutive
  iterates is at most xtol. therefore, xtol measures the relative error desired in the approximate solution.
  gtol is a nonnegative input variable. termination occurs when the cosine of the angle between fvec and
  any column of the jacobian is at most gtol in absolute value. therefore, gtol measures the orthogonality
  desired between the function vector and the columns of the jacobian.
  maxfev is a positive integer input variable. termination occurs when the number of calls to fcn is at least
  maxfev by the end of an iteration.
  epsfcn is an input variable used in determining a suitable step length for the forward-difference approximation. this
  approximation assumes that the relative errors in the functions are of the order of epsfcn. if epsfcn is less
  than the machine precision, it is assumed that the relative errors in the functions are of the order of the machine
  precision.
  diag is an array of length n. if mode = 1 (see below), diag is internally set. if mode = 2, diag
  must contain positive entries that serve as multiplicative scale factors for the variables.
  mode is an integer input variable. if mode = 1, the variables will be scaled internally. if mode = 2,
  the scaling is specified by the input diag. other values of mode are equivalent to mode = 1.
  factor is a positive input variable used in determining the initial step bound. this bound is set to the product of
  factor and the euclidean norm of diag*x if nonzero, or else to factor itself. in most cases factor should lie in the
  interval (.1,100.). 100. is a generally recommended value.
  info is an integer output variable. if the user has terminated execution, info is set to the (negative)
  value of iflag. see description of fcn. otherwise, info is set as follows.

    info = 0  improper input parameters.

    info = 1  both actual and predicted relative reductions
              in the sum of squares are at most ftol.

    info = 2  relative error between two consecutive iterates
              is at most xtol.

    info = 3  conditions for info = 1 and info = 2 both hold.

    info = 4  the cosine of the angle between fvec and any
              column of the jacobian is at most gtol in
              absolute value.

    info = 5  number of calls to fcn has reached or
              exceeded maxfev.

    info = 6  ftol is too small. no further reduction in
              the sum of squares is possible.

    info = 7  xtol is too small. no further improvement in
              the approximate solution x is possible.

    info = 8  gtol is too small. fvec is orthogonal to the
              columns of the jacobian to machine precision.

  nfev is an integer output variable set to the number of calls to f.

  user_data is optional and is intented to provide user specific data to function pointer f.

  this method calls qrfac and qrsolv.
  */
  template <class T, class Container_n, class Container_m>
  void lmdif(bool(*f)(const matrix<T, Container_n>&, matrix<T, Container_m>&, void*), matrix<T, Container_n>& x, matrix<T, Container_m>& fvec, T ftol,
    T xtol, T gtol, uint64_t maxfev, T epsfcn, matrix<T, Container_n>& diag,
    uint64_t mode, T factor, uint64_t& info, uint64_t& nfev, void* user_data = nullptr);

  /*
  levenberg-marquardt optimization.
  the method lmdif0 is intended to minimize the sum of squares of m nonlinear functions
  in n variables. the user must provide a function pointer which calculates the functions.
  the jacobian is then calculated by a forward-difference approximation.
  the parameters f, x, fvec, info, nfev, and user_data are as in method lmdif.
  this method calls lmdif.
  */
  template <class T, class Container_n, class Container_m>
  void lmdif0(bool(*f)(const matrix<T, Container_n>&, matrix<T, Container_m>&, void*), matrix<T, Container_n>& x, matrix<T, Container_m>& fvec, uint64_t& info, uint64_t& nfev, void* user_data = nullptr);

  /*
  householder reduction of a symmetric matrix a. on output, a is replaced
  by the orthogonal matrix q effecting the transformation. diagonal returns the diagonal elements
  of the tridiagonal matrix, and subdiagonal the off-diagonal elements, with subdiagonal(0)=0.
  */
  template <class T, class Container>
  void tred2(matrix<T, Container>& diagonal, matrix<T, Container>& subdiagonal, matrix<T, Container>& a);

  /*
  ql algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric,
  tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2. on
  input, diagonal contains the diagonal elements of the tridiagonal matrix. on output, diagonal returns
  the eigenvalues. the vector subdiagonal inputs the subdiagonal elements of the tridiagonal matrix,
  with subdiagonal(0) arbitrary. on output subdiagonal is destroyed.
  if the eigenvectors of a tridiagonal matrix are desired,
  the matrix a is input as the identity matrix. if the eigenvectors of a matrix
  that has been reduced by tred2 are required, then a is input as the matrix output by tred2.
  in either case, the kth column of a returns the normalized eigenvector corresponding to diagonal(k)
  */
  template <class T, class Container>
  bool tqli(matrix<T, Container>& diagonal, matrix<T, Container>& a, matrix<T, Container>& subdiagonal);

  /*
  computes the eigenvalues and eigenvectors of a symmetric matrix a.
  eigenvalues are stored in eigenvalues. eigenvectors are stored in a.
  */
  template <class T, class Container>
  bool eig_symm(matrix<T, Container>& eigenvalues, matrix<T, Container>& a);

  /*
  given a matrix a, this routine replaces it by a balanced matrix with identical
  eigenvalues. a symmetric matrix is already balanced and is unaffected by this procedure.
  */
  template <class T, class Container>
  void balanc(matrix<T, Container>& a);

  /*
  reduction to hessenberg form by the elimination method. the real, nonsymmetric matrix a
  is replaced by an upper hessenberg matrix with identical eigenvalues.
  recommended, but not required, is that this routine be preceded by balanc. on output, the
  hessenberg matrix is in elements a[i][j] with i <= j+1. elements with i > j+1 are to be
  thought of as zero, but are returned with random values.
  */
  template <class T, class Container>
  void elmhes(matrix<T, Container>& a);

  /*
  finds all eigenvalues of an upper hessenberg matrix a. on input a can be
  exactly as output from elmhes, on output it is destroyed. the real and imaginary parts
  of the eigenvalues are returned in wr and wi, respectively.
  */
  template <class T, class Container>
  bool hqr(matrix<T, Container>& a, matrix<T, Container>&  wr, matrix<T, Container>&  wi);

  /*
  finds all eigenvalues of a square matrix a. on output a is destroyed. the real and imaginary parts
  of the eigenvalues are returned in wr and wi, respectively.
  */
  template <class T, class Container>
  bool eig(matrix<T, Container>& a, matrix<T, Container>&  wr, matrix<T, Container>&  wi);

  namespace implementation_details
    {
    template <class T, class Container>
    struct resize
      {
      void operator()(Container&, uint64_t) {}
      };

    template <class T>
    struct resize<T, std::vector<T>>
      {
      void operator()(std::vector<T>& c, uint64_t size)
        {
        c.resize(size, (T)0);
        }
      };

    template <class T>
    struct get_value_type
      {
      using value_type = typename T::value_type;
      };

    template <class T>
    struct get_value_type<T*>
      {
      using value_type = T;
      };

    template <class T>
    struct get_value_type<const T*>
      {
      using value_type = T;
      };

    }

  template <class T, class Container>
  struct CommaInitializer
    {

    CommaInitializer(matrix<T, Container>* m, const T& first_value) : _m(m)
      {
      _it = _m->begin();
      *_it = first_value;
      ++_it;
      }

    CommaInitializer& operator,(const T& i)
      {
      *_it = i;
      ++_it;
      return *this;
      }

    matrix<T, Container>* _m;
    typename matrix<T, Container>::iterator _it;
    };

  template <class T, class Container>
  struct CommaConcatenator
    {
    CommaConcatenator(matrix<T, Container>* m, const matrix<T, Container>& first_matrix) : _m(m)
      {
      assert(_m->rows() == first_matrix.rows() || _m->cols() == first_matrix.cols());
      vertical = _m->rows() > first_matrix.rows();
      if (vertical)
        {
        _it = _m->begin();
        for (auto v : first_matrix)
          *_it++ = v;
        }
      else
        {
        _it = _m->begin();
        auto tmp_it = _it;
        auto tgt_it = first_matrix.begin();
        for (int r = 0; r < _m->rows() - 1; ++r)
          {
          for (int c = 0; c < first_matrix.cols(); ++c)
            {
            *tmp_it++ = *tgt_it++;
            }
          tmp_it += _m->cols() - first_matrix.cols();
          }
        for (int c = 0; c < first_matrix.cols(); ++c)
          {
          *tmp_it++ = *tgt_it++;
          }

        _it += first_matrix.cols();
        }
      }

    CommaConcatenator& operator,(const matrix<T, Container>& i)
      {
      assert(vertical ? _m->cols() == i.cols() : _m->rows() == i.rows());
      if (vertical)
        {
        for (auto v : i)
          *_it++ = v;
        }
      else
        {
        auto tmp_it = _it;
        auto tgt_it = i.begin();
        for (int r = 0; r < _m->rows() - 1; ++r)
          {
          for (int c = 0; c < i.cols(); ++c)
            {
            *tmp_it++ = *tgt_it++;
            }
          tmp_it += _m->cols() - i.cols();
          }
        for (int c = 0; c < i.cols(); ++c)
          {
          *tmp_it++ = *tgt_it++;
          }
        _it += i.cols();
        }
      return *this;
      }

    matrix<T, Container>* _m;
    typename matrix<T, Container>::iterator _it;
    bool vertical;
    };

  template <class T1, class T2>
  struct gettype
    {
    typedef T1 ty;
    };

  template <>
  struct gettype<float, double>
    {
    typedef double ty;
    };

  template <class T>
  struct OpAdd
    {
    using value_type = T;
    template <class T2, class T3>
    static inline T apply(T2 a, T3 b)
      {
      return static_cast<T>(a + b);
      }
    };

  template <class T>
  struct OpSub
    {
    using value_type = T;
    template <class T2, class T3>
    static inline T apply(T2 a, T3 b)
      {
      return static_cast<T>(a - b);
      }
    };

  template <class T>
  struct OpMul
    {
    using value_type = T;
    template <class T2, class T3>
    static inline T apply(T2 a, T3 b)
      {
      return static_cast<T>(a * b);
      }
    };

  template <class T>
  struct OpDiv
    {
    using value_type = T;
    template <class T2, class T3>
    static inline T apply(T2 a, T3 b)
      {
      return static_cast<T>(a / b);
      }
    };

  template <class T>
  struct OpDivInv
    {
    using value_type = T;
    template <class T2, class T3>
    static inline T apply(T2 a, T3 b)
      {
      return static_cast<T>(b / a);
      }
    };

  template <class T>
  struct OpNeg
    {
    using value_type = T;
    template <class T2>
    static inline T apply(T2 a)
      {
      return static_cast<T>(-a);
      }
    };

  template <class T>
  struct OpId
    {
    using value_type = T;
    template <class T2>
    static inline T apply(T2 a)
      {
      return static_cast<T>(a);
      }
    };

  ///////////////////////////////////////////////////////////////////////////////
  // Expression templates
  ///////////////////////////////////////////////////////////////////////////////

  template <class A, class B, class Op>
  class BinExprOp
    {
    public:
      using value_type = typename Op::value_type;

      BinExprOp(const A& a, const B& b, uint64_t rows, uint64_t cols, bool eval_before_assigning) : _a(a), _b(b), _rows(rows), _cols(cols),
        _evaluate_before_assigning(eval_before_assigning)
        {
        }

      value_type operator * () const
        {
        return Op::apply<typename A::value_type, typename B::value_type>(*_a, *_b);
        }

      void operator++()
        {
        ++_a; ++_b;
        }

      BinExprOp& operator += (const uint64_t offset)
        {
        _a += offset;
        _b += offset;
        return *this;
        }

      BinExprOp operator + (const uint64_t offset) const
        {
        BinExprOp tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      B _b;
      uint64_t _rows, _cols;
      bool _evaluate_before_assigning;
    };

  template <class A, class B, class Op>
  class BinExprSparseOp
    {
    public:
      using value_type = typename Op::value_type;

      BinExprSparseOp(const A& a, const B& b, uint64_t rows, uint64_t cols, bool eval_before_assigning) : _a(a), _b(b), _rows(rows), _cols(cols),
        _evaluate_before_assigning(eval_before_assigning)
        {
        _set_entry_1();
        _set_entry_2();
        }

      BinExprSparseOp(const BinExprSparseOp&) = default;

      value_type operator * () const
        {
        if (iterator_1_equal_to_iterator_2())
          return Op::apply<typename A::value_type, typename B::value_type>(*_a, *_b);
        else if (iterator_1_smaller_than_iterator_2())
          return Op::apply<typename A::value_type, typename B::value_type>(*_a, static_cast<typename B::value_type>(0));
        else
          return Op::apply<typename A::value_type, typename B::value_type>(static_cast<A::value_type>(0), *_b);
        }

      void operator++()
        {
        if (iterator_1_equal_to_iterator_2())
          {
          ++_a;
          ++_b;
          _set_entry_1();
          _set_entry_2();
          }
        else if (iterator_1_smaller_than_iterator_2())
          {
          ++_a;
          _set_entry_1();
          }
        else
          {
          ++_b;
          _set_entry_2();
          }
        }

      BinExprSparseOp end() const
        {
        BinExprSparseOp out(*this);
        out._a = out._a.end();
        out._b = out._b.end();
        return out;
        }

      BinExprSparseOp row(uint64_t idx) const
        {
        BinExprSparseOp out(*this);
        out._a = out._a.row(idx);
        out._b = out._b.row(idx);
        out._first_entry_1 = idx;
        out._first_entry_2 = idx;
        out._second_entry_1 = out._a.second_entry();
        out._second_entry_2 = out._b.second_entry();
        return out;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

      uint64_t first_entry() const
        {
        if (iterator_1_smaller_than_iterator_2())
          return _first_entry_1;
        else
          return _first_entry_2;
        }

      uint64_t second_entry() const
        {
        if (iterator_1_smaller_than_iterator_2())
          return _second_entry_1;
        else
          return _second_entry_2;
        }

      bool operator == (const BinExprSparseOp& i_other) const
        {
        return (_a == i_other._a && _b == i_other._b);
        }

      bool operator != (const BinExprSparseOp& i_other) const
        {
        return !((*this) == i_other);
        }

    private:

      void _set_entry_1()
        {
        if (_a == _a.end())
          {
          _first_entry_1 = static_cast<size_t>(-1);
          _second_entry_1 = static_cast<size_t>(-1);
          }
        else
          {
          _first_entry_1 = _a.first_entry();
          _second_entry_1 = _a.second_entry();
          }
        }

      void _set_entry_2()
        {
        if (_b == _b.end())
          {
          _first_entry_2 = static_cast<size_t>(-1);
          _second_entry_2 = static_cast<size_t>(-1);
          }
        else
          {
          _first_entry_2 = _b.first_entry();
          _second_entry_2 = _b.second_entry();
          }
        }

      bool iterator_1_equal_to_iterator_2() const
        {
        return (_first_entry_1 == _first_entry_2 &&
          _second_entry_1 == _second_entry_2);
        }

      bool iterator_1_smaller_than_iterator_2() const
        {
        if (_first_entry_1 < _first_entry_2)
          return true;
        else if (_first_entry_1 == _first_entry_2)
          {
          return (_second_entry_1 < _second_entry_2);
          }
        return false;
        }

      bool iterator_2_smaller_than_iterator_1() const
        {
        if (_first_entry_2 < _first_entry_1)
          return true;
        else if (_first_entry_1 == _first_entry_2)
          {
          return (_second_entry_2 < _second_entry_1);
          }
        return false;
        }


    private:
      A _a;
      B _b;
      uint64_t _rows, _cols;
      uint64_t _first_entry_1, _first_entry_2, _second_entry_1, _second_entry_2;
      bool _evaluate_before_assigning;
    };

  template <class A, class B>
  class MatMatMul
    {
    public:
      using value_type = typename gettype<typename ::jtk::implementation_details::get_value_type<A>::value_type, typename ::jtk::implementation_details::get_value_type<B>::value_type>::ty;

      MatMatMul(const A& a, const B& b, uint64_t rows, uint64_t cols, uint64_t mid_dim) : _a(a), _b(b), _rows(rows), _cols(cols), _mid_dim(mid_dim),
        _index(0), _evaluate_before_assigning(true)
        {
        }

      value_type operator * () const
        {
        const uint64_t r = _index / _cols;
        const uint64_t c = _index % _cols;
        auto a_it = _a + r * _mid_dim;
        auto b_it = _b + c;
        value_type res = (value_type)(*a_it) * (value_type)(*b_it);
        for (uint64_t k = 1; k < _mid_dim; ++k)
          {
          ++a_it;
          b_it += _cols;
          res += (value_type)(*a_it) * (value_type)(*b_it);
          }
        return res;
        }

      void operator++()
        {
        ++_index;
        }

      MatMatMul& operator += (const uint64_t offset)
        {
        _index += offset;
        return *this;
        }

      MatMatMul operator + (const uint64_t offset) const
        {
        MatMatMul tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      B _b;
      uint64_t _rows, _cols, _mid_dim;
      uint64_t _index;
      bool _evaluate_before_assigning;
    };

  template <class A, class B>
  class SparseMatMatMul
    {
    public:
      using value_type = typename gettype<typename ::jtk::implementation_details::get_value_type<A>::value_type, typename ::jtk::implementation_details::get_value_type<B>::value_type>::ty;

      SparseMatMatMul(const A& a, const B& b, uint64_t rows, uint64_t cols, uint64_t mid_dim) : _a(a), _b(b), _rows(rows), _cols(cols), _mid_dim(mid_dim),
        _index(0), _evaluate_before_assigning(true)
        {
        }

      value_type operator * () const
        {
        const uint64_t r = _index / _cols;
        const uint64_t c = _index % _cols;
        value_type res = (value_type)0;
        auto a_it = _a.row(r);
        for (; a_it.first_entry() == r; ++a_it)
          {
          res += static_cast<value_type>(*a_it) * static_cast<value_type>(*(_b + a_it.second_entry() * _cols + c));
          }
        return res;
        }

      void operator++()
        {
        ++_index;
        }

      SparseMatMatMul& operator += (const uint64_t offset)
        {
        _index += offset;
        return *this;
        }

      SparseMatMatMul operator + (const uint64_t offset) const
        {
        MatMatMul tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      B _b;
      uint64_t _rows, _cols, _mid_dim;
      uint64_t _index;
      bool _evaluate_before_assigning;
    };

  template <class A>
  class Transpose
    {
    public:
      using value_type = typename ::jtk::implementation_details::get_value_type<A>::value_type;

      Transpose(const A& a, uint64_t rows, uint64_t cols) : _a(a), _rows(rows), _cols(cols), _index(0), _evaluate_before_assigning(true)
        {
        }

      value_type operator * () const
        {
        const uint64_t r = _index / _cols;
        const uint64_t c = _index % _cols;
        uint64_t offset = c * _rows + r;
        return *(_a + offset);
        }

      void operator++()
        {
        ++_index;
        }

      Transpose& operator += (const uint64_t offset)
        {
        _index += offset;
        return *this;
        }

      Transpose operator + (const uint64_t offset) const
        {
        Transpose tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      uint64_t _rows, _cols;
      uint64_t _index;
      bool _evaluate_before_assigning;
    };

  template <class A>
  class Diagonal
    {
    public:
      using value_type = typename ::jtk::implementation_details::get_value_type<A>::value_type;

      Diagonal(const A& a, uint64_t rows, uint64_t cols, bool eval_before_assigning) : _a(a), _rows(rows), _cols(cols), _evaluate_before_assigning(eval_before_assigning)
        {
        _dim = std::min<uint64_t>(_rows, _cols);
        }

      value_type operator * () const
        {
        return *_a;
        }

      void operator++()
        {
        _a += (_cols + 1);
        }

      Diagonal& operator += (const uint64_t offset)
        {
        _a += (_cols + 1)*offset;
        return *this;
        }

      Diagonal operator + (const uint64_t offset) const
        {
        Diagonal tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _dim;
        }

      uint64_t cols() const
        {
        return 1;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      uint64_t _rows, _cols;
      uint64_t _dim;
      bool _evaluate_before_assigning;
    };

  template <class A>
  class DiagonalSparse
    {
    public:
      using value_type = typename ::jtk::implementation_details::get_value_type<A>::value_type;

      DiagonalSparse(const A& a, uint64_t rows, uint64_t cols, bool eval_before_assigning) : _a(a), _rows(rows), _cols(cols), _index(0), _evaluate_before_assigning(eval_before_assigning)
        {
        _dim = std::min<uint64_t>(_rows, _cols);
        }

      value_type operator * () const
        {
        auto row_it = _a.row(_index);
        while (row_it.first_entry() == _index)
          {
          if (row_it.second_entry() == _index)
            return *row_it;
          ++row_it;
          }
        return (value_type)0;
        }

      void operator++()
        {
        ++_index;
        }

      DiagonalSparse& operator += (const uint64_t offset)
        {
        _index += offset;
        return *this;
        }

      DiagonalSparse operator + (const uint64_t offset) const
        {
        Diagonal tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _dim;
        }

      uint64_t cols() const
        {
        return 1;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      uint64_t _rows, _cols;
      uint64_t _dim;
      uint64_t _index;
      bool _evaluate_before_assigning;
    };

  template <class A>
  class Block
    {
    public:
      using value_type = typename ::jtk::implementation_details::get_value_type<A>::value_type;

      Block(const A& a, uint64_t pos_r, uint64_t pos_c, uint64_t rows, uint64_t cols, uint64_t a_rows, uint64_t a_cols, bool eval_before_assigning) : _a(a), _pos_r(pos_r),
        _pos_c(pos_c), _rows(rows), _cols(cols), _a_rows(a_rows), _a_cols(a_cols), _evaluate_before_assigning(eval_before_assigning)
        {
        _current_r = _pos_r;
        _current_c = _pos_c;
        const uint64_t index = _current_r * _a_cols + _current_c;
        if (rows && cols) // block is not empty
          _a += index;
        }

      value_type operator * () const
        {
        return *_a;
        }

      void operator++()
        {
        ++_current_c;
        if (_current_c >= _pos_c + _cols)
          {
          ++_current_r;
          _current_c = _pos_c;
          _a += (_a_cols - _cols + 1);
          }
        else
          ++_a;
        }

      Block& operator += (const uint64_t offset)
        {
        uint64_t next_rows = offset / _cols;
        uint64_t remainder = offset % _cols;
        _a += next_rows * _a_cols;
        _current_r += next_rows;
        _current_c += remainder;
        if (_current_c >= _pos_c + _cols)
          {
          ++_current_r;
          _current_c -= _cols;
          _a += (_a_cols - _cols + (1 + _current_c - _pos_c));
          }
        else
          _a += remainder;
        return *this;
        }

      Block operator + (const uint64_t offset) const
        {
        Block tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      uint64_t _pos_r, _pos_c;
      uint64_t _rows, _cols;
      uint64_t _a_rows, _a_cols;
      uint64_t _current_r, _current_c;
      bool _evaluate_before_assigning;
    };

  template <class A>
  class BlockSparse
    {
    public:
      using value_type = typename ::jtk::implementation_details::get_value_type<A>::value_type;

      BlockSparse(const A& a, uint64_t pos_r, uint64_t pos_c, uint64_t rows, uint64_t cols, bool eval_before_assigning) : _a(a), _pos_r(pos_r),
        _pos_c(pos_c), _rows(rows), _cols(cols), _evaluate_before_assigning(eval_before_assigning)
        {
        _find_next_entry();
        }

      value_type operator * () const
        {
        return *_a;
        }

      void operator++()
        {
        ++_a;
        _find_next_entry();
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

      BlockSparse end() const
        {
        BlockSparse out(*this);
        out._a = out._a.end();
        return out;
        }

      BlockSparse row(uint64_t idx) const
        {
        BlockSparse out(*this);
        out._a = out._a.row(idx);
        return out;
        }

      uint64_t first_entry() const
        {
        return _a.first_entry() - _pos_r;
        }

      uint64_t second_entry() const
        {
        return _a.second_entry() - _pos_c;
        }

      bool operator == (const BlockSparse& other) const
        {
        return (_a == other._a);
        }

      bool operator != (const BlockSparse& other) const
        {
        return !((*this) == other);
        }

    private:
      void _find_next_entry()
        {
        while ((_a != _a.end()) && (_a.first_entry() < _pos_r || _a.first_entry() >= _pos_r + _rows || _a.second_entry() < _pos_c || _a.second_entry() >= _pos_c + _cols))
          {
          ++_a;
          }
        }

    private:
      A _a;
      uint64_t _pos_r, _pos_c;
      uint64_t _rows, _cols;
      bool _evaluate_before_assigning;
    };

  template <class T>
  class Constant
    {
    public:
      using value_type = T;

      Constant(T value, uint64_t rows, uint64_t cols) : _value(value), _rows(rows), _cols(cols), _evaluate_before_assigning(false)
        {
        }

      value_type operator * () const
        {
        return _value;
        }

      void operator++()
        {
        }

      Constant& operator += (const uint64_t)
        {
        return *this;
        }

      Constant operator + (const uint64_t) const
        {
        return *this;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      T _value;
      uint64_t _rows, _cols;
      bool _evaluate_before_assigning;
    };


  template <class T>
  class Identity
    {
    public:
      using value_type = T;

      Identity(uint64_t rows, uint64_t cols) : _rows(rows), _cols(cols), _index(0), _evaluate_before_assigning(false)
        {
        }

      value_type operator * () const
        {
        const uint64_t r = _index / _cols;
        const uint64_t c = _index % _cols;
        return r == c ? (value_type)1 : (value_type)0;
        }

      void operator++()
        {
        ++_index;
        }

      Identity& operator += (const uint64_t offset)
        {
        _index += offset;
        return *this;
        }

      Identity operator + (const uint64_t offset) const
        {
        Identity tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      uint64_t _rows, _cols;
      uint64_t _index;
      bool _evaluate_before_assigning;
    };

  template <class A, class Op>
  class UnExprOp
    {
    public:
      using value_type = typename Op::value_type;

      UnExprOp(const A& a, uint64_t rows, uint64_t cols, bool eval_before_assigning) : _a(a), _rows(rows), _cols(cols), _evaluate_before_assigning(eval_before_assigning)
        {
        }

      value_type operator * () const
        {
        return Op::apply(*_a);
        }

      void operator++()
        {
        ++_a;
        }

      UnExprOp& operator += (const uint64_t offset)
        {
        _a += offset;
        return *this;
        }

      UnExprOp operator + (const uint64_t offset) const
        {
        UnExprOp tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      A _a;
      uint64_t _rows, _cols;
      bool _evaluate_before_assigning;
    };

  template <class A, class Op>
  class UnExprSparseOp
    {
    public:
      using value_type = typename Op::value_type;

      UnExprSparseOp(const A& a, uint64_t rows, uint64_t cols, bool eval_before_assigning) : _a(a), _rows(rows), _cols(cols), _evaluate_before_assigning(eval_before_assigning)
        {
        }

      UnExprSparseOp(const UnExprSparseOp&) = default;

      value_type operator * () const
        {
        return Op::apply(*_a);
        }

      void operator++()
        {
        ++_a;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

      UnExprSparseOp end() const
        {
        UnExprSparseOp out(*this);
        out._a = out._a.end();
        return out;
        }

      UnExprSparseOp row(uint64_t idx) const
        {
        UnExprSparseOp out(*this);
        out._a = out._a.row(idx);
        return out;
        }

      uint64_t first_entry() const
        {
        return _a.first_entry();
        }

      uint64_t second_entry() const
        {
        return _a.second_entry();
        }

      bool operator == (const UnExprSparseOp& other) const
        {
        return (_a == other._a);
        }

      bool operator != (const UnExprSparseOp& other) const
        {
        return !((*this) == other);
        }

    private:
      A _a;
      uint64_t _rows, _cols;
      bool _evaluate_before_assigning;
    };

  template <class A, class Op>
  class BinExprSparseScalarOp
    {
    public:
      using value_type = typename Op::value_type;

      BinExprSparseScalarOp(const A& a, value_type scalar, uint64_t rows, uint64_t cols, bool evaluate_before_assigning) : _a(a), _scalar(scalar), _rows(rows), _cols(cols), _evaluate_before_assigning(evaluate_before_assigning)
        {
        }

      BinExprSparseScalarOp(const BinExprSparseScalarOp&) = default;

      value_type operator * () const
        {
        return Op::apply<typename A::value_type, value_type>(*_a, _scalar);
        }

      void operator++()
        {
        ++_a;
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

      BinExprSparseScalarOp end() const
        {
        BinExprSparseScalarOp out(*this);
        out._a = out._a.end();
        return out;
        }

      BinExprSparseScalarOp row(uint64_t idx) const
        {
        BinExprSparseScalarOp out(*this);
        out._a = out._a.row(idx);
        return out;
        }

      uint64_t first_entry() const
        {
        return _a.first_entry();
        }

      uint64_t second_entry() const
        {
        return _a.second_entry();
        }

      bool operator == (const BinExprSparseScalarOp& other) const
        {
        return (_a == other._a);
        }

      bool operator != (const BinExprSparseScalarOp& other) const
        {
        return !((*this) == other);
        }

    private:
      A _a;
      value_type _scalar;
      uint64_t _rows, _cols;
      bool _evaluate_before_assigning;
    };

  template <class ExprOp>
  class Expr
    {
    public:
      using value_type = typename ExprOp::value_type;

      Expr(const ExprOp& expr_op) : _expr_op(expr_op), _evaluate_before_assigning(expr_op.evaluate_before_assigning())
        {}

      value_type operator * () const
        {
        return *_expr_op;
        }

      void operator++()
        {
        ++_expr_op;
        }

      Expr& operator += (const uint64_t offset)
        {
        _expr_op += offset;
        return *this;
        }

      Expr operator + (const uint64_t offset) const
        {
        Expr tmp = *this;
        tmp += offset;
        return tmp;
        }

      uint64_t rows() const
        {
        return _expr_op.rows();
        }

      uint64_t cols() const
        {
        return _expr_op.cols();
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

    private:
      ExprOp _expr_op;
      bool _evaluate_before_assigning;
    };

  template <class ExprOp>
  class SparseExpr
    {
    public:
      using value_type = typename ExprOp::value_type;

      SparseExpr(const ExprOp& expr_op) : _expr_op(expr_op), _evaluate_before_assigning(expr_op.evaluate_before_assigning())
        {}

      value_type operator * () const
        {
        return *_expr_op;
        }

      void operator++()
        {
        ++_expr_op;
        }

      uint64_t rows() const
        {
        return _expr_op.rows();
        }

      uint64_t cols() const
        {
        return _expr_op.cols();
        }

      bool evaluate_before_assigning() const
        {
        return _evaluate_before_assigning;
        }

      SparseExpr end() const
        {
        SparseExpr out(*this);
        out._expr_op = out._expr_op.end();
        return out;
        }

      SparseExpr row(uint64_t idx) const
        {
        SparseExpr out(*this);
        out._expr_op = out._expr_op.row(idx);
        return out;
        }

      uint64_t first_entry() const
        {
        return _expr_op.first_entry();
        }

      uint64_t second_entry() const
        {
        return _expr_op.second_entry();
        }

      bool operator == (const SparseExpr& other) const
        {
        return (_expr_op == other._expr_op);
        }

      bool operator != (const SparseExpr& other) const
        {
        return !((*this) == other);
        }

    private:
      ExprOp _expr_op;
      bool _evaluate_before_assigning;
    };

  template <class T, class Container>
  class NoAlias
    {
    public:
      explicit NoAlias(matrix<T, Container>& m) : _matrix(m)
        {
        }

      template <class ExprOp>
      matrix<T, Container>& operator = (const Expr<ExprOp>& result)
        {
        _matrix.assign_no_alias(result);
        return _matrix;
        }

      template <class ExprOp>
      matrix<T, Container>& operator += (const Expr<ExprOp>& result)
        {
        _matrix.add_assign_no_alias(result);
        return _matrix;
        }

      template <class ExprOp>
      matrix<T, Container>& operator -= (const Expr<ExprOp>& result)
        {
        _matrix.subtract_assign_no_alias(result);
        return _matrix;
        }

    private:
      matrix<T, Container>& _matrix;
    };

  ///////////////////////////////////////////////////////////////////////////////
  // matrix class
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container = std::vector<T> >
  class matrix
    {
    public:
      using value_type = T;
      using pointer = T * ;
      using reference = T & ;
      using const_pointer = const T *;
      using const_reference = const T &;
      using iterator = typename Container::iterator;
      using const_iterator = typename Container::const_iterator;

      matrix() : _rows(0), _cols(0)
        {
        }

      matrix(uint64_t rows, uint64_t cols) : _rows(rows), _cols(cols)
        {
        ::jtk::implementation_details::resize<T, Container> r;
        r(_entries, _rows*_cols);
        }

      matrix(uint64_t rows) : _rows(rows), _cols(1)
        {
        ::jtk::implementation_details::resize<T, Container> r;
        r(_entries, _rows);
        }

      void swap(matrix& other)
        {
        std::swap(_cols, other._cols);
        std::swap(_rows, other._rows);
        std::swap(_entries, other._entries);
        }

      matrix(const matrix& other) : _rows(other._rows), _cols(other._cols), _entries(other._entries)
        {
        }

      matrix(matrix&&) = default;

      matrix& operator = (const matrix& other)
        {
        matrix temp(other);
        swap(temp);
        return *this;
        }

      matrix& operator = (matrix&&) = default;

      template <class T2, class Container2>
      matrix(const matrix<T2, Container2>& other) : _rows(other.rows()), _cols(other.cols())
        {
        resize(_rows, _cols);
        auto it = _entries.begin();
        auto other_it = other.begin();
        const auto other_it_end = other.end();
        for (; other_it != other_it_end; ++other_it)
          {
          *it = *other_it;
          ++it;
          }
        }

      template <class T2, class Container2>
      matrix& operator = (const matrix<T2, Container2>& other)
        {
        matrix temp(other);
        swap(temp);
        return *this;
        }

      template <class ExprOp>
      matrix(const Expr<ExprOp>& result)
        {
        assign(result);
        }

      template <class ExprOp>
      matrix& operator = (const Expr<ExprOp>& result)
        {
        assign(result);
        return *this;
        }

      template <class ExprOp>
      void assign(const Expr<ExprOp>& result)
        {
        if (result.evaluate_before_assigning())
          assign_alias(result);
        else
          assign_no_alias(result);
        }

      template <class ExprOp>
      void assign_alias(Expr<ExprOp> result)
        {
        matrix temp(result.rows(), result.cols());
        auto it = temp._entries.begin();
        uint64_t sz = result.rows()*result.cols();
        *it = (T)*result;
        for (uint64_t i = 1; i < sz; ++i)
          {
          ++result;
          ++it;
          *it = (T)*result;
          }
        swap(temp);
        }

      template <class ExprOp>
      void assign_no_alias(Expr<ExprOp> result)
        {
        resize(result.rows(), result.cols());
        auto it = _entries.begin();
        uint64_t sz = result.rows()*result.cols();
        *it = (T)*result;
        for (uint64_t i = 1; i < sz; ++i)
          {
          ++result;
          ++it;
          *it = (T)*result;
          }
        }

      template <class T2>
      matrix(const sparse_matrix<T2>& m)
        {
        UnExprSparseOp<typename sparse_matrix<T2>::const_iterator, OpId<T2>> expr(m.begin(), m.rows(), m.cols(), false);
        assign_no_alias(SparseExpr<UnExprSparseOp<typename sparse_matrix<T2>::const_iterator, OpId<T2>>>(expr));
        }

      template <class T2>
      matrix& operator = (const sparse_matrix<T2>& m)
        {
        UnExprSparseOp<typename sparse_matrix<T2>::const_iterator, OpId<T2>> expr(m.begin(), m.rows(), m.cols(), false);
        assign_no_alias(SparseExpr<UnExprSparseOp<typename sparse_matrix<T2>::const_iterator, OpId<T2>>>(expr));
        return *this;
        }

      template <class ExprOp>
      matrix(const SparseExpr<ExprOp>& result)
        {
        assign(result);
        }

      template <class ExprOp>
      matrix& operator = (const SparseExpr<ExprOp>& result)
        {
        assign(result);
        return *this;
        }

      template <class ExprOp>
      void assign(const SparseExpr<ExprOp>& result)
        {
        if (result.evaluate_before_assigning())
          assign_alias(result);
        else
          assign_no_alias(result);
        }

      template <class ExprOp>
      void assign_alias(SparseExpr<ExprOp> result)
        {
        matrix temp(result.rows(), result.cols());
        auto it = temp._entries.begin();
        for (uint64_t r = 0; r < rows(); ++r)
          {
          for (uint64_t c = 0; c < cols(); ++c)
            {
            if (r == result.first_entry() && c == result.second_entry())
              {
              *it = *result;
              ++result;
              }
            else
              *it = (T)0;
            ++it;
            }
          }
        swap(temp);
        }

      template <class ExprOp>
      void assign_no_alias(SparseExpr<ExprOp> result)
        {
        resize(result.rows(), result.cols());
        auto it = _entries.begin();
        for (uint64_t r = 0; r < rows(); ++r)
          {
          for (uint64_t c = 0; c < cols(); ++c)
            {
            if (r == result.first_entry() && c == result.second_entry())
              {
              *it = *result;
              ++result;
              }
            else
              *it = (T)0;
            ++it;
            }
          }
        }

      template <class T2, class Container2>
      matrix& operator += (matrix<T2, Container2>& other)
        {
        assert(other.rows() == _rows);
        assert(other.cols() == _cols);
        auto it = _entries.begin();
        for (const auto& value : other)
          {
          *it += (T)value;
          ++it;
          }
        return *this;
        }

      template <class ExprOp>
      matrix& operator += (const Expr<ExprOp>& result)
        {
        if (result.evaluate_before_assigning())
          add_assign_alias(result);
        else
          add_assign_no_alias(result);
        return *this;
        }

      template <class ExprOp>
      void add_assign_alias(Expr<ExprOp> result)
        {
        assert(result.rows() == _rows);
        assert(result.cols() == _cols);
        matrix temp(*this);
        auto it = temp._entries.begin();
        const uint64_t sz = _rows * _cols;
        *it += (T)*result;
        for (uint64_t i = 1; i < sz; ++i)
          {
          ++it;
          ++result;
          *it += (T)*result;
          }
        swap(temp);
        }

      template <class ExprOp>
      void add_assign_no_alias(Expr<ExprOp> result)
        {
        assert(result.rows() == _rows);
        assert(result.cols() == _cols);
        auto it = _entries.begin();
        const uint64_t sz = _rows * _cols;
        *it += (T)*result;
        for (uint64_t i = 1; i < sz; ++i)
          {
          ++it;
          ++result;
          *it += (T)*result;
          }
        }

      template <class T2, class Container2>
      matrix& operator -= (matrix<T2, Container2>& other)
        {
        assert(other.rows() == _rows);
        assert(other.cols() == _cols);
        auto it = _entries.begin();
        for (const auto& value : other)
          {
          *it -= (T)value;
          ++it;
          }
        return *this;
        }

      template <class ExprOp>
      matrix& operator -= (const Expr<ExprOp>& result)
        {
        if (result.evaluate_before_assigning())
          subtract_assign_alias(result);
        else
          subtract_assign_no_alias(result);
        return *this;
        }

      template <class ExprOp>
      void subtract_assign_alias(Expr<ExprOp> result)
        {
        assert(result.rows() == _rows);
        assert(result.cols() == _cols);
        matrix temp(*this);
        auto it = temp._entries.begin();
        const uint64_t sz = _rows * _cols;
        *it -= (T)*result;
        for (uint64_t i = 1; i < sz; ++i)
          {
          ++it;
          ++result;
          *it -= (T)*result;
          }
        swap(temp);
        }

      template <class ExprOp>
      void subtract_assign_no_alias(Expr<ExprOp> result)
        {
        assert(result.rows() == _rows);
        assert(result.cols() == _cols);
        auto it = _entries.begin();
        const uint64_t sz = _rows * _cols;
        *it -= (T)*result;
        for (uint64_t i = 1; i < sz; ++i)
          {
          ++it;
          ++result;
          *it -= (T)*result;
          }
        }

      template <class T2, class Container2>
      matrix& operator *= (matrix<T2, Container2>& other)
        {
        matrix tmp;
        tmp.noalias() = *this * other;
        swap(tmp);
        return *this;
        }

      template <class ExprOp>
      matrix& operator *= (const Expr<ExprOp>& result)
        {
        matrix tmp;
        tmp.noalias() = *this * result;
        swap(tmp);
        return *this;
        }

      NoAlias<T, Container> noalias()
        {
        NoAlias<T, Container> na(*this);
        return na;
        }

      bool empty() const
        {
        return _rows == 0 && _cols == 0;
        }

      uint64_t capacity() const
        {
        return _entries.max_size();
        }

      void resize(uint64_t rows, uint64_t cols)
        {
        assert(rows*cols <= capacity());
        _rows = rows;
        _cols = cols;
        ::jtk::implementation_details::resize<T, Container> r;
        r(_entries, _rows*_cols);
        }

      const_reference operator()(uint64_t r, uint64_t c) const
        {
        return _entries[r * _cols + c];
        }

      reference operator()(uint64_t r, uint64_t c)
        {
        return _entries[r * _cols + c];
        }

      const_reference operator()(uint64_t r) const
        {
        return _entries[r * _cols];
        }

      reference operator()(uint64_t r)
        {
        return _entries[r * _cols];
        }

      const_pointer operator[](uint64_t row) const
        {
        return _entries.data() + row * _cols;
        }

      pointer operator[](uint64_t row)
        {
        return _entries.data() + row * _cols;
        }

      pointer data()
        {
        return _entries.data();
        }

      const_pointer data() const
        {
        return _entries.data();
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      iterator begin()
        {
        return _entries.begin();
        }

      iterator end()
        {
        return _entries.begin() + _rows * _cols;
        }

      const_iterator cbegin() const
        {
        return _entries.cbegin();
        }

      const_iterator cend() const
        {
        return _entries.cbegin() + _rows * _cols;
        }

      const_iterator begin() const
        {
        return _entries.begin();
        }

      const_iterator end() const
        {
        return _entries.begin() + _rows * _cols;
        }

      template <class T2, class Container2>
      friend std::ostream& operator << (std::ostream&, const matrix<T2, Container2>&);

    private:
      uint64_t _rows, _cols;
      Container _entries;
    };

  ///////////////////////////////////////////////////////////////////////////////
  // sparse vector class
  ///////////////////////////////////////////////////////////////////////////////

  template <class T>
  class sparse_vector;

  template <class T>
  class sparse_vector_iterator
    {
    public:
      typedef T value_type;
      typedef T& reference;
      typedef uint64_t difference_type;
      typedef T* pointer;
      typedef std::random_access_iterator_tag iterator_category;

      sparse_vector_iterator() : _vector(nullptr) {}
      sparse_vector_iterator(sparse_vector<T>* vector) : _vector(vector), _iter(vector->_container.begin()) {}
      sparse_vector_iterator(const sparse_vector_iterator& other) : _vector(other._vector), _iter(other._iter) {}

      void swap(sparse_vector_iterator& other)
        {
        std::swap(_vector, other._vector);
        std::swap(_iter, other._iter);
        }

      sparse_vector_iterator& operator = (const sparse_vector_iterator& other)
        {
        sparse_vector_iterator temp(other);
        swap(temp);
        return *this;
        }

      sparse_vector_iterator end() const
        {
        sparse_vector_iterator output_iterator(*this);
        if (_vector)
          output_iterator._iter = _vector->_container.end();
        return output_iterator;
        }

      sparse_vector_iterator& operator++ ()
        {
        ++_iter;
        return *this;
        }

      sparse_vector_iterator operator++ (int)
        {
        sparse_vector_iterator output_iterator(*this);
        ++_iter;
        return output_iterator;
        }

      sparse_vector_iterator& operator-- ()
        {
        --_iter;
        return *this;
        }

      sparse_vector_iterator operator-- (int)
        {
        sparse_vector_iterator output_iterator(*this);
        --_iter;
        return output_iterator;
        }

      reference operator * () const
        {
        return _iter->second;
        }

      uint64_t entry() const
        {
        return _iter->first;
        }

      sparse_vector_iterator operator + (difference_type idx) const
        {
        sparse_vector_iterator output(*this);
        output._iter += idx;
        return output;
        }

      sparse_vector_iterator operator - (difference_type idx) const
        {
        sparse_vector_iterator output(*this);
        output._iter -= idx;
        return output;
        }

      sparse_vector_iterator& operator += (difference_type idx)
        {
        _iter += idx;
        return *this;
        }

      sparse_vector_iterator& operator -= (difference_type idx)
        {
        _iter -= idx;
        return *this;
        }

      difference_type operator + (const sparse_vector_iterator& other) const
        {
        return _iter + other._iter;
        }

      difference_type operator - (const sparse_vector_iterator& other) const
        {
        return _iter - other._iter;
        }

      bool operator == (const sparse_vector_iterator& other) const
        {
        return _iter == other._iter;
        }

      bool operator != (const sparse_vector_iterator& other) const
        {
        return !(*this == other);
        }

      bool operator < (const sparse_vector_iterator& other) const
        {
        if (_iter == end())
          return false;
        if (other._iter == other.end())
          return true;
        return _iter->first < other._iter->first;
        }

      bool operator <= (const sparse_vector_iterator& other) const
        {
        return (*this == other || *this < other);
        }

      bool operator > (const sparse_vector_iterator& other) const
        {
        return !(*this <= other);
        }

      bool operator >= (const sparse_vector_iterator& other) const
        {
        return !(*this < other);
        }

    private:
      sparse_vector<T>* _vector;
      typename std::vector<std::pair<uint64_t, T>>::iterator _iter;
    };

  template <class T>
  class sparse_vector_const_iterator
    {
    public:
      typedef T value_type;
      typedef const T& reference;
      typedef uint64_t difference_type;
      typedef const T* pointer;
      typedef std::random_access_iterator_tag iterator_category;

      sparse_vector_const_iterator() : _vector(nullptr) {}
      sparse_vector_const_iterator(const sparse_vector<T>* vector) : _vector(vector), _iter(vector->_container.begin()) {}
      sparse_vector_const_iterator(const sparse_vector_const_iterator& other) : _vector(other._vector), _iter(other._iter) {}

      void swap(sparse_vector_const_iterator& other)
        {
        std::swap(_vector, other._vector);
        std::swap(_iter, other._iter);
        }

      sparse_vector_const_iterator& operator = (const sparse_vector_const_iterator& other)
        {
        sparse_vector_const_iterator temp(other);
        swap(temp);
        return *this;
        }

      sparse_vector_const_iterator end() const
        {
        sparse_vector_const_iterator output_iterator(*this);
        if (_vector)
          output_iterator._iter = _vector->_container.end();
        return output_iterator;
        }

      sparse_vector_const_iterator& operator++ ()
        {
        ++_iter;
        return *this;
        }

      sparse_vector_const_iterator operator++ (int)
        {
        sparse_vector_const_iterator output_iterator(*this);
        ++_iter;
        return output_iterator;
        }

      sparse_vector_const_iterator& operator-- ()
        {
        --_iter;
        return *this;
        }

      sparse_vector_const_iterator operator-- (int)
        {
        sparse_vector_const_iterator output_iterator(*this);
        --_iter;
        return output_iterator;
        }

      reference operator * () const
        {
        return _iter->second;
        }

      uint64_t entry() const
        {
        return _iter->first;
        }

      sparse_vector_const_iterator operator + (difference_type idx) const
        {
        sparse_vector_const_iterator output(*this);
        output._iter += idx;
        return output;
        }

      sparse_vector_const_iterator operator - (difference_type idx) const
        {
        sparse_vector_const_iterator output(*this);
        output._iter -= idx;
        return output;
        }

      sparse_vector_const_iterator& operator += (difference_type idx)
        {
        _iter += idx;
        return *this;
        }

      sparse_vector_const_iterator& operator -= (difference_type idx)
        {
        _iter -= idx;
        return *this;
        }

      difference_type operator + (const sparse_vector_const_iterator& other) const
        {
        return _iter + other._iter;
        }

      difference_type operator - (const sparse_vector_const_iterator& other) const
        {
        return _iter - other._iter;
        }

      bool operator == (const sparse_vector_const_iterator& other) const
        {
        return _iter == other._iter;
        }

      bool operator != (const sparse_vector_const_iterator& other) const
        {
        return !(*this == other);
        }

      bool operator < (const sparse_vector_const_iterator& other) const
        {
        if (_iter == end())
          return false;
        if (other._iter == other.end())
          return true;
        return _iter->first < other._iter->first;
        }

      bool operator <= (const sparse_vector_const_iterator& other) const
        {
        return (*this == other || *this < other);
        }

      bool operator > (const sparse_vector_const_iterator& other) const
        {
        return !(*this <= other);
        }

      bool operator >= (const sparse_vector_const_iterator& other) const
        {
        return !(*this < other);
        }

    private:
      const sparse_vector<T>* _vector;
      typename std::vector<std::pair<uint64_t, T>>::const_iterator _iter;
    };

  template <class T>
  class sparse_vector
    {
    template <class T2> friend class sparse_vector_iterator;
    template <class T2> friend class sparse_vector_const_iterator;

    struct compare
      {
      inline bool operator () (const std::pair<uint64_t, T>& lhs, const std::pair<uint64_t, T>& rhs) const
        {
        return lhs.first < rhs.first;
        }

      inline bool operator () (const std::pair<uint64_t, T>& lhs, uint64_t value) const
        {
        return lhs.first < value;
        }

      inline bool operator () (uint64_t value, const std::pair<uint64_t, T>& rhs) const
        {
        return value < rhs.first;
        }
      };

    public:
      using value_type = T;
      using pointer = T * ;
      using reference = T & ;
      using const_pointer = const T *;
      using const_reference = const T &;
      using iterator = sparse_vector_iterator<T>;
      using const_iterator = sparse_vector_const_iterator<T>;

      sparse_vector() : _size(0), _zero((T)0) {}

      sparse_vector(uint64_t size) : _size(size), _zero((T)0) {}

      void swap(sparse_vector<T>& other)
        {
        std::swap(_container, other._container);
        std::swap(_size, other._size);
        std::swap(_zero, other._zero);
        }

      sparse_vector(const sparse_vector& other)
        {
        sparse_vector tmp(other._size);
        const_iterator iter_other = other.begin();
        const_iterator end_other = other.end();
        for (; iter_other != end_other; ++iter_other)
          {
          if (*iter_other != static_cast<T>(0))
            tmp.put(iter_other.entry()) = *iter_other;
          }
        swap(tmp);
        }

      sparse_vector(sparse_vector&&) = default;

      sparse_vector& operator = (const sparse_vector& other)
        {
        sparse_vector temp(other);
        swap(temp);
        return *this;
        }

      sparse_vector& operator = (sparse_vector&&) = default;

      template <class T2>
      sparse_vector(const sparse_vector<T2>& other)
        {
        sparse_vector tmp(other._size);
        const_iterator iter_other = other.begin();
        const_iterator end_other = other.end();
        for (; iter_other != end_other; ++iter_other)
          {
          if (*iter_other != static_cast<T2>(0))
            tmp.put(iter_other.entry()) = static_cast<T>(*iter_other);
          }
        swap(tmp);
        }

      template <class T2>
      sparse_vector& operator = (const sparse_vector<T2>& other)
        {
        sparse_vector temp(other);
        swap(temp);
        return *this;
        }

      reference put(uint64_t idx)
        {
        auto iter = std::lower_bound(_container.begin(), _container.end(), idx, compare());
        if (iter != _container.end())
          {
          if (iter->first == idx)
            {
            return iter->second;
            }
          else
            {
            iter = _container.insert(iter, std::pair<uint64_t, T>(idx, static_cast<T>(0)));
            return iter->second;
            }
          }
        else
          {
          _container.push_back(std::pair<uint64_t, T>(idx, static_cast<T>(0)));
          return _container.back().second;
          }
        }

      const_reference get(uint64_t idx) const
        {
        const auto iter = std::lower_bound(_container.begin(), _container.end(), idx, compare());
        return (iter != _container.end() && iter->first == idx ? iter->second : _zero);
        }

      iterator begin()
        {
        return iterator(this);
        }

      iterator end()
        {
        return iterator(this).end();
        }

      const_iterator begin() const
        {
        return const_iterator(this);
        }

      const_iterator end() const
        {
        return const_iterator(this).end();
        }

      void resize(uint64_t size)
        {
        _size = size;
        }

      void operator += (const sparse_vector<T>& other)
        {
        assert(other.size() <= size());
        const_iterator iter = other.begin();
        const_iterator iter_end = other.end();
        for (; iter != iter_end; ++iter)
          put(iter.entry()) += *iter;
        }

      void operator -= (const sparse_vector<T>& other)
        {
        assert(other.size() <= size());
        const_iterator iter = other.begin();
        const_iterator iter_end = other.end();
        for (; iter != iter_end; ++iter)
          put(iter.entry()) -= *iter;
        }

      void operator *= (value_type coeff)
        {
        iterator iter = begin();
        iterator iter_end = end();
        for (; iter != iter_end; ++iter)
          *iter *= coeff;
        }

      void operator /= (value_type coeff)
        {
        iterator iter = begin();
        iterator iter_end = end();
        for (; iter != iter_end; ++iter)
          *iter /= coeff;
        }

      uint64_t size() const
        {
        return _size;
        }

      uint64_t entries_stored() const
        {
        return (uint64_t)_container.size();
        }

      void clear()
        {
        _container.clear();
        }

      bool operator == (const sparse_vector& other) const
        {
        return _size == other._size && _container == other._container;
        }

      bool operator != (const sparse_vector& other) const
        {
        return !(*this == other);
        }

    private:
      std::vector<std::pair<uint64_t, T>> _container;
      uint64_t _size;
      T _zero;
    };

  ///////////////////////////////////////////////////////////////////////////////
  // sparse_matrix class
  ///////////////////////////////////////////////////////////////////////////////

  template <class T>
  class sparse_matrix_iterator
    {
    public:
      typedef T value_type;
      typedef T& reference;
      typedef uint64_t difference_type;
      typedef T* pointer;
      typedef std::forward_iterator_tag iterator_category;

      sparse_matrix_iterator() : _matrix(nullptr) {}
      sparse_matrix_iterator(sparse_matrix<T>* matrix) : _matrix(matrix), _row(0), _iter()
        {
        if (_matrix && !_matrix->empty())
          {
          _iter = _matrix->_entries[0].begin();
          _update_iterator_forward();
          }
        }
      sparse_matrix_iterator(const sparse_matrix_iterator& other) : _matrix(other._matrix), _row(other._row), _iter(other._iter) {}

      void swap(sparse_matrix_iterator& other)
        {
        std::swap(_matrix, other._matrix);
        std::swap(_row, other._row);
        std::swap(_iter, other._iter);
        }

      sparse_matrix_iterator& operator = (const sparse_matrix_iterator& other)
        {
        sparse_matrix_iterator temp(other);
        swap(temp);
        return *this;
        }

      sparse_matrix_iterator end() const
        {
        sparse_matrix_iterator output_iterator(*this);
        if (_matrix)
          {
          output_iterator._row = _matrix->rows();
          output_iterator._iter = _matrix->_entries[output_iterator._row - 1].end();
          }
        return output_iterator;
        }

      sparse_matrix_iterator row(uint64_t idx) const
        {
        sparse_matrix_iterator output_iterator(*this);
        output_iterator._row = idx;
        output_iterator._iter = _matrix->_entries[idx].begin();
        output_iterator._update_iterator_forward();
        return output_iterator;
        }

      sparse_matrix_iterator& operator++ ()
        {
        ++_iter;
        _update_iterator_forward();
        return *this;
        }

      sparse_matrix_iterator operator++ (int)
        {
        sparse_matrix_iterator output_iterator(*this);
        ++_iter;
        _update_iterator_forward();
        return output_iterator;
        }

      reference operator * () const
        {
        return *_iter;
        }

      uint64_t first_entry() const
        {
        return _row;
        }

      uint64_t second_entry() const
        {
        return _iter.entry();
        }

      bool operator == (const sparse_matrix_iterator& other) const
        {
        return _row == other._row && _iter == other._iter;
        }

      bool operator != (const sparse_matrix_iterator& other) const
        {
        return !(*this == other);
        }

      bool operator < (const sparse_matrix_iterator& other) const
        {
        if (_row == other._row)
          return _iter < other._iter;
        return _row < other._row;
        }

      bool operator <= (const sparse_matrix_iterator& other) const
        {
        return (*this == other || *this < other);
        }

      bool operator > (const sparse_matrix_iterator& other) const
        {
        return !(*this <= other);
        }

      bool operator >= (const sparse_matrix_iterator& other) const
        {
        return !(*this < other);
        }

    private:
      void _update_iterator_forward()
        {
        while (_iter == _iter.end() && ++_row < _matrix->rows())
          _iter = _matrix->_entries[_row].begin();
        }

    private:
      sparse_matrix<T>* _matrix;
      uint64_t _row;
      typename sparse_vector<T>::iterator _iter;
    };

  template <class T>
  class sparse_matrix_const_iterator
    {
    public:
      typedef T value_type;
      typedef const T& reference;
      typedef uint64_t difference_type;
      typedef const T* pointer;
      typedef std::forward_iterator_tag iterator_category;

      sparse_matrix_const_iterator() : _matrix(nullptr) {}
      sparse_matrix_const_iterator(const sparse_matrix<T>* matrix) : _matrix(matrix), _row(0), _iter()
        {
        if (_matrix && !_matrix->empty())
          {
          _iter = _matrix->_entries[0].begin();
          _update_iterator_forward();
          }
        }
      sparse_matrix_const_iterator(const sparse_matrix_const_iterator& other) : _matrix(other._matrix), _row(other._row), _iter(other._iter) {}

      void swap(sparse_matrix_const_iterator& other)
        {
        std::swap(_matrix, other._matrix);
        std::swap(_row, other._row);
        std::swap(_iter, other._iter);
        }

      sparse_matrix_const_iterator& operator = (const sparse_matrix_const_iterator& other)
        {
        sparse_matrix_const_iterator temp(other);
        swap(temp);
        return *this;
        }

      sparse_matrix_const_iterator end() const
        {
        sparse_matrix_const_iterator output_iterator(*this);
        if (_matrix)
          {
          output_iterator._row = _matrix->rows();
          output_iterator._iter = _matrix->_entries[output_iterator._row - 1].end();
          }
        return output_iterator;
        }

      sparse_matrix_const_iterator row(uint64_t idx) const
        {
        sparse_matrix_const_iterator output_iterator(*this);
        output_iterator._row = idx;
        output_iterator._iter = _matrix->_entries[idx].begin();
        output_iterator._update_iterator_forward();
        return output_iterator;
        }

      sparse_matrix_const_iterator& operator++ ()
        {
        ++_iter;
        _update_iterator_forward();
        return *this;
        }

      sparse_matrix_const_iterator operator++ (int)
        {
        sparse_matrix_const_iterator output_iterator(*this);
        ++_iter;
        _update_iterator_forward();
        return output_iterator;
        }

      reference operator * () const
        {
        return *_iter;
        }

      uint64_t first_entry() const
        {
        return _row;
        }

      uint64_t second_entry() const
        {
        return _iter.entry();
        }

      bool operator == (const sparse_matrix_const_iterator& other) const
        {
        return _row == other._row && _iter == other._iter;
        }

      bool operator != (const sparse_matrix_const_iterator& other) const
        {
        return !(*this == other);
        }

      bool operator < (const sparse_matrix_const_iterator& other) const
        {
        if (_row == other._row)
          return _iter < other._iter;
        return _row < other._row;
        }

      bool operator <= (const sparse_matrix_const_iterator& other) const
        {
        return (*this == other || *this < other);
        }

      bool operator > (const sparse_matrix_const_iterator& other) const
        {
        return !(*this <= other);
        }

      bool operator >= (const sparse_matrix_const_iterator& other) const
        {
        return !(*this < other);
        }

    private:
      void _update_iterator_forward()
        {
        while (_iter == _iter.end() && ++_row < _matrix->rows())
          _iter = _matrix->_entries[_row].begin();
        }

    private:
      const sparse_matrix<T>* _matrix;
      uint64_t _row;
      typename sparse_vector<T>::const_iterator _iter;
    };

  template <class T>
  class sparse_matrix
    {
    template <class T2> friend class sparse_matrix_iterator;
    template <class T2> friend class sparse_matrix_const_iterator;

    public:
      using value_type = T;
      using pointer = T * ;
      using reference = T & ;
      using const_pointer = const T *;
      using const_reference = const T &;
      using iterator = sparse_matrix_iterator<T>;
      using const_iterator = sparse_matrix_const_iterator<T>;

      sparse_matrix() : _rows(0), _cols(0)
        {
        }

      sparse_matrix(uint64_t rows, uint64_t cols) : _rows(rows), _cols(cols), _entries(rows, sparse_vector<T>(cols))
        {
        }

      void swap(sparse_matrix& other)
        {
        std::swap(_cols, other._cols);
        std::swap(_rows, other._rows);
        std::swap(_entries, other._entries);
        }

      sparse_matrix(const sparse_matrix& other) : _rows(other._rows), _cols(other._cols), _entries(other._entries)
        {
        }

      sparse_matrix(sparse_matrix&&) = default;

      sparse_matrix& operator = (const sparse_matrix& other)
        {
        sparse_matrix temp(other);
        swap(temp);
        return *this;
        }

      sparse_matrix& operator = (sparse_matrix&&) = default;

      template <class T2>
      sparse_matrix(const sparse_matrix<T2>& other) : _rows(other.rows()), _cols(other.cols())
        {
        resize(_rows, _cols);
        auto it = _entries.begin();
        auto other_it = other.begin();
        const auto other_it_end = other.end();
        for (; other_it != other_it_end; ++other_it)
          {
          *it = *other_it;
          ++it;
          }
        }

      template <class T2>
      sparse_matrix& operator = (const sparse_matrix<T2>& other)
        {
        sparse_matrix temp(other);
        swap(temp);
        return *this;
        }

      template <class ExprOp>
      sparse_matrix(const SparseExpr<ExprOp>& result)
        {
        assign(result);
        }

      template <class ExprOp>
      sparse_matrix& operator = (const SparseExpr<ExprOp>& result)
        {
        assign(result);
        return *this;
        }

      template <class ExprOp>
      void assign(const SparseExpr<ExprOp>& result)
        {
        if (result.evaluate_before_assigning())
          assign_alias(result);
        else
          assign_no_alias(result);
        }

      template <class ExprOp>
      void assign_alias(SparseExpr<ExprOp> result)
        {
        sparse_matrix temp(result.rows(), result.cols());
        auto iter_end = result.end();
        while (result != iter_end)
          {
          auto value = *result;
          if (value)
            temp.put(result.first_entry(), result.second_entry()) = value;
          ++result;
          }
        swap(temp);
        }

      template <class ExprOp>
      void assign_no_alias(SparseExpr<ExprOp> result)
        {
        resize(result.rows(), result.cols());
        for (auto& row : _entries)
          row.clear();
        auto iter_end = result.end();
        while (result != iter_end)
          {
          auto value = *result;
          if (value)
            put(result.first_entry(), result.second_entry()) = value;
          ++result;
          }
        }

      template <class T2, class TContainer>
      sparse_matrix(const matrix<T2, TContainer>& m)
        {
        UnExprOp<typename matrix<T2, TContainer>::const_iterator, OpId<T2>> expr(m.begin(), m.rows(), m.cols(), false);
        assign_no_alias(Expr<UnExprOp< typename matrix<T2, TContainer>::const_iterator, OpId<T2>>>(expr));
        }

      template <class T2, class TContainer>
      sparse_matrix& operator = (const matrix<T2, TContainer>& m)
        {
        UnExprOp<typename matrix<T2, TContainer>::const_iterator, OpId<T2>> expr(m.begin(), m.rows(), m.cols(), false);
        assign_no_alias(Expr<UnExprOp< typename matrix<T2, TContainer>::const_iterator, OpId<T2>>>(expr));
        return *this;
        }

      template <class ExprOp>
      sparse_matrix(const Expr<ExprOp>& result)
        {
        assign(result);
        }

      template <class ExprOp>
      sparse_matrix& operator = (const Expr<ExprOp>& result)
        {
        assign(result);
        return *this;
        }

      template <class ExprOp>
      void assign(const Expr<ExprOp>& result)
        {
        if (result.evaluate_before_assigning())
          assign_alias(result);
        else
          assign_no_alias(result);
        }

      template <class ExprOp>
      void assign_alias(Expr<ExprOp> result)
        {
        sparse_matrix temp(result.rows(), result.cols());
        for (uint64_t r = 0; r < result.rows(); ++r)
          {
          for (uint64_t c = 0; c < result.cols(); ++c)
            {
            T value = (T)*result;
            if (value)
              put(r, c) = value;
            ++result;
            }
          }
        swap(temp);
        }

      template <class ExprOp>
      void assign_no_alias(Expr<ExprOp> result)
        {
        resize(result.rows(), result.cols());
        for (auto& row : _entries)
          row.clear();
        for (uint64_t r = 0; r < result.rows(); ++r)
          {
          for (uint64_t c = 0; c < result.cols(); ++c)
            {
            T value = (T)*result;
            if (value)
              put(r, c) = value;
            ++result;
            }
          }
        }

      bool empty() const
        {
        return _rows == 0 && _cols == 0;
        }

      uint64_t entries_stored() const
        {
        uint64_t output = 0;
        for (const auto& row : _entries)
          output += row.entries_stored();
        return output;
        }

      void resize(uint64_t rows, uint64_t cols)
        {
        _rows = rows;
        _cols = cols;
        _entries.resize(rows);
        for (auto& row : _entries)
          row.resize(cols);
        }

      const_reference operator()(uint64_t r, uint64_t c) const
        {
        return _entries[r].get(c);
        }

      const_reference get(uint64_t r, uint64_t c) const
        {
        return _entries[r].get(c);
        }

      reference put(uint64_t r, uint64_t c)
        {
        return _entries[r].put(c);
        }

      uint64_t rows() const
        {
        return _rows;
        }

      uint64_t cols() const
        {
        return _cols;
        }

      iterator begin()
        {
        return iterator(this);
        }

      iterator end()
        {
        return iterator(this).end();
        }

      const_iterator cbegin() const
        {
        return const_iterator(this);
        }

      const_iterator cend() const
        {
        return const_iterator(this).end();
        }

      const_iterator begin() const
        {
        return const_iterator(this);
        }

      const_iterator end() const
        {
        return const_iterator(this).end();
        }

      sparse_vector<T>& row(uint64_t idx)
        {
        return _entries[idx];
        }

      const sparse_vector<T>& row(uint64_t idx) const
        {
        return _entries[idx];
        }

    private:
      uint64_t _rows, _cols;
      std::vector<sparse_vector<T>> _entries;
    };

  ///////////////////////////////////////////////////////////////////////////////
  // Comparison operators
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2, class Container2>
  bool operator == (const matrix<T, Container>& left, const matrix<T2, Container2>& right)
    {
    if (left.rows() != right.rows())
      return false;
    if (left.cols() != right.cols())
      return false;
    auto it = left.begin();
    auto it2 = right.begin();
    const auto it_end = left.end();
    for (; it != it_end; ++it, ++it2)
      {
      if (*it != *it2)
        return false;
      }
    return true;
    }

  template <class T, class Container, class T2, class Container2>
  bool operator != (const matrix<T, Container>& left, const matrix<T2, Container2>& right)
    {
    return !(left == right);
    }

  template <class T, class Container, class ExprOp>
  bool operator == (const matrix<T, Container>& left, Expr<ExprOp> right)
    {
    if (left.rows() != right.rows())
      return false;
    if (left.cols() != right.cols())
      return false;
    auto it = left.begin();
    const auto it_end = left.end();
    for (; it != it_end; ++it, ++right)
      {
      if (*it != *right)
        return false;
      }
    return true;
    }

  template <class T, class Container, class ExprOp>
  bool operator != (const matrix<T, Container>& left, const Expr<ExprOp>& right)
    {
    return !(left == right);
    }

  template <class T, class Container, class ExprOp>
  bool operator == (const Expr<ExprOp>& left, const matrix<T, Container>& right)
    {
    return right == left;
    }

  template <class T, class Container, class ExprOp>
  bool operator != (const Expr<ExprOp>& left, const matrix<T, Container>& right)
    {
    return !(right == left);
    }

  template <class ExprOp1, class ExprOp2>
  bool operator == (Expr<ExprOp1> left, Expr<ExprOp2> right)
    {
    if (left.rows() != right.rows())
      return false;
    if (left.cols() != right.cols())
      return false;
    if (*left != *right)
      return false;
    const uint64_t sz = left.rows()*left.cols();
    for (uint64_t i = 1; i < sz; ++i)
      {
      ++left;
      ++right;
      if (*left != *right)
        return false;
      }
    return true;
    }

  template <class ExprOp1, class ExprOp2>
  bool operator != (const Expr<ExprOp1>& left, const Expr<ExprOp2>& right)
    {
    return !(left == right);
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Stream operators
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container>
  std::ostream& operator << (std::ostream& os, const matrix<T, Container>& m)
    {
    if (m.empty())
      {
      return os;
      }
    auto p = os.precision();
    auto ff = os.flags();
    int width = 0;
    auto it = m.begin();
    for (uint64_t r = 0; r < m.rows(); ++r)
      {
      for (uint64_t c = 0; c < m.cols(); ++c)
        {
        std::stringstream sstr;
        sstr << std::setprecision(p);
        sstr.flags(ff);
        sstr << *it;
        width = std::max<int>(width, (int)sstr.str().length());
        ++it;
        }
      }

    it = m.begin();
    for (uint64_t r = 0; r < m.rows(); ++r)
      {
      os.width(width);
      os << *it;
      ++it;
      for (uint64_t c = 1; c < m.cols(); ++c)
        {
        os << " ";
        os.width(width);
        os << *it;
        ++it;
        }
      os << std::endl;
      }
    return os;
    }

  template <class T, class Container, class T2>
  CommaInitializer<T, Container> operator << (matrix<T, Container>& m, const T2& value)
    {
    CommaInitializer<T, Container> ci(&m, (T)value);
    return ci;
    }

  template <class T, class Container>
  CommaConcatenator<T, Container> operator << (matrix<T, Container>& m, const matrix<T, Container>& matrix_to_concat)
    {
    CommaConcatenator<T, Container> ci(&m, matrix_to_concat);
    return ci;
    }

  template <class T, class Container, class ExprOp>
  CommaConcatenator<T, Container> operator << (matrix<T, Container>& m, const Expr<ExprOp>& matrix_to_concat)
    {
    CommaConcatenator<T, Container> ci(&m, matrix_to_concat);
    return ci;
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Adding matrices
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2, class Container2 >
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    typename matrix<T2, Container2>::const_iterator,
    OpAdd<typename gettype<T, T2>::ty> > > operator + (const matrix<T, Container>& a, const matrix<T2, Container2>& b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      typename matrix<T2, Container2>::const_iterator,
      OpAdd<typename gettype<T, T2>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a.begin(), b.begin(), a.rows(), a.cols(), false));
    }

  template <class T, class Container, class ExprOp >
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    Expr<ExprOp>,
    OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > > operator + (const matrix<T, Container>& a, const Expr<ExprOp>& b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      Expr<ExprOp>,
      OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a.begin(), b, a.rows(), a.cols(), b.evaluate_before_assigning()));
    }

  template <class T, class Container, class ExprOp >
  Expr<
    BinExprOp<
    Expr<ExprOp>,
    typename matrix<T, Container>::const_iterator,
    OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > > operator + (const Expr<ExprOp>& a, const matrix<T, Container>& b)
    {
    typedef BinExprOp<
      Expr<ExprOp>,
      typename matrix<T, Container>::const_iterator,
      OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a, b.begin(), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp1, class ExprOp2 >
  Expr<
    BinExprOp<
    Expr<ExprOp1>,
    Expr<ExprOp2>,
    OpAdd<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > > operator + (const Expr<ExprOp1>& a, const Expr<ExprOp2>& b)
    {
    typedef BinExprOp<
      Expr<ExprOp1>,
      Expr<ExprOp2>,
      OpAdd<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a, b, a.rows(), a.cols(), a.evaluate_before_assigning() || b.evaluate_before_assigning()));
    }


  template <class T, class T2>
  SparseExpr<
    BinExprSparseOp<
    typename sparse_matrix<T>::const_iterator,
    typename sparse_matrix<T2>::const_iterator,
    OpAdd<typename gettype<T, T2>::ty> > > operator + (const sparse_matrix<T>& a, const sparse_matrix<T2>& b)
    {
    typedef BinExprSparseOp<
      typename sparse_matrix<T>::const_iterator,
      typename sparse_matrix<T2>::const_iterator,
      OpAdd<typename gettype<T, T2>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a.begin(), b.begin(), a.rows(), a.cols(), false));
    }

  template <class T, class ExprOp >
  SparseExpr<
    BinExprSparseOp<
    typename sparse_matrix<T>::const_iterator,
    SparseExpr<ExprOp>,
    OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > > operator + (const sparse_matrix<T>& a, const SparseExpr<ExprOp>& b)
    {
    typedef BinExprSparseOp<
      typename sparse_matrix<T>::const_iterator,
      SparseExpr<ExprOp>,
      OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a.begin(), b, a.rows(), a.cols(), b.evaluate_before_assigning()));
    }

  template <class T, class ExprOp >
  SparseExpr<
    BinExprSparseOp<
    SparseExpr<ExprOp>,
    typename sparse_matrix<T>::const_iterator,
    OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > > operator + (const SparseExpr<ExprOp>& a, const sparse_matrix<T>& b)
    {
    typedef BinExprSparseOp<
      SparseExpr<ExprOp>,
      typename sparse_matrix<T>::const_iterator,
      OpAdd<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a, b.begin(), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp1, class ExprOp2 >
  SparseExpr<
    BinExprSparseOp<
    SparseExpr<ExprOp1>,
    SparseExpr<ExprOp2>,
    OpAdd<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > > operator + (const SparseExpr<ExprOp1>& a, const SparseExpr<ExprOp2>& b)
    {
    typedef BinExprSparseOp<
      SparseExpr<ExprOp1>,
      SparseExpr<ExprOp2>,
      OpAdd<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a, b, a.rows(), a.cols(), a.evaluate_before_assigning() || b.evaluate_before_assigning()));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Subtracting matrices
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2, class Container2 >
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    typename matrix<T2, Container2>::const_iterator,
    OpSub<typename gettype<T, T2>::ty> > > operator - (const matrix<T, Container>& a, const matrix<T2, Container2>& b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      typename matrix<T2, Container2>::const_iterator,
      OpSub<typename gettype<T, T2>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a.begin(), b.begin(), a.rows(), a.cols(), false));
    }

  template <class T, class Container, class ExprOp >
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    Expr<ExprOp>,
    OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > > operator - (const matrix<T, Container>& a, const Expr<ExprOp>& b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      Expr<ExprOp>,
      OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a.begin(), b, a.rows(), a.cols(), b.evaluate_before_assigning()));
    }

  template <class T, class Container, class ExprOp >
  Expr<
    BinExprOp<
    Expr<ExprOp>,
    typename matrix<T, Container>::const_iterator,
    OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > > operator - (const Expr<ExprOp>& a, const matrix<T, Container>& b)
    {
    typedef BinExprOp<
      Expr<ExprOp>,
      typename matrix<T, Container>::const_iterator,
      OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a, b.begin(), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp1, class ExprOp2 >
  Expr<
    BinExprOp<
    Expr<ExprOp1>,
    Expr<ExprOp2>,
    OpSub<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > > operator - (const Expr<ExprOp1>& a, const Expr<ExprOp2>& b)
    {
    typedef BinExprOp<
      Expr<ExprOp1>,
      Expr<ExprOp2>,
      OpSub<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > ExprT;
    assert(a.rows()*a.cols() == b.rows()*b.cols());
    return Expr<ExprT>(ExprT(a, b, a.rows(), a.cols(), a.evaluate_before_assigning() || b.evaluate_before_assigning()));
    }


  template <class T, class T2>
  SparseExpr<
    BinExprSparseOp<
    typename sparse_matrix<T>::const_iterator,
    typename sparse_matrix<T2>::const_iterator,
    OpSub<typename gettype<T, T2>::ty> > > operator - (const sparse_matrix<T>& a, const sparse_matrix<T2>& b)
    {
    typedef BinExprSparseOp<
      typename sparse_matrix<T>::const_iterator,
      typename sparse_matrix<T2>::const_iterator,
      OpSub<typename gettype<T, T2>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a.begin(), b.begin(), a.rows(), a.cols(), false));
    }

  template <class T, class ExprOp >
  SparseExpr<
    BinExprSparseOp<
    typename sparse_matrix<T>::const_iterator,
    SparseExpr<ExprOp>,
    OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > > operator - (const sparse_matrix<T>& a, const SparseExpr<ExprOp>& b)
    {
    typedef BinExprSparseOp<
      typename sparse_matrix<T>::const_iterator,
      SparseExpr<ExprOp>,
      OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a.begin(), b, a.rows(), a.cols(), b.evaluate_before_assigning()));
    }

  template <class T, class ExprOp >
  SparseExpr<
    BinExprSparseOp<
    SparseExpr<ExprOp>,
    typename sparse_matrix<T>::const_iterator,
    OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > > operator - (const SparseExpr<ExprOp>& a, const sparse_matrix<T>& b)
    {
    typedef BinExprSparseOp<
      SparseExpr<ExprOp>,
      typename sparse_matrix<T>::const_iterator,
      OpSub<typename gettype<T, typename ExprOp::value_type>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a, b.begin(), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp1, class ExprOp2 >
  SparseExpr<
    BinExprSparseOp<
    SparseExpr<ExprOp1>,
    SparseExpr<ExprOp2>,
    OpSub<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > > operator - (const SparseExpr<ExprOp1>& a, const SparseExpr<ExprOp2>& b)
    {
    typedef BinExprSparseOp<
      SparseExpr<ExprOp1>,
      SparseExpr<ExprOp2>,
      OpSub<typename gettype<typename ExprOp1::value_type, typename ExprOp2::value_type>::ty> > ExprT;
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    return SparseExpr<ExprT>(ExprT(a, b, a.rows(), a.cols(), a.evaluate_before_assigning() || b.evaluate_before_assigning()));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Negating matrices
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container >
  Expr<
    UnExprOp<
    typename matrix<T, Container>::const_iterator,
    OpNeg<T> > > operator - (const matrix<T, Container>& a)
    {
    typedef UnExprOp<
      typename matrix<T, Container>::const_iterator,
      OpNeg<T> > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), a.rows(), a.cols(), false));
    }

  template <class ExprOp>
  Expr<
    UnExprOp<
    Expr<ExprOp>,
    OpNeg<typename ExprOp::value_type> > > operator - (const Expr<ExprOp>& a)
    {
    typedef UnExprOp<
      Expr<ExprOp>,
      OpNeg<typename ExprOp::value_type> > ExprT;
    return Expr<ExprT>(ExprT(a, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class T >
  SparseExpr<
    UnExprSparseOp<
    typename sparse_matrix<T>::const_iterator,
    OpNeg<T> > > operator - (const sparse_matrix<T>& a)
    {
    typedef UnExprSparseOp<
      typename sparse_matrix<T>::const_iterator,
      OpNeg<T> > ExprT;
    return SparseExpr<ExprT>(ExprT(a.begin(), a.rows(), a.cols(), false));
    }

  template <class ExprOp>
  SparseExpr<
    UnExprSparseOp<
    SparseExpr<ExprOp>,
    OpNeg<typename ExprOp::value_type> > > operator - (const SparseExpr<ExprOp>& a)
    {
    typedef UnExprSparseOp<
      SparseExpr<ExprOp>,
      OpNeg<typename ExprOp::value_type> > ExprT;
    return SparseExpr<ExprT>(ExprT(a, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Scalar multiplication
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2>
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    Constant<T>,
    OpMul<typename gettype<T, T2>::ty> > > operator * (const matrix<T, Container>& a, T2 b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      Constant<T>,
      OpMul<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), Constant<T>((T)b, a.rows(), a.cols()), a.rows(), a.cols(), false));
    }

  template <class T, class Container, class T2 >
  Expr<
    BinExprOp<
    Constant<T>,
    typename matrix<T, Container>::const_iterator,
    OpMul<typename gettype<T, T2>::ty> > > operator * (T2 a, const matrix<T, Container>& b)
    {
    typedef BinExprOp<
      Constant<T>,
      typename matrix<T, Container>::const_iterator,
      OpMul<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<T>((T)a, b.rows(), b.cols()), b.begin(), b.rows(), b.cols(), false));
    }

  template <class ExprOp, class T2>
  Expr<
    BinExprOp<
    Expr<ExprOp>,
    Constant<typename ExprOp::value_type>,
    OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator * (const Expr<ExprOp>& a, T2 b)
    {
    typedef BinExprOp<
      Expr<ExprOp>,
      Constant<typename ExprOp::value_type>,
      OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a, Constant<typename ExprOp::value_type>((typename ExprOp::value_type)b, a.rows(), a.cols()), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp, class T2 >
  Expr<
    BinExprOp<
    Constant<typename ExprOp::value_type>,
    Expr<ExprOp>,
    OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator * (T2 a, const Expr<ExprOp>& b)
    {
    typedef BinExprOp<
      Constant<typename ExprOp::value_type>,
      Expr<ExprOp>,
      OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<typename ExprOp::value_type>((typename ExprOp::value_type)a, b.rows(), b.cols()), b, b.rows(), b.cols(), b.evaluate_before_assigning()));
    }


  template <class T, class T2>
  SparseExpr <
    BinExprSparseScalarOp<
    typename sparse_matrix<T>::const_iterator,
    OpMul<typename gettype<T, T2>::ty> > > operator * (const sparse_matrix<T>& a, T2 b)
    {
    typedef BinExprSparseScalarOp<
      typename sparse_matrix<T>::const_iterator,
      OpMul<typename gettype<T, T2>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a.begin(), (T)b, a.rows(), a.cols(), false));
    }

  template <class T, class T2>
  SparseExpr <
    BinExprSparseScalarOp<
    typename sparse_matrix<T>::const_iterator,
    OpMul<typename gettype<T, T2>::ty> > > operator * (T2 b, const sparse_matrix<T>& a)
    {
    typedef BinExprSparseScalarOp<
      typename sparse_matrix<T>::const_iterator,
      OpMul<typename gettype<T, T2>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a.begin(), (T)b, a.rows(), a.cols(), false));
    }

  template <class ExprOp, class T2>
  SparseExpr<
    BinExprSparseScalarOp<
    SparseExpr<ExprOp>,
    OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator * (const SparseExpr<ExprOp>& a, T2 b)
    {
    typedef BinExprSparseScalarOp<
      SparseExpr<ExprOp>,
      OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a, (typename ExprOp::value_type)b, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp, class T2>
  SparseExpr<
    BinExprSparseScalarOp<
    SparseExpr<ExprOp>,
    OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator * (T2 b, const SparseExpr<ExprOp>& a)
    {
    typedef BinExprSparseScalarOp<
      SparseExpr<ExprOp>,
      OpMul<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a, (typename ExprOp::value_type)b, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Scalar division
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2 >
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    Constant<T>,
    OpDiv<typename gettype<T, T2>::ty> > > operator / (const matrix<T, Container>& a, T2 b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      Constant<T>,
      OpDiv<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), Constant<T>((T)b, a.rows(), a.cols()), a.rows(), a.cols(), false));
    }

  template <class T, class Container, class T2 >
  Expr<
    BinExprOp<
    Constant<T>,
    typename matrix<T, Container>::const_iterator,
    OpDiv<typename gettype<T, T2>::ty> > > operator / (T2 a, const matrix<T, Container>& b)
    {
    typedef BinExprOp<
      Constant<T>,
      typename matrix<T, Container>::const_iterator,
      OpDiv<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<T>((T)a, b.rows(), b.cols()), b.begin(), b.rows(), b.cols(), false));
    }

  template <class ExprOp, class T2>
  Expr<
    BinExprOp<
    Expr<ExprOp>,
    Constant<typename ExprOp::value_type>,
    OpDiv<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator / (const Expr<ExprOp>& a, T2 b)
    {
    typedef BinExprOp<
      Expr<ExprOp>,
      Constant<typename ExprOp::value_type>,
      OpDiv<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a, Constant<typename ExprOp::value_type>((typename ExprOp::value_type)b, a.rows(), a.cols()), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp, class T2 >
  Expr<
    BinExprOp<
    Constant<typename ExprOp::value_type>,
    Expr<ExprOp>,
    OpDiv<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator / (T2 a, const Expr<ExprOp>& b)
    {
    typedef BinExprOp<
      Constant<typename ExprOp::value_type>,
      Expr<ExprOp>,
      OpDiv<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<typename ExprOp::value_type>((typename ExprOp::value_type)a, b.rows(), b.cols()), b, b.rows(), b.cols(), b.evaluate_before_assigning()));
    }




  template <class T, class T2>
  SparseExpr <
    BinExprSparseScalarOp<
    typename sparse_matrix<T>::const_iterator,
    OpDiv<typename gettype<T, T2>::ty> > > operator / (const sparse_matrix<T>& a, T2 b)
    {
    typedef BinExprSparseScalarOp<
      typename sparse_matrix<T>::const_iterator,
      OpDiv<typename gettype<T, T2>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a.begin(), (T)b, a.rows(), a.cols(), false));
    }

  template <class T, class T2>
  SparseExpr <
    BinExprSparseScalarOp<
    typename sparse_matrix<T>::const_iterator,
    OpDivInv<typename gettype<T, T2>::ty> > > operator / (T2 b, const sparse_matrix<T>& a)
    {
    typedef BinExprSparseScalarOp<
      typename sparse_matrix<T>::const_iterator,
      OpDivInv<typename gettype<T, T2>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a.begin(), (T)b, a.rows(), a.cols(), false));
    }

  template <class ExprOp, class T2>
  SparseExpr<
    BinExprSparseScalarOp<
    SparseExpr<ExprOp>,
    OpDiv<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator / (const SparseExpr<ExprOp>& a, T2 b)
    {
    typedef BinExprSparseScalarOp<
      SparseExpr<ExprOp>,
      OpDiv<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a, (typename ExprOp::value_type)b, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp, class T2>
  SparseExpr<
    BinExprSparseScalarOp<
    SparseExpr<ExprOp>,
    OpDivInv<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator / (T2 b, const SparseExpr<ExprOp>& a)
    {
    typedef BinExprSparseScalarOp<
      SparseExpr<ExprOp>,
      OpDivInv<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return SparseExpr<ExprT>(ExprT(a, (typename ExprOp::value_type)b, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Scalar addition
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2 >
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    Constant<T>,
    OpAdd<typename gettype<T, T2>::ty> > > operator + (const matrix<T, Container>& a, T2 b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      Constant<T>,
      OpAdd<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), Constant<T>((T)b, a.rows(), a.cols()), a.rows(), a.cols(), false));
    }

  template <class T, class Container, class T2 >
  Expr<
    BinExprOp<
    Constant<T>,
    typename matrix<T, Container>::const_iterator,
    OpAdd<typename gettype<T, T2>::ty> > > operator + (T2 a, const matrix<T, Container>& b)
    {
    typedef BinExprOp<
      Constant<T>,
      typename matrix<T, Container>::const_iterator,
      OpAdd<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<T>((T)a, b.rows(), b.cols()), b.begin(), b.rows(), b.cols(), false));
    }

  template <class ExprOp, class T2>
  Expr<
    BinExprOp<
    Expr<ExprOp>,
    Constant<typename ExprOp::value_type>,
    OpAdd<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator + (const Expr<ExprOp>& a, T2 b)
    {
    typedef BinExprOp<
      Expr<ExprOp>,
      Constant<typename ExprOp::value_type>,
      OpAdd<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a, Constant<typename ExprOp::value_type>((typename ExprOp::value_type)b, a.rows(), a.cols()), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp, class T2 >
  Expr<
    BinExprOp<
    Constant<typename ExprOp::value_type>,
    Expr<ExprOp>,
    OpAdd<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator + (T2 a, const Expr<ExprOp>& b)
    {
    typedef BinExprOp<
      Constant<typename ExprOp::value_type>,
      Expr<ExprOp>,
      OpAdd<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<typename ExprOp::value_type>((typename ExprOp::value_type)a, b.rows(), b.cols()), b, b.rows(), b.cols(), b.evaluate_before_assigning()));
    }


  ///////////////////////////////////////////////////////////////////////////////
  // Scalar subtraction
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2 >
  Expr<
    BinExprOp<
    typename matrix<T, Container>::const_iterator,
    Constant<T>,
    OpSub<typename gettype<T, T2>::ty> > > operator - (const matrix<T, Container>& a, T2 b)
    {
    typedef BinExprOp<
      typename matrix<T, Container>::const_iterator,
      Constant<T>,
      OpSub<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), Constant<T>((T)b, a.rows(), a.cols()), a.rows(), a.cols(), false));
    }

  template <class T, class Container, class T2 >
  Expr<
    BinExprOp<
    Constant<T>,
    typename matrix<T, Container>::const_iterator,
    OpSub<typename gettype<T, T2>::ty> > > operator - (T2 a, const matrix<T, Container>& b)
    {
    typedef BinExprOp<
      Constant<T>,
      typename matrix<T, Container>::const_iterator,
      OpSub<typename gettype<T, T2>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<T>((T)a, b.rows(), b.cols()), b.begin(), b.rows(), b.cols(), false));
    }

  template <class ExprOp, class T2>
  Expr<
    BinExprOp<
    Expr<ExprOp>,
    Constant<typename ExprOp::value_type>,
    OpSub<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator - (const Expr<ExprOp>& a, T2 b)
    {
    typedef BinExprOp<
      Expr<ExprOp>,
      Constant<typename ExprOp::value_type>,
      OpSub<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(a, Constant<typename ExprOp::value_type>((typename ExprOp::value_type)b, a.rows(), a.cols()), a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class ExprOp, class T2 >
  Expr<
    BinExprOp<
    Constant<typename ExprOp::value_type>,
    Expr<ExprOp>,
    OpSub<typename gettype<T2, typename ExprOp::value_type>::ty> > > operator - (T2 a, const Expr<ExprOp>& b)
    {
    typedef BinExprOp<
      Constant<typename ExprOp::value_type>,
      Expr<ExprOp>,
      OpSub<typename gettype<T2, typename ExprOp::value_type>::ty> > ExprT;
    return Expr<ExprT>(ExprT(Constant<typename ExprOp::value_type>((typename ExprOp::value_type)a, b.rows(), b.cols()), b, b.rows(), b.cols(), b.evaluate_before_assigning()));
    }


  ///////////////////////////////////////////////////////////////////////////////
  // Matrix multiplication
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class T2, class Container2 >
  Expr<
    MatMatMul<
    typename matrix<T, Container>::const_iterator,
    typename matrix<T2, Container2>::const_iterator > > operator * (const matrix<T, Container>& a, const matrix<T2, Container2>& b)
    {
    typedef MatMatMul<
      typename matrix<T, Container>::const_iterator,
      typename matrix<T2, Container2>::const_iterator > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a.begin(), b.begin(), a.rows(), b.cols(), a.cols()));
    }

  template <class T, class Container, class ExprOp >
  Expr<
    MatMatMul<
    typename matrix<T, Container>::const_iterator,
    Expr<ExprOp> > > operator * (const matrix<T, Container>& a, const Expr<ExprOp>& b)
    {
    typedef MatMatMul<
      typename matrix<T, Container>::const_iterator,
      Expr<ExprOp> > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a.begin(), b, a.rows(), b.cols(), a.cols()));
    }

  template <class T, class Container, class ExprOp >
  Expr<
    MatMatMul<
    Expr<ExprOp>,
    typename matrix<T, Container>::const_iterator > > operator * (const Expr<ExprOp>& a, const matrix<T, Container>& b)
    {
    typedef MatMatMul<
      Expr<ExprOp>,
      typename matrix<T, Container>::const_iterator > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a, b.begin(), a.rows(), b.cols(), a.cols()));
    }

  template <class ExprOp1, class ExprOp2 >
  Expr<
    MatMatMul<
    Expr<ExprOp1>,
    Expr<ExprOp2> > > operator * (const Expr<ExprOp1>& a, const Expr<ExprOp2>& b)
    {
    typedef MatMatMul<
      Expr<ExprOp1>,
      Expr<ExprOp2> > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a, b, a.rows(), b.cols(), a.cols()));
    }

  template <class T, class T2>
  T dot(const sparse_vector<T>& a, const sparse_vector<T2>& b)
    {
    double out = 0.0;
    const auto a_end = a.end();
    const auto b_end = b.end();
    auto a_it = a.begin();
    auto b_it = b.begin();
    while (a_it != a_end && b_it != b_end)
      {
      if (a_it.entry() == b_it.entry())
        {
        out += (double)(*a_it)*(double)(*b_it);
        ++a_it; ++b_it;
        }
      else if (a_it.entry() < b_it.entry())
        ++a_it;
      else
        ++b_it;
      }
    return (T)out;
    }

  template <class T, class T2>
  sparse_matrix<T> operator * (const sparse_matrix<T>& a, const sparse_matrix<T2>& b)
    {
    sparse_matrix<T> out(a.rows(), b.cols());
    for (uint64_t i = 0; i < a.rows(); ++i)
      {
      auto it = a.row(i).begin();
      auto it_end = it.end();
      for (; it != it_end; ++it)
        {
        uint64_t k = it.entry();
        auto it2 = b.row(k).begin();
        auto it2_end = it2.end();
        for (; it2 != it2_end; ++it2)
          {
          uint64_t j = it2.entry();
          out.put(i, j) += *it*(T)(*it2);
          }
        }
      }
    return out;
    }

  template <class ExprOp1, class T>
  sparse_matrix<typename ExprOp1::value_type> operator * (const SparseExpr<ExprOp1>& a, const sparse_matrix<T>& b)
    {
    sparse_matrix<typename ExprOp1::value_type> out(a.rows(), b.cols());
    sparse_matrix<typename ExprOp1::value_type> a1 = a;
    for (uint64_t i = 0; i < a1.rows(); ++i)
      {
      auto it = a1.row(i).begin();
      auto it_end = it.end();
      for (; it != it_end; ++it)
        {
        uint64_t k = it.entry();
        auto it2 = b.row(k).begin();
        auto it2_end = it2.end();
        for (; it2 != it2_end; ++it2)
          {
          uint64_t j = it2.entry();
          out.put(i, j) += *it*(typename ExprOp1::value_type)(*it2);
          }
        }
      }
    return out;
    }

  template <class ExprOp1, class T>
  sparse_matrix<T> operator * (const sparse_matrix<T>& a, const SparseExpr<ExprOp1>& b)
    {
    sparse_matrix<T> out(a.rows(), b.cols());
    sparse_matrix<typename ExprOp1::value_type> bb = b;
    for (uint64_t i = 0; i < a.rows(); ++i)
      {
      auto it = a.row(i).begin();
      auto it_end = it.end();
      for (; it != it_end; ++it)
        {
        uint64_t k = it.entry();
        auto it2 = bb.row(k).begin();
        auto it2_end = it2.end();
        for (; it2 != it2_end; ++it2)
          {
          uint64_t j = it2.entry();
          out.put(i, j) += *it*(T)(*it2);
          }
        }
      }
    return out;
    }

  template <class ExprOp1, class ExprOp2>
  sparse_matrix<typename ExprOp1::value_type> operator * (const SparseExpr<ExprOp1>& a, const SparseExpr<ExprOp2>& b)
    {
    sparse_matrix<typename ExprOp1::value_type> out(a.rows(), b.cols());
    sparse_matrix<typename ExprOp1::value_type> a1 = a;
    sparse_matrix<typename ExprOp2::value_type> bb = b;
    for (uint64_t i = 0; i < a1.rows(); ++i)
      {
      auto it = a1.row(i).begin();
      auto it_end = it.end();
      for (; it != it_end; ++it)
        {
        uint64_t k = it.entry();
        auto it2 = bb.row(k).begin();
        auto it2_end = it2.end();
        for (; it2 != it2_end; ++it2)
          {
          uint64_t j = it2.entry();
          out.put(i, j) += *it*(typename ExprOp1::value_type)(*it2);
          }
        }
      }
    return out;
    }

  /*
  template <class T, class T2, class Container2>
  matrix<T, Container2> operator *  (const sparse_matrix<T>& a, const matrix<T2, Container2>& b)
    {
    matrix<T, Container2> out(a.rows(), b.cols());
    for (uint64_t i = 0; i < a.rows(); ++i)
      {
      for (uint64_t j = 0; j < b.cols(); ++j)
        {
        T val = (T)0;
        auto it = a.row(i).begin();
        auto it_end = it.end();
        for (; it != it_end; ++it)
          {
          val += *it * b(it.entry(), j);
          }
        out(i, j) = val;
        }
      }
    return out;
    }
  */

  template <class T, class T2, class Container2>
  Expr<
    SparseMatMatMul<
    typename sparse_matrix<T>::const_iterator,
    typename matrix<T2, Container2>::const_iterator > > operator * (const sparse_matrix<T>& a, const matrix<T2, Container2>& b)
    {
    typedef SparseMatMatMul<
      typename sparse_matrix<T>::const_iterator,
      typename matrix<T2, Container2>::const_iterator > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a.begin(), b.begin(), a.rows(), b.cols(), a.cols()));
    }


  template <class ExprOp, class T2, class Container2>
  Expr<
    SparseMatMatMul<
    SparseExpr<ExprOp>,
    typename matrix<T2, Container2>::const_iterator > > operator * (const SparseExpr<ExprOp>& a, const matrix<T2, Container2>& b)
    {
    typedef SparseMatMatMul<
      SparseExpr<ExprOp>,
      typename matrix<T2, Container2>::const_iterator > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a, b.begin(), a.rows(), b.cols(), a.cols()));
    }

  template <class T, class ExprOp>
  Expr<
    SparseMatMatMul<
    typename sparse_matrix<T>::const_iterator,
    Expr<ExprOp> > > operator * (const sparse_matrix<T>& a, const Expr<ExprOp>& b)
    {
    typedef SparseMatMatMul<
      typename sparse_matrix<T>::const_iterator,
      Expr<ExprOp> > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a.begin(), b, a.rows(), b.cols(), a.cols()));
    }

  template <class ExprOp, class ExprOp2>
  Expr<
    SparseMatMatMul<
    SparseExpr<ExprOp>,
    Expr<ExprOp2>  > > operator * (const SparseExpr<ExprOp>& a, const Expr<ExprOp2>& b)
    {
    typedef SparseMatMatMul<
      SparseExpr<ExprOp>,
      Expr<ExprOp2> > ExprT;
    assert(a.cols() == b.rows());
    return Expr<ExprT>(ExprT(a, b, a.rows(), b.cols(), a.cols()));
    }

  template <class T, class Container, class T2, class Container2>
  T dot(const matrix<T, Container>& a, const matrix<T2, Container2>& b)
    {
    double out = 0.0;
    auto it = a.begin();
    auto it_end = a.end();
    auto it2 = b.begin();
    for (; it != it_end; ++it, ++it2)
      out += (double)(*it) * (double)(*it2);
    return (T)out;
    }

  template <class Container, class Container2>
  float dot(const matrix<float, Container>& a, const matrix<float, Container2>& b)
    {
    /*
    const float* ita = a.data();
    const float* itb = b.data();
    uint64_t len = a.rows()*a.cols();
    uint64_t len4 = len / 4;
    __m128d sum = _mm_setzero_pd();
    for (uint64_t i = 0; i < len4; ++i)
      {
      __m128 v1 = _mm_loadu_ps(ita + (i << 2));
      __m128 v2 = _mm_loadu_ps(itb + (i << 2));
      __m128 d = _mm_mul_ps(v1, v2);
      __m128d d1 = _mm_cvtps_pd(d);
      __m128d d2 = _mm_cvtps_pd(_mm_shuffle_ps(d, d, _MM_SHUFFLE(1, 0, 3, 2)));
      sum = _mm_add_pd(sum, d1);
      sum = _mm_add_pd(sum, d2);
      }
    double buffer[2];
    _mm_storeu_pd(buffer, sum);
    double totalsum = buffer[0] + buffer[1];
    for (uint64_t i = len4 * 4; i < len; ++i)
      {
      float v1 = *(ita + i);
      float v2 = *(itb + i);
      float d = v1 * v2;
      totalsum += (double)d;
      }
    return (float)totalsum;
    */
    const float* ita = a.data();
    const float* itb = b.data();
    uint64_t len = a.rows()*a.cols();
    uint64_t len4 = len / 4;
    double sum = 0.0;
    for (uint64_t i = 0; i < len4; ++i)
      {
      __m128 v1 = _mm_loadu_ps(ita + (i << 2));
      __m128 v2 = _mm_loadu_ps(itb + (i << 2));
      __m128 d = _mm_dp_ps(v1, v2, 0xf1);
      sum += (double)_mm_cvtss_f32(d);
      }
    for (uint64_t i = len4 * 4; i < len; ++i)
      {
      float v1 = *(ita + i);
      float v2 = *(itb + i);
      float d = v1 * v2;
      sum += (double)d;
      }
    return (float)sum;
    }

  template <class Container, class Container2>
  double dot(const matrix<double, Container>& a, const matrix<double, Container2>& b)
    {    
    const double* ita = a.data();
    const double* itb = b.data();
    uint64_t len = a.rows()*a.cols();
    uint64_t len2 = len / 2;
    __m128d sum = _mm_setzero_pd();
    for (uint64_t i = 0; i < len2; ++i)
      {
      __m128d v1 = _mm_loadu_pd(ita + (i << 1));
      __m128d v2 = _mm_loadu_pd(itb + (i << 1));
      __m128d d = _mm_mul_pd(v1, v2);     
      sum = _mm_add_pd(sum, d);      
      }
    double buffer[2];
    _mm_storeu_pd(buffer, sum);
    double totalsum = buffer[0] + buffer[1];
    for (uint64_t i = len2 * 2; i < len; ++i)
      {
      double v1 = *(ita + i);
      double v2 = *(itb + i);
      double d = v1 * v2;
      totalsum += (double)d;
      }
    return (float)totalsum;    
    }
  ///////////////////////////////////////////////////////////////////////////////
  // Transposing matrices
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container >
  Expr<
    Transpose<
    typename matrix<T, Container>::const_iterator > > transpose(const matrix<T, Container>& a)
    {
    typedef Transpose< typename matrix<T, Container>::const_iterator > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), a.cols(), a.rows()));
    }

  template <class ExprOp>
  Expr<
    Transpose<
    Expr<ExprOp> > > transpose(const Expr<ExprOp>& a)
    {
    typedef Transpose< Expr<ExprOp> > ExprT;
    return Expr<ExprT>(ExprT(a, a.cols(), a.rows()));
    }

  template <class T>
  sparse_matrix<T> transpose(const sparse_matrix<T>& a)
    {
    sparse_matrix<T> out(a.cols(), a.rows());
    const auto iter_end = a.end();
    auto iter = a.begin();
    while (iter != iter_end)
      {
      out.put(iter.second_entry(), iter.first_entry()) = *iter;
      ++iter;
      }
    return out;
    }

  template <class ExprOp>
  sparse_matrix<typename ExprOp::value_type> transpose(const SparseExpr<ExprOp>& a)
    {
    sparse_matrix<typename ExprOp::value_type> out(a.cols(), a.rows());
    const auto iter_end = a.end();
    auto iter = a;
    while (iter != iter_end)
      {
      typename ExprOp::value_type value = *iter;
      if (value)
        out.put(iter.second_entry(), iter.first_entry()) = value;
      ++iter;
      }
    return out;
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Diagonal of matrices
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container >
  Expr<
    Diagonal<
    typename matrix<T, Container>::const_iterator > > diagonal(const matrix<T, Container>& a)
    {
    typedef Diagonal< typename matrix<T, Container>::const_iterator > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), a.rows(), a.cols(), false));
    }

  template <class ExprOp>
  Expr<
    Diagonal<
    Expr<ExprOp> > > diagonal(const Expr<ExprOp>& a)
    {
    typedef Diagonal< Expr<ExprOp> > ExprT;
    return Expr<ExprT>(ExprT(a, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class T>
  Expr<
    DiagonalSparse<
    typename sparse_matrix<T>::const_iterator > > diagonal(const sparse_matrix<T>& a)
    {
    typedef DiagonalSparse< typename sparse_matrix<T>::const_iterator > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), a.rows(), a.cols(), false));
    }

  template <class ExprOp>
  Expr<
    DiagonalSparse<
    SparseExpr<ExprOp> > > diagonal(const SparseExpr<ExprOp>& a)
    {
    typedef DiagonalSparse< SparseExpr<ExprOp> > ExprT;
    return Expr<ExprT>(ExprT(a, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Block matrices
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container >
  Expr<
    Block<
    typename matrix<T, Container>::const_iterator > > block(const matrix<T, Container>& a, uint64_t i, uint64_t j, uint64_t rows, uint64_t cols)
    {
    typedef Block< typename matrix<T, Container>::const_iterator > ExprT;
    return Expr<ExprT>(ExprT(a.begin(), i, j, rows, cols, a.rows(), a.cols(), false));
    }

  template <class ExprOp>
  Expr<
    Block<
    Expr<ExprOp> > > block(const Expr<ExprOp>& a, uint64_t i, uint64_t j, uint64_t rows, uint64_t cols)
    {
    typedef Block< Expr<ExprOp> > ExprT;
    return Expr<ExprT>(ExprT(a, i, j, rows, cols, a.rows(), a.cols(), a.evaluate_before_assigning()));
    }

  template <class T >
  SparseExpr<
    BlockSparse<
    typename sparse_matrix<T>::const_iterator > > block(const sparse_matrix<T>& a, uint64_t i, uint64_t j, uint64_t rows, uint64_t cols)
    {
    typedef BlockSparse< typename sparse_matrix<T>::const_iterator > ExprT;
    return SparseExpr<ExprT>(ExprT(a.begin(), i, j, rows, cols, false));
    }

  template <class ExprOp>
  SparseExpr<
    BlockSparse<
    SparseExpr<ExprOp> > > block(const SparseExpr<ExprOp>& a, uint64_t i, uint64_t j, uint64_t rows, uint64_t cols)
    {
    typedef BlockSparse< SparseExpr<ExprOp> > ExprT;
    return SparseExpr<ExprT>(ExprT(a, i, j, rows, cols, a.evaluate_before_assigning()));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Matrix initializers (zeros, ones, identity)
  ///////////////////////////////////////////////////////////////////////////////

  inline Expr<Constant<double> > zeros(uint64_t rows, uint64_t cols)
    {
    return Expr<Constant<double> >(Constant<double>(0.0, rows, cols));
    }

  template <class T>
  inline Expr<Constant<T> > zeros(uint64_t rows, uint64_t cols)
    {
    return Expr<Constant<T> >(Constant<T>((T)0, rows, cols));
    }

  inline Expr<Constant<double> > ones(uint64_t rows, uint64_t cols)
    {
    return Expr<Constant<double> >(Constant<double>(1.0, rows, cols));
    }

  template <class T>
  inline Expr<Constant<T> > ones(uint64_t rows, uint64_t cols)
    {
    return Expr<Constant<T> >(Constant<T>((T)1, rows, cols));
    }

  inline Expr<Identity<double> > identity(uint64_t rows, uint64_t cols)
    {
    return Expr<Identity<double> >(Identity<double>(rows, cols));
    }

  template <class T>
  inline Expr<Identity<T> > identity(uint64_t rows, uint64_t cols)
    {
    return Expr<Identity<T> >(Identity<T>(rows, cols));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Norm
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container >
  T norm_sqr(const matrix<T, Container>& a)
    {
    double sum = 0.0;
    for (const auto& value : a)
      sum += value * value;
    return (T)sum;
    }

  template <class Container>
  float norm_sqr(const matrix<float, Container>& a)
    {
    const float* ita = a.data();
    uint64_t len = a.rows()*a.cols();
    uint64_t len4 = len / 4;
    double sum = 0.0;
    for (uint64_t i = 0; i < len4; ++i)
      {
      __m128 v1 = _mm_loadu_ps(ita + (i << 2));
      __m128 d = _mm_dp_ps(v1, v1, 0xf1);
      sum += (double)_mm_cvtss_f32(d);
      }
    for (uint64_t i = len4 * 4; i < len; ++i)
      {
      float v1 = *(ita + i);
      float d = v1 * v1;
      sum += (double)d;
      }
    return (float)sum;
    }

  template <class Container>
  double norm_sqr(const matrix<double, Container>& a)
    {
    const double* ita = a.data();
    uint64_t len = a.rows()*a.cols();
    uint64_t len2 = len / 2;
    __m128d sum = _mm_setzero_pd();
    for (uint64_t i = 0; i < len2; ++i)
      {
      __m128d v1 = _mm_loadu_pd(ita + (i << 1));
      __m128d d = _mm_mul_pd(v1, v1);
      sum = _mm_add_pd(sum, d);
      }
    double buffer[2];
    _mm_storeu_pd(buffer, sum);
    double totalsum = buffer[0] + buffer[1];
    for (uint64_t i = len2 * 2; i < len; ++i)
      {
      double v1 = *(ita + i);
      double d = v1 * v1;
      totalsum += (double)d;
      }
    return (float)totalsum;
    }

  template <class ExprOp>
  double norm_sqr(Expr<ExprOp> expr)
    {
    uint64_t sz = expr.rows()*expr.cols();
    if (sz == 0)
      return 0.0;
    auto value = *expr;
    double sum = (double)value*value;
    for (uint64_t i = 1; i < sz; ++i)
      {
      ++expr;
      value = *expr;
      sum += (double)value*value;
      }
    return sum;
    }

  template <class T >
  T norm_sqr(const sparse_matrix<T>& a)
    {
    double sum = 0.0;
    for (const auto& value : a)
      sum += value * value;
    return (T)sum;
    }

  template <class ExprOp>
  double norm_sqr(SparseExpr<ExprOp> expr)
    {
    double sum = 0.0;
    auto expr_end = expr.end();
    while (expr != expr_end)
      {
      auto value = *expr;
      sum += (double)value*value;
      ++expr;
      }
    return sum;
    }

  template <class T, class Container >
  T norm(const matrix<T, Container>& a)
    {
    return std::sqrt(norm_sqr(a));
    }

  template <class ExprOp>
  double norm(Expr<ExprOp> expr)
    {
    return std::sqrt(norm_sqr(expr));
    }

  template <class T >
  T norm(const sparse_matrix<T>& a)
    {
    return std::sqrt(norm_sqr(a));
    }

  template <class ExprOp>
  double norm(SparseExpr<ExprOp> expr)
    {
    return std::sqrt(norm_sqr(expr));
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Trace
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container >
  T trace(const matrix<T, Container>& a)
    {
    auto it = diagonal(a);
    uint64_t sz = std::min(a.rows(), a.cols());
    double sum = (double)*it;
    for (uint64_t i = 1; i < sz; ++i)
      {
      ++it;
      sum += (double)*it;
      }
    return (T)sum;
    }

  template <class ExprOp>
  double trace(Expr<ExprOp> expr)
    {
    auto it = diagonal(expr);
    uint64_t sz = std::min(expr.rows(), expr.cols());
    double sum = *it;
    for (uint64_t i = 1; i < sz; ++i)
      {
      ++it;
      sum += (double)*it;
      }
    return sum;
    }

  template <class T >
  T trace(const sparse_matrix<T>& a)
    {
    auto it = diagonal(a);
    uint64_t sz = std::min(a.rows(), a.cols());
    double sum = (double)*it;
    for (uint64_t i = 1; i < sz; ++i)
      {
      ++it;
      sum += (double)*it;
      }
    return (T)sum;
    }

  template <class ExprOp>
  double trace(SparseExpr<ExprOp> expr)
    {
    auto it = diagonal(expr);
    uint64_t sz = std::min(expr.rows(), expr.cols());
    double sum = *it;
    for (uint64_t i = 1; i < sz; ++i)
      {
      ++it;
      sum += (double)*it;
      }
    return sum;
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Singular value decomposition
  ///////////////////////////////////////////////////////////////////////////////

  namespace svd_details
    {
    template <class T>
    T sign(T a, T b)
      {
      if (b >= (T)0)
        return std::abs(a);
      else
        return -std::abs(a);
      }

    template <class Matrix_mxn>
    void swap_cols(int m, Matrix_mxn& A, size_t i, size_t j)
      {
      for (int row = 0; row < m; ++row)
        std::swap(A[row][i], A[row][j]);
      }

    template <class T>
    T pythag(T a, T b)
      /* compute (a2 + b2)^1/2 without destructive underflow or overflow */
      {
      T absa, absb;
      absa = std::abs(a);
      absb = std::abs(b);
      if (absa > absb) return absa * std::sqrt((T)1 + (absb / absa)*(absb / absa));
      else return (absb == (T)0 ? (T)0 : absb * std::sqrt((T)1 + (absa / absb)*(absa / absb)));
      }

    template <class T, class Matrix_mxn, class Matrix_nxn, class Vector_n>
    bool _svd(int m, int n, Matrix_mxn& a, Vector_n& w, Matrix_nxn& v, Vector_n& rv1)
      {
      int flag, i, its, j, jj, k, l = 0, nm = 0;
      T anorm, c, f, g, h, s, scale, x, y, z;
      g = scale = anorm = (T)0; /* Householder reduction to bidiagonal form */
      for (i = 0; i < n; ++i)
        {
        l = i + 1;
        rv1(i) = scale * g;
        g = s = scale = (T)0;
        if (i < m)
          {
          for (k = i; k < m; ++k)
            scale += std::abs(a[k][i]);
          if (scale)
            {
            for (k = i; k < m; ++k)
              {
              a[k][i] /= scale;
              s += a[k][i] * a[k][i];
              }
            f = a[i][i];
            g = -sign(std::sqrt(s), f);
            h = f * g - s;
            a[i][i] = f - g;
            for (j = l; j < n; ++j)
              {
              for (s = (T)0, k = i; k < m; ++k)
                s += a[k][i] * a[k][j];
              f = s / h;
              for (k = i; k < m; ++k)
                a[k][j] += f * a[k][i];
              }
            for (k = i; k < m; ++k)
              a[k][i] *= scale;
            }
          }
        w(i) = scale * g;
        g = s = scale = (T)0;
        if (i < m && i != n - 1) {
          for (k = l; k < n; ++k)
            scale += std::abs(a[i][k]);
          if (scale)
            {
            for (k = l; k < n; ++k)
              {
              a[i][k] /= scale;
              s += a[i][k] * a[i][k];
              }
            f = a[i][l];
            g = -sign(std::sqrt(s), f);
            h = f * g - s;
            a[i][l] = f - g;
            for (k = l; k < n; ++k)
              rv1(k) = a[i][k] / h;
            for (j = l; j < m; ++j)
              {
              for (s = (T)0, k = l; k < n; ++k)
                s += a[j][k] * a[i][k];
              for (k = l; k < n; ++k)
                a[j][k] += s * rv1(k);
              }
            for (k = l; k < n; ++k)
              a[i][k] *= scale;
            }
          }
        anorm = std::max<T>(anorm, (std::abs(w(i)) + std::abs(rv1(i))));
        }
      for (i = n - 1; i >= 0; --i) { /* Accumulation of right-hand transformations. */
        if (i < n - 1)
          {
          if (g)
            {
            for (j = l; j < n; ++j) /* Double division to avoid possible underflow. */
              v[j][i] = (a[i][j] / a[i][l]) / g;
            for (j = l; j < n; ++j)
              {
              for (s = (T)0, k = l; k < n; ++k)
                s += a[i][k] * v[k][j];
              for (k = l; k < n; ++k)
                v[k][j] += s * v[k][i];
              }
            }
          for (j = l; j < n; ++j)
            v[i][j] = v[j][i] = (T)0;
          }
        v[i][i] = (T)1;
        g = rv1(i);
        l = i;
        }
      for (i = std::min<int>(m - 1, n - 1); i >= 0; --i)
        { /* Accumulation of left-hand transformations. */
        l = i + 1;
        g = w(i);
        for (j = l; j < n; ++j)
          a[i][j] = (T)0;
        if (g)
          {
          g = (T)1 / g;
          for (j = l; j < n; ++j)
            {
            for (s = (T)0, k = l; k < m; ++k)
              s += a[k][i] * a[k][j];
            f = (s / a[i][i])*g;
            for (k = i; k < m; ++k)
              a[k][j] += f * a[k][i];
            }
          for (j = i; j < m; ++j)
            a[j][i] *= g;
          }
        else for (j = i; j < m; ++j)
          a[j][i] = (T)0;
        ++a[i][i];
        }
      for (k = n - 1; k >= 0; --k) { /* Diagonalization of the bidiagonal form. */
        for (its = 0; its < 30; ++its)
          {
          flag = 1;
          for (l = k; l >= 0; --l)
            { /* Test for splitting. */
            nm = l - 1; /* Note that rv1(0] is always zero. */
            if ((T)(std::abs(rv1(l)) + anorm) == anorm)
              {
              flag = 0;
              break;
              }
            if ((T)(std::abs(w(nm)) + anorm) == anorm) break;
            }
          if (flag) {
            c = (T)0; /* Cancellation of rv1(l), if l > 0. */
            s = (T)1;
            for (i = l; i <= k; ++i)
              {
              f = s * rv1(i);
              rv1(i) = c * rv1(i);
              if ((T)(std::abs(f) + anorm) == anorm) break;
              g = w(i);
              h = pythag(f, g);
              w(i) = h;
              h = (T)1 / h;
              c = g * h;
              s = -f * h;
              for (j = 0; j < m; ++j)
                {
                y = a[j][nm];
                z = a[j][i];
                a[j][nm] = y * c + z * s;
                a[j][i] = z * c - y * s;
                }
              }
            }
          z = w(k);
          if (l == k)
            { /* Convergence. */
            if (z < (T)0)
              { /* Singular value is made nonnegative. */
              w(k) = -z;
              for (j = 0; j < n; ++j) v[j][k] = -v[j][k];
              }
            break;
            }
          if (its >= 30)
            {
            return false;
            }
          x = w(l); /* Shift from bottom 2-by-2 minor. */
          nm = k - 1;
          y = w(nm);
          g = rv1(nm);
          h = rv1(k);
          f = ((y - z)*(y + z) + (g - h)*(g + h)) / ((T)2 * h*y);
          g = pythag(f, (T)1);
          f = ((x - z)*(x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
          c = s = (T)1; /* Next QR transformation: */
          for (j = l; j <= nm; ++j)
            {
            i = j + 1;
            g = rv1(i);
            y = w(i);
            h = s * g;
            g = c * g;
            z = pythag(f, h);
            rv1(j) = z;
            c = f / z;
            s = h / z;
            f = x * c + g * s;
            g = g * c - x * s;
            h = y * s;
            y *= c;
            for (jj = 0; jj < n; ++jj)
              {
              x = v[jj][j];
              z = v[jj][i];
              v[jj][j] = x * c + z * s;
              v[jj][i] = z * c - x * s;
              }
            z = pythag(f, h);
            w(j) = z; /* Rotation can be arbitrary if z = 0. */
            if (z) {
              z = (T)1 / z;
              c = f * z;
              s = h * z;
              }
            f = c * g + s * y;
            x = c * y - s * g;
            for (jj = 0; jj < m; ++jj)
              {
              y = a[jj][j];
              z = a[jj][i];
              a[jj][j] = y * c + z * s;
              a[jj][i] = z * c - y * s;
              }
            }
          rv1(l) = (T)0;
          rv1(k) = f;
          w(k) = x;
          }
        }
      return true;
      }

    template <class T, class Matrix_mxn, class Matrix_nxn, class Vector_n>
    bool _svd_sorted(int m, int n, Matrix_mxn& a, Vector_n& w, Matrix_nxn& v, Vector_n& rv1)
      {
      bool res = _svd<T, Matrix_mxn, Matrix_nxn, Vector_n>(m, n, a, w, v, rv1);
      if (res)
        {
        std::vector<std::pair<T, int>> w_idx_sorted;
        w_idx_sorted.reserve(n);
        for (int i = 0; i < n; ++i)
          w_idx_sorted.push_back(std::make_pair(w(i), i));
        std::sort(w_idx_sorted.begin(), w_idx_sorted.end(), [](const std::pair<T, int>& p0, const std::pair<T, int>& p1) {return p0.first > p1.first; });
        for (int i = 0; i < n; ++i)
          {
          int j = w_idx_sorted[i].second;
          while (j < i)
            j = w_idx_sorted[j].second;
          if (j == i)
            continue;
          std::swap(w(i), w(j));
          swap_cols(m, a, i, j);
          swap_cols(n, v, i, j);
          }
        }
      return res;
      }
    }

  template <class T, class Container_mn, class Container_n, class Container_nn>
  bool svd(matrix<T, Container_mn>& a, matrix<T, Container_n>& w, matrix<T, Container_nn>& v)
    {
    using namespace ::jtk::svd_details;
    auto m = a.rows();
    auto n = a.cols();
    w.resize(n, 1);
    v.resize(n, n);
    matrix<T, Container_n> rv1(n);
    bool res = _svd_sorted<T>((int)m, (int)n, a, w, v, rv1);
    return res;
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Pseudo inverse
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container>
  int pseudo_inverse(matrix<T, Container>& inv, matrix<T, Container>& a, T tol)
    {
    using namespace ::jtk::svd_details;
    int p = 0;
    matrix<T, Container> w, v;
    bool res = svd(a, w, v);
    if (res)
      {
      auto m = a.rows();
      auto n = a.cols();
      inv.resize(n, m);
      while ((p < n) && (std::abs(w(p)) > tol))
        ++p;
      for (int r = 0; r < n; ++r)
        {
        for (int c = 0; c < m; ++c)
          {
          T val = (T)0;
          for (int k = 0; k < p; ++k)
            val += v[r][k] * a[c][k] / w(k);
          inv[r][c] = val;
          }
        }
      }
    return p;
    }

  ///////////////////////////////////////////////////////////////////////////////
  // linear Least Squares problem solver using singular value Decomposition
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container, class Container2>
  void lsd(matrix<T, Container>& A, matrix<T, Container2>& b, T tol)
    {
    matrix<T, Container2> sigma;
    matrix<T, Container> V;
    svd(A, sigma, V);
    matrix < T, Container2 > res;
    res.noalias() = transpose(A)*b;
    for (uint64_t i = 0; i < res.rows(); ++i)
      {
      if (std::abs(sigma(i)) < tol)
        res(i) = (T)0;
      else
        res(i) /= sigma(i);
      }
    b.noalias() = V * res;
    }

  ///////////////////////////////////////////////////////////////////////////////
  // LU decomposition
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container>
  void ludcmp(matrix<T, Container>& A, std::vector<uint64_t>& permutations, T& d)
    {
    assert(A.rows() == A.cols());
    uint64_t n = A.rows();
    uint64_t i, imax = 0, j, k;
    T big, dum, sum, temp;
    std::vector<T> vv;
    vv.resize(n);
    permutations.resize(n);
    d = 1;
    for (i = 0; i < n; ++i)
      {
      big = 0;
      for (j = 0; j < n; ++j)
        {
        temp = std::abs(A(i, j));
        if (temp > big)
          big = temp;
        }
      if (big == 0)
        {
        return; // singular matrix
        }
      vv[i] = 1 / big;
      }
    for (j = 0; j < n; ++j)
      {
      for (i = 0; i < j; ++i)
        {
        sum = A(i, j);
        for (k = 0; k < i; ++k)
          sum -= A(i, k)*A(k, j);
        A(i, j) = sum;
        }
      big = 0;
      for (i = j; i < n; ++i)
        {
        sum = A(i, j);
        for (k = 0; k < j; ++k)
          sum -= A(i, k)*A(k, j);
        A(i, j) = sum;
        if ((dum = vv[i] * std::abs(sum)) >= big)
          {
          big = dum;
          imax = i;
          }
        }
      if (j != imax)
        {
        for (k = 0; k < n; ++k)
          {
          dum = A(imax, k);
          A(imax, k) = A(j, k);
          A(j, k) = dum;
          }
        d = -d;
        vv[imax] = vv[j];
        }
      permutations[j] = imax;
      if (A(j, j) == 0)
        A(j, j) = std::numeric_limits<T>::epsilon();
      if (j != (n - 1))
        {
        dum = 1 / A(j, j);
        for (i = j + 1; i < n; ++i)
          A(i, j) *= dum;
        }
      }
    }

  template <class T, class Container, class Container2>
  void lubksb(matrix<T, Container2>& b, const matrix<T, Container>& A, const std::vector<uint64_t>& permutations)
    {
    assert(A.rows() == A.cols());
    uint64_t i, ii = 0xffffffffffffffff, ip, j;
    T sum;
    uint64_t n = A.rows();

    for (i = 0; i < n; ++i)
      {
      ip = permutations[i];
      sum = b(ip);
      b(ip) = b(i);
      if (ii < 0xffffffffffffffff)
        for (j = ii; j < i; ++j)
          sum -= A(i, j)*b(j);
      else if (sum)
        ii = i;
      b(i) = sum;
      }
    for (i = n; i > 0; --i)
      {
      sum = b(i - 1);
      for (j = i; j < n; ++j)
        sum -= A(i - 1, j)*b(j);
      b(i - 1) = sum / A(i - 1, i - 1);
      }
    }

  template <class T, class Container, class Container2>
  void solve(matrix<T, Container>& A, matrix<T, Container2>& b)
    {
    std::vector<uint64_t> permutations;
    T d;
    ludcmp(A, permutations, d);
    lubksb(b, A, permutations);
    }

  template <class T, class Container>
  void invert(matrix<T, Container>& Ainv, const matrix<T, Container>& A)
    {
    matrix<T, Container> m(A);
    Ainv.resize(A.rows(), A.cols());
    std::vector<uint64_t> permutations;
    T d;
    ludcmp(m, permutations, d);
    for (uint64_t j = 0; j < A.cols(); ++j)
      {
      matrix<T, Container> col = zeros(A.rows(), 1);
      col(j) = 1;
      lubksb(col, m, permutations);
      for (uint64_t i = 0; i < A.rows(); ++i)
        Ainv(i, j) = col(i);
      }
    }

  template <class T, class Container>
  T determinant(const matrix<T, Container>& A)
    {
    matrix<T, Container> m(A);
    std::vector<uint64_t> permutations;
    T d;
    ludcmp(m, permutations, d);
    for (uint64_t j = 0; j < A.cols(); ++j)
      d *= m(j, j);
    return d;
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Cholesky decomposition
  ///////////////////////////////////////////////////////////////////////////////

  namespace cholesky_details
    {
    // Divides a row by its diagonal element
    template <class T, class Container>
    void _rdiv(matrix<T, Container>& A, uint64_t r)
      {
      T* it = A[r] + r;
      T denom = static_cast<T>(std::sqrt(*it));
      *it /= denom;
      for (uint64_t i = r + 1; i < A.cols(); ++i)
        {
        ++it;
        *it /= denom;
        }
      }

    // Replaces row r1 by
    // row r1 = row r1 - A(r2,r1)/A(r2,r2)*row r2
    template <class T, class Container>
    void _rmod(matrix<T, Container>& A, uint64_t r1, uint64_t r2)
      {
      T factor = A(r2, r1) / A(r2, r2);
      if (factor == (T)0)
        return;
      T* it1 = A[r1] + r1;
      T* it2 = A[r2] + r1;
      for (uint64_t i = r1; i < A.cols(); ++i)
        {
        *it1 -= factor * (*it2);
        ++it1;
        ++it2;
        }
      }
    }

  template <class T, class Container>
  void cholesky(matrix<T, Container>& A)
    {
    assert(A.rows() == A.cols());
    using namespace ::jtk::cholesky_details;
    for (uint64_t k = 0; k < A.rows(); ++k)
      {
      for (uint64_t j = k + 1; j < A.rows(); ++j)
        {
        _rmod(A, j, k);
        }
      _rdiv(A, k);
      }
    }

  template <class T, class Container, class Container2>
  void back_substitution(matrix<T, Container2>& x, const matrix<T, Container2>& b, const matrix<T, Container>& A)
    {
    assert(A.rows() == A.cols());
    x.resize(A.rows(), 1);
    for (uint64_t j = A.rows(); j > 0; --j)
      {
      uint64_t ind = j - 1;
      T val = b(ind);
      auto it = A[ind] + ind;
      const auto it_end = A[ind] + A.cols();
      T diag = *it;
      ++it;
      auto xit = x[ind + 1];
      for (; it != it_end; ++it, ++xit)
        {
        val -= *xit*(*it);
        }
      val /= diag;
      x(ind) = val;
      }
    }

  template <class T, class Container, class Container2>
  void forward_substitution(matrix<T, Container2>& x, const matrix<T, Container2>& b, const matrix<T, Container>& A)
    {
    assert(A.rows() == A.cols());
    x.resize(A.rows(), 1);
    for (uint64_t j = 0; j < A.rows(); ++j)
      {
      T val = b(j);
      T diag = A(j, j);
      auto it = A[j];
      const auto it_end = A[j] + j;
      auto xit = x[0];
      for (; it != it_end; ++it, ++xit)
        {
        val -= *xit*(*it);
        }
      val /= diag;
      x(j) = val;
      }
    }

  template <class T, class Container, class Container2>
  void solve_cholesky(matrix<T, Container>& A, matrix<T, Container2>& b)
    {
    cholesky(A);
    matrix<T, Container2> y;
    matrix<T, Container> At;
    At.noalias() = transpose(A);
    forward_substitution(y, b, At);
    back_substitution(b, y, A);
    }

  ///////////////////////////////////////////////////////////////////////////////
  // QR decomposition
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container>
  void qrdcmp(matrix<T, Container>& a, matrix<T, Container>& qt, std::vector<uint64_t>& permutations, bool pivot)
    {
    uint64_t m(a.rows());
    uint64_t n(a.cols());
    uint64_t mn = m > n ? n : m;
    matrix<T, Container> rdiag;
    matrix<T, Container> acnorm;
    qrfac(a, pivot, permutations, rdiag, acnorm);
    qt = identity(m, m);
    for (uint64_t k = 0; k < mn; ++k)
      {
      if (a(k, k) != 0.0)
        {
        for (uint64_t j = 0; j < m; ++j)
          {
          T sum = (T)0;
          for (uint64_t i = k; i < m; ++i)
            sum += a(i, k) * qt(i, j);
          sum /= a(k, k);
          for (uint64_t i = k; i < m; ++i)
            qt(i, j) -= sum * a(i, k);
          }
        }
      }

    for (uint64_t i = 0; i < mn; ++i)
      {
      a(i, i) = rdiag(i);
      for (uint64_t j = 0; j < i; ++j)
        a(i, j) = (T)0;
      }
    for (uint64_t i = mn; i < m; ++i)
      {
      for (uint64_t j = 0; j < n; ++j)
        a(i, j) = (T)0;
      }
    }

  template <class T, class Container, class Container2>
  void solve_qr(matrix<T, Container>& a, matrix<T, Container2>& b)
    {
    std::vector<uint64_t> permutations;
    matrix<T, Container2> rdiag;
    matrix<T, Container2> acnorm;
    qrfac(a, true, permutations, rdiag, acnorm);
    uint64_t m = a.rows();
    uint64_t n = a.cols();
    matrix<T, Container2> qtb(n, 1);
    for (uint64_t j = 0; j < n; ++j)
      {
      if (a[j][j] != (T)0)
        {
        T sum = (T)0;
        for (uint64_t i = j; i < m; ++i)
          sum += a[i][j] * b(i);
        T temp = -sum / a[j][j];
        for (uint64_t i = j; i < m; ++i)
          b(i) += a[i][j] * temp;
        }
      a[j][j] = rdiag(j);
      qtb(j) = b(j);
      }
    matrix<T, Container2> sdiag;
    matrix<T, Container2> diag = zeros(n, 1);
    qrsolv(a, permutations, diag, qtb, b, sdiag);
    }

  template <class T, class Container, class Container2>
  void qrfac(matrix<T, Container>& a, bool pivot, std::vector<uint64_t>& ipvt,
    matrix<T, Container2>& rdiag, matrix<T, Container2>& acnorm)
    {
    uint64_t m = a.rows();
    uint64_t n = a.cols();

    if (pivot)
      ipvt.resize(n);
    rdiag.resize(n, 1);
    acnorm.resize(n, 1);
    matrix<T, Container> wa(n, 1);

    uint64_t i, j, jp1, k, kmax, minmn;
    T ajnorm, epsmch, sum, temp;

    /* get machine precision */
    epsmch = std::numeric_limits<T>::epsilon();
    /* compute the initial column norms and initialize several arrays */
    for (j = 0; j < n; j++)
      {
      acnorm(j) = (T)norm(block(a, 0, j, m, 1));
      rdiag(j) = acnorm(j);
      wa(j) = rdiag(j);
      if (pivot)
        ipvt[j] = j;
      }
    /* reduce a to r with householder transformations */
    minmn = (m < n) ? m : n;
    for (j = 0; j < minmn; ++j)
      {
      if (pivot)
        {
        /* bring column with largest norm into the pivot position */
        kmax = j;
        for (k = j; k < n; ++k)
          if (rdiag(k) > rdiag(kmax))
            kmax = k;
        if (kmax != j)
          {
          for (i = 0; i < m; ++i)
            {
            temp = a[i][j];
            a[i][j] = a[i][kmax];
            a[i][kmax] = temp;
            }
          rdiag(kmax) = rdiag(j);
          wa(kmax) = wa(j);
          k = ipvt[j];
          ipvt[j] = ipvt[kmax];
          ipvt[kmax] = k;
          }
        }
      /* compute the householder transformation */
      ajnorm = (T)norm(block(a, j, j, m - j, 1));
      if (ajnorm != (T)0)
        {
        if (a[j][j] < (T)0)
          ajnorm = -ajnorm;
        for (i = j; i < m; ++i)
          a[i][j] /= ajnorm;
        a[j][j] += 1.0;
        jp1 = j + 1;
        if (n > jp1)
          {
          for (k = jp1; k < n; ++k)
            {
            sum = (T)0;
            for (i = j; i < m; ++i)
              sum += a[i][j] * a[i][k];
            temp = sum / a[j][j];
            for (i = j; i < m; ++i)
              a[i][k] -= temp * a[i][j];
            if (!pivot || !rdiag(k))
              continue;
            temp = a[j][k] / rdiag(k);
            rdiag(k) *= std::sqrt(std::max<T>(0, 1 - temp * temp));
            if ((T)0.5 * (rdiag(k) * rdiag(k) / (wa(k) * wa(k))) > epsmch)
              continue;
            rdiag(k) = (T)norm(block(a, jp1, k, m - jp1, 1));
            wa(k) = rdiag(k);
            }
          }
        }
      rdiag(j) = -ajnorm;
      }
    }

  template <class T, class Container, class Container2>
  void qrsolv(matrix<T, Container>& r, const std::vector<uint64_t>& ipvt, const matrix<T, Container2>& diag,
    const matrix<T, Container2>& qtb, matrix<T, Container2>& x, matrix<T, Container2>& sdiag)
    {
    uint64_t i, j, jp1, k, kp1, l, nsing;
    T qtbpj, sum, temp, dsin, dcos, dtan, dcotan;

    uint64_t n = r.cols();
    matrix<T, Container2> wa(n, 1);
    x.resize(n, 1);
    sdiag.resize(n, 1);

    for (j = 0; j < n; ++j)
      {
      for (i = j; i < n; ++i)
        r[i][j] = r[j][i];
      x(j) = r[j][j];
      wa(j) = qtb(j);
      }
    for (j = 0; j < n; ++j)
      {
      l = ipvt[j];
      if (diag(l) != (T)0)
        {
        for (k = j; k < n; ++k)
          sdiag(k) = (T)0;
        sdiag(j) = diag(l);
        qtbpj = 0.0;
        for (k = j; k < n; ++k)
          {
          if (sdiag(k) == (T)0)
            continue;
          if (std::abs(r[k][k]) < std::abs(sdiag(k)))
            {
            dcotan = r[k][k] / sdiag(k);
            dsin = (T)1 / std::sqrt((T)1 + dcotan * dcotan);
            dcos = dsin * dcotan;
            }
          else
            {
            dtan = sdiag(k) / r[k][k];
            dcos = (T)1 / std::sqrt((T)1 + dtan * dtan);
            dsin = dcos * dtan;
            }
          r[k][k] = dcos * r[k][k] + dsin * sdiag(k);
          temp = dcos * wa(k) + dsin * qtbpj;
          qtbpj = -dsin * wa(k) + dcos * qtbpj;
          wa(k) = temp;
          kp1 = k + 1;
          if (n <= kp1)
            continue;
          for (i = kp1; i < n; ++i)
            {
            temp = dcos * r[i][k] + dsin * sdiag(i);
            sdiag(i) = -dsin * r[i][k] + dcos * sdiag(i);
            r[i][k] = temp;
            }
          }
        }
      sdiag(j) = r[j][j];
      r[j][j] = x(j);
      }
    nsing = n;
    for (j = 0; j < n; ++j)
      {
      if ((sdiag(j) == (T)0) && (nsing == n))
        nsing = j;
      if (nsing < n)
        wa(j) = (T)0;
      }
    if (nsing >= 1)
      {
      for (k = 0; k < nsing; ++k)
        {
        j = nsing - k - 1;
        sum = (T)0;
        jp1 = j + 1;
        if (nsing > jp1)
          {
          for (i = jp1; i < nsing; ++i)
            sum += r[i][j] * wa(i);
          }
        wa(j) = (wa(j) - sum) / sdiag(j);
        }
      }
    for (j = 0; j < n; ++j)
      {
      l = ipvt[j];
      x(l) = wa(j);
      }
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Levenberg-Marquardt optimization
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container_mn, class Container_n, class Container_m>
  bool fdjac(bool(*f)(const matrix<T, Container_n>&, matrix<T, Container_m>&, void*), matrix<T, Container_n>& x, matrix<T, Container_m>& fvec, matrix<T, Container_mn>& fjac, T epsfcn, void* user_data)
    {
    uint64_t m = fvec.rows();
    uint64_t n = x.rows();
    matrix<T, Container_m> wa(m, 1);
    uint64_t i, j;
    T eps, epsmch, h, temp;

    fjac.resize(m, n);

    epsmch = (epsfcn > std::numeric_limits<T>::epsilon()) ? epsfcn : std::numeric_limits<T>::epsilon();
    eps = std::sqrt(epsmch);

    for (j = 0; j < n; ++j)
      {
      temp = x(j);
      if (temp == (T)0)
        h = eps;
      else
        h = eps * std::abs(temp);
      x(j) = temp + h;
      bool res = f(x, wa, user_data);
      if (!res)
        return false;
      x(j) = temp;
      for (i = 0; i < m; ++i)
        fjac[i][j] = (wa(i) - fvec(i)) / h;
      }
    return true;
    }

  namespace lm_details
    {

    template <class T>
    struct lm_limits
      {
      static T minimum()
        {
        return 0;
        }
      };

    template <>
    struct lm_limits<double>
      {
      static double minimum()
        {
        return 2.225073858507201e-308;
        }
      };

    template <>
    struct lm_limits<float>
      {
      static float minimum()
        {
        return 1.1754942107e-38f;
        }
      };

    template <class T, class Container_mn, class Container_n>
    void lmpar(matrix<T, Container_mn>& r, const std::vector<uint64_t>& ipvt, const matrix<T, Container_n>& diag, const matrix<T, Container_n>& qtb,
      T delta, T& par, matrix<T, Container_n>& x, matrix<T, Container_n>& sdiag)
      {
      uint64_t n = r.cols();
      uint64_t i, iter, j, jp1, k, l, nsing;
      T dxnorm, dwarf, fp, gnorm, parc, parl, paru;
      T sum, temp;

      matrix<T, Container_n> wa1(n, 1), wa2(n, 1);

      dwarf = lm_limits<T>::minimum();
      nsing = n;
      for (j = 0; j < n; ++j)
        {
        wa1(j) = qtb(j);
        if ((r[j][j] == (T)0) && (nsing == n))
          nsing = j;
        if (nsing < n)
          wa1(j) = (T)0;
        }
      if (nsing >= 1)
        {
        for (k = 0; k < nsing; ++k)
          {
          j = nsing - k - 1;
          wa1(j) /= r[j][j];
          temp = wa1(j);
          if (j < 1)
            continue;
          for (i = 0; i < j; ++i)
            wa1(i) -= r[i][j] * temp;
          }
        }
      for (j = 0; j < n; ++j)
        {
        l = ipvt[j];
        x(l) = wa1(j);
        }
      iter = 0;
      for (j = 0; j < n; ++j)
        wa2(j) = diag(j) * x(j);
      dxnorm = norm(wa2);
      fp = dxnorm - delta;
      if (fp <= (T)0.1*delta)
        {
        if (iter == 0)
          par = (T)0;
        return;
        }
      parl = (T)0;
      if (nsing >= n)
        {
        for (j = 0; j < n; ++j)
          {
          l = ipvt[j];
          wa1(j) = diag(l) * wa2(l) / dxnorm;
          }
        for (j = 0; j < n; ++j)
          {
          sum = (T)0;
          if (j >= 1)
            {
            for (i = 0; i < j; ++i)
              sum += r[i][j] * wa1(i);
            }
          wa1(j) = (wa1(j) - sum) / r[j][j];
          }
        temp = norm(wa1);
        parl = ((fp / delta) / temp) / temp;
        }
      for (j = 0; j < n; ++j)
        {
        sum = (T)0;
        for (i = 0; i <= j; ++i)
          sum += r[i][j] * qtb(i);
        l = ipvt[j];
        wa1(j) = sum / diag(l);
        }
      gnorm = norm(wa1);
      paru = gnorm / delta;
      if (paru == (T)0)
        paru = dwarf / std::min<T>(delta, (T)0.1);
      par = std::max<T>(par, parl);
      par = std::min<T>(par, paru);
      if (par == (T)0)
        par = gnorm / dxnorm;
      for (;;)
        {
        iter++;
        if (par == (T)0)
          par = std::max<T>(dwarf, (T)0.001 * paru);
        temp = std::sqrt(par);
        for (j = 0; j < n; ++j)
          wa1(j) = temp * diag(j);
        qrsolv(r, ipvt, wa1, qtb, x, sdiag);
        for (j = 0; j < n; ++j)
          wa2(j) = diag(j) * x(j);
        dxnorm = norm(wa2);
        temp = fp;
        fp = dxnorm - delta;

        if ((std::abs(fp) <= (T)0.1*delta) || ((parl == (T)0) && (fp <= temp) && (temp > (T)0)) || iter == 10)
          {
          if (iter == 0)
            par = (T)0;
          return;
          }
        for (j = 0; j < n; ++j)
          {
          l = ipvt[j];
          wa1(j) = diag(l) * wa2(l) / dxnorm;
          }
        for (j = 0; j < n; ++j)
          {
          wa1(j) /= sdiag(j);
          temp = wa1(j);
          jp1 = j + 1;
          if (jp1 < n)
            for (i = jp1; i < n; ++i)
              wa1(i) -= r[i][j] * temp;
          }
        temp = norm(wa1);
        parc = ((fp / delta) / temp) / temp;
        if (fp > (T)0)
          parl = std::max<T>(parl, par);
        if (fp < (T)0)
          paru = std::min<T>(paru, par);
        par = std::max<T>(parl, par + parc);
        }
      }
    }

  template <class T, class Container_n, class Container_m>
  void lmdif(bool(*f)(const matrix<T, Container_n>&, matrix<T, Container_m>&, void*), matrix<T, Container_n>& x, matrix<T, Container_m>& fvec, T ftol,
    T xtol, T gtol, uint64_t maxfev, T epsfcn, matrix<T, Container_n>& diag,
    uint64_t mode, T factor, uint64_t& info, uint64_t& nfev, void* user_data)
    {
    uint64_t n = x.rows();

    uint64_t i, iter, j, l, iflag;
    T actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm;
    T par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm;

    /* initialize */
    epsmch = std::numeric_limits<T>::epsilon();
    delta = (T)0;
    xnorm = (T)0;
    temp = (T)0;
    info = 0;
    nfev = 0;

    /* check for input parameter errors */
    if ((n <= 0) || (maxfev <= 0)
      || (factor <= 0)) return;
    if (mode == 2)
      {
      for (j = 0; j < n; ++j)
        if (diag[j] <= 0)
          return;
      }

    iflag = 1;
    nfev = 1;
    /* evaluate the function at the starting point and calculate its norm */
    if (!f(x, fvec, user_data))
      {
      info = 9; // break by user
      return;
      }

    uint64_t m = fvec.rows();

    if (m < n)
      {
      info = 0;
      return;
      }

    std::vector<uint64_t> ipvt;
    matrix<T, Container_n> qtf(n, 1), wa1(n, 1), wa2(n, 1), wa3(n, 1);
    matrix<T, Container_m> wa4(m, 1);
    matrix<T> fjac;

    fnorm = norm(fvec);

    /* initialize levenberg-marquardt counters */
    par = 0;
    iter = 1;

    /* outer loop */
    for (;;)
      {
      /* calculate jacobian matrix */
      iflag = 2;
      bool res = fdjac(f, x, fvec, fjac, epsfcn, user_data);
      nfev += n;
      if (!res)
        {
        info = 9; // break by user
        return;
        }
      res = f(x, fvec, user_data);
      /* compute the qr factorization of the jacobian */
      qrfac(fjac, true, ipvt, wa1, wa2);
      if (iter == 1)
        {
        if (mode != 2)
          {
          for (j = 0; j < n; ++j)
            {
            diag(j) = wa2(j);
            if (wa2(j) == (T)0)
              diag(j) = (T)1;
            }
          }
        for (j = 0; j < n; ++j)
          wa3(j) = diag(j) * x(j);
        xnorm = norm(wa3);
        delta = factor * xnorm;
        if (delta == 0)
          delta = factor;
        }
      for (i = 0; i < m; ++i)
        wa4(i) = fvec(i);
      for (j = 0; j < n; ++j)
        {
        if (fjac[j][j] != (T)0)
          {
          sum = (T)0;
          for (i = j; i < m; ++i)
            sum += fjac[i][j] * wa4(i);
          temp = -sum / fjac[j][j];
          for (i = j; i < m; ++i)
            wa4(i) += fjac[i][j] * temp;
          }
        fjac[j][j] = wa1(j);
        qtf(j) = wa4(j);
        }
      /* compute the norm of the scaled gradient */
      gnorm = (T)0;
      if (fnorm != (T)0)
        {
        for (j = 0; j < n; ++j)
          {
          l = ipvt[j];
          if (wa2(l) == (T)0)
            continue;
          sum = (T)0;
          for (i = l; i <= j; ++i)
            sum += fjac[i][j] * qtf(i) / fnorm;
          gnorm = std::max<T>(gnorm, std::abs(sum / wa2(l)));
          }
        }
      /* test for convergence of the gradient norm */
      if (gnorm <= gtol)
        info = 4;
      if (info != 0)
        {
        info = iflag;
        return;
        }
      /* rescale if necessary */
      if (mode != 2)
        {
        for (j = 0; j < n; ++j)
          diag(j) = std::max<T>(diag(j), wa2(j));
        }
      /* beginning of inner loop */
      do {
        /* determine the levenberg-marquardt parameter */
        ::jtk::lm_details::lmpar(fjac, ipvt, diag, qtf, delta, par, wa1, wa2);
        for (j = 0; j < n; ++j)
          {
          wa1(j) = -wa1(j);
          wa2(j) = x(j) + wa1(j);
          wa3(j) = diag(j) * wa1(j);
          }
        pnorm = norm(wa3);
        if (iter == 1)
          delta = std::min<T>(delta, pnorm);
        iflag = 1;
        res = f(wa2, wa4, user_data);
        ++(nfev);
        if (!res)
          {
          info = 9; // break by user
          return;
          }
        fnorm1 = norm(wa4);
        actred = (T)-1;
        if ((T)0.1 * fnorm1 < fnorm)
          actred = (T)1 - (fnorm1*fnorm1 / (fnorm*fnorm));
        for (j = 0; j < n; j++)
          {
          wa3(j) = (T)0;
          l = ipvt[j];
          temp = wa1(l);
          for (i = 0; i <= j; ++i)
            wa3(i) += fjac[i][j] * temp;
          }
        temp1 = norm(wa3) / fnorm;
        temp2 = std::sqrt(par) * pnorm / fnorm;
        prered = temp1 * temp1 + temp2 * temp2 / (T)0.5;
        dirder = -(temp1*temp1 + temp2 * temp2);
        ratio = (T)0;
        if (prered != (T)0)
          ratio = actred / prered;
        if (ratio <= (T)0.25)
          {
          if (actred > (T)0) temp = (T)0.5;
          if (actred < (T)0) temp = (T)0.5*dirder / (dirder + (T)0.5*actred);
          delta = temp * std::min<T>(delta, pnorm / (T)0.1);
          par /= temp;
          }
        else {
          if ((par == (T)0) || (ratio >= (T)0.75))
            {
            delta = pnorm / (T)0.5;
            par *= (T)0.5;
            }
          }
        if (ratio >= (T)0.0001)
          {
          for (j = 0; j < n; ++j)
            {
            x(j) = wa2(j);
            wa2(j) = diag(j) * x(j);
            }
          for (i = 0; i < m; ++i)
            fvec(i) = wa4(i);
          xnorm = norm(wa2);
          fnorm = fnorm1;
          iter++;
          }
        if ((std::abs(actred) <= ftol) && (prered <= ftol) &&
          ((T)0.5*ratio <= (T)1))
          info = 1;
        if (delta <= xtol * xnorm)
          info = 2;
        if ((std::abs(actred) <= ftol) && (prered <= ftol) &&
          ((T)0.5*ratio <= (T)1) && (info == 2))
          info = 3;
        if (nfev >= maxfev)
          info = 5;
        if ((std::abs(actred) <= epsmch) && (prered <= epsmch) &&
          ((T)0.5*ratio <= (T)1))
          info = 6;
        if (delta <= epsmch * xnorm)
          info = 7;
        if (gnorm <= epsmch)
          info = 8;
        if (info != 0)
          {
          info = iflag;
          return;
          }
        } while (ratio <= (T)0.0001);

      }
    }

  template <class T, class Container_n, class Container_m>
  void lmdif0(bool(*f)(const matrix<T, Container_n>&, matrix<T, Container_m>&, void*), matrix<T, Container_n>& x, matrix<T, Container_m>& fvec, uint64_t& info, uint64_t& nfev, void* user_data)
    {
    uint64_t j, maxfev;
    uint64_t mode;
    T ftol, xtol, gtol, epsfcn, factor;

    uint64_t n = x.rows();
    matrix<T, Container_n> diag(n, 1);

    /* Set convergence tolerances */
    ftol = std::numeric_limits<T>::epsilon() * 10000;
    xtol = std::numeric_limits<T>::epsilon() * 10000;
    gtol = (T)0;

    maxfev = n * 1000;
    epsfcn = std::numeric_limits<T>::epsilon();

    mode = 2;
    factor = (T)100;

    for (j = 0; j < n; ++j)
      diag(j) = (T)1.0001;

    lmdif(f, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor, info, nfev, user_data);
    }

  ///////////////////////////////////////////////////////////////////////////////
  // Eigenvalue computation
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container>
  void tred2(matrix<T, Container>& diagonal, matrix<T, Container>& subdiagonal, matrix<T, Container>& a)
    {
    size_t l, k, j, i;
    T scale, hh, h, g, f;

    size_t n = a.cols();
    diagonal.resize(n, 1);
    subdiagonal.resize(n, 1);
    for (i = n - 1; i > 0; --i)
      {
      l = i - 1;
      h = 0.0;
      scale = 0.0;
      if (l > 0)
        {
        for (k = 0; k <= l; ++k)
          scale += std::abs(a(i, k));
        if (scale == 0.0)
          subdiagonal(i) = a(i, l);
        else
          {
          for (k = 0; k <= l; ++k)
            {
            a(i, k) /= scale;
            h += a(i, k)*a(i, k);
            }
          f = a(i, l);
          g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
          subdiagonal(i) = scale * g;
          h -= f * g;
          a(i, l) = f - g;
          f = 0.0;
          for (j = 0; j <= l; ++j)
            {
            a(j, i) = a(i, j) / h;
            g = 0.0;
            for (k = 0; k <= j; ++k)
              g += a(j, k)*a(i, k);
            for (k = j + 1; k <= l; ++k)
              g += a(k, j)*a(i, k);
            subdiagonal(j) = g / h;
            f += subdiagonal(j)*a(i, j);
            } // for j
          hh = f / (h + h);
          for (j = 0; j <= l; ++j)
            {
            f = a(i, j);
            g = subdiagonal(j) - hh * f;
            subdiagonal(j) = g;
            for (k = 0; k <= j; ++k)
              {
              a(j, k) -= (f*subdiagonal(k) + g * a(i, k));
              } // for k
            } // for j
          }
        }
      else
        {
        subdiagonal(i) = a(i, l);
        }
      diagonal(i) = h;
      } // for i
    diagonal(0) = 0.0;
    subdiagonal(0) = 0.0;
    for (i = 0; i < n; ++i)
      {
      l = i - 1;
      if (diagonal(i))
        {
        if (i > 0)
          for (j = 0; j <= l; ++j)
            {
            g = 0.0;
            for (k = 0; k <= l; ++k)
              g += a(i, k)*a(k, j);
            for (k = 0; k <= l; ++k)
              a(k, j) -= g * a(k, i);
            } // for j
        }
      diagonal(i) = a(i, i);
      a(i, i) = 1.0;
      if (i > 0)
        for (j = 0; j <= l; ++j)
          {
          a(j, i) = 0.0;
          a(i, j) = 0.0;
          }
      } // for i
    }

  template <class T, class Container>
  bool tqli(matrix<T, Container>& diagonal, matrix<T, Container>& a, matrix<T, Container>& subdiagonal)
    {
    size_t n = diagonal.rows();
    size_t m, iterations, l, i, k;
    T s, r, p, g, f, dd, c, b;
    for (size_t count = 1; count < n; ++count)
      {
      subdiagonal(count - 1) = subdiagonal(count);
      }
    subdiagonal(n - 1) = 0.0;
    for (l = 1; l <= n; ++l)
      {
      iterations = 0;
      do
        {
        for (m = l; m <= n - 1; ++m)
          {
          dd = std::abs(diagonal(m - 1)) + std::abs(diagonal(m));
          if (std::abs(subdiagonal(m - 1)) < 1e-8)
            break;
          } // for m
        if (m != l)
          {
          if (iterations++ == 100)
            return false; // to many iterations					
          g = (diagonal(l) - diagonal(l - 1)) / (2.0*subdiagonal(l - 1));
          r = sqrt(g*g + 1.0);
          g = diagonal(m - 1) - diagonal(l - 1) + subdiagonal(l - 1) / (g + svd_details::sign(r, g));
          s = 1.0;
          c = 1.0;
          p = 0.0;
          for (i = m - 1; i >= l; --i)
            {
            f = s * subdiagonal(i - 1);
            b = c * subdiagonal(i - 1);
            r = sqrt(f*f + g * g);
            subdiagonal(i) = r;
            if (r == 0.0)
              {
              diagonal(i) -= p;
              subdiagonal(m - 1) = 0.0;
              break;
              }
            s = f / r;
            c = g / r;
            g = diagonal(i) - p;
            r = (diagonal(i - 1) - g)*s + 2.0*c*b;
            p = s * r;
            diagonal(i) = g + p;
            g = c * r - b;
            for (k = 1; k <= n; ++k)
              {
              f = a(k - 1, i);
              a(k - 1, i) = s * a(k - 1, i - 1) + c * f;
              a(k - 1, i - 1) = c * a(k - 1, i - 1) - s * f;
              }
            }
          if (r == 0.0 && i >= l)
            continue;
          diagonal(l - 1) -= p;
          subdiagonal(l - 1) = g;
          subdiagonal(m - 1) = 0.0;
          } // if (m != l)
        } while (m != l); // do-loop
      } // for l
    return true;
    }

  template <class T, class Container>
  bool eig_symm(matrix<T, Container>& eigenvalues, matrix<T, Container>& a)
    {
    matrix<T, Container> subdiag;
    tred2(eigenvalues, subdiag, a);
    return tqli(eigenvalues, a, subdiag);
    }

  template <class T, class Container>
  void balanc(matrix<T, Container>& a)
    {
    assert(a.rows() == a.cols());
    constexpr double RADIX = 2.0;
    uint64_t last, j, i;
    uint64_t n = a.rows();
    double s, r, g, f, c, sqrdx;
    sqrdx = RADIX * RADIX;
    last = 0;
    while (last == 0)
      {
      last = 1;
      for (i = 0; i < n; ++i)
        {
        //Calculate row and column norms.
        r = c = 0.0;
        for (j = 0; j < n; ++j)
          if (j != i)
            {
            c += fabs(a[j][i]);
            r += fabs(a[i][j]);
            }
        if (c && r)
          {
          //If both are nonzero,
          g = r / RADIX;
          f = 1.0;
          s = c + r;
          while (c < g)
            {
            //find the integer power of the machine radix that
            f *= RADIX; //comes closest to balancing the matrix.
            c *= sqrdx;
            }
          g = r * RADIX;
          while (c > g)
            {
            f /= RADIX;
            c /= sqrdx;
            }
          if ((c + r) / f < 0.95*s)
            {
            last = 0;
            g = 1.0 / f;
            for (j = 0; j < n; ++j)
              a[i][j] *= g; //Apply similarity transformation
            for (j = 0; j < n; ++j)
              a[j][i] *= f;
            }
          }
        }
      }
    }

  template <class T, class Container>
  void elmhes(matrix<T, Container>& a)
    {
    assert(a.rows() == a.cols());
    uint64_t n = a.rows();
    uint64_t m, j, i;
    T y, x;
    for (m = 1; m < n - 1; ++m)
      {
      x = (T)0;
      i = m;
      for (j = m; j < n; ++j)
        {
        //Find the pivot.
        if (std::abs(a[j][m - 1]) > std::abs(x))
          {
          x = a[j][m - 1];
          i = j;
          }
        }
      if (i != m)
        {
        //Interchange rows and columns.
        for (j = m - 1; j < n; ++j)
          std::swap(a[i][j], a[m][j]);
        for (j = 0; j < n; ++j)
          std::swap(a[j][i], a[j][m]);
        }
      if (x)
        {
        //Carry out the elimination.
        for (i = m + 1; i < n; ++i)
          {
          if ((y = a[i][m - 1]) != (T)0)
            {
            y /= x;
            a[i][m - 1] = y;
            for (j = m; j < n; ++j)
              a[i][j] -= y * a[m][j];
            for (j = 0; j < n; ++j)
              a[j][m] += y * a[j][i];
            }
          }
        }
      }
    }

  template <class T, class Container>
  bool hqr(matrix<T, Container>& a, matrix<T, Container>&  wr, matrix<T, Container>&  wi)
    {
    assert(a.rows() == a.cols());
    uint64_t n = a.rows();
    wr.resize(n, 1);
    wi.resize(n, 1);
    int nn, m, l, k, j, its, i, mmin;
    double z, y, x, w, v, u, t, s, r, q, p, anorm;
    p = q = r = (T)0;
    anorm = (T)0; //Compute matrix norm for possible use in locating single small subdiagonal element.
    for (i = 0; i < n; ++i)
      for (j = std::max<int>(i - 1, 0); j < n; ++j)
        anorm += std::abs(a[i][j]);
    nn = (int)n;
    t = 0.0; //Gets changed only by an exceptional shift.
    while (nn >= 1)
      {
      //Begin search for next eigenvalue.
      its = 0;
      do {
        for (l = nn - 1; l >= 1; --l)
          {
          //Begin iteration : look for single small subdiagonal element.
          s = fabs(a[l - 1][l - 1]) + fabs(a[l][l]);
          if (s == (T)0)
            s = anorm;
          if ((double)(std::abs(a[l][l - 1]) + s) == s)
            break;
          }
        x = a[nn - 1][nn - 1];
        if (l == nn - 1)
          {
          //One root found.
          wr(nn - 1) = x + t;
          wi(nn - 1) = (T)0;
          --nn;
          }
        else
          {
          y = a[nn - 2][nn - 2];
          w = a[nn - 1][nn - 2] * a[nn - 2][nn - 1];
          if (l == (nn - 2))
            {
            //Two roots found...
            p = (T)0.5*(y - x);
            q = p * p + w;
            z = std::sqrt(std::abs(q));
            x += t;
            if (q >= (T)0)
              {
              //...a real pair.
              z = p + svd_details::sign(z, p);
              wr(nn - 2) = wr(nn - 1) = x + z;
              if (z) wr(nn - 1) = x - w / z;
              wi(nn - 2) = wi(nn - 1) = (T)0;
              }
            else
              {
              //...a complex pair.
              wr(nn - 2) = wr(nn - 1) = x + p;
              wi(nn - 2) = -(wi(nn - 1) = z);
              }
            nn -= 2;
            }
          else
            {
            //No roots found.Continue iteration.
            if (its == 30)
              return false;
            if (its == 10 || its == 20)
              {
              //Form exceptional shift.
              t += x;
              for (i = 0; i < nn; ++i)
                a[i][i] -= x;
              s = std::abs(a[nn - 1][nn - 2]) + std::abs(a[nn - 2][nn - 3]);
              y = x = (T)0.75*s;
              w = (T)-0.4375*s*s;
              }
            ++its;
            for (m = (nn - 3); m >= l; --m)
              {
              //Form shift and then look for
              //  2 consecutive small subdiagonal elements.
              z = a[m][m];
              r = x - z;
              s = y - z;
              p = (r*s - w) / a[m + 1][m] + a[m][m + 1];
              q = a[m + 1][m + 1] - z - r - s;
              r = a[m + 2][m + 1];
              s = std::abs(p) + std::abs(q) + std::abs(r); //Scale to prevent overflow or underflow.
              p /= s;
              q /= s;
              r /= s;
              if (m == l)
                break;
              u = std::abs(a[m][m - 1])*(fabs(q) + fabs(r));
              v = std::abs(p)*(std::abs(a[m - 1][m - 1]) + std::abs(z) + std::abs(a[m + 1][m + 1]));
              if ((T)(u + v) == v) break;
              }
            for (i = m + 2; i < nn; ++i)
              {
              a[i][i - 2] = (T)0;
              if (i != (m + 2)) a[i][i - 3] = (T)0;
              }
            for (k = m; k < nn - 1; ++k)
              {
              // Double QR step on rows l to nn and columns m to nn.
              if (k != m)
                {
                p = a[k][k - 1]; //Begin setup of Householder vector
                q = a[k + 1][k - 1];
                r = (T)0;
                if (k != (nn - 2)) r = a[k + 2][k - 1];
                if ((x = std::abs(p) + std::abs(q) + std::abs(r)) != (T)0) {
                  p /= x; //Scale to prevent overflow or underflow.
                  q /= x;
                  r /= x;
                  }
                }
              if ((s = svd_details::sign(std::sqrt(p*p + q * q + r * r), p)) != (T)0)
                {
                if (k == m) {
                  if (l != m)
                    a[k][k - 1] = -a[k][k - 1];
                  }
                else
                  a[k][k - 1] = -s * x;
                p += s;
                x = p / s;
                y = q / s;
                z = r / s;
                q /= p;
                r /= p;
                for (j = k; j < nn; ++j)
                  {
                  //Row modification.
                  p = a[k][j] + q * a[k + 1][j];
                  if (k != (nn - 2))
                    {
                    p += r * a[k + 2][j];
                    a[k + 2][j] -= p * z;
                    }
                  a[k + 1][j] -= p * y;
                  a[k][j] -= p * x;
                  }
                mmin = nn < k + 4 ? nn : k + 4;
                for (i = l; i < mmin; ++i) {
                  //Column modification.
                  p = x * a[i][k] + y * a[i][k + 1];
                  if (k != (nn - 2))
                    {
                    p += z * a[i][k + 2];
                    a[i][k + 2] -= p * r;
                    }
                  a[i][k + 1] -= p * q;
                  a[i][k] -= p;
                  }
                }
              }
            }
          }
        } while (l < nn - 2);
      }
    return true;
    }

  template <class T, class Container>
  bool eig(matrix<T, Container>& a, matrix<T, Container>&  wr, matrix<T, Container>&  wi)
    {
    assert(a.rows() == a.cols());
    balanc(a);
    elmhes(a);
    return hqr(a, wr, wi);
    }

  ///////////////////////////////////////////////////////////////////////////////
  // iterative methods
  ///////////////////////////////////////////////////////////////////////////////

  template <class T, class Container>
  void sparse_matrix_vector_multiply(matrix<T, Container>& out, const sparse_matrix<T>& a, const matrix<T, Container>& b)
    {
    assert(b.cols() == 1);
    out.resize(a.rows(), 1);
    for (uint64_t i = 0; i < a.rows(); ++i)
      {
      T val = (T)0;
      auto it = a.row(i).begin();
      const auto it_end = it.end();
      for (; it != it_end; ++it)
        {
        val += *it * b(it.entry());
        }
      out(i) = val;
      }
    }

  template <class T, class Container>
  void sparse_symmetric_matrix_vector_multiply(matrix<T, Container>& out, const sparse_matrix<T>& a, const matrix<T, Container>& b)
    {
    assert(b.cols() == 1);
    out.resize(a.rows(), 1);
    for (uint64_t i = 0; i < a.rows(); ++i)
      {
      T val = (T)0;
      const T bi = b(i);
      auto it = a.row(i).begin();
      const auto it_end = it.end();

      for (; (it != it_end) && (it.entry() < i); ++it)
        {
        const auto Aij = *it;
        val += Aij * b(it.entry());
        out(it.entry()) += Aij * bi;
        }
      if ((it != it_end) && (it.entry() == i))
        {
        val += (*it)*bi;
        }
      out(i) = val;
      }
    }

  template <class T>
  class DiagonalPreconditioner
    {
    public:
      explicit DiagonalPreconditioner(const sparse_matrix<T>& A)
        {
        invD = diagonal(A);
        for (auto& v : invD)
          v = (v==0) ? 1 : (1 / v);
        }

      template <class T2, class Container>
      Expr<
        BinExprOp<
        typename matrix<T>::const_iterator,
        typename matrix<T2, Container>::const_iterator,
        OpMul<typename gettype<T, T2>::ty> > >
      solve(const matrix<T2, Container>& residual) const
        {
        typedef BinExprOp<
          typename matrix<T>::const_iterator,
          typename matrix<T2, Container>::const_iterator,
          OpMul<typename gettype<T, T2>::ty> > ExprT;
        return Expr<ExprT>(ExprT(invD.begin(), residual.begin(), invD.rows(), 1, false));
        }

    private:
      matrix<T> invD;
    };

  template <class T, class Container, class TPreconditioner>
  void preconditioned_conjugate_gradient(matrix<T, Container>& out,
    T& residu,
    uint64_t& iterations,
    const sparse_matrix<T>& A,
    const TPreconditioner& P,
    const matrix<T, Container>& b,
    const matrix<T, Container>& x0,
    T tolerance)
    {
    out = x0;
    matrix<T, Container> r = b - A * out;
    T rhsnorm = norm_sqr(b);
    if (rhsnorm == (T)0)
      {
      out = zeros(b.rows(), 1);
      iterations = 0;
      residu = (T)0;
      return;
      }
    const T threshold = tolerance * tolerance * rhsnorm;
    T residualnorm = norm_sqr(r);
    if (residualnorm < threshold)
      {
      iterations = 0;
      residu = residualnorm / rhsnorm;
      return;
      }

    matrix<T, Container> p = P.solve(r);

    matrix<T, Container> z, tmp;
    T absnew = dot(r, p);

    double time1 = 0;
    double time2 = 0;

    timer t1, t2;

    uint64_t i = 0;
    for (; i < A.rows(); ++i)
      {
      t1.start();
      //tmp.noalias() = A * p;
      //sparse_matrix_vector_multiply(tmp, A, p);
      sparse_symmetric_matrix_vector_multiply(tmp, A, p);
      double multime = t1.time_elapsed();
      //printf("A*x, time: %f\n", multime);
      time1 += multime;

      t2.start();
      T alpha = absnew / dot(p, tmp);
      out += alpha * p;
      r -= alpha * tmp;
      residualnorm = norm_sqr(r);
      if (residualnorm < threshold)
        break;
      z = P.solve(r);
      T absold = absnew;
      absnew = dot(r, z);
      T beta = absnew / absold;
      p = z + beta * p;
      time2 += t2.time_elapsed();
      }
    printf("time1: %f\n", time1);
    printf("time2: %f\n", time2);
    iterations = i + 1;
    residu = residualnorm / rhsnorm;
    }

  template <class T, class Container>
  void conjugate_gradient(matrix<T, Container>& out,
    T& residu,
    uint64_t& iterations,
    const sparse_matrix<T>& A,
    const matrix<T, Container>& b,
    const matrix<T, Container>& x0,
    T tolerance)
    {
    DiagonalPreconditioner P(A);
    preconditioned_conjugate_gradient(out, residu, iterations, A, P, b, x0, tolerance);
    }

  template <class T, class Container, class TPreconditioner>
  void bipcgstab(matrix<T, Container>& out,
    T& residu,
    uint64_t& iterations,
    const sparse_matrix<T>& A,
    const TPreconditioner& P,
    const matrix<T, Container>& b,
    const matrix<T, Container>& x0,
    T tolerance)
    {
    matrix<T, Container> r, r_, p, s, phat, shat, t, v;
    matrix<T, Container> rho_1(1), rho_2(1), alpha(1), omega(1);
    T beta = 0;
    out = x0;
    T normb = norm(b);
    r = b - A * out;
    r_ = r;

    if (normb == (T)0)
      normb = (T)1;

    T resid = norm(r) / normb;
    if (resid <= tolerance)
      {
      residu = resid;
      iterations = 0;
      return;
      }

    for (uint64_t i = 0; i < A.rows(); ++i)
      {
      rho_1 = transpose(r_) * r;
      if (rho_1(0) == (T)0)
        {
        residu = norm(r) / normb;
        iterations = i + 1;
        return;
        }
      if (i == 0)
        p = r;
      else
        {
        beta = (rho_1(0) / rho_2(0))*(alpha(0) / omega(0));
        p = r + beta * (p - omega(0)*v);
        }
      phat = P.solve(p);
      v = A * phat;
      alpha = rho_1(0) / (transpose(r_)*v);
      s = r - alpha(0)*v;
      resid = norm(s) / normb;
      if (resid < tolerance)
        {
        out += alpha(0)*phat;
        residu = resid;
        iterations = i + 1;
        return;
        }
      shat = P.solve(s);
      t = A * shat;
      omega = (transpose(t)*s) / matrix<T, Container>(transpose(t)*t)(0);
      out += alpha(0)*phat + omega(0)*shat;
      r = s - omega(0)*t;
      rho_2 = rho_1;
      resid = norm(r) / normb;
      if (resid < tolerance)
        {
        residu = resid;
        iterations = i + 1;
        return;
        }
      if (omega(0) == 0)
        {
        residu = norm(r) / normb;
        iterations = i + 1;
        return;
        }
      }
    residu = resid;
    iterations = A.rows();
    }

  template <class T, class Container>
  void bicgstab(matrix<T, Container>& out,
    T& residu,
    uint64_t& iterations,
    const sparse_matrix<T>& A,
    const matrix<T, Container>& b,
    const matrix<T, Container>& x0,
    T tolerance)
    {
    DiagonalPreconditioner P(A);
    bipcgstab(out, residu, iterations, A, P, b, x0, tolerance);
    }

  ///////////////////////////////////////////////////////////////////////////////
  // typedefs
  ///////////////////////////////////////////////////////////////////////////////

  typedef matrix<double> mat;
  typedef matrix<float> matf;
  typedef sparse_matrix<double> smat;
  typedef sparse_matrix<float> smatf;

  typedef matrix<double, std::array<double, 1>> mat1;
  typedef matrix<double, std::array<double, 2>> mat2;
  typedef matrix<double, std::array<double, 3>> mat3;
  typedef matrix<double, std::array<double, 4>> mat4;
  typedef matrix<double, std::array<double, 5>> mat5;
  typedef matrix<double, std::array<double, 6>> mat6;
  typedef matrix<double, std::array<double, 7>> mat7;
  typedef matrix<double, std::array<double, 8>> mat8;
  typedef matrix<double, std::array<double, 9>> mat9;
  typedef matrix<double, std::array<double, 10>> mat10;
  typedef matrix<double, std::array<double, 11>> mat11;
  typedef matrix<double, std::array<double, 12>> mat12;
  typedef matrix<double, std::array<double, 13>> mat13;
  typedef matrix<double, std::array<double, 14>> mat14;
  typedef matrix<double, std::array<double, 15>> mat15;
  typedef matrix<double, std::array<double, 16>> mat16;
  typedef matrix<double, std::array<double, 17>> mat17;
  typedef matrix<double, std::array<double, 18>> mat18;
  typedef matrix<double, std::array<double, 19>> mat19;
  typedef matrix<double, std::array<double, 20>> mat20;
  typedef matrix<double, std::array<double, 21>> mat21;
  typedef matrix<double, std::array<double, 22>> mat22;
  typedef matrix<double, std::array<double, 23>> mat23;
  typedef matrix<double, std::array<double, 24>> mat24;
  typedef matrix<double, std::array<double, 25>> mat25;
  typedef matrix<double, std::array<double, 26>> mat26;
  typedef matrix<double, std::array<double, 27>> mat27;
  typedef matrix<double, std::array<double, 28>> mat28;
  typedef matrix<double, std::array<double, 29>> mat29;
  typedef matrix<double, std::array<double, 30>> mat30;
  typedef matrix<double, std::array<double, 31>> mat31;
  typedef matrix<double, std::array<double, 32>> mat32;
  typedef matrix<double, std::array<double, 33>> mat33;
  typedef matrix<double, std::array<double, 34>> mat34;
  typedef matrix<double, std::array<double, 35>> mat35;
  typedef matrix<double, std::array<double, 36>> mat36;
  typedef matrix<double, std::array<double, 37>> mat37;
  typedef matrix<double, std::array<double, 38>> mat38;
  typedef matrix<double, std::array<double, 39>> mat39;
  typedef matrix<double, std::array<double, 40>> mat40;
  typedef matrix<double, std::array<double, 41>> mat41;
  typedef matrix<double, std::array<double, 42>> mat42;
  typedef matrix<double, std::array<double, 43>> mat43;
  typedef matrix<double, std::array<double, 44>> mat44;
  typedef matrix<double, std::array<double, 45>> mat45;
  typedef matrix<double, std::array<double, 46>> mat46;
  typedef matrix<double, std::array<double, 47>> mat47;
  typedef matrix<double, std::array<double, 48>> mat48;
  typedef matrix<double, std::array<double, 49>> mat49;
  typedef matrix<double, std::array<double, 50>> mat50;
  typedef matrix<double, std::array<double, 51>> mat51;
  typedef matrix<double, std::array<double, 52>> mat52;
  typedef matrix<double, std::array<double, 53>> mat53;
  typedef matrix<double, std::array<double, 54>> mat54;
  typedef matrix<double, std::array<double, 55>> mat55;
  typedef matrix<double, std::array<double, 56>> mat56;
  typedef matrix<double, std::array<double, 57>> mat57;
  typedef matrix<double, std::array<double, 58>> mat58;
  typedef matrix<double, std::array<double, 59>> mat59;
  typedef matrix<double, std::array<double, 60>> mat60;
  typedef matrix<double, std::array<double, 61>> mat61;
  typedef matrix<double, std::array<double, 62>> mat62;
  typedef matrix<double, std::array<double, 63>> mat63;
  typedef matrix<double, std::array<double, 64>> mat64;
  typedef matrix<double, std::array<double, 65>> mat65;
  typedef matrix<double, std::array<double, 66>> mat66;
  typedef matrix<double, std::array<double, 67>> mat67;
  typedef matrix<double, std::array<double, 68>> mat68;
  typedef matrix<double, std::array<double, 69>> mat69;
  typedef matrix<double, std::array<double, 70>> mat70;
  typedef matrix<double, std::array<double, 71>> mat71;
  typedef matrix<double, std::array<double, 72>> mat72;
  typedef matrix<double, std::array<double, 73>> mat73;
  typedef matrix<double, std::array<double, 74>> mat74;
  typedef matrix<double, std::array<double, 75>> mat75;
  typedef matrix<double, std::array<double, 76>> mat76;
  typedef matrix<double, std::array<double, 77>> mat77;
  typedef matrix<double, std::array<double, 78>> mat78;
  typedef matrix<double, std::array<double, 79>> mat79;
  typedef matrix<double, std::array<double, 81>> mat81;
  typedef matrix<double, std::array<double, 82>> mat82;
  typedef matrix<double, std::array<double, 83>> mat83;
  typedef matrix<double, std::array<double, 84>> mat84;
  typedef matrix<double, std::array<double, 85>> mat85;
  typedef matrix<double, std::array<double, 86>> mat86;
  typedef matrix<double, std::array<double, 87>> mat87;
  typedef matrix<double, std::array<double, 88>> mat88;
  typedef matrix<double, std::array<double, 89>> mat89;
  typedef matrix<double, std::array<double, 91>> mat91;
  typedef matrix<double, std::array<double, 92>> mat92;
  typedef matrix<double, std::array<double, 93>> mat93;
  typedef matrix<double, std::array<double, 94>> mat94;
  typedef matrix<double, std::array<double, 95>> mat95;
  typedef matrix<double, std::array<double, 96>> mat96;
  typedef matrix<double, std::array<double, 97>> mat97;
  typedef matrix<double, std::array<double, 98>> mat98;
  typedef matrix<double, std::array<double, 99>> mat99;
  typedef matrix<double, std::array<double, 100>> mat100;

  typedef matrix<float, std::array<float, 1>> matf1;
  typedef matrix<float, std::array<float, 2>> matf2;
  typedef matrix<float, std::array<float, 3>> matf3;
  typedef matrix<float, std::array<float, 4>> matf4;
  typedef matrix<float, std::array<float, 5>> matf5;
  typedef matrix<float, std::array<float, 6>> matf6;
  typedef matrix<float, std::array<float, 7>> matf7;
  typedef matrix<float, std::array<float, 8>> matf8;
  typedef matrix<float, std::array<float, 9>> matf9;
  typedef matrix<float, std::array<float, 10>> matf10;
  typedef matrix<float, std::array<float, 11>> matf11;
  typedef matrix<float, std::array<float, 12>> matf12;
  typedef matrix<float, std::array<float, 13>> matf13;
  typedef matrix<float, std::array<float, 14>> matf14;
  typedef matrix<float, std::array<float, 15>> matf15;
  typedef matrix<float, std::array<float, 16>> matf16;
  typedef matrix<float, std::array<float, 17>> matf17;
  typedef matrix<float, std::array<float, 18>> matf18;
  typedef matrix<float, std::array<float, 19>> matf19;
  typedef matrix<float, std::array<float, 20>> matf20;
  typedef matrix<float, std::array<float, 21>> matf21;
  typedef matrix<float, std::array<float, 22>> matf22;
  typedef matrix<float, std::array<float, 23>> matf23;
  typedef matrix<float, std::array<float, 24>> matf24;
  typedef matrix<float, std::array<float, 25>> matf25;
  typedef matrix<float, std::array<float, 26>> matf26;
  typedef matrix<float, std::array<float, 27>> matf27;
  typedef matrix<float, std::array<float, 28>> matf28;
  typedef matrix<float, std::array<float, 29>> matf29;
  typedef matrix<float, std::array<float, 30>> matf30;
  typedef matrix<float, std::array<float, 31>> matf31;
  typedef matrix<float, std::array<float, 32>> matf32;
  typedef matrix<float, std::array<float, 33>> matf33;
  typedef matrix<float, std::array<float, 34>> matf34;
  typedef matrix<float, std::array<float, 35>> matf35;
  typedef matrix<float, std::array<float, 36>> matf36;
  typedef matrix<float, std::array<float, 37>> matf37;
  typedef matrix<float, std::array<float, 38>> matf38;
  typedef matrix<float, std::array<float, 39>> matf39;
  typedef matrix<float, std::array<float, 40>> matf40;
  typedef matrix<float, std::array<float, 41>> matf41;
  typedef matrix<float, std::array<float, 42>> matf42;
  typedef matrix<float, std::array<float, 43>> matf43;
  typedef matrix<float, std::array<float, 44>> matf44;
  typedef matrix<float, std::array<float, 45>> matf45;
  typedef matrix<float, std::array<float, 46>> matf46;
  typedef matrix<float, std::array<float, 47>> matf47;
  typedef matrix<float, std::array<float, 48>> matf48;
  typedef matrix<float, std::array<float, 49>> matf49;
  typedef matrix<float, std::array<float, 50>> matf50;
  typedef matrix<float, std::array<float, 51>> matf51;
  typedef matrix<float, std::array<float, 52>> matf52;
  typedef matrix<float, std::array<float, 53>> matf53;
  typedef matrix<float, std::array<float, 54>> matf54;
  typedef matrix<float, std::array<float, 55>> matf55;
  typedef matrix<float, std::array<float, 56>> matf56;
  typedef matrix<float, std::array<float, 57>> matf57;
  typedef matrix<float, std::array<float, 58>> matf58;
  typedef matrix<float, std::array<float, 59>> matf59;
  typedef matrix<float, std::array<float, 60>> matf60;
  typedef matrix<float, std::array<float, 61>> matf61;
  typedef matrix<float, std::array<float, 62>> matf62;
  typedef matrix<float, std::array<float, 63>> matf63;
  typedef matrix<float, std::array<float, 64>> matf64;
  typedef matrix<float, std::array<float, 65>> matf65;
  typedef matrix<float, std::array<float, 66>> matf66;
  typedef matrix<float, std::array<float, 67>> matf67;
  typedef matrix<float, std::array<float, 68>> matf68;
  typedef matrix<float, std::array<float, 69>> matf69;
  typedef matrix<float, std::array<float, 70>> matf70;
  typedef matrix<float, std::array<float, 71>> matf71;
  typedef matrix<float, std::array<float, 72>> matf72;
  typedef matrix<float, std::array<float, 73>> matf73;
  typedef matrix<float, std::array<float, 74>> matf74;
  typedef matrix<float, std::array<float, 75>> matf75;
  typedef matrix<float, std::array<float, 76>> matf76;
  typedef matrix<float, std::array<float, 77>> matf77;
  typedef matrix<float, std::array<float, 78>> matf78;
  typedef matrix<float, std::array<float, 79>> matf79;
  typedef matrix<float, std::array<float, 81>> matf81;
  typedef matrix<float, std::array<float, 82>> matf82;
  typedef matrix<float, std::array<float, 83>> matf83;
  typedef matrix<float, std::array<float, 84>> matf84;
  typedef matrix<float, std::array<float, 85>> matf85;
  typedef matrix<float, std::array<float, 86>> matf86;
  typedef matrix<float, std::array<float, 87>> matf87;
  typedef matrix<float, std::array<float, 88>> matf88;
  typedef matrix<float, std::array<float, 89>> matf89;
  typedef matrix<float, std::array<float, 91>> matf91;
  typedef matrix<float, std::array<float, 92>> matf92;
  typedef matrix<float, std::array<float, 93>> matf93;
  typedef matrix<float, std::array<float, 94>> matf94;
  typedef matrix<float, std::array<float, 95>> matf95;
  typedef matrix<float, std::array<float, 96>> matf96;
  typedef matrix<float, std::array<float, 97>> matf97;
  typedef matrix<float, std::array<float, 98>> matf98;
  typedef matrix<float, std::array<float, 99>> matf99;
  typedef matrix<float, std::array<float, 100>> matf100;

  }