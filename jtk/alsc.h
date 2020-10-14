#pragma once

#include <stdint.h>
#include <cmath>
#include <stdio.h>

namespace jtk {

/*
Adaptive Least Squrares Correlation: A powerful image matching technique. Armin Gruen. S. Afr. J. of Photogrammetry,
Remote Sensing and Cartography 14 (3), 1985: 175-187.
*/

void print_patch(double* patch, int m, int n);
/*
resample_patch

double* patch: allocated memory of size mxn of type double. Will be filled with resampled data.
int m, int n: size of the patch
double a,b,c,d,e,f : coefficients for affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y
const uint8_t* rhimage: the image where we sample from
int stride: the stride for the image (typically the width of the image)

Note that there is no check for the affine transformation to go outside the borders of the rhimage.
*/
void resample_patch(double* patch, int m, int n,
  double a, double b, double c, double d, double e, double f,
  const uint8_t* rhimage, int stride);

/*
interpolate

Performs bilinear interpolation.
const uint8_t* rhimage: the input image.
int stride: the stride for the image (typically the width of the image)
double x,y: the location where to apply the bilinear interpolation.

returns the interpolated value as double
*/
double interpolate(const uint8_t* rhimage, int stride, double x, double y);

/*
fill_design_matrix

Fills the design matrix A. Actually we fill the transposed design matrix.
double* At: The transposed design matrix that will be filled. This pointer should point to allocated memory of mxnx8 double values.
const double* patch: The input patch that we use to fill At. Typically this patch is generated with the method resample_patch.
double rg: current radiometric gain.
int m, n: the patch is of size mxn. The design matrix A has mxn rows (or At has mxn columns).
*/
void fill_design_matrix(double* A,
  const double* patch, int m, int n);

/*
fill_patch_difference

Fills vector l as the difference between patch_f and patch_g.
double* l: Memory of size m.n of type double.
const double* patch_f: patch of size mxn
const double* patch_g: patch of size mxn
int m,n: dimensions of the patches.
*/
void fill_patch_difference(double* l, const double* patch_f, const double* patch_g, double rs, double rg, int m, int n);

double compute_residual_l2_norm_sqr(const double* patch_f, const double* patch_g, double rs, double rg, int m, int n);

/*
compute_At_A

Computes the matrix product At * A.
double* AtA: memory allocated of size 8x8 of type double. Will contain At * A.
const double* At: transposed design matrix, filled with method fill_design_matrix
int m,n: dimensions of the patches. At is a matrix of size 8 x (mxn).
*/
void compute_At_A(double* AtA, const double* A, int m, int n);

/*
compute_At_l

Computes the matrix vector product At * l.
double* Atl: memory allocated of size 8 of type double. Will contain the matrix vector product At*l.
const double* At: transposed design matrix, filled with method fill_design_matrix
const double* l: patch difference vector, filled with method fill_patch_difference
int m,n: dimensions of the patches. At is a matrix of size 8 x (mxn). l is a vector of size mxn.
*/
void compute_At_l(double* Atl, const double* A, const double* l, int m, int n);

/*
cholesky

Computes the Cholesky factorization of matrix AtA, i.e. AtA = U**T * U.
double* AtA: Points to memory of size 8x8 of type double. The memory was filled with the method compute_At_A.
AtA will be replaced by an upper triangular matrix U, such that AtA = U**T * U.
*/
int cholesky(double* AtA);

/*
solve_alsc_cholesky

Solves the system AtA x = Atl.
double* U: the upper triangular matrix obtained from the Cholesky factorization of matrix AtA, where AtA was obtained with the
method compute_At_A, and subsequently fed to the method cholesky to obtain the upper triangular matrix of the Cholesky
factorization AtA = U**T * U. U is a matrix of size 8x8 of type double.
double* Atl: vector of size 8 of type double obtained with the method compute_At_l. Upon exit, Atl contains the solution x
of the system AtA x = Atl.
*/
int solve_alsc_cholesky(double* U, double* Atl);

/*
update_coefficients

Updates the coefficients a,b,c,d,e,f for the affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y by the
solution of the system AtA x = Atl.
double &a,&b,&c,&d,&e,&f: coefficients for affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y that will be updated.
double &rg: radiometric gain
const double* x: Solution of the system AtA x = Atl, computed with method solve_alsc_cholesky.
*/
bool update_coefficients(double& a, double& b, double& c, double& d, double& e, double& f, double& rs, double& rg, const double* x);


/*
terminate

Termination test when to stop the ALSC iteration. The iterations finish when the changes in the affine transformation coefficients
a,b,c,d,e,f fall below a threshold.
const double* x: Solution of the system AtA x = Atl, computed with method solve_alsc_cholesky.
double tolerance[6]: threshold value for the six affine transformation coefficients.
*/
bool terminate(const double* x, double rg, double tolerance[6]);

/*
inverse_from_cholesky

double* U: the upper triangular matrix obtained from the Cholesky factorization of matrix AtA, where AtA was obtained with the
method compute_At_A, and subsequently fed to the method cholesky to obtain the upper triangular matrix of the Cholesky
factorization AtA = U**T * U. U is a matrix of size 8x8 of type double.
Upon exit, U is replaced by (AtA)**(-1).
*/
int inverse_from_cholesky(double* U);

/*
get_coefficients_unskewed_patch

Computes the coefficients a,b,c,d,e,f such that the method resample_patch generates an unskewed patch around point (x,y) of size mxn.
double &a,&b,&c,&d,&e,&f: coefficients for affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y
double &rg: radiometric gain.
int m,n: dimensions of the patch.
double x,y: center location of the patch.
*/
void get_coefficients_unskewed_patch(double& a, double& b, double& c, double& d, double& e, double& f, double& rs, double& rg, int m, int n, double x, double y);


void get_shift_coefficients(double& a, double& d, double b, double c, double e, double f, int m, int n, double x, double y);

/*
get_rectangle_around_patch
Given the coefficients a,b,c,d,e,f, the corners (x1,y1),(x2,y2),(x3,y3),(x4,y4) of the patch obtained via method resample_patch are
computed. Useful for drawing debug output.
double &x1,&y1,&x2,&y2,&x3,&y3,&x4,&y4: output coordinates of the four corners of the patch.
double a,b,c,d,e,f: coefficients for affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y
int m,n: dimensions of the patch.
*/
void get_rectangle_around_patch(double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, double& x4, double& y4, double a, double b, double c, double d, double e, double f, int m, int n);

/*
get_residual

Computes the residual v = A*x - l, where A is the design matrix, l is the patch difference, and x is the solution of the
system AtA x = Atl.
const double* At: transposed design matrix, filled with method fill_design_matrix.
const double* x: Solution of the system AtA x = Atl, obtained from method solve_alsc_cholesky.
double* l: patch difference vector, filled with method fill_patch_difference. Upon exit, it contains the residual v = A*x - l.
int m,n: dimensions of the patches.
*/
void get_residual(const double* A, const double* x, double* l, int m, int n);

/*
get_variance_factor

Returns the variance factor sigma**2 = 1/r * v**t * v, where v is the residual, and r = mxn - 8.
const double* residual: The residual, obtained with mtehod get_residual.
int m,n: dimensions of the patches.
Returns a double value representing the variance factor 1/r * v**t * v.
*/
double get_variance_factor(const double* residual, int m, int n);

/*
compute_standard_deviations

Computes the standard deviations of the coefficients for the affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y and
for the radiometric shift.
double &sigma_a, &sigma_b, &sigma_c, &sigma_d, &sigma_e, &sigma_f, &sigma_radiometric_shift : output values that are computed.
double variance_factor: Computed with the method get_variance_factor.
const double* AtA_inv: Inverse of matrix AtA. Compute with method inverse_from_cholesky.
*/
void compute_standard_deviations(
  double& sigma_a, double& sigma_b, double& sigma_c,
  double& sigma_d, double& sigma_e, double& sigma_f, double& sigma_radiometric_shift, double& sigma_radiometric_gain,
  double variance_factor, const double* AtA_inv);

/*
compute_standard_deviations

Computes the standard deviations of the more interesting shift coefficients a and d for the
affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y.
double &sigma_a, &sigma_d: output values that are computed.
double variance_factor: Computed with the method get_variance_factor.
const double* AtA_inv: Inverse of matrix AtA. Compute with method inverse_from_cholesky.
*/
void compute_standard_deviations(double& sigma_a, double& sigma_d, double variance_factor, const double* AtA_inv);

/*
*/
double quality_of_match(double variance_factor, double* AtA_inv);
/*
compute_patch_center

Computes the center of the patch.
double &x, &y: output coordinates representing the center of the patch.
int m, int n: size of the patch
double a,b,c,d,e,f : coefficients for affine transformation x_ = a + b*x + c*y, y_ = d + e*x + f*y
*/
void compute_patch_center(double& x, double& y, double a, double b, double c, double d, double e, double f, int m, int n);



bool alsc_single_iteration(double& a, double& b, double& c, double& d, double& e, double& f, double& r_shift, double& r_gain,
  const double* patch_f, double* patch_g,
  double* At,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  const uint8_t* rhimage, int stride);

void get_extrema_patch(double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, double& x4, double& y4, double a, double b, double c, double d, double e, double f, int m, int n);

struct alsc_result
  {
  double right_x, right_y;
  double a, b, c, d, e, f;
  double r_shift;
  double r_gain;
  double sigma_a;
  double sigma_b;
  double sigma_c;
  double sigma_d;
  double sigma_e;
  double sigma_f;
  double sigma_rs;
  double sigma_rg;
  double variance_factor;
  int iterations;
  bool success;
  double quality;
  double residu;
  };

enum alsc_output_variable
  {
  RIGHT_X,
  RIGHT_Y,
  A,
  B,
  C,
  D,
  E,
  F,
  RS,
  RG,
  SIGMA_X,
  SIGMA_Y,
  VARIANCE_FACTOR,
  QUALITY,
  RESIDU,
  ITERATIONS
  };

bool alsc(
  double* output,
  double left_x, double left_y,
  double right_x, double right_y,
  const uint8_t* lhimage,
  const uint8_t* rhimage,
  double sin_angle,
  double cos_angle,
  int w,
  int h,
  double* patch_f,
  double* patch_g,
  double* A,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  double tol[6],
  int max_iter);

alsc_result alsc(double right_x, double right_y,
  const uint8_t* rhimage,
  int w,
  int h,
  const double* patch_f,
  double* patch_g,
  double* A,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  double tol[6],
  int max_iter);

alsc_result find_match_on_line(
  double x, double y,
  const double* line,
  const uint8_t* rhimage,
  int w,
  int h,
  const double* patch_f, double* patch_g,
  double* A,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  int step_size,
  int max_disparity,
  double tol[6],
  int max_iter
);

////////////////////////////////////
// Implementation
////////////////////////////////////


#include <cmath>
#include <stdio.h>

inline void print_patch(double* patch, int m, int n)
  {
  for (int y = -1; y <= m; ++y)
    {
    for (int x = -1; x <= n; ++x)
      {
      printf("%.2f ", *patch++);
      }
    printf("\n");
    }
  }

inline void print_mat(double* mat, int m, int n)
  {
  for (int r = 0; r < m; ++r)
    {
    for (int c = 0; c < n; ++c)
      {
      printf("%.2f ", mat[c*m + r]);
      }
    printf("\n");
    }
  }

inline void resample_patch(double* patch, int m, int n,
  double a, double b, double c, double d, double e, double f,
  const uint8_t* rhimage, int stride)
  {
  double x_from, y_from;
  for (int y = -1; y <= m; ++y)
    {
    for (int x = -1; x <= n; ++x)
      {
      x_from = a + b * x + c * y;
      y_from = d + e * x + f * y;
      //printf("a:%2.f  b:%2.f c:%2.f u:%2.f  v:%2.f\n", a, b, c, x_from, y_from);
      *patch++ = interpolate(rhimage, stride, x_from, y_from);
      }
    }
  }

inline void compute_patch_center(double& x, double& y, double a, double b, double c, double d, double e, double f, int m, int n)
  {
  double x_ = (n - 1.0) / 2.0;
  double y_ = (m - 1.0) / 2.0;
  x = a + b * x_ + c * y_;
  y = d + e * x_ + f * y_;
  }

inline double interpolate(const uint8_t* rhimage, int stride, double x, double y)
  {
  int x_ = int(std::floor(x));
  int y_ = int(std::floor(y));
  double fract_x = x - double(x_);
  double fract_y = y - double(y_);
  const uint8_t* p_im = rhimage + y_ * stride + x_;
  double topleft = *p_im;
  ++p_im;
  double topright = *p_im;
  p_im += stride - 1;
  double bottomleft = *p_im;
  ++p_im;
  double bottomright = *p_im;
  double temp1 = topleft + fract_x * (topright - topleft);
  double temp2 = bottomleft + fract_x * (bottomright - bottomleft);
  double result = temp1 + fract_y * (temp2 - temp1);
  return result;
  }

inline void fill_design_matrix(double* A,
  const double* patch, int m, int n)
  {
  int mn = m * n;
  double gx, gy, g;
  double* A0 = A;
  double* A1 = A + mn;
  double* A2 = A + 2 * mn;
  double* A3 = A + 3 * mn;
  double* A4 = A + 4 * mn;
  double* A5 = A + 5 * mn;
  double* A6 = A + 6 * mn;
  double* A7 = A + 7 * mn;
  for (int y = 1; y <= m; ++y)
    {
    for (int x = 1; x <= n; ++x)
      {
      g = double(*(patch + y * (n + 2) + x));
      gx = (double(*(patch + y * (n + 2) + x + 1)) - double(*(patch + y * (n + 2) + x - 1)))*0.5;
      gy = (double(*(patch + (y + 1)*(n + 2) + x)) - double(*(patch + (y - 1)*(n + 2) + x)))*0.5;
      *A0++ = gx;
      *A1++ = (x - 1)*gx;
      *A2++ = (y - 1)*gx;
      *A3++ = gy;
      *A4++ = (x - 1)*gy;
      *A5++ = (y - 1)*gy;
      *A6++ = 1.0;
      *A7++ = g;
      }
    }
  }

inline void fill_patch_difference(double* l, const double* patch_f, const double* patch_g, double rs, double rg, int m, int n)
  {
  double f, g;
  for (int y = 1; y <= m; ++y)
    {
    for (int x = 1; x <= n; ++x)
      {
      f = *(patch_f + y * (n + 2) + x);
      g = *(patch_g + y * (n + 2) + x);
      *l++ = f - ((1.0 + rg)*g + rs);
      }
    }
  }

namespace alsc_details
  {

  inline double sqr(double a)
    {
    return a * a;
    }

  inline void dswap(double& a, double& b)
    {
    double c = a;
    a = b;
    b = c;
    }
  }

inline double compute_residual_l2_norm_sqr(const double* patch_f, const double* patch_g, double rs, double rg, int m, int n)
  {
  using namespace alsc_details;
  double f, g;
  double residu = 0.0;
  for (int y = 1; y <= m; ++y)
    {
    for (int x = 1; x <= n; ++x)
      {
      f = *(patch_f + y * (n + 2) + x);
      g = *(patch_g + y * (n + 2) + x);
      residu += alsc_details::sqr(f - ((1.0 + rg)*g + rs));
      }
    }
  return residu;
  }

inline void compute_At_A(double* AtA, const double* A, int m, int n)
  {
  int mn = m * n;

  for (int c = 0; c < 8; ++c)
    {
    for (int r = 0; r <= c; ++r)
      {
      const double* A_it_1 = A + mn * c;
      const double* A_it_2 = A + mn * r;
      double val = 0.0;
      for (int k = 0; k < mn; ++k)
        {
        val += *A_it_1++ * *A_it_2++;
        }
      AtA[c * 8 + r] = val;
      }
    }
  }

inline void compute_At_l(double* Atl, const double* A, const double* l, int m, int n)
  {
  int mn = m * n;
  for (int i = 0; i < 8; ++i)
    {
    double val = 0.0;
    const double* lit = l;
    for (int k = 0; k < mn; ++k)
      {
      val += (*A++) * (*lit++);
      }
    *Atl++ = val;
    }
  }

inline void get_coefficients_unskewed_patch(double& a, double& b, double& c, double& d, double& e, double& f, double& rs, double& rg, int m, int n, double x, double y)
  {
  b = 1.0;
  c = 0.0;
  e = 0.0;
  f = 1.0;
  a = x - (n - 1.0) / 2.0;
  d = y - (m - 1.0) / 2.0;
  rg = 0.0;
  rs = 0.0;
  }

inline void get_shift_coefficients(double& a, double& d, double b, double c, double e, double f, int m, int n, double x, double y)
  {
  a = x - b * (n - 1.0) / 2.0 - c * (m - 1.0) / 2.0;
  d = y - e * (n - 1.0) / 2.0 - f * (m - 1.0) / 2.0;
  }

inline void get_rectangle_around_patch(double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, double& x4, double& y4, double a, double b, double c, double d, double e, double f, int m, int n)
  {
  x1 = a;
  y1 = d;
  x2 = b * (n - 1) + a;
  y2 = e * (n - 1) + d;
  x3 = b * (n - 1) + c * (m - 1) + a;
  y3 = e * (n - 1) + f * (m - 1) + d;
  x4 = c * (m - 1) + a;
  y4 = f * (m - 1) + d;
  }

inline void get_residual(const double* A, const double* x, double* l, int m, int n)
  {
  int mn = m * n;
  const double* A0 = A;
  const double* A1 = A + mn;
  const double* A2 = A + 2 * mn;
  const double* A3 = A + 3 * mn;
  const double* A4 = A + 4 * mn;
  const double* A5 = A + 5 * mn;
  const double* A6 = A + 6 * mn;
  const double* A7 = A + 7 * mn;
  for (int i = 0; i < mn; ++i)
    {
    *l++ -= (*A0++) * x[0] + (*A1++) * x[1] + (*A2++) * x[2] + (*A3++) * x[3] + (*A4++) * x[4] + (*A5++) * x[5] + (*A6++) * x[6] + (*A7++) * x[7];
    }
  }

inline double get_variance_factor(const double* residual, int m, int n)
  {
  double sigma = 0.0;
  int mn = m * n;
  for (int i = 0; i < mn; ++i)
    sigma += residual[i] * residual[i];
  double r = mn - 8;
  return sigma / r;
  }

inline int cholesky(double* AtA)
  {
  int n = 8;
  for (int i = 0; i < n; ++i)
    {
    for (int j = 0; j <= i; ++j)
      {
      double s = 0;
      for (int k = 0; k < j; ++k)
        s += AtA[i * n + k] * AtA[j * n + k];
      AtA[i * n + j] = (i == j) ?
        std::sqrt(AtA[i * n + i] - s) :
        (1.0 / AtA[j * n + j] * (AtA[i * n + j] - s));
      }
    }
  return 0;
  }

inline int solve_alsc_cholesky(double* U, double* Atl)
  {
  int n = 8;
  for (int i = 0; i < n; ++i)
    {
    double sum = Atl[i];
    for (int k = 0; k < i; ++k)
      sum -= U[i*n + k] * Atl[k];
    Atl[i] = sum / U[i*n + i];
    }
  for (int i = n - 1; i >= 0; --i)
    {
    double sum = Atl[i];
    for (int k = i + 1; k < n; ++k)
      {
      sum -= U[k*n + i] * Atl[k];
      }
    Atl[i] = sum / U[i*n + i];
    }
  return 0;
  }

inline bool update_coefficients(double& a, double& b, double& c, double& d, double& e, double& f, double& rs, double& rg, const double* x)
  {
  for (int i = 0; i < 8; ++i)
    if (x[i] != x[i])
      return false;
  rs += x[6];
  rg += x[7];
  const double denom = 1.0 + rg;
  if (fabs(denom) < 1e-12)
    return false;
  a += x[0] / denom;
  b += x[1] / denom;
  c += x[2] / denom;
  d += x[3] / denom;
  e += x[4] / denom;
  f += x[5] / denom;
  return true;
  }

inline bool terminate(const double* x, double rg, double tolerance[6])
  {
  const double denom = 1.0 + rg;
  return (fabs(x[0] / denom) < tolerance[0] &&
    fabs(x[1] / denom) < tolerance[1] &&
    fabs(x[2] / denom) < tolerance[2] &&
    fabs(x[3] / denom) < tolerance[3] &&
    fabs(x[4] / denom) < tolerance[4] &&
    fabs(x[5] / denom) < tolerance[5]);
  }

inline int inverse_from_cholesky(double* U)
  {
  double A[64];
  std::memcpy((void*)A, (void*)U, sizeof(double) * 64);
  int n = 8;
  for (int j = 0; j < n; ++j)
    {
    for (int i = 0; i < n; ++i)
      {
      double sum = i == j ? 1.0 : 0.0;
      for (int k = 0; k < i; ++k)
        sum -= A[i*n + k] * U[j*n + k];
      U[j*n + i] = sum / A[i*n + i];
      }
    }
  for (int j = 0; j < n; ++j)
    {
    for (int i = n - 1; i >= 0; --i)
      {
      double sum = U[j*n + i];
      for (int k = i + 1; k < n; ++k)
        {
        sum -= A[k*n + i] * U[j*n + k];
        }
      U[j*n + i] = sum / A[i*n + i];
      }
    }
  return 0;
  }

inline void compute_standard_deviations(
  double& sigma_a, double& sigma_b, double& sigma_c,
  double& sigma_d, double& sigma_e, double& sigma_f, double& sigma_radiometric_shift, double& sigma_radiometric_gain,
  double variance_factor, const double* AtA_inv)
  {
  sigma_a = std::sqrt(variance_factor * AtA_inv[0]);
  sigma_b = std::sqrt(variance_factor * AtA_inv[1 * 8 + 1]);
  sigma_c = std::sqrt(variance_factor * AtA_inv[2 * 8 + 2]);
  sigma_d = std::sqrt(variance_factor * AtA_inv[3 * 8 + 3]);
  sigma_e = std::sqrt(variance_factor * AtA_inv[4 * 8 + 4]);
  sigma_f = std::sqrt(variance_factor * AtA_inv[5 * 8 + 5]);
  sigma_radiometric_shift = std::sqrt(variance_factor * AtA_inv[6 * 8 + 6]);
  sigma_radiometric_gain = std::sqrt(variance_factor * AtA_inv[7 * 8 + 7]);
  }

inline void compute_standard_deviations(double& sigma_a, double& sigma_d, double variance_factor, const double* AtA_inv)
  {
  sigma_a = std::sqrt(variance_factor * AtA_inv[0]);
  sigma_d = std::sqrt(variance_factor * AtA_inv[3 * 8 + 3]);
  }

inline double quality_of_match(double variance_factor, double* AtA_inv)
  {
  double a = variance_factor * AtA_inv[0];
  double b = variance_factor * AtA_inv[3 * 8 + 3];
  double c = variance_factor * AtA_inv[3 * 8]; // take upper part of AtA_inv, because lower part is invalidated
  double apb = a + b;
  double D = apb * apb - 4.0*(a*b - c * c);
  double eig0 = std::abs(apb + std::sqrt(D)) / 2.0;
  double eig1 = std::abs(apb - std::sqrt(D)) / 2.0;
  return eig0 > eig1 ? eig0 : eig1;
  }


inline void get_extrema_patch(double& x1, double& y1, double& x2, double& y2, double& x3, double& y3, double& x4, double& y4, double a, double b, double c, double d, double e, double f, int m, int n)
  {
  x1 = a - b - c;
  y1 = d - e - f;
  x2 = a + b * (n)-c;
  y2 = d + e * (n)-f;
  x3 = a + b * (n)+c * (m);
  y3 = d + e * (n)+f * (m);
  x4 = a - b + c * (m);
  y4 = d - e + f * (m);
  }

inline bool alsc_single_iteration(double& a, double& b, double& c, double& d, double& e, double& f, double& r_shift, double& r_gain,
  const double* patch_f, double* patch_g,
  double* A,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  const uint8_t* rhimage, int stride)
  {
  resample_patch(patch_g, m, n, a, b, c, d, e, f, rhimage, stride);
  //print_patch(patch_g, m, n);
  fill_design_matrix(A, patch_g, m, n);
  //print_mat(A, m*n, 8);
  fill_patch_difference(l, patch_f, patch_g, r_shift, r_gain, m, n);
  //print_mat(l, m*n, 1);
  compute_At_A(AtA, A, m, n);
  //print_mat(AtA, 8, 8);
  compute_At_l(Atl, A, l, m, n);
  //print_mat(Atl, 8, 1);
  int res = cholesky(AtA);
  if (res)
    return false;
  res = solve_alsc_cholesky(AtA, Atl);
  if (res)
    return false;
  return update_coefficients(a, b, c, d, e, f, r_shift, r_gain, Atl);
  }

inline bool in_safety_zone(double x, double y, int w, int h, int m, int n)
  {
  if (x < m)
    return false;
  if (y < n)
    return false;
  if (x >= (w - m))
    return false;
  if (y >= (h - n))
    return false;
  return true;
  }

inline bool alsc(
  double* output,
  double left_x, double left_y,
  double right_x, double right_y,
  const uint8_t* lhimage,
  const uint8_t* rhimage,
  double sin_angle,
  double cos_angle,
  int w,
  int h,
  double* patch_f,
  double* patch_g,
  double* A,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  double tol[6],
  int max_iter)
  {
  if (!in_safety_zone(left_x, left_y, w, h, m, n))
    return false;
  if (!in_safety_zone(right_x, right_y, w, h, m, n))
    return false;
  double a, d;
  get_shift_coefficients(a, d, cos_angle, -sin_angle, sin_angle, cos_angle, m, n, left_x, left_y);
  resample_patch(patch_f, m, n, a, cos_angle, -sin_angle, d, sin_angle, cos_angle, lhimage, w);

  output[alsc_output_variable::B] = 1.0;
  output[alsc_output_variable::C] = 0.0;
  output[alsc_output_variable::E] = 0.0;
  output[alsc_output_variable::F] = 1.0;
  output[alsc_output_variable::RG] = 0.0;
  output[alsc_output_variable::RS] = 0.0;
  output[alsc_output_variable::RESIDU] = -1.0;
  bool success = false;
  bool error = false;
  get_shift_coefficients(output[alsc_output_variable::A], output[alsc_output_variable::D], output[alsc_output_variable::B], output[alsc_output_variable::C], output[alsc_output_variable::E], output[alsc_output_variable::F], m, n, right_x, right_y);
  int iterations;
  for (iterations = 0; iterations < max_iter; ++iterations)
    {
    if (!error)
      error = !alsc_single_iteration(output[alsc_output_variable::A], output[alsc_output_variable::B], output[alsc_output_variable::C], output[alsc_output_variable::D], output[alsc_output_variable::E], output[alsc_output_variable::F], output[alsc_output_variable::RS], output[alsc_output_variable::RG],
        patch_f, patch_g, A, l, AtA, Atl, m, n, rhimage, w);
    double x1, y1, x2, y2, x3, y3, x4, y4;
    get_extrema_patch(x1, y1, x2, y2, x3, y3, x4, y4, output[alsc_output_variable::A], output[alsc_output_variable::B], output[alsc_output_variable::C], output[alsc_output_variable::D], output[alsc_output_variable::E], output[alsc_output_variable::F], m, n);
    if (x1 < 0.0 || y1 < 0.0 || x2 < 0.0 || y2 < 0.0 || x3 < 0.0 || y3 < 0.0 || x4 < 0.0 || y4 < 0.0 ||
      x1 > w - 1 || y1 > h - 1 || x2 > w - 1 || y2 > h - 1 || x3 > w - 1 || y3 > h - 1 || x4 > w - 1 || y4 > h - 1)
      {
      error = true;
      success = false;
      }
    if (!error && terminate(Atl, output[RG], tol))
      {
      ++iterations;
      success = true;
      break;
      }
    }
  if (!error)
    {
    output[ITERATIONS] = double(iterations);
    inverse_from_cholesky(AtA);
    get_residual(A, Atl, l, m, n);
    output[VARIANCE_FACTOR] = get_variance_factor(l, m, n);
    compute_standard_deviations(output[alsc_output_variable::SIGMA_X], output[alsc_output_variable::SIGMA_Y], output[alsc_output_variable::VARIANCE_FACTOR], AtA);

    compute_patch_center(output[alsc_output_variable::RIGHT_X], output[alsc_output_variable::RIGHT_Y], output[alsc_output_variable::A], output[alsc_output_variable::B], output[alsc_output_variable::C], output[alsc_output_variable::D], output[alsc_output_variable::E], output[alsc_output_variable::F], m, n);
    output[QUALITY] = quality_of_match(output[alsc_output_variable::VARIANCE_FACTOR], AtA);
    resample_patch(patch_g, m, n, output[alsc_output_variable::A], output[alsc_output_variable::B], output[alsc_output_variable::C], output[alsc_output_variable::D], output[alsc_output_variable::E], output[alsc_output_variable::F], rhimage, w);
    output[RESIDU] = std::sqrt(compute_residual_l2_norm_sqr(patch_f, patch_g, output[alsc_output_variable::RS], output[alsc_output_variable::RG], m, n)) / float(m*n);
    }
  else
    success = false;
  return success;
  }

inline alsc_result alsc(double right_x, double right_y,
  const uint8_t* rhimage,
  int w,
  int h,
  const double* patch_f,
  double* patch_g,
  double* A,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  double tol[6],
  int max_iter)
  {
  alsc_result result;
  result.success = false;
  result.b = 1.0;
  result.c = 0.0;
  result.e = 0.0;
  result.f = 1.0;
  result.r_gain = 0.0;
  result.r_shift = 0.0;
  result.residu = -1;
  get_shift_coefficients(result.a, result.d, result.b, result.c, result.e, result.f, m, n, right_x, right_y);
  bool error = false;
  for (result.iterations = 0; result.iterations < max_iter; ++result.iterations)
    {
    if (!error)
      error = !alsc_single_iteration(result.a, result.b, result.c, result.d, result.e, result.f, result.r_shift, result.r_gain,
        patch_f, patch_g, A, l, AtA, Atl, m, n, rhimage, w);

    double x1, y1, x2, y2, x3, y3, x4, y4;
    get_extrema_patch(x1, y1, x2, y2, x3, y3, x4, y4, result.a, result.b, result.c, result.d, result.e, result.f, m, n);
    if (x1 < 0.0 || y1 < 0.0 || x2 < 0.0 || y2 < 0.0 || x3 < 0.0 || y3 < 0.0 || x4 < 0.0 || y4 < 0.0 ||
      x1 > w - 1 || y1 > h - 1 || x2 > w - 1 || y2 > h - 1 || x3 > w - 1 || y3 > h - 1 || x4 > w - 1 || y4 > h - 1
      || x1 != x1 || y1 != y1 || x2 != x2 || y2 != y2 || x3 != x3 || y3 != y3 || x4 != x4 || y4 != y4)
      {
      error = true;
      result.success = false;
      }
    if (terminate(Atl, result.r_gain, tol))
      {
      ++result.iterations;
      result.success = true;
      break;
      }
    }

  if (!error)
    {
    inverse_from_cholesky(AtA);
    get_residual(A, Atl, l, m, n);
    result.variance_factor = get_variance_factor(l, m, n);
    compute_standard_deviations(result.sigma_a, result.sigma_b, result.sigma_c, result.sigma_d, result.sigma_e, result.sigma_f, result.sigma_rs, result.sigma_rg, result.variance_factor, AtA);

    compute_patch_center(result.right_x, result.right_y, result.a, result.b, result.c, result.d, result.e, result.f, m, n);
    result.quality = quality_of_match(result.variance_factor, AtA);
    resample_patch(patch_g, m, n, result.a, result.b, result.c, result.d, result.e, result.f, rhimage, w);
    result.residu = std::sqrt(compute_residual_l2_norm_sqr(patch_f, patch_g, result.r_shift, result.r_gain, m, n)) / double(m*n);
    }
  else
    result.success = false;
  return result;
  }


inline alsc_result find_match_on_line(
  double x, double y,
  const double* line,
  const uint8_t* rhimage,
  int w,
  int h,
  const double* patch_f, double* patch_g,
  double* A,
  double* l,
  double* AtA,
  double* Atl,
  int m, int n,
  int step_size,
  int max_disparity,
  double tol[6],
  int max_iter
)
  {
  using namespace alsc_details;
  double x1, y1, x2, y2;
  if (fabs(line[1]) < fabs(line[0]))
    {
    y1 = 0.0;
    x1 = (-line[2] - y1 * line[1]) / line[0];
    y2 = h - 1;
    x2 = (-line[2] - y2 * line[1]) / line[0];
    }
  else
    {
    x1 = 0.0;
    y1 = (-line[2] - x1 * line[0]) / line[1];
    x2 = w - 1;
    y2 = (-line[2] - x2 * line[0]) / line[1];
    }
  alsc_result result;
  result.success = false;
  result.residu = 1000.0;
  result.variance_factor = 1000.0;
  const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
  if (steep)
    {
    if (y1 > y2)
      {
      dswap(x1, x2);
      dswap(y1, y2);
      }
    }
  else
    {
    if (x1 > x2)
      {
      dswap(x1, x2);
      dswap(y1, y2);
      }
    }
  double px, py;
  int max_it = (max_disparity / step_size) * 2;
  for (int i = 0; i < max_it; ++i)
    {
    int sgn = (i % 2) ? 1 : -1;
    if (steep)
      {
      py = y + sgn * (i / 2)*step_size;
      px = (py - y1)*(x2 - x1) / (y2 - y1) + x1;
      }
    else
      {
      px = x + sgn * (i / 2)*step_size;
      py = (px - x1)*(y2 - y1) / (x2 - x1) + y1;
      }
    if (px > n+1  && px < w - n-2 && py > m+1  && py < h - m-2)
      {
      auto res = alsc(px, py, rhimage, w, h, patch_f, patch_g, A, l, AtA, Atl, m, n, tol, max_iter);
      if (res.success && res.residu < result.residu)
        {
        result = res;
        }
      }
    }
  return result;
  }


}
