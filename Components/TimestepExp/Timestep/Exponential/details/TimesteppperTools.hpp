/**
 * @file ITimestepper.hpp
 * @brief Implementation of base for the templated (coupled) equation
 * timestepper
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_DETAILS_TIMESTEPPERTOOLS_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_DETAILS_TIMESTEPPERTOOLS_HPP

// System includes
//
#include <Eigen/Dense>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "View/View.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace details {
/**
 * @brief Compute z = y
 */
template <typename TData> void computeSet(TData& z, const TData& y);

/**
 * @brief Compute z = y
 */
template <typename TData> void computeSet(TData& y, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow);

/**
 * @brief Compute y = x
 */
template <typename TData> void computeSet(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& y, const TData& x, const std::size_t startRow);

/**
 * @brief Compute y = x
 */
void computeSet(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& y, const DecoupledZMatrix& x, const std::size_t startRow);

/**
 * @brief Compute y = x
 */
template <typename TData> void computeSet(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& y, const TData& x, const std::size_t startRow);

/**
 * @brief Compute z = y
 */
void computeSet(DecoupledZMatrix& z, const DecoupledZMatrix& y);

/**
 * @brief Compute z = a*y
 */
template <typename TData>
void computeSet(TData& z, const MHDFloat a, const TData& y);

/**
 * @brief Compute z = a*y
 */
void computeSet(DecoupledZMatrix& z, const MHDFloat a,
   const DecoupledZMatrix& y);

/**
 * @brief Flatten data y_2n = x_n.re, y_2n+1 = x_n.im
 */
template <typename TData> void flatten2Real(Matrix& y, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow, const std::size_t col);

/**
 * @brief Unflatten data y_n = x_2n + x_2n+1 * j
 */
void unflatten2Complex(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& y, const Matrix& x, const std::size_t startRow, const std::size_t col);

/**
 * @brief Compute y = a*x + y
 */
void computeAXPY(Matrix& y, const MHDFloat a, const Matrix& x);

/**
 * @brief Compute y = a*x + y
 */
void flattenAXPY(Matrix& y, const MHDFloat a, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow, const std::size_t col);

/**
 * @brief Apply correction from tuple
 */
void addCorrection(Matrix& rVal, const std::vector<std::tuple<MHDComplex,int,int>>& corr, const int rows, const int cols, const std::size_t startRow, const std::size_t col);

/**
 * @brief Add (scaled) real part of decoupled storage to sparse matrix
 */
void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat);
//
//
//
//
//

template <typename TData> inline void computeSet(TData& y, const TData& x)
{
   y = x;
}

inline void computeSet(DecoupledZMatrix& y, const DecoupledZMatrix& x)
{
   y.real() = x.real();

   y.imag() = x.imag();
}

template <typename TData> void computeSet(TData& y, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow)
{
   assert(static_cast<std::size_t>(y.cols()) == x.dims()[1]);
   assert(static_cast<std::size_t>(y.rows()) >= x.dims()[0] + startRow);

   for(std::size_t j = 0;  j < x.dims()[1]; j++)
   {
      for(std::size_t i = 0;  i < x.dims()[0]; i++)
      {
         y(i + startRow,j) = x(i,j);
      }
   }
}

inline void computeSet(DecoupledZMatrix& y, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow)
{
   assert(static_cast<std::size_t>(y.real().cols()) == x.dims()[1]);
   assert(static_cast<std::size_t>(y.imag().cols()) == x.dims()[1]);
   assert(static_cast<std::size_t>(y.real().rows()) >= x.dims()[0] + startRow);
   assert(static_cast<std::size_t>(y.imag().rows()) >= x.dims()[0] + startRow);

   for(std::size_t j = 0;  j < x.dims()[1]; j++)
   {
      for(std::size_t i = 0;  i < x.dims()[0]; i++)
      {
         std::size_t i_ = i + startRow;
         y.real()(i_,j) = x(i,j).real();
         y.imag()(i_,j) = x(i,j).imag();
      }
   }
}

template <typename TData> void computeSet(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& y, const TData& x, const std::size_t startRow)
{
   assert(static_cast<std::size_t>(x.cols()) == y.dims()[1]);
   assert(static_cast<std::size_t>(x.rows()) >= y.dims()[0] + startRow);

   for(std::size_t j = 0;  j < y.dims()[1]; j++)
   {
      for(std::size_t i = 0;  i < y.dims()[0]; i++)
      {
         y(i,j) = x(i + startRow,j);
      }
   }
}

inline void computeSet(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& y, const DecoupledZMatrix& x, const std::size_t startRow)
{
   assert(static_cast<std::size_t>(x.real().cols()) == y.dims()[1]);
   assert(static_cast<std::size_t>(x.imag().cols()) == y.dims()[1]);
   assert(static_cast<std::size_t>(x.real().rows()) >= y.dims()[0] + startRow);
   assert(static_cast<std::size_t>(x.imag().rows()) >= y.dims()[0]  + startRow);

   for(std::size_t j = 0;  j < y.dims()[1]; j++)
   {
      for(std::size_t i = 0;  i < y.dims()[0]; i++)
      {
         std::size_t i_ = i + startRow;
         y(i,j) = MHDComplex(x.real()(i_,j), x.imag()(i_,j));
      }
   }
}

template <typename TData>
inline void computeSet(TData& y, const MHDFloat a, const TData& x)
{
   y = a * x;
}

inline void computeSet(DecoupledZMatrix& y, const MHDFloat a,
   const DecoupledZMatrix& x)
{
   y.real() = a * x.real();

   y.imag() = a * x.imag();
}

inline void flatten2Real(Matrix& y, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow, const std::size_t col)
{
   assert(static_cast<std::size_t>(y.cols()) > col);
   assert(static_cast<std::size_t>(y.rows()) >= 2*x.dims()[0]*x.dims()[1] + startRow);

   size_t kRe = startRow;
   size_t kIm = kRe + x.dims()[1]*x.dims()[0];
   for(std::size_t j = 0;  j < x.dims()[1]; j++)
   {
      for(std::size_t i = 0;  i < x.dims()[0]; i++)
      {
         const MHDComplex& z = x(i,j);
         y(kRe, col) = z.real();
         y(kIm, col) = z.imag();
         kRe++;
         kIm++;
      }
   }
}

inline void unflatten2Complex(View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& y, const Matrix& x, const std::size_t startRow, const std::size_t col)
{
   assert(static_cast<std::size_t>(x.cols()) > col);
   assert(static_cast<std::size_t>(x.rows()) >= 2*y.dims()[0]*y.dims()[1] + startRow);

   size_t kRe = startRow;
   size_t kIm = kRe + y.dims()[1]*y.dims()[0];
   for(std::size_t j = 0;  j < y.dims()[1]; j++)
   {
      for(std::size_t i = 0;  i < y.dims()[0]; i++)
      {
         MHDComplex z(x(kRe,col), x(kIm,col));
         y(i,j) = z;
         kRe++;
         kIm++;
      }
   }
}

inline void computeAXPY(Matrix& y, const MHDFloat a, const Matrix& x)
{
   assert(x.rows() == y.rows());
   assert(x.cols() == y.cols());

   if(a != 0)
   {
      if(a == 1.0)
      {
         y += x;
      }
      else
      {
         y += a*x;
      }
   }
}

inline void flattenAXPY(Matrix& y, const MHDFloat a, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow,  const std::size_t col)
{
   assert(static_cast<std::size_t>(y.cols()) > col);
   assert(static_cast<std::size_t>(y.rows()) >= 2*x.dims()[0]*x.dims()[1] + startRow);

   if(a != 0)
   {
      size_t kRe = startRow;
      size_t kIm = kRe + x.dims()[1]*x.dims()[0];
      if(a == 1.0)
      {
         for(std::size_t j = 0;  j < x.dims()[1]; j++)
         {
            for(std::size_t i = 0;  i < x.dims()[0]; i++)
            {
               const MHDComplex& z = x(i,j);
               y(kRe, col) += z.real();
               y(kIm, col) += z.imag();
               kRe++;
               kIm++;
            }
         }
      }
      else
      {
         for(std::size_t j = 0;  j < x.dims()[1]; j++)
         {
            for(std::size_t i = 0;  i < x.dims()[0]; i++)
            {
               const MHDComplex& z = a*x(i,j);
               y(kRe, col) += z.real();
               y(kIm, col) += z.imag();
               kRe++;
               kIm++;
            }
         }
      }
   }
}

inline void addCorrection(Matrix& rVal, const std::vector<std::tuple<MHDComplex,int,int>>& corr, const int rows, const int cols, const std::size_t startRow, const std::size_t col)
{
   for(auto&& c: corr)
   {
      auto&& val = std::get<0>(c);
      auto&& i = std::get<1>(c);
      auto&& j = std::get<2>(c);
      std::size_t kRe = startRow + (i + j*rows);
      std::size_t kIm = kRe + rows*cols;

      rVal(kRe, col) += val.real();
      rVal(kIm, col) += val.imag();
   }
}

inline void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat)
{
   assert(decMat.real().rows() > 0);
   assert(decMat.real().cols() > 0);
   assert(decMat.imag().size() == 0 || decMat.imag().nonZeros() == 0);

   if(c != 1.0)
   {
      mat += c*decMat.real();
   } else
   {
      mat += decMat.real();
   }
}

} // namespace details
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_DETAILS_TIMESTEPPERTOOLS_HPP
