/**
 * @file ITimestepper.hpp
 * @brief Implementation of base for the templated (coupled) equation
 * timestepper
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_DETAILS_TIMESTEPPERTOOLS_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_DETAILS_TIMESTEPPERTOOLS_HPP

// System includes
//
#include <Eigen/Dense>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "View/View.hpp"

#include <iostream>
namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

namespace Views {

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
 * @brief Compute z = A*y
 */
template <typename TOperator, typename TData>
void computeMV(TData& z, const TOperator& mat, const TData& y);

/**
 * @brief Compute z = A*y
 */
void computeMV(DecoupledZMatrix& z, const SparseMatrix& mat,
   const DecoupledZMatrix& y);

/**
 * @brief Compute z = a*x + b*y + z
 */
template <typename TData>
void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x,
   const MHDFloat b, const TData& y);

/**
 * @brief Compute z = a*x + b*y + z
 */
void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a,
   const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y);

/**
 * @brief Compute y = a*M*x + y
 */
template <typename TOperator, typename TData>
void computeAMXPY(TData& y, const TOperator& mat, const MHDFloat a,
   const TData& x);

/**
 * @brief Compute y = a*M*x + y
 */
void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x);

/**
 * @brief Compute z = a*M*x + y + b*z
 */
template <typename TData>
void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const TData& y, const MHDFloat b);

/**
 * @brief Compute z = a*M*x + y + b*z
 */
void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y,
   const MHDFloat b);

/**
 * @brief Compute z = a*M*x + b*y + z
 */
template <typename TData>
void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y);

/**
 * @brief Compute z = a*M*x + b*y + z
 */
void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y);

/**
 * @brief Compute z = a*M*x + b*y + M*z
 */
template <typename TData>
void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y);

/**
 * @brief Compute z = a*M*x + b*y + M*z
 */
void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y);

/**
 * @brief Compute y = a*x + y
 */
template <typename TData>
void computeAXPY(TData& y, const MHDFloat a, const TData& x);

/**
 * @brief Compute y = a*x + y
 */
template <typename TData>
void computeAXPY(TData& y, const MHDFloat a, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow);

/**
 * @brief Compute y = a*x + y
 */
void computeAXPY(DecoupledZMatrix& y, const MHDFloat a,
   const DecoupledZMatrix& x);

/**
 * @brief Compute y = a*x + y
 */
void computeAXPY(DecoupledZMatrix& y, const MHDFloat a, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow);

/**
 * @brief Compute y = x + a*y
 */
template <typename TData>
void computeXPAY(TData& y, const TData& x, const MHDFloat a);

/**
 * @brief Compute y = x + a*y
 */
void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x,
   const MHDFloat a);

/**
 * @brief Compute error
 */
template <typename TData>
void computeErrorFromDiff(MHDFloat& err, const TData& diff, const TData& ref);

/**
 * @brief Compute error
 */
void computeErrorFromDiff(MHDFloat& err, const DecoupledZMatrix& diff,
   const DecoupledZMatrix& ref);

/**
 * @brief Compute error
 */
template <typename TData>
void computeError(MHDFloat& err, const TData& x, const TData& y);

/**
 * @brief Compute error
 */
void computeError(MHDFloat& err, const DecoupledZMatrix& x,
   const DecoupledZMatrix& y);

/**
 * @brief Initialize influence matrix kernel and boundary condition
 *
 * @param kernel     3 parts, initialize boundary condition and boundary value
 * @param val        Boundary value of inhomogeneous problem
 * @param bc         Boundary conditions
 */
template <typename TData>
void initInfluence(TData& kernel, const DecoupledZSparse& val,
   const DecoupledZSparse& bc);

/**
 * @brief Initialize influence matrix kernel and boundary condition
 *
 * @param kernel     3 parts, initialize boundary condition and boundary value
 * @param val        Boundary value of inhomogeneous problem
 * @param bc         Boundary conditions
 */
void initInfluence(DecoupledZMatrix& kernel, const DecoupledZSparse& val,
   const DecoupledZSparse& bc);

/**
 * @brief Store influence matrix kernel and boundary condition
 *
 * @param reg  Kernel storage with 3 parts: kernel, boundary value and boundary
 * condition
 * @param x    Kernel solution to store
 */
template <typename TData> void computeSetInfluence(TData& reg, const TData& x);

/**
 * @brief Store influence matrix kernel and boundary condition
 *
 * @param reg  Kernel storage with 3 parts: kernel, boundary value and boundary
 * condition
 * @param x    Kernel solution to store
 */
void computeSetInfluence(DecoupledZMatrix& reg, const DecoupledZMatrix& x);

/**
 * @brief Compute correction from influence matrix (Green's function)
 *
 * @param y Split solution to correct
 * @param x 3 parts of influence matrix solution
 */
template <typename TData>
void computeInfluenceCorrection(TData& y, const TData& x);

/**
 * @brief Compute correction from influence matrix (Green's function)
 *
 * @param y Split solution to correct
 * @param x 3 parts of influence matrix solution
 */
void computeInfluenceCorrection(DecoupledZMatrix& y, const DecoupledZMatrix& x);

/**
 * @brief Add (scaled) real part of decoupled storage to sparse matrix
 */
void addOperators(SparseMatrix& mat, const MHDFloat c, const DecoupledZSparse& decMat);

/**
 * @brief Add (scaled) decoupled storage to sparse matrix
 */
void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat);

/**
 * @brief Apply correction to decoupled storage from sparse matrix
 */
void addCorrection(DecoupledZMatrix& rVal, const SparseMatrixZ& corr);

/**
 * @brief Apply correction
 */
template <typename TData, typename TCorr> void addCorrection(TData& rVal, const TCorr& corr);

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

template <typename TData>
inline void computeAXPY(TData& y, const MHDFloat a, const TData& x)
{
   if (a != 0.0)
   {
      y += a * x;
   }
}

inline void computeAXPY(DecoupledZMatrix& y, const MHDFloat a,
   const DecoupledZMatrix& x)
{
   if (a != 0.0)
   {
      y.real() += a * x.real();

      y.imag() += a * x.imag();
   }
}

template <typename TData>
void computeAXPY(TData& y, const MHDFloat a, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow)
{
   assert(static_cast<std::size_t>(y.cols()) == x.dims()[1]);
   assert(static_cast<std::size_t>(y.rows()) >= x.dims()[0] + startRow);

   if(a != 0)
   {
      if(a == 1.0)
      {
         for(std::size_t j = 0;  j < x.dims()[1]; j++)
         {
            for(std::size_t i = 0;  i < x.dims()[0]; i++)
            {
               y(i+startRow,j) += x(i,j);
            }
         }
      }
      else
      {
         for(std::size_t j = 0;  j < x.dims()[1]; j++)
         {
            for(std::size_t i = 0;  i < x.dims()[0]; i++)
            {
               y(i+startRow,j) += a*x(i,j);
            }
         }
      }
   }
}

inline void computeAXPY(DecoupledZMatrix& y, const MHDFloat a, const View::View<MHDComplex, View::Attributes<View::DimLevelType<View::dense_t, View::dense_t>>>& x, const std::size_t startRow)
{
   assert(static_cast<std::size_t>(y.real().cols()) == x.dims()[1]);
   assert(static_cast<std::size_t>(y.imag().cols()) == x.dims()[1]);
   assert(static_cast<std::size_t>(y.real().rows()) >= x.dims()[0] + startRow);
   assert(static_cast<std::size_t>(y.imag().rows()) >= x.dims()[0] + startRow);

   if(a != 0)
   {
      if(a == 1.0)
      {
         for(std::size_t j = 0;  j < x.dims()[1]; j++)
         {
            for(std::size_t i = 0;  i < x.dims()[0]; i++)
            {
               std::size_t i_ = i + startRow;
               y.real()(i_,j) += x(i,j).real();
               y.imag()(i_,j) += x(i,j).imag();
            }
         }
      }
      else
      {
         for(std::size_t j = 0;  j < x.dims()[1]; j++)
         {
            for(std::size_t i = 0;  i < x.dims()[0]; i++)
            {
               std::size_t i_ = i + startRow;
               y.real()(i_,j) += a*x(i,j).real();
               y.imag()(i_,j) += a*x(i,j).imag();
            }
         }
      }
   }
}

template <typename TData>
inline void computeXPAY(TData& y, const TData& x, const MHDFloat a)
{
   if (a != 0.0)
   {
      y = x + a * y;
   }
   else
   {
      computeSet<TData>(y, x);
   }
}

inline void computeXPAY(DecoupledZMatrix& y, const DecoupledZMatrix& x,
   const MHDFloat a)
{
   if (a != 0.0)
   {
      y.real() = x.real() + a * y.real();

      y.imag() = x.imag() + a * y.imag();
   }
   else
   {
      computeSet(y, x);
   }
}

template <typename TData>
inline void computeXPAYPBZ(TData& z, const TData& x, const MHDFloat a,
   const TData& y, const MHDFloat b)
{
   if (a == 0.0)
   {
      z = x + b * z;
   }
   else if (b == 0.0)
   {
      z = x + a * y;
   }
   else
   {
      z = x + a * y + b * z;
   }
}

inline void computeXPAYPBZ(DecoupledZMatrix& z, const DecoupledZMatrix& x,
   const MHDFloat a, const DecoupledZMatrix& y, const MHDFloat b)
{
   if (a == 0.0)
   {
      z.real() = x.real() + b * z.real();

      z.imag() = x.imag() + b * z.imag();
   }
   else if (b == 0.0)
   {
      z.real() = x.real() + a * y.real();

      z.imag() = x.imag() + a * y.imag();
   }
   else
   {
      z.real() = x.real() + a * y.real() + b * z.real();

      z.imag() = x.imag() + a * y.imag() + b * z.imag();
   }
}

template <typename TData>
inline void computeAXPBYPZ(TData& z, const MHDFloat a, const TData& x,
   const MHDFloat b, const TData& y)
{
   if (a == 0.0)
   {
      z += b * y;
   }
   else if (b == 0.0)
   {
      z += a * x;
   }
   else
   {
      z += a * x + b * y;
   }
}

inline void computeAXPBYPZ(DecoupledZMatrix& z, const MHDFloat a,
   const DecoupledZMatrix& x, const MHDFloat b, const DecoupledZMatrix& y)
{
   if (a == 0.0)
   {
      z.real() += b * y.real();

      z.imag() += b * y.imag();
   }
   else if (b == 0.0)
   {
      z.real() += a * x.real();

      z.imag() += a * x.imag();
   }
   else
   {
      z.real() += a * x.real() + b * y.real();

      z.imag() += a * x.imag() + b * y.imag();
   }
}

template <typename TOperator, typename TData>
void computeAMXPY(TData& y, const TOperator& mat, const MHDFloat a,
   const TData& x)
{
   if (a != 0.0)
   {
      y += mat * (a * x);
   }
}

inline void computeAMXPY(DecoupledZMatrix& y, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x)
{
   if (a != 0.0)
   {
      y.real() += mat * (a * x.real());

      y.imag() += mat * (a * x.imag());
   }
}

template <typename TData>
void computeAMXPYPBZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const TData& y, const MHDFloat b)
{
   if (a == 0.0)
   {
      z = y + b * z;
   }
   else
   {
      z = mat * (a * x) + y + b * z;
   }
}

inline void computeAMXPYPBZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const DecoupledZMatrix& y,
   const MHDFloat b)
{
   if (a == 0.0)
   {
      z.real() = y.real() + b * z.real();

      z.imag() = y.imag() + b * z.imag();
   }
   else
   {
      z.real() = mat * (a * x.real()) + y.real() + b * z.real();

      z.imag() = mat * (a * x.imag()) + y.imag() + b * z.imag();
   }
}

template <typename TData>
void computeAMXPBYPZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y)
{
   if (a == 0.0)
   {
      z += b * y;
   }
   else if (b == 0.0)
   {
      z += mat * (a * x);
   }
   else
   {
      z += mat * (a * x) + b * y;
   }
}

inline void computeAMXPBYPZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y)
{
   if (a == 0.0)
   {
      z.real() += b * y.real();

      z.imag() += b * y.imag();
   }
   else if (b == 0.0)
   {
      z.real() += mat * (a * x.real());

      z.imag() += mat * (a * x.imag());
   }
   else
   {
      z.real() += mat * (a * x.real()) + b * y.real();

      z.imag() += mat * (a * x.imag()) + b * y.imag();
   }
}

template <typename TData>
void computeAMXPBYPMZ(TData& z, const SparseMatrix& mat, const MHDFloat a,
   const TData& x, const MHDFloat b, const TData& y)
{
   if (a == 0.0)
   {
      z = b * y + mat * z;
   }
   else if (b == 0.0)
   {
      z = mat * (a * x + z);
   }
   else
   {
      z = mat * (a * x + z) + b * y;
   }
}

inline void computeAMXPBYPMZ(DecoupledZMatrix& z, const SparseMatrix& mat,
   const MHDFloat a, const DecoupledZMatrix& x, const MHDFloat b,
   const DecoupledZMatrix& y)
{
   if (a == 0.0)
   {
      z.real() = b * y.real() + mat * z.real();

      z.imag() = b * y.imag() + mat * z.imag();
   }
   else if (b == 0.0)
   {
      z.real() = mat * (a * x.real() + z.real());

      z.imag() = mat * (a * x.imag() + z.imag());
   }
   else
   {
      z.real() = mat * (a * x.real() + z.real()) + b * y.real();

      z.imag() = mat * (a * x.imag() + z.imag()) + b * y.imag();
   }
}

template <typename TOperator, typename TData>
inline void computeMV(TData& y, const TOperator& A, const TData& x)
{
   y = A * x;
}

inline void computeMV(DecoupledZMatrix& y, const SparseMatrix& A,
   const DecoupledZMatrix& x)
{
   y.real() = A * x.real();

   y.imag() = A * x.imag();
}

template <typename TData>
inline void computeErrorFromDiff(MHDFloat& err, const TData& diff,
   const TData& ref)
{
   err = std::max(err,
      (diff.array() / (1.0 + ref.array().abs())).abs().maxCoeff());
}

inline void computeErrorFromDiff(MHDFloat& err, const DecoupledZMatrix& diff,
   const DecoupledZMatrix& ref)
{
   err = std::max(err, (diff.real().array() / (1.0 + ref.real().array().abs()))
                          .abs()
                          .maxCoeff());

   err = std::max(err, (diff.imag().array() / (1.0 + ref.imag().array().abs()))
                          .abs()
                          .maxCoeff());
}

template <typename TData>
inline void computeError(MHDFloat& err, const TData& x, const TData& y)
{
   err = std::max(err,
      ((x.array() - y.array()) / (1.0 + x.array().abs())).abs().maxCoeff());
}

inline void computeError(MHDFloat& err, const DecoupledZMatrix& x,
   const DecoupledZMatrix& y)
{
   err = std::max(err,
      ((x.real().array() - y.real().array()) / (1.0 + x.real().array().abs()))
         .abs()
         .maxCoeff());

   err = std::max(err,
      ((x.imag().array() - y.imag().array()) / (1.0 + x.imag().array().abs()))
         .abs()
         .maxCoeff());
}

template <typename TData>
inline void initInfluence(TData& y, const DecoupledZSparse& val,
   const DecoupledZSparse& bc)
{
   throw std::logic_error("Not yet implemented");
}

inline void initInfluence(DecoupledZMatrix& y, const DecoupledZSparse& val,
   const DecoupledZSparse& bc)
{
   assert(bc.real().rows() >= val.real().cols());

   // Real value
   int cols = val.real().cols();
   y.real().resize(val.real().rows(), 3 * cols);

   for (int i = 0; i < cols; i++)
   {
      y.real().col(3 * i + 1) = bc.real().row(i).transpose();
      y.real().col(3 * i + 2) = val.real().col(i);
   }

   // Imaginary value
   cols = val.imag().cols();
   y.imag().resize(val.imag().rows(), 3 * cols);

   for (int i = 0; i < cols; i++)
   {
      y.imag().col(3 * i + 1) = bc.real().row(i).transpose();
      y.imag().col(3 * i + 2) = val.imag().col(i);
   }
}

template <typename TData>
inline void computeSetInfluence(TData& reg, const TData& x)
{
   assert(reg.cols() == x.cols() * 3);

   for (int i = 0; i < x.cols(); i++)
   {
      reg.col(3 * i) = x.col(i);
   }
}

inline void computeSetInfluence(DecoupledZMatrix& reg,
   const DecoupledZMatrix& x)
{
   assert(reg.real().cols() == 3 * x.real().cols());
   assert(reg.imag().cols() == 3 * x.imag().cols());

   for (int i = 0; i < x.real().cols(); i++)
   {
      reg.real().col(3 * i) = x.real().col(i);
   }

   for (int i = 0; i < x.imag().cols(); i++)
   {
      reg.imag().col(3 * i) = x.imag().col(i);
   }
}

template <typename TData>
void computeInfluenceCorrection(TData& y, const TData& x)
{
   assert(x.cols() % 3 == 0);

   // compute influence matrix
   int nBC = x.cols() / 3;
   TData mat(nBC, nBC);
   for (int j = 0; j < nBC; j++)
   {
      const auto bc = x.col(3 * j + 1).transpose();
      for (int k = 0; k < nBC; k++)
      {
         mat(j, k) = (bc * x.col(3 * k)).value();
      }
   }
   mat = mat.inverse();

   TData bcVal(nBC, 1);
   for (int j = 0; j < y.cols(); j++)
   {
      for (int k = 0; k < nBC; k++)
      {
         bcVal(k, 0) = (x.col(3 * k + 1).transpose() * y.col(j)).value();
      }
      bcVal = mat * bcVal;
      for (int k = 0; k < nBC; k++)
      {
         y.col(j) -= bcVal(k) * x.col(3 * k);
      }
   }
}

inline void computeInfluenceCorrection(DecoupledZMatrix& y,
   const DecoupledZMatrix& x)
{
   assert(x.real().cols() % 3 == 0);
   assert(x.imag().cols() % 3 == 0);

   // compute influence matrix
   int nBC = x.real().cols() / 3;
   Matrix mat(nBC, nBC);
   for (int j = 0; j < nBC; j++)
   {
      const auto bc = x.real().col(3 * j + 1).transpose();
      for (int k = 0; k < nBC; k++)
      {
         mat(j, k) = bc * x.real().col(3 * k);
      }
   }
   mat = mat.inverse();

   Matrix bcVal(nBC, 1);
   for (int j = 0; j < y.real().cols(); j++)
   {
      for (int k = 0; k < nBC; k++)
      {
         auto bc = x.real().col(3 * k + 1).transpose();
         bcVal(k, 0) = bc * y.real().col(j);
      }
      bcVal = mat * bcVal;
      for (int k = 0; k < nBC; k++)
      {
         y.real().col(j) -= bcVal(k) * x.real().col(3 * k);
      }
   }

   for (int j = 0; j < y.imag().cols(); j++)
   {
      for (int k = 0; k < nBC; k++)
      {
         auto bc = x.imag().col(3 * k + 1).transpose();
         bcVal(k, 0) = bc * y.imag().col(j);
      }
      bcVal = mat * bcVal;
      for (int k = 0; k < nBC; k++)
      {
         y.imag().col(j) -= bcVal(k) * x.imag().col(3 * k);
      }
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

inline void addOperators(SparseMatrixZ& mat, const MHDFloat c, const DecoupledZSparse& decMat)
{
   assert(decMat.real().rows() > 0);
   assert(decMat.real().cols() > 0);
   assert(decMat.imag().rows() > 0);
   assert(decMat.imag().cols() > 0);
   assert(decMat.real().rows() == decMat.imag().rows());
   assert(decMat.real().cols() == decMat.imag().cols());

   if(c != 1.0)
   {
      mat += c*decMat.real().cast<MHDComplex>() + c*Math::cI*decMat.imag();
   } else
   {
      mat += decMat.real().cast<MHDComplex>() + Math::cI*decMat.imag();
   }
}

inline void addCorrection(DecoupledZMatrix& rVal, const SparseMatrixZ& corr)
{
   assert(rVal.real().rows() > 0);
   assert(rVal.real().cols() > 0);
   assert(rVal.imag().rows() > 0);
   assert(rVal.imag().cols() > 0);
   assert(rVal.real().rows() == rVal.imag().rows());
   assert(rVal.real().cols() == rVal.imag().cols());

   rVal.real() += corr.real();
   rVal.imag() += corr.imag();
}

template <typename TData, typename TCorr> inline void addCorrection(TData& rVal, const TCorr& corr)
{
   rVal += corr;
}

} // namespace details
} // namespace Views
} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_VIEWS_DETAILS_TIMESTEPPERTOOLS_HPP
