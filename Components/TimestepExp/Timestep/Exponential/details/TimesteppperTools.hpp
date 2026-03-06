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

} // namespace details
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_DETAILS_TIMESTEPPERTOOLS_HPP
