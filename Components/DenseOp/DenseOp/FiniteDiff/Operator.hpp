/**
 * @file Operator.hpp
 * @brief Implementation of the generic interface for a dense finite difference operator
 */

#ifndef QUICC_DENSEOP_FINITEDIFF_OPERATOR_HPP
#define QUICC_DENSEOP_FINITEDIFF_OPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "DenseOp/FiniteDiff/IBaseOperator.hpp"
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {
namespace DenseOp {
namespace FiniteDiff {

/// @brief Wrapper for generic Finite Differences operator
/// @tparam TFdBuilder FD builder
template <class TFdBuilder> class Operator : public IBaseOperator
{
public:
   /// @brief Pass-by-value FD builder ctor
   /// @param fdBuilder to be stored and used
   Operator(TFdBuilder fdBuilder) : mFdBuilder(fdBuilder){};

   /// @brief default ctor
   Operator() = default;

   /// @brief dtor
   ~Operator() = default;

   /// @brief populate op matrix
   /// @param op
   /// @param grid
   /// @param weights
   /// @param l
   void compute(Eigen::Ref<Matrix> op, const Internal::Array& grid,
      const std::uint32_t l) final;

private:
   TFdBuilder mFdBuilder;
};

template <class TFdBuilder>
void Operator<TFdBuilder>::compute(Eigen::Ref<Matrix> op,
   const Internal::Array& grid, const std::uint32_t l)
{
   mFdBuilder.compute(op, static_cast<int>(l), grid);
}


} // namespace FiniteDiff
} // namespace DenseOp
} // namespace QuICC

#endif // define QUICC_DENSEOP_FINITEDIFF_OPERATOR_HPP
