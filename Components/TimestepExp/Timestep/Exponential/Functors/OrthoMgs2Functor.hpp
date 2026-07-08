/**
 * @file OrthoMgs2Functor.hpp
 * @brief Modified Gram-Schmidt orthogonalization with re-orthogonalization
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOMGS2FUNCTOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOMGS2FUNCTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

namespace Functors {

/**
 * @brief Modified Gram-Schmidt with re-orthogonalization
 */
class OrthoMgs2Functor
{
public:
   /// @brief ctor
   OrthoMgs2Functor(const int p);

   /// @brief dtor
   ~OrthoMgs2Functor() = default;

   /**
    * Orthogonalize
    */
   double operator()(Matrix& matV, Matrix& matH, const int j, const int n);

   /**
    * @brief Gram-schmidt order
    */
   int p() const;

private:
   /**
    * @brief Length of incomplete orthogonalization
    */
   const int mcP;

};

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOMGS2FUNCTOR_HPP
