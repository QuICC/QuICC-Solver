/**
 * @file OrthoCgs2Functor.hpp
 * @brief Classical Gram-Schmidt orthogonalization with re-orthogonalization
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGS2FUNCTOR_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGS2FUNCTOR_HPP

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
 * @brief Classical Gram-Schmidt with re-orthogonalization
 */
class OrthoCgs2Functor
{
public:
   /// @brief ctor
   OrthoCgs2Functor(const int p);

   /// @brief dtor
   ~OrthoCgs2Functor() = default;

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

#endif // QUICC_TIMESTEP_EXPONENTIAL_FUNCTORS_ORTHOCGS2FUNCTOR_HPP
