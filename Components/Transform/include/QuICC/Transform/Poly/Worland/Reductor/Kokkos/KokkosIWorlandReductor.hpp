/**
 * @file IWorlandReductor.hpp
 * @brief Interface for a parallel Worland based reduction operator
 */

#ifndef QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_IWORLANDREDUCTOR_HPP
#define QUICC_TRANSFORM_POLY_WORLAND_REDUCTOR_KOKKOS_IWORLANDREDUCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandReductor.hpp"
#include "QuICC/Transform/Poly/KokkosIOperatorTypes.hpp"
#include "QuICC/Transform/Poly/KokkosIOperatorGemmUtils.hpp"
#ifdef QUICC_USE_KOKKOS_CUDA
#include "QuICC/Transform/Poly/CudaIOperatorTypes.hpp"
#include "QuICC/Transform/Poly/CudaIOperatorGemmUtils.hpp"
#endif


namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

/**
 * @brief Interface for a Worland based energy operator
 */
class KokkosIWorlandReductor : public IWorlandReductor
{
public:
   using OpMatrixL = typename KokkosIOperatorTypes::OpMatrixL;
   /**
    * @brief Constructor
    */
   KokkosIWorlandReductor();

   /**
    * @brief Destructor
    */
   virtual ~KokkosIWorlandReductor();

protected:
   /**
    * @brief Storage for the operators
    */
   mutable OpMatrixL vmOps;

   /**
    * @brief Storage for the quadrature grid
    */
   mutable Internal::Array mGrid;

   /**
    * @brief Storage for the quadrature weights
    */
   mutable Internal::Array mWeights;
};

} // namespace Reductor
} // namespace Worland
} // namespace Poly
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_POLY_WORLAND_IWORLANDREDUCTOR_HPP
