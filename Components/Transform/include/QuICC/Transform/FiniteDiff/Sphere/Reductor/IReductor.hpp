/**
 * @file IReductor.hpp
 * @brief Interface for a Finite Difference based reduction operator
 */

#ifndef QUICC_TRANSFORM_FINITEDIFF_SPHERE_REDUCTOR_IREDUCTOR_HPP
#define QUICC_TRANSFORM_FINITEDIFF_SPHERE_REDUCTOR_IREDUCTOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/IOperator.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

namespace Sphere {

namespace Reductor {

   /**
    * @brief Interface for a Finite Differences based energy operator
    */
   class IReductor: public IOperator
   {
      public:
         /**
          * @brief Constructor
          */
         IReductor();

         /**
          * @brief Destructor
          */
         virtual ~IReductor() = default;

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;

      protected:
         /**
          * @brief Storage for the operators
          */
         mutable std::vector<Matrix>  mOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable Internal::Array  mGrid;

      private:
   };

} // namespace Reductor
} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FINITEDIFF_SPHERE_REDUCTOR_IREDUCTOR_HPP
