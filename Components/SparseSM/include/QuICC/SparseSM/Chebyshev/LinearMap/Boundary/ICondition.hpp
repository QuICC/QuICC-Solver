/**
 * @file ICondition.hpp
 * @brief Implementation of the generic interface for a boundary condition for the Cheyshev sparse operator based on a linear map y = ax + b, x = [-1, 1] (natural chebyshev grid)
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_ICONDITION_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_ICONDITION_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   /**
    * @brief Implementation of the generic interface for a boundary condition for the Cheyshev sparse operator based on a linear map y = ax + b, x = [-1, 1] (natural chebyshev grid)
    */
   class ICondition
   {
      public:
         /// Typedef for scalar
         typedef Internal::MHDFloat Scalar_t;

         /// Typedef for coefficient array
         typedef Internal::ACoeff ACoeff_t;

         /// Position enum
         enum class Position { BOTTOM = -1, TOP = 1 };

         /**
          * @brief Constructor
          *
          * @param lower   Lower bound of y
          * @param upper   Upper bound of y
          * @param position   Boundary position
          */
         ICondition(const Scalar_t lower, const Scalar_t upper, const Position pos);

         /**
          * @brief Destructor
          */
         virtual ~ICondition() = default;

      protected:
         /**
          * @brief Get position
          */
         Position position() const;

         /**
          * @brief Get mapping a coefficient from y = ax + b
          */
         Scalar_t a() const;

         /**
          * @brief Get mapping b coefficient from y = ax + b
          */
         Scalar_t b() const;

         /**
          * Chebyshev normalization factor
          */
         Scalar_t c() const;

      private:
         /**
          * @brief Compute mapping
          *
          * @param lower   Lower bound
          * @param upper   Upper bound
          */
         void setBounds(const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Boundary position
          */
         Position mPosition;

         /**
          * @brief a coefficienct of y = ax + b
          */
         Scalar_t mA;

         /**
          * @brief b coefficienct of y = ax + b
          */
         Scalar_t mB;
   };

} // Boundary
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_BOUNDARY_ICONDITION_HPP
