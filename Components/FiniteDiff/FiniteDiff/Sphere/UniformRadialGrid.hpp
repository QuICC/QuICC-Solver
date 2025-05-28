/**
 * @file UniformRadialGrid.hpp
 * @brief Implementation of a uniform radial grid
 */

#ifndef QUICC_FINITEDIFF_SPHERE_UNIFORMRADIALGRID_HPP
#define QUICC_FINITEDIFF_SPHERE_UNIFORMRADIALGRID_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   /**
    * @brief Implementation of a uniforma radial grid
    */
   class UniformRadialGrid
   {
      public:
         /**
          * @brief Constructor
          */
         UniformRadialGrid() = default;

         /**
          * @brief Destructor
          */
         ~UniformRadialGrid() = default;

         /**
          * @brief Compute grid for r [0, 1]
          */
         void computeGrid(Internal::Array& igrid, const int size);

      private:

   };

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_UNIFORMRADIALGRID_HPP
