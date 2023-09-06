/** 
 * @file ICondition.hpp
 * @brief Interface to generic Worland boundary condition
 */

#ifndef QUICC_SPARSESM_WORLAND_BOUNDARY_ICONDITION_HPP
#define QUICC_SPARSESM_WORLAND_BOUNDARY_ICONDITION_HPP

// System includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/SparseSM/Worland/IDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   /**
    * @brief Interface to generic Worland boundary condition
    */ 
   class ICondition: public IDiags
   {
      public:
         /**
          * @brief Constructor
          *
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          * @param q       Truncation q (only consider rows - q equations)
          */
         ICondition(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q);

         /**
          * @brief Destructor
          */
         virtual ~ICondition() = default;
         
      protected:

      private:
   };

} // Boundary
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_BOUNDARY_ICONDITION_HPP
