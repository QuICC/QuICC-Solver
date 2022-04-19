/** 
 * @file EnergyD1R1.hpp
 * @brief Implementation of the Worland based D R energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGYD1R1_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGYD1R1_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Worland/Reductor/IWorlandEnergy.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Reductor {

   /**
    * @brief Implementation of the Worland based D R energy operator
    */ 
   class EnergyD1R1: public IWorlandEnergy
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyD1R1();

         /**
          * @brief Destructor
          */
         virtual ~EnergyD1R1();
         
      protected:

      private:
         /**
          * @brief Make operator
          */
         virtual void makeOperator(Matrix& op, Matrix& eop, const internal::Array& igrid, const internal::Array& iweights, const Array& eweights, const int i) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGYD1R1_HPP
