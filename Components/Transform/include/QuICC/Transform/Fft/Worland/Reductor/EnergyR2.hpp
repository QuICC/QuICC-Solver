/** 
 * @file EnergyR2.hpp
 * @brief Implementation of the Worland based R^2 energy operator
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGYR2_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGYR2_HPP

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
    * @brief Implementation of the Worland based R^2 energy operator
    */ 
   class EnergyR2: public IWorlandEnergy
   {
      public:
         /**
          * @brief Constructor
          */
         EnergyR2();

         /**
          * @brief Destructor
          */
         virtual ~EnergyR2();
         
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

#endif // QUICC_TRANSFORM_FFT_WORLAND_REDUCTOR_ENERGYR2_HPP
