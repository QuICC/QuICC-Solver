/**
 * @file PLl.hpp
 * @brief Implementation of the associated Legendre based LL integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_PLl_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_PLl_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/PP.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Implementation of the associated Legendre based P integrator
    */
   template<typename OpTypes = PIALegendreOperatorTypes>
   class PLl: public PP<OpTypes>
   {
      public:
        using OpArray = typename OpTypes::OpArray;
        using OpMatrix = typename OpTypes::OpMatrix;
        using OpMatrixR = typename OpTypes::OpMatrixR;
        using OpMatrixCR = typename OpTypes::OpMatrixCR;

        using OpVectorI = typename OpTypes::OpVectorI;
        using OpMatrixI = typename OpTypes::OpMatrixI;
        using OpMatrixLZ = typename OpTypes::OpMatrixLZ;
        using OpMatrixL = typename OpTypes::OpMatrixL;

        using DataType = typename OpTypes::DataType;

        /**
         * @brief Constructor
         */
        PLl();

        /**
         * @brief Destructor
         */
        virtual ~PLl();

      private:
        virtual void applyUnitOperator(const OpMatrixLZ &rOut,
           const OpMatrixLZ &in, const OpVectorI &scan,
           const int totalOpsCols) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const;

        /**
         * @brief l(l+1) scaling factors
         */
        virtual void initSpecial() const override;

        /**
         * @brief Storage for l(l+1) factors
         */
        mutable Array mLl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_P_HPP
