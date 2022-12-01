/**
 * @file IALegendreIntegrator.hpp
 * @brief Interface for a associated Legendre based integrator
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_PIALEGENDREINTEGRATOR_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_PIALEGENDREINTEGRATOR_HPP

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
#include "QuICC/Transform/Poly/ALegendre/Integrator/IALegendreIntegrator.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

namespace Integrator {

   /**
    * @brief Interface for a associated Legendre based integrator
    */
    template<typename OpTypes = PIALegendreOperatorTypes>
   class PIALegendreIntegrator: public IALegendreIntegrator<OpTypes>
   /* class PIALegendreIntegrator: public IALegendreOperator<OpTypes> */
   {
      public:
         using OpArray = typename OpTypes::OpArray;
         using OpMatrix = typename OpTypes::OpMatrix;
         using OpMatrixZ = typename OpTypes::OpMatrixZ;
         using OpMatrixR = typename OpTypes::OpMatrixR;
         using OpMatrixCR = typename OpTypes::OpMatrixCR;
         using DataType = typename OpTypes::DataType;

         using OpVectorI = typename OpTypes::OpVectorI;
         using OpMatrixI = typename OpTypes::OpMatrixI;
         using OpMatrixLZ = typename OpTypes::OpMatrixLZ;
         using OpMatrixL = typename OpTypes::OpMatrixL;

         /**
          * @brief Constructor
          */
         PIALegendreIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~PIALegendreIntegrator();

         /**
          * @brief Rows of output data
          */
         virtual int outRows() const override;

         /**
          * @brief Columns of output data
          */
         virtual int outCols() const override;

         /**
          * @brief Get the memory requirements
          */
         /* MHDFloat requiredStorage() const; */

      protected:
         /**
          * @brief Storage for the operators
          */
         mutable OpMatrixL vmOps;

         /**
          * @brief Storage for the quadrature grid
          */
         mutable OpArray  mGrid;

         /**
          * @brief Storage for the quadrature weights
          */
         mutable OpArray  mWeights;

         //Apply operator using the entire OpMatrix in parallel instead  of block by block

         virtual void applyUnitOperator(const OpMatrixLZ &rOut,
            const OpMatrixLZ &in, const OpVectorI &scan,
            const int totalOpsCols) const = 0;

      private:
         /**
          * @brief Initialise the operators
          */
         virtual void initOperators(const OpArray& igrid, const OpArray& iweights) const;

         /**
          * @brief Make operator
          */
         virtual void makeOperator(OpMatrix& op, const OpArray& igrid, const OpArray& iweights, const int i) const = 0;

         /**
          * @brief Compute projection
          *
          * @param rOut Output values
          * @param in   Input values
          */
         void applyOperators(OpMatrixZ &rOut, const OpMatrixZ &in) const;

         /**
          * @brief Apply ith operator
          */
         virtual void applyOperator(OpMatrixR rOut, const int i, const OpMatrixCR& in) const = 0;

         /**
          * @brief Allow specific initialisation for operators
          */
         virtual void initSpecial() const;

   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_INTEGRATOR_IALEGENDREINTEGRATOR_HPP
