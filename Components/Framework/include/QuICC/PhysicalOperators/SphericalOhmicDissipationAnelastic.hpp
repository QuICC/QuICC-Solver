/**
 * @file SphericalOhmicDissipationAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALOHMICDISSIPATIONANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALOHMICDISSIPATIONANELASTIC_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"


//#include "QuICC/DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalOhmicDissipationAnelastic
   {
      public:
         /**
          * @brief Add Coriolis term to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, 
                         const Resolution& res, 
                         const Array& r, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, 
                         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF,
                         const QuICC::Equations::EquationParameters &eqParams,
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalOhmicDissipationAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalOhmicDissipationAnelastic() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALOHMICDISSIPATIONANELASTIC_HPP
