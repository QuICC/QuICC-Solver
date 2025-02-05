/**
 * @file SphericalSelfAdvectionAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

#ifndef QUICC_PHYSICAL_SPHERICALSELFADVECTIONANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALSELFADVECTIONANELASTIC_HPP

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

//#include "QuICC/DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical coriolis term
    */
   class SphericalSelfAdvectionAnelastic
   {
      public:
         /**
          * @brief Set S to Coriolis term
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, 
                         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF,
                         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, 
                         const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis term to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, 
                         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF,
                         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, 
                         const MHDFloat c = 1.0);

         /**
          * @brief Substract Coriolis term from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                         const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w, 
                         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF,
                         std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, 
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalSelfAdvectionAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalSelfAdvectionAnelastic() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALSELFADVECTIONANELASTIC_HPP
