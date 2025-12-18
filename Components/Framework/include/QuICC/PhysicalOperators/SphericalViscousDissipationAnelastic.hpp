/**
 * @file SphericalViscousDissipationAnelastic.hpp
 * @brief Implementation of the spherical coriolis term
 */

 #ifndef QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP
 #define QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP

 // System includes
 //
 
 // Project includes
 //
 #include "Types/Typedefs.hpp"
 #include "QuICC/Enums/FieldIds.hpp"
 #include "QuICC/VectorFields/VectorField.hpp"
 #include "QuICC/Resolutions/Resolution.hpp"
 #include "QuICC/ScalarFields/ScalarField.hpp"
 #include "QuICC/Equations/IVectorEquation.hpp"
 
#include "DenseSM/IGenericProfile.hpp"
 
 namespace QuICC {
 
 namespace Physical {
 
    /**
     * @brief Implementation of the spherical coriolis term
     */
    class SphericalViscousDissipationAnelastic
    {
       public:
          /**
           * @brief Set S to viscous dissipation term
           */
          static void set(Framework::Selector::PhysicalScalarField &rS, 
                          const Resolution& res, 
                          const Array& r, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                          const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv, 
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, 
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);
 
          /**
           * @brief Add viscous dissipation term to S
           */
          static void add(Framework::Selector::PhysicalScalarField &rS, 
                          const Resolution& res, 
                          const Array& r, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                          const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv, 
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, 
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);
 
          /**
           * @brief Substract viscous dissipation term from S
           */
          static void sub(Framework::Selector::PhysicalScalarField &rS, 
                          const Resolution& res, 
                          const Array& r, 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                          const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv, 
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, 
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);

         static void test(Framework::Selector::PhysicalScalarField &rS, 
                          const Resolution& res, 
                          const Array& r, 
                          const Array& thGrid,    // Add theta grid
                          const Array& phGrid,    // Add phi grid 
                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, 
                          const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv, 
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF,
                          std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, 
                          const QuICC::Equations::EquationParameters &eqParams,
                          const MHDFloat c = 1.0);
 
       protected:
 
       private:
          /**
           * @brief Empty constructor
           */
          SphericalViscousDissipationAnelastic() = default;
 
          /**
           * @brief Empty destructor
           */
          ~SphericalViscousDissipationAnelastic() = default;

          /**
           * @brief helper function to calculate viscous dissipation slice
           */
          static Eigen::Matrix<MHDFloat, 
                               Eigen::Dynamic, 
                               Eigen::Dynamic> computeViscousSlice(const int iR,
                                                                     const int iR_,
                                                                     const MHDFloat c,
                                                                     const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                                                                           FieldComponents::Physical::Id>& v,
                                                                     const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, 
                                                                           FieldComponents::Physical::Id>& Dv,
                                                                     const MHDFloat nu,
                                                                     const MHDFloat T,
                                                                     const MHDFloat Rho,
                                                                     const MHDFloat dLogRho);

    };
 }
 }
 
 #endif // QUICC_PHYSICAL_SPHERICALVISCOUSDISSIPATIONANELASTIC_HPP
 