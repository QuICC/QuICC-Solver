/**
 * @file AnnulusExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector states in an annulus
 */

#ifndef QUICC_EQUATIONS_ANNULUSEXACTVECTORSTATE_HPP
#define QUICC_EQUATIONS_ANNULUSEXACTVECTORSTATE_HPP

// Configuration includes
//

// System includes
//
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Generator/States/AnnulusExactStateIds.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact vector state in an annulus
    */
   class AnnulusExactVectorState: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         AnnulusExactVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AnnulusExactVectorState();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const;

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDVariant sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set the state type id
          */
         void setStateType(const FieldComponents::Physical::Id compId, const AnnulusExactStateIds::Id id);

         /**
          * @brief Set the options for the solution states
          *
          * @param a1   Amplitude of the first direction
          * @param k1   Wave number of the first direction
          * @param a2   Amplitude of the second direction
          * @param k2   Wave number of the second direction
          * @param a3   Amplitude of the second direction
          * @param k3   Wave number of the second direction
          */
         void setModeOptions(const FieldComponents::Physical::Id compId, const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents();

      private:
         /**
          * @brief Type of the state to generate
          */
         std::map<FieldComponents::Physical::Id,AnnulusExactStateIds::Id> mTypeId;

         /**
          * @brief Amplitude of the state
          */
         std::map<FieldComponents::Physical::Id,Array> mModeA;

         /**
          * @brief Mode number of the state (wave number of polynomial order)
          */
         std::map<FieldComponents::Physical::Id,Array> mModeK;
   };

   /// Typedef for a shared AnnulusExactVectorState
   typedef std::shared_ptr<AnnulusExactVectorState> SharedAnnulusExactVectorState;

}
}

#endif // QUICC_EQUATIONS_ANNULUSEXACTVECTORSTATE_HPP
