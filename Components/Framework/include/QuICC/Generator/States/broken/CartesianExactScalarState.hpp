/**
 * @file CartesianExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a cartesian geometries
 */

#ifndef QUICC_EQUATIONS_CARTESIANEXACTSCALARSTATE_HPP
#define QUICC_EQUATIONS_CARTESIANEXACTSCALARSTATE_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Generator/States/CartesianExactStateIds.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact scalar state in cartesian geometries
    */
   class CartesianExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         CartesianExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CartesianExactScalarState();

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
         void setStateType(const CartesianExactStateIds::Id id);

         /**
          * @brief Set the options for the solution states in 2D
          *
          * @param a1   Amplitude of the first direction
          * @param k1   Wave number of the first direction
          * @param a2   Amplitude of the second direction
          * @param k2   Wave number of the second direction
          */
         void setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2);

         /**
          * @brief Set the options for the solution states in 3D
          *
          * @param a1   Amplitude of the first direction
          * @param k1   Wave number of the first direction
          * @param a2   Amplitude of the second direction
          * @param k2   Wave number of the second direction
          * @param a3   Amplitude of the second direction
          * @param k3   Wave number of the second direction
          */
         void setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Build grids
          */
         void buildGrid(Array& g1D, Array& g2D, Array& g3D) const;

         /**
          * @brief Type of the state to generate
          */
         CartesianExactStateIds::Id mTypeId;

         /**
          * @brief Amplitude of the state
          */
         Array mModeA;

         /**
          * @brief Mode number of the state (wave number of polynomial order)
          */
         Array mModeK;
   };

   /// Typedef for a shared CartesianExactScalarState
   typedef std::shared_ptr<CartesianExactScalarState> SharedCartesianExactScalarState;

}
}

#endif // QUICC_EQUATIONS_CARTESIANEXACTSCALARSTATE_HPP
