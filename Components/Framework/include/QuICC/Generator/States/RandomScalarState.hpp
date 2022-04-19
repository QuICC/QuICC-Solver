/**
 * @file RandomScalarState.hpp
 * @brief Implementation of the equation to generate a random scalar state
 */

#ifndef QUICC_EQUATIONS_RANDOMSCALARSTATE_HPP
#define QUICC_EQUATIONS_RANDOMSCALARSTATE_HPP

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
#include "QuICC/Framework/Selector/ScalarField.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Generator/States/RandomStateData.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate a random scalar state
    */
   class RandomScalarState: public IScalarEquation
   {
      public:
         /// Typedef for SpecialId
         typedef RandomStateData::SpecialId SpecialId;

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         RandomScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RandomScalarState();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false);

         /**
          * @brief Compute the random state as a source term
          *
          * @param compId  ID of the spectral component
          * @param i      Fast direction
          * @param j      Middle direction
          * @param k      Slow direction
          */
         virtual MHDVariant sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId special = RandomStateData::NOTHING);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId special = RandomStateData::NOTHING);

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
          * @brief tools and data for random data generation
          */
         RandomStateData mRandom;
   };

   /// Typedef for a shared RandomScalarState
   typedef std::shared_ptr<RandomScalarState> SharedRandomScalarState;

}
}

#endif // QUICC_EQUATIONS_RANDOMSCALARSTATE_HPP
