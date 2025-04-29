/**
 * @file RandomVectorState.hpp
 * @brief Implementation of the equation to generate a random vector state
 */

#ifndef QUICC_EQUATIONS_RANDOMVECTORSTATE_HPP
#define QUICC_EQUATIONS_RANDOMVECTORSTATE_HPP

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
#include "QuICC/Generator/States/RandomStateData.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate a random vector state
    */
   class RandomVectorState: public IVectorEquation
   {
      public:
         /// Typedef for SpecialId
         typedef RandomStateData::SpecialId SpecialId;

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         RandomVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RandomVectorState();

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
         void setSpectrum(const FieldComponents::Spectral::Id comp, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId special = RandomStateData::NOTHING);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const FieldComponents::Spectral::Id compId, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId special = RandomStateData::NOTHING);

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
          * @brief tools and data for random data generation
          */
         RandomStateData mRandom;
   };

   /// Typedef for a shared RandomVectorState
   typedef std::shared_ptr<RandomVectorState> SharedRandomVectorState;

}
}

#endif // QUICC_EQUATIONS_RANDOMVECTORSTATE_HPP
