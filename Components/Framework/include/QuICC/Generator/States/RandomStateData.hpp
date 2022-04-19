/**
 * @file RandomStateData.hpp
 * @brief Generic data storage for random state generator
 */

#ifndef QUICC_EQUATIONS_RANDOMSTATEDATA_HPP
#define QUICC_EQUATIONS_RANDOMSTATEDATA_HPP

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
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Generic data storage for random state generator
    */
   class RandomStateData
   {
      public:
         enum SpecialId {
            ONLYMEAN, ///< Only mean
            NOMEAN, ///< No mean
            NOTHING, ///< Nothing special
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         RandomStateData();

         /**
          * @brief Simple empty destructor
          */
         virtual ~RandomStateData();

         /**
          * @brief Set Shared resolution
          */
         void setResolution(SharedResolution spRes);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const FieldComponents::Spectral::Id comp, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const SpecialId special = NOTHING);

         /**
          * @brief Set spectrum variables
          */
         void setSpectrum(const FieldComponents::Spectral::Id comp, const MHDFloat min, const MHDFloat max, const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D, const SpecialId special = NOTHING);

         /**
          * Generate Random value
          */
         void makeRandom(MHDFloat& val, const int i1D, const int i3D, const int i2D, const FieldComponents::Spectral::Id comp, const unsigned int seed = 1) const;

         /**
          * Generate Random value
          */
         void makeRandom(MHDComplex& val, const int i1D, const int i3D, const int i2D, const FieldComponents::Spectral::Id comp) const;

         /**
          * @brief Compute scaling ratio
          */
         MHDFloat scalingRatio(const FieldComponents::Spectral::Id comp, const int i, const int dim, const int n) const;

         /**
          * @brief Compute the random state as a source term
          *
          * @param compId  ID of the spectral component
          * @param i      Fast direction
          * @param j      Middle direction
          * @param k      Slow direction
          */
         MHDVariant sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

      protected:
         /**
          * @brief Get resolution
          */
         const Resolution& res() const;

      private:
         /**
          * @brief Minimum value in coefficient range
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mMin;

         /**
          * @brief Maximum value in coefficient range
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mMax;

         /**
          * @brief Ratio between first and last coefficient in first direction
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mRatio1D;

         /**
          * @brief Ratio between first and last coefficient in second direction
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mRatio2D;

         /**
          * @brief Ratio between first and last coefficient in third direction
          */
         std::map<FieldComponents::Spectral::Id,MHDFloat> mRatio3D;

         /**
          * @brief Shared resolution
          */
         SharedResolution mspRes;

         /**
          * @brief Starting seed used for random generator
          */
         int mStartSeed;

         /**
          * @brief ID flag for special setup
          */
         SpecialId mSpecial;
   };

}
}

#endif // QUICC_EQUATIONS_RANDOMSTATEDATA_HPP
