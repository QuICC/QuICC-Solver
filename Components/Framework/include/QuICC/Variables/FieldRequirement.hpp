/**
 * @file FieldRequirement.hpp
 * @brief Implementation of a class to store requirements for a field
 */

#ifndef QUICC_FIELDREQUIREMENT_HPP
#define QUICC_FIELDREQUIREMENT_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Tools/ComponentAlias.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a class to store requirements for a field
    */
   class FieldRequirement
   {
      public:
         /**
          * @brief General constructor with default requirements
          */
         FieldRequirement(const bool isScalar, const Tools::ComponentAlias<FieldComponents::Spectral::Id> spec, const Tools::ComponentAlias<FieldComponents::Physical::Id> phys);

         /**
          * @brief Destructor
          */
         ~FieldRequirement();

         /**
          * brief Set field to scalar
          */
         void enableScalar();

         /**
          * brief Set field to vector
          */
         void enableVector();

         /**
          * brief Enable spectral space field values
          */
         void enableSpectral();

         /**
          * brief Enable physical space field values
          */
         void enablePhysical();

         /**
          * brief Enable physical space gradient values
          */
         void enableGradient();

         /**
          * brief Enable physical space curl values
          */
         void enableCurl();

         /**
          * brief enable physical space second gradient values
          */
         void enableGradient2();

         /**
          * @brief Check if it is scalar field
          */
         bool isScalar() const;

         /**
          * @brief Get spectral requirement
          */
         bool needSpectral() const;

         /**
          * @brief Get physical field requirement
          */
         bool needPhysical() const;

         /**
          * @brief Get physical gradient requirement
          */
         bool needPhysicalGradient() const;

         /**
          * @brief Get physical curl requirement
          */
         bool needPhysicalCurl() const;

         /**
          * @brief Get physical 2nd order gradient requirement
          */
         bool needPhysicalGradient2() const;

         /**
          * @brief Get the physical field components requirements
          */
         const ArrayB& physicalComps() const;

         /**
          * @brief Get the physical gradient components requirements
          */
         const ArrayB& gradientComps(const FieldComponents::Spectral::Id id) const;

         /**
          * @brief Get the physical curl components requirements
          */
         const ArrayB& curlComps() const;

         /**
          * @brief Get the physical 2nd order gradient components requirements
          */
         const MatrixB& gradient2Comps(const FieldComponents::Spectral::Id id) const;

         /**
          * @brief Get vector of spectral field components
          */
         const std::vector<FieldComponents::Spectral::Id>&  spectralIds() const;

         /**
          * @brief Get vector of physical field components
          */
         const std::vector<FieldComponents::Physical::Id>& physicalIds() const;

         /**
          * @brief Get map for physical components to field requirements
          */
         std::map<FieldComponents::Physical::Id,bool> mapPhysicalComps() const;

         /**
          * @brief Get map for gradient components to field requirements
          */
         std::map<FieldComponents::Physical::Id,bool> mapGradientComps(const FieldComponents::Spectral::Id id) const;

         /**
          * @brief Get map for curl components to field requirements
          */
         std::map<FieldComponents::Physical::Id,bool> mapCurlComps() const;

         /**
          * @brief Get map for 2nd order gradient components to field requirements
          */
         std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool> mapGradient2Comps(const FieldComponents::Spectral::Id id) const;

         /**
          * @brief Update the physical component requirements
          */
         void updatePhysical(const ArrayB& comps);

         /**
          * @brief Update the gradient component requirements
          */
         void updateGradient(const std::map<FieldComponents::Spectral::Id,ArrayB>& comps);

         /**
          * @brief Update the curl component requirements
          */
         void updateCurl(const ArrayB& comps);

         /**
          * @brief Update the 2nd order gradient component requirements
          */
         void updateGradient2(const std::map<FieldComponents::Spectral::Id,MatrixB>& comps);

         /**
          * @brief Merge information from other requirements
          *
          * @param req Requirements to merge
          */
         void merge(const FieldRequirement& req);

      protected:

      private:
         /**
          * @brief Initialise defaut spectral and physical IDs
          */
         void initDefaultIds(const Tools::ComponentAlias<FieldComponents::Spectral::Id> spec, const Tools::ComponentAlias<FieldComponents::Physical::Id> phys);

         /**
          * @brief Is field a scalar field?
          */
         bool  mIsScalar;

         /**
          * @brief Storage for spectral storage requirements
          */
         bool  mNeedSpectral;

         /**
          * @brief Storage for physical storage requirements
          */
         bool  mNeedPhysical;

         /**
          * @brief Storage for physical gradient storage requirements
          */
         bool  mNeedGradient;

         /**
          * @brief Storage for physical curl storage requirements
          */
         bool  mNeedCurl;

         /**
          * @brief Storage for physical 2nd order gradient storage requirements
          */
         bool  mNeedGradient2;

         /**
          * @brief List of spectral field components
          */
         std::vector<FieldComponents::Spectral::Id>  mSpectralIds;

         /**
          * @brief List of physical field components
          */
         std::vector<FieldComponents::Physical::Id>  mPhysicalIds;

         /**
          * @brief Detailed requirements for physical field components
          */
         ArrayB mPhysicalComps;

         /**
          * @brief Detailed requirements for gradient field components
          */
         std::map<FieldComponents::Spectral::Id,ArrayB> mGradientComps;

         /**
          * @brief Detailed requirements for curl field components
          */
         ArrayB mCurlComps;

         /**
          * @brief Detailed requirements for 2nd order gradient field components
          */
         std::map<FieldComponents::Spectral::Id,MatrixB> mGradient2Comps;
   };

}

#endif // QUICC_FIELDREQUIREMENT_HPP
