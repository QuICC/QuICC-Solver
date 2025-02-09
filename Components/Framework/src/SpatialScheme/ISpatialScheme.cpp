/** 
 * @file ISpatialScheme.cpp
 * @brief Source of the interface for a spatial scheme
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Transform/Setup/Default.hpp"

namespace QuICC {

namespace SpatialScheme {

   ISpatialScheme::ISpatialScheme(const VectorFormulation::Id formulation, const GridPurpose::Id purpose, const int dimension, const std::size_t id, const std::string tag, const std::string formatted)
      : mFormulation(formulation), mPurpose(purpose), mDimension(dimension), mId(id), mTag(tag), mFormatted(formatted)
   {
      // Initialize default implementation options
      std::vector<std::size_t> def = {Transform::Setup::Default::id()};
      for(int i = 0; i < dimension; i++)
      {
         this->mImplType.emplace(i, def);
      }
   }

   VectorFormulation::Id ISpatialScheme::formulation() const
   {
      return this->mFormulation;
   }

   GridPurpose::Id ISpatialScheme::purpose() const
   {
      return this->mPurpose;
   }

   int ISpatialScheme::dimension() const
   {
      return this->mDimension;
   }

   std::size_t ISpatialScheme::id() const
   {
      return this->mId;
   }

   std::string ISpatialScheme::tag() const
   {
      return this->mTag;
   }

   std::string ISpatialScheme::formatted() const
   {
      return this->mFormatted;
   }

   const Tools::ComponentAlias<FieldComponents::Physical::Id>& ISpatialScheme::physical() const
   {
      return this->mPhys;
   }

   const Tools::ComponentAlias<FieldComponents::Spectral::Id>& ISpatialScheme::spectral() const
   {
      return this->mSpec;
   }

   void ISpatialScheme::setImplementation(const std::map<std::size_t,std::vector<std::size_t>>& type)
   {
      this->mImplType = type;
   }

   bool ISpatialScheme::has(const Feature f) const
   {
      return this->mFeatures.count(f);
   }

   void ISpatialScheme::enable(const Feature f)
   {
      this->mFeatures.insert(f);
   }

   void ISpatialScheme::enable(const std::set<Feature>& f)
   {
      this->mFeatures.insert(f.begin(), f.end());
   }

   void ISpatialScheme::disable(const Feature f)
   {
      if(this->has(f))
      {
         this->mFeatures.erase(f);
      }
   }

   const std::set<Feature>& ISpatialScheme::features() const
   {
      return this->mFeatures;
   }

   void ISpatialScheme::useCustomMesher(std::shared_ptr<IMesher> m)
   {
      this->mspCustomMesher = m;
   }

} // SpatialScheme
} // QuICC
