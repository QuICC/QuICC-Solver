/**
 * @file FieldRequirement.cpp
 * @brief Source of the variable requirements
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Variables/FieldRequirement.hpp"

// Project includes
//

namespace QuICC {

   FieldRequirement::FieldRequirement(const bool isScalar, const Tools::ComponentAlias<FieldComponents::Spectral::Id> spec, const Tools::ComponentAlias<FieldComponents::Physical::Id> phys)
      : mIsScalar(isScalar), mNeedSpectral(false), mNeedPhysical(false), mNeedGradient(false), mNeedCurl(false), mNeedGradient2(false), mPhysicalComps(3), mGradientComps(), mCurlComps(3), mGradient2Comps()
   {
      // Init default physical and spectral IDs
      this->initDefaultIds(spec, phys);

      // Set default physical
      this->mPhysicalComps.setConstant(this->mNeedPhysical);

      // Set default curl needs
      this->mCurlComps.setConstant(this->mNeedCurl);

      // Set default gradient needs
      ArrayB arr(3);
      arr.setConstant(this->mNeedGradient);
      for(auto id: this->mSpectralIds)
      {
         this->mGradientComps.insert(std::make_pair(id, arr));
      }

      // Set default 2nd order gradient needs
      MatrixB mat = MatrixB::Zero(3,3);
      mat.triangularView<Eigen::Upper>().setConstant(this->mNeedGradient2);
      for(auto id: this->mSpectralIds)
      {
         this->mGradient2Comps.insert(std::make_pair(id, mat));
      }
   }

   FieldRequirement::~FieldRequirement()
   {
   }

   void FieldRequirement::initDefaultIds(const Tools::ComponentAlias<FieldComponents::Spectral::Id> spec, const Tools::ComponentAlias<FieldComponents::Physical::Id> phys)
   {
      // Set scalar spectral component
      if(this->mIsScalar)
      {
         this->mSpectralIds.push_back(FieldComponents::Spectral::SCALAR);

      // Create default spectral components
      } else
      {
         for(auto it = spec.cbegin(); it != spec.cend(); ++it)
         {
            if(*it != FieldComponents::Spectral::NOTUSED)
            {
               this->mSpectralIds.push_back(*it);
            }
         }
      }

      // Create default physical components
      for(auto it = phys.cbegin(); it != phys.cend(); ++it)
      {
         if(*it != FieldComponents::Physical::NOTUSED)
         {
            this->mPhysicalIds.push_back(*it);
         }
      }
   }

   void FieldRequirement::enableSpectral()
   {
      this->mNeedSpectral = true;
   }

   void FieldRequirement::enablePhysical()
   {
      this->mNeedPhysical = true;

      // Enable all components
      this->mPhysicalComps.setConstant(this->mNeedPhysical);
   }

   void FieldRequirement::enableGradient()
   {
      this->mNeedGradient = true;

      // Enable all components
      for(auto& c: this->mGradientComps)
      {
         c.second.setConstant(this->mNeedGradient);
      }
   }

   void FieldRequirement::enableCurl()
   {
      if(this->isScalar())
      {
         throw std::logic_error("Tried to setup curl calculation for a scalar!");
      }

      this->mNeedCurl = true;

      // Enable all components
      this->mCurlComps.setConstant(this->mNeedCurl);
   }

   void FieldRequirement::enableGradient2()
   {
      this->mNeedGradient2 = true;

      // Enable all components
      for(auto& c: this->mGradient2Comps)
      {
         c.second.triangularView<Eigen::Upper>().setConstant(this->mNeedGradient2);
      }
   }

   bool FieldRequirement::isScalar() const
   {
      return this->mIsScalar;
   }

   bool FieldRequirement::needSpectral() const
   {
      return this->mNeedSpectral;
   }

   bool FieldRequirement::needPhysical() const
   {
      return this->mNeedPhysical;
   }

   bool FieldRequirement::needPhysicalGradient() const
   {
      return this->mNeedGradient;
   }

   bool FieldRequirement::needPhysicalCurl() const
   {
      return this->mNeedCurl;
   }

   bool FieldRequirement::needPhysicalGradient2() const
   {
      return this->mNeedGradient2;
   }

   const ArrayB& FieldRequirement::physicalComps() const
   {
      return this->mPhysicalComps;
   }

   const ArrayB& FieldRequirement::gradientComps(const FieldComponents::Spectral::Id id) const
   {
      assert(this->mGradientComps.count(id));

      return this->mGradientComps.find(id)->second;
   }

   const ArrayB& FieldRequirement::curlComps() const
   {
      return this->mCurlComps;
   }

   const MatrixB& FieldRequirement::gradient2Comps(const FieldComponents::Spectral::Id id) const
   {
      assert(this->mGradient2Comps.count(id));

      return this->mGradient2Comps.find(id)->second;
   }

   std::map<FieldComponents::Physical::Id,bool> FieldRequirement::mapPhysicalComps() const
   {
      std::map<FieldComponents::Physical::Id,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         comps.insert(std::make_pair(this->mPhysicalIds.at(i), this->mPhysicalComps(i)));
      }

      return comps;
   }

   std::map<FieldComponents::Physical::Id,bool> FieldRequirement::mapGradientComps(const FieldComponents::Spectral::Id id) const
   {
      std::map<FieldComponents::Physical::Id,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         comps.insert(std::make_pair(this->mPhysicalIds.at(i), this->mGradientComps.find(id)->second(i)));
      }

      return comps;
   }

   std::map<FieldComponents::Physical::Id,bool> FieldRequirement::mapCurlComps() const
   {
      std::map<FieldComponents::Physical::Id,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         comps.insert(std::make_pair(this->mPhysicalIds.at(i), this->mCurlComps(i)));
      }

      return comps;
   }

   std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool> FieldRequirement::mapGradient2Comps(const FieldComponents::Spectral::Id id) const
   {
      std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool> comps;

      for(unsigned int i = 0; i < this->mPhysicalIds.size(); i++)
      {
         for(unsigned int j = 0; j < this->mPhysicalIds.size(); j++)
         {
            comps.insert(std::make_pair(std::make_pair(this->mPhysicalIds.at(i),this->mPhysicalIds.at(j)), this->mGradient2Comps.find(id)->second(i,j)));
         }
      }

      return comps;
   }

   const std::vector<FieldComponents::Physical::Id>& FieldRequirement::physicalIds() const
   {
      return this->mPhysicalIds;
   }

   const std::vector<FieldComponents::Spectral::Id>& FieldRequirement::spectralIds() const
   {
      return this->mSpectralIds;
   }

   void FieldRequirement::updatePhysical(const ArrayB& comps)
   {
      this->mPhysicalComps = comps;

      this->mNeedPhysical = this->mPhysicalComps.any();
   }

   void FieldRequirement::updateGradient(const std::map<FieldComponents::Spectral::Id,ArrayB>& comps)
   {
      this->mGradientComps = comps;

      // Update gradient need
      for(auto it = this->mGradientComps.cbegin(); it != this->mGradientComps.cend(); ++it)
      {
         this->mNeedGradient = this->mNeedGradient || it->second.any();
      }
   }

   void FieldRequirement::updateCurl(const ArrayB& comps)
   {
      this->mCurlComps = comps;

      this->mNeedCurl = this->mCurlComps.any();
   }

   void FieldRequirement::updateGradient2(const std::map<FieldComponents::Spectral::Id,MatrixB>& comps)
   {
      this->mGradient2Comps = comps;

      // Update 2nd order gradient need
      for(auto it = this->mGradient2Comps.begin(); it != this->mGradient2Comps.end(); ++it)
      {
         // Clear strictly lower triangular part
         it->second.triangularView<Eigen::StrictlyLower>().setZero();

         this->mNeedGradient2 = this->mNeedGradient2 || it->second.any();
      }
   }

   void FieldRequirement::merge(const FieldRequirement& req)
   {
      // Assert for same type
      assert(this->mIsScalar == req.isScalar());

      // Do OR operation on spectral requirement
      this->mNeedSpectral = this->mNeedSpectral || req.needSpectral();

      // Do OR operation on physical requirement
      this->mNeedPhysical = this->mNeedPhysical || req.needPhysical();

      // Do OR operation on physical gradient requirement
      this->mNeedGradient = this->mNeedGradient || req.needPhysicalGradient();

      // Do OR operation on physical curl requirement
      this->mNeedCurl = this->mNeedCurl || req.needPhysicalCurl();

      // Do OR operation on physical 2nd order gradient requirement
      this->mNeedGradient2 = this->mNeedGradient2 || req.needPhysicalGradient2();

      // Do OR operation of physical components requirement
      this->mPhysicalComps = this->mPhysicalComps || req.physicalComps();

      // Do OR operation of physical gradient components requirement
      if(req.needPhysicalGradient())
      {
         for(auto it = this->mSpectralIds.cbegin(); it != this->mSpectralIds.cend(); ++it)
         {
            this->mGradientComps.find(*it)->second = this->mGradientComps.find(*it)->second || req.gradientComps(*it);
         }
      }

      // Do OR operation of physical curl components requirement
      if(req.needPhysicalCurl())
      {
         this->mCurlComps = this->mCurlComps || req.curlComps();
      }

      // Do OR operation of physical 2nd order gradient components requirement
      if(req.needPhysicalGradient2())
      {
         for(auto it = this->mSpectralIds.cbegin(); it != this->mSpectralIds.cend(); ++it)
         {
            this->mGradient2Comps.find(*it)->second = this->mGradient2Comps.find(*it)->second || req.gradient2Comps(*it);
         }
      }
   }

}
