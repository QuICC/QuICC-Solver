/**
 * @file VariableDomain.hpp
 * @brief Implementation of the variable's domain
 */

#ifndef QUICC_DATATYPES_VARIABLEDOMAIN_HPP
#define QUICC_DATATYPES_VARIABLEDOMAIN_HPP

// Configuration includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation for a full sphere field
    *
    * \tparam TVariable Type of the variable
    * \tparam DOMAINS   Number of different domains
    */
   template <typename TVariable, int DOMAINS> class VariableDomain
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spRes Resolution information
          */
         VariableDomain(SharedResolution spRes);

         /**
          * @brief Destructor
          */
         virtual ~VariableDomain();

         /**
          * @brief Get variable of requested domain
          *
          * @param i Index of the domain
          */
         const TVariable&  dom(const int i) const;

         /**
          * @brief Set Physical variable in full sphere
          *
          * @param i Index of the domain
          */
         TVariable&  rDom(const int i);

         /**
          * @brief initialise to zeros
          */
         void setZeros();

         /**
          * @brief Initialise the spectral values storage
          */
         void initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps);

         /**
          * @brief Initialise the physical values storage
          */
         void initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical gradient values storage
          */
         void initPhysicalGradient(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical curl values storage
          */
         void initPhysicalCurl(const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical 2nd order gradient values storage
          */
         void initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& comps);

         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;

         /**
         * @brief Profile the memory requirements
         */
         void profileStorage() const;

      protected:
         /**
          * @brief Storage for the OC (or full sphere) part
          */
         std::vector<TVariable>   mDomains;

      private:
   };

   template <typename TVariable, int DOMAINS> inline const TVariable& VariableDomain<TVariable,DOMAINS>::dom(const int i) const
   {
      // Assert for valid index
      assert(this->mDomains.size() > static_cast<size_t>(i));

      return this->mDomains.at(i);
   }

   template <typename TVariable, int DOMAINS> inline TVariable& VariableDomain<TVariable,DOMAINS>::rDom(const int i)
   {
      // Assert for valid index
      assert(i >= 0);
      assert(this->mDomains.size() > static_cast<size_t>(i));

      return this->mDomains.at(i);
   }

   template <typename TVariable, int DOMAINS> VariableDomain<TVariable,DOMAINS>::VariableDomain(SharedResolution spRes)
   {
      for(int i = 0; i < DOMAINS; i++)
      {
         this->mDomains.push_back(TVariable(spRes));
      }
   }

   template <typename TVariable, int DOMAINS> VariableDomain<TVariable,DOMAINS>::~VariableDomain()
   {
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::setZeros()
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).setZeros();
      }
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps)
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).initSpectral(comps);
      }
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).initPhysical(comps);
      }
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::initPhysicalGradient(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).initPhysicalGradient(id, comps);
      }
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::initPhysicalCurl(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).initPhysicalCurl(comps);
      }
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& comps)
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).initPhysicalGradient2(id, comps);
      }
   }

   template <typename TVariable, int DOMAINS> MHDFloat  VariableDomain<TVariable,DOMAINS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         mem += this->mDomains.at(i).requiredStorage();
      }
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::profileStorage() const
   {
#ifdef QUICC_STORAGEPROFILE
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).profileStorage();
      }
#endif // QUICC_STORAGEPROFILE
   }

}
}

#endif // QUICC_DATATYPES_VARIABLEDOMAIN_HPP
