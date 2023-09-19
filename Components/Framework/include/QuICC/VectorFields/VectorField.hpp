/**
 * @file VectorField.hpp
 * @brief Implementation of a generic vector field
 */

#ifndef QUICC_DATATYPES_VECTORFIELD_HPP
#define QUICC_DATATYPES_VECTORFIELD_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation of a generic vector field
    */
   template <typename TScalar, typename TType> class VectorField
   {
      public:
         /**
          * @brief Constructor with simplified interface
          *
          * @param spSetup Shared setup object for the scalar fields
          */
         VectorField(typename TScalar::SharedSetupType spSetup, const typename std::map<TType,bool>& comps);

         /**
          * @brief Destructor
          */
         virtual ~VectorField();

         /**
          * @brief Get field component
          *
          * @param i Index of the component
          */
         const TScalar& comp(const TType id) const;

         /**
          * @brief Set field component
          *
          * @param i Index of the component
          */
         TScalar& rComp(const TType id);

         /**
          * @brief Set field components to zero
          */
         void setZeros();

         /**
          * @brief Get field components
          */
         const std::map<TType,TScalar>& data() const;

         /**
          * @brief Get map of enabled components
          */
         const std::map<TType,bool>& enabled() const;

         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;

         /**
          * @brief Set internal storage field data
          *
          * \warning This routine should only be used in exceptional cases. Use setData, addData, subData when you can!
         */
         std::map<TType,TScalar>& rData();

      protected:
         /**
          * @brief Storage for the field components
          */
         std::map<TType,TScalar> mComponents;

      private:
         /**
          * @brief Map of enabled components
          */
         typename std::map<TType,bool> mEnabled;
   };

   template <typename TScalar, typename TType> inline const TScalar& VectorField<TScalar,TType>::comp(const TType id) const
   {
      // Assert that index is valid
      assert(this->mComponents.count(id) == 1);

      return this->mComponents.find(id)->second;
   }

   template <typename TScalar, typename TType> inline TScalar& VectorField<TScalar,TType>::rComp(const TType id)
   {
      // Assert that index is valid
      assert(this->mComponents.count(id) == 1);

      return this->mComponents.find(id)->second;
   }

   template <typename TScalar, typename TType> inline const std::map<TType,bool>& VectorField<TScalar,TType>::enabled() const
   {
      return this->mEnabled;
   }

   template <typename TScalar, typename TType> inline const std::map<TType,TScalar>& VectorField<TScalar,TType>::data() const
   {
      return this->mComponents;
   }

   template <typename TScalar, typename TType> inline std::map<TType,TScalar>& VectorField<TScalar,TType>::rData()
   {
      return this->mComponents;
   }

   template <typename TScalar, typename TType> VectorField<TScalar,TType>::VectorField(typename TScalar::SharedSetupType spSetup, const typename std::map<TType,bool>& comps)
   {
      // Initialise the components
      for(auto it = comps.cbegin(); it != comps.cend(); ++it)
      {
         if(it->second)
         {
            this->mComponents.insert(std::make_pair(it->first, TScalar(spSetup)));
         }
      }

      // Store map of enabled components
      this->mEnabled = comps;
   }

   template <typename TScalar, typename TType> VectorField<TScalar, TType>::~VectorField()
   {
   }

   template <typename TScalar, typename TType> void VectorField<TScalar,TType>::setZeros()
   {
      // Initialise the components
      for(auto it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         it->second.setZeros();
      }
   }

   template <typename TScalar, typename TType> MHDFloat VectorField<TScalar,TType>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      for(auto it = this->mComponents.cbegin(); it != this->mComponents.cend(); it++)
      {
         mem += it->second.requiredStorage();
      }
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}

#endif // QUICC_DATATYPES_VECTORFIELD_HPP
