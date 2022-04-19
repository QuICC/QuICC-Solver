/**
 * @file SlicedScalarField.hpp
 * @brief Base for a  scalar field with sliced data layout
 *
 *  \mhdBug Needs corrected implementation
 */

#ifndef QUICC_DATATYPES_SLICEDSCALARFIELD_HPP
#define QUICC_DATATYPES_SLICEDSCALARFIELD_HPP

// Debug includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//
#include <cassert>
#include <vector>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Base for a  scalar field with sliced data layout
    */
   template <typename TData> class SlicedScalarField
   {
      public:
         /// Typedef for the coefficient type
         typedef ScalarFieldSetup SetupType;

         /// Typedef for the coefficient type
         typedef TData PointType;

         /// Typedef for the storage type
         typedef typename std::vector<Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic> > StorageType;

         /// Typedef for the profiles of the storage type
         typedef typename Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic>::ColXpr  ProfileType;

         /// Typedef for the slices of the storage type
         typedef typename Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic>  SliceType;

         /**
          * @brief Constructor
          */
         explicit SlicedScalarField(std::shared_ptr<ScalarFieldSetup> spSetup);

         /**
          * @brief Destructor
          */
         ~SlicedScalarField();

         /**
          * @brief Get a point value of the field
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         PointType point(const int i, const int j = 0, const int k = 0) const;

         /**
          * @brief Get a 1D profile of the field
          *
          * A profile is defined as all value for fixed indexes in the second and third dimension.
          *
          * @param j Index of the profile
          * @param k Index of the slice
          */
         ProfileType profile(const int j, const int k = 0) const;

         /**
          * @brief Get a 2D slice of the field
          *
          * A slice is the matrix of values for a fixed index in the third dimension
          *
          * @param k Index of the slice
          */
         SliceType slice(const int k) const;

         /**
          * @brief Set a point value of the field
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         void setPoint(const PointType pt, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Set a profile of the field
          *
          * @param pf   Profile values
          * @param j    Index of the profile
          * @param k    Index of the slice
          */
         template <typename Derived> void setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k = 0);

         /**
          * @brief Set a 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void setSlice(const Eigen::MatrixBase<Derived>& sl, const int k);

         /**
          * @brief Get internal storage field data
          */
         const StorageType& data() const;

         /**
          * @brief Set internal storage field data
          */
         template <typename Derived> void setData(const Eigen::MatrixBase<Derived>& field);

         /**
          * @brief Set the complete field to zero
          */
         void setZeros();

         /**
          * @brief Rescale the complete field by a real coefficient
          *
          * @param scale Scaling factor
          */
         void rescale(const MHDFloat scale);

     #ifdef QUICC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE

      protected:

      private:
         /**
          * @brief Setup object for the scalar field
          */
         std::shared_ptr<ScalarFieldSetup> mspSetup;

         /**
          * @brief Field values
          */
         std::shared_ptr<StorageType>  mspField;
   };

   template <typename TData> inline typename SlicedScalarField<TData>::PointType SlicedScalarField<TData>::point(const int i, const int j, const int k) const
   {
      return this->mspField->at(k)(i,j);
   }

   template <typename TData> inline typename SlicedScalarField<TData>::ProfileType SlicedScalarField<TData>::profile(const int j, const int k) const
   {
      return this->mspField->at(k).col(j);
   }

   template <typename TData> inline typename SlicedScalarField<TData>::SliceType SlicedScalarField<TData>::slice(const int k) const
   {
      return this->mField->at(k);
   }

   template <typename TData> inline const typename SlicedScalarField<TData>::StorageType& SlicedScalarField<TData>::data() const
   {
      return *this->mspField;
   }

   template <typename TData> void SlicedScalarField<TData>::setPoint(const SlicedScalarField<TData>::PointType pt, const int i, const int j, const int k)
   {
      this->mspField->at(k)(i,j) = pt;
   }

   template <typename TData> template<typename Derived> void SlicedScalarField<TData>::setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      this->mspField->at(k).col(j) = pf;
   }

   template <typename TData> template<typename Derived> void SlicedScalarField<TData>::setSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      this->mspField->at(k) = sl;
   }

   template <typename TData> template<typename Derived> void SlicedScalarField<TData>::setData(const Eigen::MatrixBase<Derived>& field)
   {
      *this->mspField = field;
   }

   template <typename TData> SlicedScalarField<TData>::SlicedScalarField(std::shared_ptr<ScalarFieldSetup> spSetup)
      : mspSetup(spSetup)
   {
   }

   template <typename TData> SlicedScalarField<TData>::~SlicedScalarField()
   {
   }

   template <typename TData> void SlicedScalarField<TData>::setZeros()
   {
      // Iterate over all slices
      for(auto it = this->mspField->begin(); it != this->mspField->end(); ++it)
      {
         it->setConstant(0.0);
      }
   }

   template <typename TData> void SlicedScalarField<TData>::rescale(const MHDFloat scale)
   {
      // Iterate over all slices
      for(auto it = this->mspField->begin(); it != this->mspField->end(); ++it)
      {
         (*it) *= scale;
      }
   }

#ifdef QUICC_STORAGEPROFILE
   template <typename TData> MHDFloat SlicedScalarField<TData>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      for(auto it = this->mspField->cbegin(); it != this->mspField->cend(); ++it)
      {
         mem += it->size();
      }

      return static_cast<MHDFloat>(Debug::MemorySize<TData>::BYTES)*mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // QUICC_DATATYPES_SLICEDSCALARFIELD_HPP
