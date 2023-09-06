/**
 * @file FlatScalarField.hpp
 * @brief Base for a  scalar field with flat data layout
 */

#ifndef QUICC_DATATYPES_FLATSCALARFIELD_HPP
#define QUICC_DATATYPES_FLATSCALARFIELD_HPP

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
    * @brief Base for a  scalar field with flat data layout
    */
   template <typename TData> class FlatScalarField
   {
      public:
         /// Typedef for the coefficient type
         typedef ScalarFieldSetup SetupType;

         /// Typedef for the coefficient type
         typedef std::shared_ptr<ScalarFieldSetup> SharedSetupType;

         /// Typedef for the coefficient type
         typedef TData PointType;

         /// Typedef for the storage type
         typedef typename Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic> StorageType;

         /// Typedef for the profiles of the storage type
         typedef typename StorageType::ColXpr  ProfileType;

         /// Typedef for the slices of the storage type
         typedef typename Eigen::Block<StorageType>  SliceType;

         /**
          * @brief Constructor
          */
         explicit FlatScalarField(std::shared_ptr<ScalarFieldSetup> spSetup);

         /**
          * @brief Copy constructor
          */
         FlatScalarField(const FlatScalarField<TData>& other);

         /**
          * @brief Destructor
          */
         ~FlatScalarField();

         /**
          * @brief Get a point value of the field
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         PointType point(const int i, const int j = 0, const int k = 0) const;

         /**
          * @brief Get a point value of the field
          *
          * @param coord Generic coordinate
          */
         PointType point(const std::vector<int>& coord) const;

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
          * @brief Set a point value of the field from variant
          *
          * @param i Index of the point
          * @param j Index of the profile
          * @param k Index of the slice
          */
         void setPoint(const MHDVariant pt, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Set a profile of the field
          *
          * @param pf   Profile values
          * @param j    Index of the profile
          * @param k    Index of the slice
          */
         template <typename Derived> void setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k = 0);

         /**
          * @brief Add a profile of the field
          *
          * @param pf   Profile values
          * @param j    Index of the profile
          * @param k    Index of the slice
          */
         template <typename Derived> void addProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k = 0);

         /**
          * @brief Substract a profile of the field
          *
          * @param pf   Profile values
          * @param j    Index of the profile
          * @param k    Index of the slice
          */
         template <typename Derived> void subProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k = 0);

         /**
          * @brief Set a 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void setSlice(const Eigen::MatrixBase<Derived>& sl, const int k);

         /**
          * @brief Multiply a 2D slice of the field by a matrix
          *
          * @param mat  Left matrix multiplication
          * @param k    Index of the slice
          */
         template <typename Derived> void multSlice(const Eigen::SparseMatrixBase<Derived>& mat, const int k);

         /**
          * @brief Add to 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void addSlice(const Eigen::MatrixBase<Derived>& sl, const int k);

         /**
          * @brief Substract from 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void subSlice(const Eigen::MatrixBase<Derived>& sl, const int k);

         /**
          * @brief Set the top rows of a 2D slice of the field
          *
          * Use this to adapt to differentenlty dealiased data
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void setSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows);

         /**
          * @brief Multiply the top rows of a 2D slice of the field by a matrix
          *
          * Use this to adapt to differentenlty dealiased data
          *
          * @param mat  Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void multSlice(const Eigen::SparseMatrixBase<Derived>& mat, const int k, const int rows);

         /**
          * @brief Add to the top rows of a 2D slice of the field
          *
          * Use this to adapt to differentenlty dealiased data
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void addSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows);

         /**
          * @brief Substract from the top rows of a 2D slice of the field
          *
          * Use this to adapt to differentenlty dealiased data
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         template <typename Derived> void subSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows);

         /**
          * @brief Get internal storage field data pointer
          */
         const TData* data(const int k) const;

         /**
          * @brief Get internal storage field data
          */
         const StorageType& data() const;

         /**
          * @brief Set internal storage field data
          */
         template <typename Derived> void setData(const Eigen::MatrixBase<Derived>& field);

         /**
          * @brief Add to internal storage field data
          */
         template <typename Derived> void addData(const Eigen::MatrixBase<Derived>& field);

         /**
          * @brief Substract from internal storage field data
          */
         template <typename Derived> void subData(const Eigen::MatrixBase<Derived>& field);

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

         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;

         /**
          * @brief Get the number of slices
          */
         int nSlice() const;

         /**
          * @brief Set internal storage field data pointer
          *
          * \warning This routine should only be used in exceptional cases. Use setData, addData, subData when you can!
          */
         TData* rData(const int k);

         /**
          * @brief Set internal storage field data
          *
          * \warning This routine should only be used in exceptional cases. Use setData, addData, subData when you can!
          */
         StorageType& rData();

         /**
          * @brief Set internal storage point data
          *
          * \warning This routine should only be used in exceptional cases. Use setPoint!
          */
         PointType& rPoint(const int i, const int j = 0, const int k = 0);

         /**
          * @brief Set internal storage point data
          *
          * @param coord Generic coordinate
          *
          * \warning This routine should only be used in exceptional cases. Use setPoint!
          */
         PointType& rPoint(const std::vector<int>& coord);

         /**
          * @brief To unify interface with vector field return self as comp
          *
          * @param i Index of the component
          */
         template <typename TType> const FlatScalarField<TData>& comp(const TType id) const;

         /**
          * @brief Set field component
          *
          * @param i Index of the component
          */
         template <typename TType> FlatScalarField<TData>& rComp(const TType id);

      protected:

      private:
         /**
          * @brief Setup object for the scalar field
          */
         std::shared_ptr<ScalarFieldSetup> mspSetup;

         /**
          * @brief Field values in shared flat storage
          */
         std::shared_ptr<StorageType>  mspField;
   };

   template <typename TData> inline typename FlatScalarField<TData>::PointType FlatScalarField<TData>::point(const int i, const int j, const int k) const
   {
      return (*this->mspField)(i,this->mspSetup->colIdx(j,k));
   }

   template <typename TData> inline typename FlatScalarField<TData>::PointType FlatScalarField<TData>::point(const std::vector<int>& coord) const
   {
      switch(coord.size())
      {
         case 3:
            return this->point(coord.at(0), coord.at(1), coord.at(2));
         case 2:
            return this->point(coord.at(0), coord.at(1));
         case 1:
            return this->point(coord.at(0));
         default:
            throw std::logic_error("Coordinate for point has more than 3 dimensions");
      }
   }

   template <typename TData> inline typename FlatScalarField<TData>::ProfileType FlatScalarField<TData>::profile(const int j, const int k) const
   {
      return this->mspField->col(this->mspSetup->colIdx(j,k));
   }

   template <typename TData> inline typename FlatScalarField<TData>::SliceType FlatScalarField<TData>::slice(const int k) const
   {
      return this->mspField->block(0,this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k));
   }

   template <typename TData> inline const TData* FlatScalarField<TData>::data(const int k) const
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      return this->mspField->data() + this->mspSetup->dataRows()*this->mspSetup->blockIdx(k);
   }

   template <typename TData> inline const typename FlatScalarField<TData>::StorageType& FlatScalarField<TData>::data() const
   {
      // Assert for positive sizes
      assert(this->mspField);

      return *this->mspField;
   }

   template <typename TData> inline typename FlatScalarField<TData>::StorageType& FlatScalarField<TData>::rData()
   {
      // Assert for positive sizes
      assert(this->mspField);

      return *this->mspField;
   }

   template <typename TData> inline TData* FlatScalarField<TData>::rData(const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      return this->mspField->data() + this->mspSetup->dataRows()*this->mspSetup->blockIdx(k);
   }

   template <typename TData> void FlatScalarField<TData>::setPoint(const FlatScalarField<TData>::PointType pt, const int i, const int j, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      (*this->mspField)(i,this->mspSetup->colIdx(j,k)) = pt;
   }

   template <typename TData> void FlatScalarField<TData>::setPoint(const MHDVariant pt, const int i, const int j, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      (*this->mspField)(i,this->mspSetup->colIdx(j,k)) = std::get<TData>(pt);
   }

   template <typename TData> inline typename FlatScalarField<TData>::PointType& FlatScalarField<TData>::rPoint(const int i, const int j, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      return (*this->mspField)(i,this->mspSetup->colIdx(j,k));
   }

   template <typename TData> inline typename FlatScalarField<TData>::PointType& FlatScalarField<TData>::rPoint(const std::vector<int>& coord)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      switch(coord.size())
      {
         case 3:
            return this->rPoint(coord.at(0), coord.at(1), coord.at(2));
         case 2:
            return this->rPoint(coord.at(0), coord.at(1));
         case 1:
            return this->rPoint(coord.at(0));
         default:
            throw std::logic_error("Coordinate for point has more than 3 dimensions");
      }
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      this->mspField->col(this->mspSetup->colIdx(j,k)) = pf;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::addProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      this->mspField->col(this->mspSetup->colIdx(j,k)) += pf;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::subProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      this->mspField->col(this->mspSetup->colIdx(j,k)) -= pf;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::setSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      this->mspField->block(0, this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k)) = sl;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::multSlice(const Eigen::SparseMatrixBase<Derived>& mat, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      // Assert for compatible matrix
      assert(mat.rows() == mat.cols());
      assert(mat.cols() == this->mspSetup->blockRows(k));

      this->mspField->block(0, this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k)) *= mat;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::addSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      this->mspField->block(0, this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k)) += sl;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::subSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      this->mspField->block(0, this->mspSetup->blockIdx(k), this->mspSetup->blockRows(k), this->mspSetup->blockCols(k)) -= sl;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::setSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      // Should not be called as a replacement of the shorter version
      assert(rows < this->mspSetup->blockRows(k));

      this->mspField->block(0, this->mspSetup->blockIdx(k), rows, this->mspSetup->blockCols(k)) = sl;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::multSlice(const Eigen::SparseMatrixBase<Derived>& mat, const int k, const int rows)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      // Assert for compatible matrix
      assert(mat.rows() == mat.cols());
      assert(mat.cols() == rows);

      // Should not be called as a replacement of the shorter version
      assert(rows < this->mspSetup->blockRows(k));

      this->mspField->block(0, this->mspSetup->blockIdx(k), rows, this->mspSetup->blockCols(k)) *= mat;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::addSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      // Should not be called as a replacement of the shorter version
      assert(rows < this->mspSetup->blockRows(k));

      this->mspField->block(0, this->mspSetup->blockIdx(k), rows, this->mspSetup->blockCols(k)) += sl;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::subSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      // Should not be called as a replacement of the shorter version
      assert(rows < this->mspSetup->blockRows(k));

      this->mspField->block(0, this->mspSetup->blockIdx(k), rows, this->mspSetup->blockCols(k)) -= sl;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::setData(const Eigen::MatrixBase<Derived>& field)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() == field.rows());
      assert(this->mspField->cols() == field.cols());

      // Assert for positive sizes
      this->mspField->block(0, 0, this->mspSetup->dataRows(), this->mspSetup->dataCols()) = field;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::addData(const Eigen::MatrixBase<Derived>& field)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() == field.rows());
      assert(this->mspField->cols() == field.cols());

      this->mspField->block(0, 0, this->mspSetup->dataRows(), this->mspSetup->dataCols()) += field;
   }

   template <typename TData> template<typename Derived> void FlatScalarField<TData>::subData(const Eigen::MatrixBase<Derived>& field)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() == field.rows());
      assert(this->mspField->cols() == field.cols());

      this->mspField->block(0, 0, this->mspSetup->dataRows(), this->mspSetup->dataCols()) -= field;
   }

   template <typename TData> FlatScalarField<TData>::FlatScalarField(std::shared_ptr<ScalarFieldSetup> spSetup)
      : mspSetup(spSetup), mspField(new StorageType(spSetup->dataRows(), spSetup->dataCols()))
   {
   }

   template <typename TData> FlatScalarField<TData>::FlatScalarField(const FlatScalarField<TData>& other)
      : mspSetup(other.mspSetup), mspField(new StorageType(other.mspField->rows(),other.mspField->cols()))
   {
      *this->mspField = *(other.mspField);
   }

   template <typename TData> FlatScalarField<TData>::~FlatScalarField()
   {
   }

   template <typename TData> void FlatScalarField<TData>::setZeros()
   {
      this->mspField->setConstant(0.0);
   }

   template <typename TData> void FlatScalarField<TData>::rescale(const MHDFloat scale)
   {
      // Assert for positive sizes
      assert(this->mspField->rows() > 0);
      assert(this->mspField->cols() > 0);

      *this->mspField *= scale;
   }

   template <typename TData> int FlatScalarField<TData>::nSlice() const
   {
      return this->mspSetup->nBlock();
   }

   template <typename TData> template <typename TType> inline const FlatScalarField<TData>& FlatScalarField<TData>::comp(const TType id) const
   {
      assert(TType::SCALAR == id);

      return *this;
   }

   template <typename TData> template <typename TType> inline FlatScalarField<TData>& FlatScalarField<TData>::rComp(const TType id)
   {
      assert(TType::SCALAR == id);

      return *this;
   }

   template <typename TData> MHDFloat FlatScalarField<TData>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<TData>::BYTES*this->mspField->size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}

#endif // QUICC_DATATYPES_FLATSCALARFIELD_HPP
