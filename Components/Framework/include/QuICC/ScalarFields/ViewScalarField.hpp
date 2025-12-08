/**
 * @file ViewScalarField.hpp
 * @brief Base for a  scalar field with View data layout
 */

#ifndef QUICC_DATATYPES_VIEWSCALARFIELD_HPP
#define QUICC_DATATYPES_VIEWSCALARFIELD_HPP

// System includes
//
#include <cassert>
#include <cstdint>
#include <vector>
#include <memory>

// Project includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Types/Typedefs.hpp"
#include "View/View.hpp"
#include "QuICC/ScalarFields/ScalarFieldSetup.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Base for a  scalar field with View data layout
    */
   template <typename TData> class ViewScalarField
   {
      public:
         /// Typedef for the coefficient type
         typedef ScalarFieldSetup SetupType;

         /// Typedef for the coefficient type
         typedef std::shared_ptr<ScalarFieldSetup> SharedSetupType;

         /// Typedef for the coefficient type
         typedef TData PointType;

         /// Typedef for the storage type
         using EigenMatrixMapType = Eigen::Map<Eigen::Matrix<PointType, Eigen::Dynamic, Eigen::Dynamic>>;

         /// Typedef for return type of Eigen profile
         using EigenVectorMapType = Eigen::Map<Eigen::Vector<PointType, Eigen::Dynamic>>;

         /// Typedef for the storage type
         using ViewStorageType = View::View<PointType, View::DCCSC3D>;

         /// Typedef for the profiles of the storage type
         using ViewProfileType = View::View<PointType, View::dense1D>;

         /// Typedef for the slices of the storage type
         using ViewSliceType = View::View<PointType, View::dense2D>;

         /**
          * @brief Constructor
          */
         explicit ViewScalarField(std::shared_ptr<ScalarFieldSetup> spSetup, std::shared_ptr<Memory::memory_resource> mem);

         /**
          * @brief Constructor
          */
         explicit ViewScalarField(std::shared_ptr<ScalarFieldSetup> spSetup);

         /**
          * @brief Copy constructor
          */
         ViewScalarField(const ViewScalarField<TData>& other);

         /**
          * @brief Destructor
          */
         ~ViewScalarField() = default;

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
         EigenVectorMapType profile(const int j, const int k = 0) const;

         /**
          * @brief Get a 2D slice of the field
          *
          * A slice is the matrix of values for a fixed index in the third dimension
          *
          * @param k Index of the slice
          */
         EigenMatrixMapType slice(const int k) const;

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
         EigenMatrixMapType data() const;

         /**
          * @brief Set internal storage field data
          */
         template <typename Derived> void setData(const Eigen::MatrixBase<Derived>& field);

         /**
          * @brief Set internal storage field data with flipped sign
          */
         template <typename Derived> void setNegData(const Eigen::MatrixBase<Derived>& field);

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
          * @brief Set the complete field to constant
          */
         void setConstant(const PointType c);

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
         EigenMatrixMapType rData();

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
         template <typename TType> const ViewScalarField<TData>& comp(const TType id) const;

         /**
          * @brief Set field component
          *
          * @param i Index of the component
          */
         template <typename TType> ViewScalarField<TData>& rComp(const TType id);

         /**
          * @brief Get a 1D profile of the field
          *
          * A profile is defined as all value for fixed indexes in the second and third dimension.
          *
          * @param j Index of the profile
          * @param k Index of the slice
          */
         ViewProfileType profileView(const int j, const int k = 0) const;

         /**
          * @brief Get a 2D slice of the field
          *
          * A slice is the matrix of values for a fixed index in the third dimension
          *
          * @param k Index of the slice
          */
         ViewSliceType sliceView(const int k) const;

         /**
          * @brief Set a profile of the field
          *
          * @param pf   Profile values
          * @param j    Index of the profile
          * @param k    Index of the slice
          */
         void setProfile(const ViewProfileType& pf, const int j, const int k = 0);

         /**
          * @brief Set a 2D slice of the field
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         void setSlice(const ViewSliceType& sl, const int k);

         /**
          * @brief Set the top rows of a 2D slice of the field
          *
          * Use this to adapt to differentenlty dealiased data
          *
          * @param sl   Slice values
          * @param k    Index of the slice
          */
         void setSlice(const ViewSliceType& sl, const int k, const int rows);

         /**
          * @brief Get full data view
          */
         const ViewStorageType& dataView() const;

         /**
          * @brief Set full data view
          */
         ViewStorageType& rDataView();

         /**
          * @brief Get full data global view
          */
         const ViewStorageType& globalView() const;

         /**
          * @brief Set full data global view
          */
         ViewStorageType& rGlobalView();

         /**
          * @brief Set full field data
          */
         void setData(const ViewStorageType& field);

      protected:

      private:
         /**
          * @brief Set metadata for local view (compressed)
          */
         void setLocalMetadata(std::shared_ptr<ScalarFieldSetup> spSetup);

         /**
          * @brief Set metadata for global view
          */
         void setGlobalMetadata(std::shared_ptr<ScalarFieldSetup> spSetup);

         /**
          * @brief Memory resources
          */
         std::shared_ptr<Memory::memory_resource> mMem;

         /**
          * @brief Data view
          */
         ViewStorageType mView;

         /**
          * @brief Global Data view
          */
         ViewStorageType mGlobalView;

         /**
          * @brief Data storage
          */
         std::shared_ptr<Memory::MemBlock<PointType>> mspData;

         /**
          * @brief Indices storage
          */
         std::shared_ptr<Memory::MemBlock<typename ViewStorageType::IndexType>> mspIndices;

         /**
          * @brief Pointer storage
          */
         std::shared_ptr<Memory::MemBlock<typename ViewStorageType::IndexType>> mspPointers;

         /**
          * @brief Global pointer storage
          */
         std::shared_ptr<Memory::MemBlock<typename ViewStorageType::IndexType>> mspGlobalPointers;

         /**
          * @brief Global indices storage
          */
         std::shared_ptr<Memory::MemBlock<typename ViewStorageType::IndexType>> mspGlobalIndices;
   };

   template <typename TData> inline typename ViewScalarField<TData>::PointType ViewScalarField<TData>::point(const int i, const int j, const int k) const
   {
      return this->mView(i,j,k);
   }

   template <typename TData> inline typename ViewScalarField<TData>::PointType ViewScalarField<TData>::point(const std::vector<int>& coord) const
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

   template <typename TData> inline TData* ViewScalarField<TData>::rData(const int k)
   {
      auto idx = this->mView.pointers()[1][k] * this->mView.lds();
      return this->mView.data() + idx;
   }

   template <typename TData> void ViewScalarField<TData>::setPoint(const ViewScalarField<TData>::PointType pt, const int i, const int j, const int k)
   {
      this->mView(i,j,k) = pt;
   }

   template <typename TData> void ViewScalarField<TData>::setPoint(const MHDVariant pt, const int i, const int j, const int k)
   {
      this->mView(i,j,k) = std::get<TData>(pt);
   }

   template <typename TData> inline typename ViewScalarField<TData>::PointType& ViewScalarField<TData>::rPoint(const int i, const int j, const int k)
   {
      auto idx = (this->mView.pointers()[1][k] + j) * this->mView.lds() + i;
      return *(this->mView.data() + idx);
   }

   template <typename TData> inline typename ViewScalarField<TData>::PointType& ViewScalarField<TData>::rPoint(const std::vector<int>& coord)
   {
      // Assert for positive sizes
      assert(this->mView.dims()[0] > 0);
      assert(this->mView.pointers()[1][this->mView.dims()[2]] > 0);

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

   template <typename TData> inline typename ViewScalarField<TData>::ViewProfileType ViewScalarField<TData>::profileView(const int j, const int k) const
   {
      std::uint32_t idx = (this->mView.pointers()[1][k] + j)*this->mView.lds();
      View::ViewBase<PointType> pView(this->mView.data() + idx, this->mView.lds());
      std::array<typename ViewProfileType::IndexType, 1> dims = {this->mView.dims()[0]};
      ViewProfileType profile(pView, dims);
      return profile;
   }

   template <typename TData> inline typename ViewScalarField<TData>::ViewSliceType ViewScalarField<TData>::sliceView(const int k) const
   {
      std::uint32_t idx = (this->mView.pointers()[1][k])*this->mView.lds();
      std::uint32_t cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      View::ViewBase<PointType> pView(this->mView.data() + idx, cols*this->mView.lds());
      std::array<typename ViewSliceType::IndexType, 2> dims = {this->mView.dims()[0], cols};
      ViewSliceType slice(pView, {dims});
      return slice;
   }

   template <typename TData> inline const TData* ViewScalarField<TData>::data(const int k) const
   {
      auto idx = this->mView.pointers()[1][k] * this->mView.lds();
      return this->mView.data() + idx;
   }

   template <typename TData> inline const typename ViewScalarField<TData>::ViewStorageType& ViewScalarField<TData>::dataView() const
   {
      return this->mView;
   }

   template <typename TData> inline typename ViewScalarField<TData>::ViewStorageType& ViewScalarField<TData>::rDataView()
   {
      return this->mView;
   }

   template <typename TData> inline const typename ViewScalarField<TData>::ViewStorageType& ViewScalarField<TData>::globalView() const
   {
      return this->mGlobalView;
   }

   template <typename TData> inline typename ViewScalarField<TData>::ViewStorageType& ViewScalarField<TData>::rGlobalView()
   {
      return this->mGlobalView;
   }

   template <typename TData> void ViewScalarField<TData>::setProfile(const ViewProfileType& pf, const int j, const int k)
   {
      assert(this->mView.dims()[0] == pf.dims()[0]);

      for(std::uint32_t i = 0; i < this->mView.dims()[0]; i++)
      {
         this->mView(i,j,k) = pf(i);
      }
   }

   template <typename TData> void ViewScalarField<TData>::setSlice(const ViewSliceType& sl, const int k)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(this->mView.dims()[0] == sl.dims()[0]);
      assert(cols == sl.dims()[1]);

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView(i,j,k) = sl(i, j);
         }
      }
   }

   template <typename TData> void ViewScalarField<TData>::setSlice(const ViewSliceType& sl, const int k, const int rows)
   {
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(static_cast<std::uint32_t>(rows) <= this->mView.dims()[0]);
      assert(static_cast<std::uint32_t>(rows) <= sl.dims()[0]);
      assert(cols == sl.dims()[1]);

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < static_cast<std::uint32_t>(rows); i++)
         {
            this->mView(i,j,k) = sl(i, j);
         }
      }
   }

   template <typename TData> void ViewScalarField<TData>::setData(const ViewStorageType& field)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][this->mView.dims()[2]];
      assert(field.dims()[0] == rows);
      assert(field.pointers()[1][this->mView.dims()[2]] == cols);

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView.data()[i + j*this->mView.lds()] = field.data()[i + j*field.lds()];
         }
      }
   }

   template <typename TData> inline typename ViewScalarField<TData>::EigenVectorMapType ViewScalarField<TData>::profile(const int j, const int k) const
   {
      std::uint32_t idx = (this->mView.pointers()[1][k] + j)*this->mView.lds();
      EigenVectorMapType mat(this->mView.data() + idx, this->mView.lds(), 1);
      return mat;
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::setProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      assert(this->mView.dims()[0] == pf.rows());

      for(std::uint32_t i = 0; i < this->mView.dims()[0]; i++)
      {
         this->mView(i,j,k) = pf(i);
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::addProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      assert(this->mView.dims()[0] == pf.rows());

      for(std::uint32_t i = 0; i < this->mView.dims()[0]; i++)
      {
         this->mView(i,j,k) += pf(i);
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::subProfile(const Eigen::MatrixBase<Derived>& pf, const int j, const int k)
   {
      assert(this->mView.dims()[0] == pf.rows());

      for(std::uint32_t i = 0; i < this->mView.dims()[0]; i++)
      {
         this->mView(i,j,k) -= pf(i);
      }
   }

   template <typename TData> inline typename ViewScalarField<TData>::EigenMatrixMapType ViewScalarField<TData>::slice(const int k) const
   {
      std::uint32_t idx = (this->mView.pointers()[1][k])*this->mView.lds();
      std::uint32_t cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      EigenMatrixMapType mat(this->mView.data() + idx, this->mView.lds(), cols);
      return mat;
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::setSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(this->mView.dims()[0] == sl.rows());
      assert(cols == sl.cols());

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView(i,j,k) = sl(i, j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::addSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(this->mView.dims()[0] == sl.rows());
      assert(cols == sl.cols());

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView(i,j,k) += sl(i, j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::subSlice(const Eigen::MatrixBase<Derived>& sl, const int k)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(this->mView.dims()[0] == sl.rows());
      assert(cols == sl.cols());

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView(i,j,k) -= sl(i, j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::setSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows)
   {
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(static_cast<std::uint32_t>(rows) <= this->mView.dims()[0]);
      assert(static_cast<std::uint32_t>(rows) <= sl.rows());
      assert(cols == sl.cols());

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < static_cast<std::uint32_t>(rows); i++)
         {
            this->mView(i,j,k) = sl(i, j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::addSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows)
   {
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(static_cast<std::uint32_t>(rows) <= this->mView.dims()[0]);
      assert(static_cast<std::uint32_t>(rows) <= sl.rows());
      assert(cols == sl.cols());

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < static_cast<std::uint32_t>(rows); i++)
         {
            this->mView(i,j,k) += sl(i, j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::subSlice(const Eigen::MatrixBase<Derived>& sl, const int k, const int rows)
   {
      auto cols = this->mView.pointers()[1][k+1] - this->mView.pointers()[1][k];
      assert(static_cast<std::uint32_t>(rows) <= this->mView.dims()[0]);
      assert(static_cast<std::uint32_t>(rows) <= sl.rows());
      assert(cols == sl.cols());

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < static_cast<std::uint32_t>(rows); i++)
         {
            this->mView(i,j,k) -= sl(i, j);
         }
      }
   }

   template <typename TData> inline typename ViewScalarField<TData>::EigenMatrixMapType ViewScalarField<TData>::data() const
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][this->mView.dims()[2]];
      EigenMatrixMapType  mat(this->mView.data(), rows, cols);
      return mat;
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::setData(const Eigen::MatrixBase<Derived>& field)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][this->mView.dims()[2]];
      assert(field.rows() == rows);
      assert(field.cols() == cols);

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView.data()[i + j*this->mView.lds()] = field(i,j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::setNegData(const Eigen::MatrixBase<Derived>& field)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][this->mView.dims()[2]];
      assert(field.rows() == rows);
      assert(field.cols() == cols);

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView.data()[i + j*this->mView.lds()] = -field(i,j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::addData(const Eigen::MatrixBase<Derived>& field)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][this->mView.dims()[2]];
      assert(field.rows() == rows);
      assert(field.cols() == cols);

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView.data()[i + j*this->mView.lds()] += field(i,j);
         }
      }
   }

   template <typename TData> template <typename Derived> void ViewScalarField<TData>::subData(const Eigen::MatrixBase<Derived>& field)
   {
      auto rows = this->mView.dims()[0];
      auto cols = this->mView.pointers()[1][this->mView.dims()[2]];
      assert(field.rows() == rows);
      assert(field.cols() == cols);

      for(std::uint32_t j = 0; j < cols; j++)
      {
         for(std::uint32_t i = 0; i < rows; i++)
         {
            this->mView.data()[i + j*this->mView.lds()] -= field(i,j);
         }
      }
   }

   template <typename TData> ViewScalarField<TData>::ViewScalarField(std::shared_ptr<ScalarFieldSetup> spSetup)
      : ViewScalarField(spSetup, spSetup->mem())
   {
   }

   template <typename TData> ViewScalarField<TData>::ViewScalarField(std::shared_ptr<ScalarFieldSetup> spSetup, std::shared_ptr<Memory::memory_resource> mem)
      : mMem(mem)
   {
      this->setLocalMetadata(spSetup);
      this->setGlobalMetadata(spSetup);
   }

   template <typename TData> ViewScalarField<TData>::ViewScalarField(const ViewScalarField<TData>& other)
   {
      this->mMem = other.mMem;
      this->mView = other.mView;
      this->mspData  = other.mspData;
      this->mspIndices = other.mspIndices;
      this->mspPointers = other.mspPointers;
      this->mGlobalView = other.mGlobalView;
      this->mspGlobalPointers = other.mspGlobalPointers;
   }

   template <typename TData> void ViewScalarField<TData>::setZeros()
   {
      auto ptr = this->mView.data();

      for(std::uint32_t i  = 0; i < this->mView.size(); i++)
      {
         *ptr = 0;
         ptr++;
      }
   }

   template <typename TData> void ViewScalarField<TData>::setConstant(const PointType c)
   {
      auto ptr = this->mView.data();

      for(std::uint32_t i  = 0; i < this->mView.size(); i++)
      {
         *ptr = c;
         ptr++;
      }
   }

   template <typename TData> void ViewScalarField<TData>::rescale(const MHDFloat scale)
   {
      auto ptr = this->mView.data();

      for(std::uint32_t i  = 0; i < this->mView.size(); i++)
      {
         *ptr *= scale;
         ptr++;
      }
   }

   template <typename TData> inline typename ViewScalarField<TData>::EigenMatrixMapType ViewScalarField<TData>::rData()
   {
      std::uint32_t cols = this->mView.pointers()[1][this->mView.dims()[2]];
      EigenMatrixMapType mat(this->mView.data(), this->mView.lds(), cols);
      return mat;
   }

   template <typename TData> int ViewScalarField<TData>::nSlice() const
   {
      return this->mView.pointers()[1].size() - 1;
   }

   template <typename TData> template <typename TType> inline const ViewScalarField<TData>& ViewScalarField<TData>::comp(const TType id) const
   {
      assert(TType::SCALAR == id);

      return *this;
   }

   template <typename TData> template <typename TType> inline ViewScalarField<TData>& ViewScalarField<TData>::rComp(const TType id)
   {
      assert(TType::SCALAR == id);

      return *this;
   }

   template <typename TData> void ViewScalarField<TData>::setLocalMetadata(std::shared_ptr<ScalarFieldSetup> spSetup)
   {
      auto meta = *spSetup->viewMeta();

      // Get dimensions
      std::uint32_t n1D = static_cast<std::uint32_t>(spSetup->dataRows());
      std::uint32_t n2D = 0;
      for(int k = 0; k < spSetup->nBlock(); k++)
      {
         n2D = std::max(n2D, static_cast<std::uint32_t>(spSetup->blockCols(k)));
      }
      std::uint32_t n3D = static_cast<std::uint32_t>(spSetup->nBlock());

      std::array<std::uint32_t, 3> dimensions{n1D, n2D, n3D};
      std::uint32_t dataSize = n1D*static_cast<std::uint32_t>(spSetup->dataCols());

      // Alloc storage
      this->mspData = std::make_shared<Memory::MemBlock<PointType>>(dataSize, this->mMem.get());
      this->mspPointers = std::make_shared<Memory::MemBlock<typename ViewStorageType::IndexType>>(n3D + 1, this->mMem.get());
      this->mspIndices = std::make_shared<Memory::MemBlock<typename ViewStorageType::IndexType>>(meta.idx2D.size(), this->mMem.get());

      // Set view
      View::ViewBase<typename ViewStorageType::IndexType> pointers[this->mView.rank()];
      View::ViewBase<typename ViewStorageType::IndexType> indices[this->mView.rank()];
      pointers[1] =
         View::ViewBase<typename ViewStorageType::IndexType>(this->mspPointers->data(), this->mspPointers->size());
      indices[1] =
         View::ViewBase<typename ViewStorageType::IndexType>(this->mspIndices->data(), this->mspIndices->size());
      this->mView = ViewStorageType(this->mspData->data(), this->mspData->size(), dimensions.data(), pointers, indices);

      // Set pointers and indices
      std::uint32_t ii = 1;
      pointers[1][0] = 0;
      for(std::uint32_t i = 1; i < meta.ptr2D.size(); i++)
      {
         if(meta.ptr2D.at(i) > meta.ptr2D.at(i-1))
         {
            assert(ii < pointers[1].size());
            pointers[1][ii] = meta.ptr2D.at(i);
            ii++;
         }
      }
      ii = 0;
      for(std::uint32_t i = 0; i < pointers[1].size() - 1; i++)
      {
         for(std::uint32_t j = 0; j < pointers[1][i+1] - pointers[1][i]; j++)
         {
            indices[1][ii] = j;
            ii++;
         }
      }
   }

   template <typename TData> void ViewScalarField<TData>::setGlobalMetadata(std::shared_ptr<ScalarFieldSetup> spSetup)
   {
      auto meta = *spSetup->viewMeta();

      // Get dimensions
      std::uint32_t n1D = meta.global1D;
      std::uint32_t n2D = meta.global2D;
      std::uint32_t n3D = meta.global3D;

      std::array<std::uint32_t, 3> dimensions{n1D, n2D, n3D};
      std::uint32_t dataSize = static_cast<std::uint32_t>(spSetup->dataRows())*static_cast<std::uint32_t>(spSetup->dataCols());

      // Alloc storage
      assert(this->mspData->size() == dataSize);
      this->mspGlobalPointers = std::make_shared<Memory::MemBlock<typename ViewStorageType::IndexType>>(n3D + 1, this->mMem.get());
      this->mspGlobalIndices = std::make_shared<Memory::MemBlock<typename ViewStorageType::IndexType>>(meta.idx2D.size(), this->mMem.get());

      // Set view
      View::ViewBase<typename ViewStorageType::IndexType> pointers[this->mView.rank()];
      View::ViewBase<typename ViewStorageType::IndexType> indices[this->mView.rank()];
      pointers[1] =
         View::ViewBase<typename ViewStorageType::IndexType>(this->mspGlobalPointers->data(), this->mspGlobalPointers->size());
      indices[1] =
         View::ViewBase<typename ViewStorageType::IndexType>(this->mspGlobalIndices->data(), this->mspGlobalIndices->size());
      this->mGlobalView = ViewStorageType(this->mspData->data(), this->mspData->size(), dimensions.data(), pointers, indices);

      // Set pointers and indices
      assert(meta.ptr2D.size() == pointers[1].size());
      for(std::uint32_t i = 0; i < meta.ptr2D.size(); i++)
      {
         pointers[1][i] = meta.ptr2D.at(i);
      }
      std::copy(meta.idx2D.begin(), meta.idx2D.end(), indices[1].data());
   }

   template <typename TData> MHDFloat ViewScalarField<TData>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<TData>::BYTES*this->mspField->size());
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}

#endif // QUICC_DATATYPES_VIEWSCALARFIELD_HPP
