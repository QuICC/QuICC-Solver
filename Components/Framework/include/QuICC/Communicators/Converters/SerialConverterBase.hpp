/**
 * @file SerialConverterBase.hpp
 * @brief Implementation of the serial data converter base
 */

#ifndef QUICC_PARALLEL_SERIALCONVERTERBASE_HPP
#define QUICC_PARALLEL_SERIALCONVERTERBASE_HPP

// Configuration includes
//

// System includes
//
#include <cassert>
#include <type_traits>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/StorageProviders/DynamicPairProvider.hpp"
#include "QuICC/Communicators/Converters/IConverter.hpp"

namespace QuICC {

namespace Parallel {

   namespace internal
   {
      /**
       * @brief Reorder 3D data
       */
      template<typename T, Dimensions::Data::Id TId> void reorder3D(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv);

      /**
       * @brief Reorder 2D data
       */
      template<typename T, Dimensions::Data::Id TId> void reorder2D(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv);

      /**
       * @brief Reorder 1D data
       */
      template<typename T, Dimensions::Data::Id TId> void reorder1D(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv);
   }

   /**
    * @brief Implementation of the serial data converter base.
    */
   class SerialConverterBase: public IConverter
   {
      public:
         /// Typedef for forward datatype
         typedef IConverter::RealFwdData RealFwdData;

         /// Typedef for forward datatype
         typedef IConverter::ComplexFwdData ComplexFwdData;

         /// Typedef for backward datatype
         typedef IConverter::RealBwdData RealBwdData;

         /// Typedef for backward datatype
         typedef IConverter::ComplexBwdData ComplexBwdData;

         /**
          * @brief Constructor
          */
         SerialConverterBase();

         /**
          * @brief Destructor
          */
         virtual ~SerialConverterBase();

         /**
          * @brief Convert real data from Fwd to Bwd
          */
         virtual void convertFwd(const RealFwdData &in, DynamicPairProvider &storage) override;

         /**
          * @brief Convert complex data from Fwd to Bwd
          */
         virtual void convertFwd(const ComplexFwdData &in, DynamicPairProvider &storage) override;

         /**
          * @brief Convert real data from Bwd to Fwd
          */
         virtual void convertBwd(const RealBwdData &in, DynamicPairProvider &storage) override;

         /**
          * @brief Convert complex data from Bwd to Fwd
          */
         virtual void convertBwd(const ComplexBwdData &in, DynamicPairProvider &storage) override;

      protected:
         /**
          * @brief Templated convert data from Fwd to Bwd
          */
         template <typename T, Dimensions::Data::Id TDataID> void processFwdImpl(const T& in, DynamicPairProvider &storage, const IIndexConv& outConv, const IIndexConv& inConv);

         /**
          * @brief Templated convert data from Fwd to Bwd
          */
         template <typename T, Dimensions::Data::Id TDataID> void processBwdImpl(const T& in, DynamicPairProvider &storage, const IIndexConv& outConv, const IIndexConv& inConv);

         /**
          * @brief Templated reorder data
          */
         template <typename T, Dimensions::Data::Id TDataId> void reorderImpl(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv);

         /**
          * @brief Store the local transform resolution of the first transform
          */
         SharedCTransformResolution mspTRes;
      private:
   };

   template <typename T, Dimensions::Data::Id TDataId>
      void SerialConverterBase::reorderImpl(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv)
   {
      // 3D case
      if(this->mDimensions == 3)
      {
         internal::reorder3D<T,TDataId>(pOut, in, tRes, outConv, inConv);
      }
      // 2D case
      else if(this->mDimensions == 2)
      {
         internal::reorder2D<T,TDataId>(pOut, in, tRes, outConv, inConv);
      }
      // 1D case
      else if(this->mDimensions == 1)
      {
         internal::reorder1D<T,TDataId>(pOut, in, tRes, outConv, inConv);
      }
      else
      {
         throw std::logic_error("Dimension of converter was not set properly!");
      }
   }

   template <typename T, Dimensions::Data::Id TDataId>
      void SerialConverterBase::processFwdImpl(const T& in, DynamicPairProvider& storage, const IIndexConv& outConv, const IIndexConv& inConv)
   {
      // Get storage for output value
      T *pOut;
      storage.provideFwd(pOut);

      this->reorderImpl<T,TDataId>(pOut, in, *this->mspTRes, outConv, inConv);

      // Hold the output data
      storage.holdFwd(*pOut);
   }

   template <typename T, Dimensions::Data::Id TDataId>
      void SerialConverterBase::processBwdImpl(const T& in, DynamicPairProvider& storage, const IIndexConv& outConv, const IIndexConv& inConv)
   {
      // Get storage for output value
      T *pOut;
      storage.provideBwd(pOut);

      this->reorderImpl<T,TDataId>(pOut, in, *this->mspTRes, outConv, inConv);

      // Hold the output data
      storage.holdBwd(*pOut);
   }

   namespace internal
   {
      template<typename T, Dimensions::Data::Id TId> void reorder3D(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv)
      {
         // Loop over slowest direction of output
         for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
         {
            int idxK = tRes.template idx<Dimensions::Data::DAT3D>(k);
            // Loop over slow direction of output
            for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int idxJ = tRes.template idx<Dimensions::Data::DAT2D>(j,k);
               // Loop over fast direction of output
               for(int i = 0; i < tRes.template dim<TId>(k); i++)
               {
                  int idxI = tRes.template idx<TId>(i,k);
                  pOut->rPoint(outConv.i(i,j,k,idxI,idxJ,idxK),outConv.j(i,j,k,idxI,idxJ,idxK),outConv.k(i,j,k,idxI,idxJ,idxK)) = in.point(inConv.i(i,j,k,idxI,idxJ,idxK),inConv.j(i,j,k,idxI,idxJ,idxK),inConv.k(i,j,k,idxI,idxJ,idxK));
               }
            }
         }
      }

      template<typename T, Dimensions::Data::Id TId> void reorder2D(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv)
      {
         // Loop over slow direction of output
         for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(); j++)
         {
            int idxJ = tRes.template idx<Dimensions::Data::DAT2D>(j);
            // Loop over fast direction of output
            for(int i = 0; i < tRes.template dim<TId>(); i++)
            {
               int idxI = tRes.template idx<TId>(i);
               pOut->rPoint(outConv.i(i,j,idxI,idxJ),outConv.j(i,j,idxI,idxJ)) = in.point(inConv.i(i,j,idxI,idxJ),inConv.j(i,j,idxI,idxJ));
            }
         }
      }

      template<typename T, Dimensions::Data::Id TId> void reorder1D(T* pOut, const T& in, const TransformResolution& tRes, const IIndexConv& outConv, const IIndexConv& inConv)
      {
         //
         // No work is required in 1D
         //
      }
   }
}
}

#endif // QUICC_PARALLEL_SERIALCONVERTERBASE_HPP
