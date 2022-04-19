/**
 * @file SerialConverterBase.hpp
 * @brief Implementation of the serial data converter base
 */

#ifndef QUICC_PARALLEL_SERIALCONVERTERBASE_HPP
#define QUICC_PARALLEL_SERIALCONVERTERBASE_HPP

// Debug includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"

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
         template <typename T> void convertFwdImpl(const T& in, DynamicPairProvider &storage);

         /**
          * @brief Templated convert data from Bwd to Fwd
          */
         template <typename T> void convertBwdImpl(const T& in, DynamicPairProvider &storage);

         /**
          * @brief Get point data from Bwd scalar (might involve modification of indexes)
          *
          * @param in   Input data
          * @param i    First index of Bwd extracted from Fwd
          * @param j    Second index of Bwd extracted from Fwd
          * @param k    Third index of Bwd extracted from Fwd
          */
         template <typename T> typename T::PointType bwdPoint(const T& in, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Set point data from Bwd scalar (might involve modification of indexes)
          *
          * @param rOut Output data
          * @param i    First index of Bwd extracted from Fwd
          * @param j    Second index of Bwd extracted from Fwd
          * @param k    Third index of Bwd extracted from Fwd
          */
         template <typename T> typename T::PointType& rBwdPoint(T& rOut, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Store the local transform resolution of the first transform
          */
         SharedCTransformResolution mspTRes;
      private:
   };

   template <typename T>
      void SerialConverterBase::convertBwdImpl(const T& in, DynamicPairProvider& storage)
   {
      // Get storage for output value
      T *pOut;
      storage.provideFwd(pOut);

      // 3D case
      if(this->mDimensions == 3)
      {
         // Loop over slowest direction of output
         for(int k = 0; k < this->mspTRes->template dim<Dimensions::Data::DAT3D>(); k++)
         {
            // Loop over slow direction of output
            for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(k); j++)
            {
               // Loop over fast direction of output
               for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(k); i++)
               {
                  pOut->rPoint(i,j,k) = this->bwdPoint(in, i, j, k);
               }
            }
         }

      // 2D case
      } else if(this->mDimensions == 2)
      {
         // Loop over slow direction of output
         for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(); j++)
         {
            // Loop over fast direction of output
            for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(); i++)
            {
               pOut->rPoint(i,j) = this->bwdPoint(in, i, j);
            }
         }

      // 1D case
      } else if(this->mDimensions == 1)
      {
         //
         // No work is required in 1D
         //
         //
      } else
      {
         throw std::logic_error("Dimension of converter was not set properly!");
      }

      // Hold the output data
      storage.holdFwd(*pOut);
   }

   template <typename T>
      void SerialConverterBase::convertFwdImpl(const T& in, DynamicPairProvider& storage)
   {
      // Get storage for output value
      T *pOut;
      storage.provideBwd(pOut);

      // 3D case
      if(this->mDimensions == 3)
      {
         // Loop over slowest direction of input
         for(int k = 0; k < this->mspTRes->template dim<Dimensions::Data::DAT3D>(); k++)
         {
            // Loop over slow direction of input
            for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(k); j++)
            {
               // Loop over fast direction of input
               for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(k); i++)
               {
                  this->rBwdPoint(*pOut, i,j,k) = in.point(i,j,k);
               }
            }
         }

      // 2D case
      } else if(this->mDimensions == 2)
      {
         // Loop over slow direction of output
         for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(); j++)
         {
            for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(); i++)
            {
               this->rBwdPoint(*pOut, i,j) = in.point(i,j);
            }
         }

      // 1D case
      } else if(this->mDimensions == 1)
      {
         //
         // No work is required in 1D
         //
      } else
      {
         throw std::logic_error("Dimension of converter was not set properly!");
      }

      // Hold the output data
      storage.holdBwd(*pOut);
   }

   template <typename T>
      typename T::PointType SerialConverterBase::bwdPoint(const T& in, const int i, const int j, const int k)
   {
      #ifdef QUICC_MPI
         if(this->mDimensions == 3)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i,k);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j,k);
            int idxK = this->mspTRes->template idx<Dimensions::Data::DAT3D>(k);
            return in.point(this->idxConv().i(i,j,k,idxI,idxJ,idxK),this->idxConv().j(i,j,k,idxI,idxJ,idxK),this->idxConv().k(i,j,k,idxI,idxJ,idxK));

         } else if(this->mDimensions == 2)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j);
            return in.point(this->idxConv().i(i,j,idxI,idxJ),this->idxConv().j(i,j,idxI,idxJ));

         } else
         {
            return in.point(this->idxConv().i(i));
         }
      #else
         if(this->mDimensions == 3)
         {
            return in.point(this->idxConv().iS(i,j,k),this->idxConv().jS(i,j,k),this->idxConv().kS(i,j,k));

         } else if(this->mDimensions == 2)
         {
            return in.point(this->idxConv().iS(i,j),this->idxConv().jS(i,j));

         } else
         {
            return in.point(this->idxConv().i(i));
         }
      #endif //QUICC_MPI
   }

   template <typename T> typename T::PointType& SerialConverterBase::rBwdPoint(T& rOut, const int i, const int j, const int k)
   {
      #ifdef QUICC_MPI
         if(this->mDimensions == 3)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i,k);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j,k);
            int idxK = this->mspTRes->template idx<Dimensions::Data::DAT3D>(k);

            return rOut.rPoint(this->idxConv().i(i,j,k,idxI,idxJ,idxK),this->idxConv().j(i,j,k,idxI,idxJ,idxK),this->idxConv().k(i,j,k,idxI,idxJ,idxK));

         } else if(this->mDimensions == 2)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j);
            return rOut.rPoint(this->idxConv().i(i,j,idxI,idxJ),this->idxConv().j(i,j,idxI,idxJ));

         } else
         {
            return rOut.rPoint(this->idxConv().i(i));
         }
      #else
         if(this->mDimensions == 3)
         {
            return rOut.rPoint(this->idxConv().iS(i,j,k),this->idxConv().jS(i,j,k),this->idxConv().kS(i,j,k));

         } else if(this->mDimensions == 2)
         {
            return rOut.rPoint(this->idxConv().iS(i,j),this->idxConv().jS(i,j));

         } else
         {
            return rOut.rPoint(this->idxConv().i(i));
         }
      #endif //QUICC_MPI
   }

}
}

#endif // QUICC_PARALLEL_SERIALCONVERTERBASE_HPP
