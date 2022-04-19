/**
 * @file PMIndexConv.hpp
 * @brief Implementation of the index converter for plus-minus FFT frequency ordering
 */

#ifndef QUICC_PARALLEL_PMINDEXCONV_HPP
#define QUICC_PARALLEL_PMINDEXCONV_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Communicators/Converters/IIndexConv.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the index converter for plus-minus FFT frequency ordering
    */
   class PMIndexConv: public IIndexConv
   {
      public:
         /**
          * @brief Constructor
          */
         PMIndexConv();

         /**
          * @brief Destructor
          */
         virtual ~PMIndexConv();

         /*
          * @brief Initialize index conersions
          *
          * @param res   resolution
          * @param id    Forward dimension index ID
          */
         virtual void init(const Resolution& res, const Dimensions::Transform::Id id) override;

         /**
          * @brief Compute shift due to central padding
          *
          * @param idx  Reference index
          * @param k    Array index of third dimension
          */
         int centralPadding(const int idx, const int k) const override;

         /**
          * @brief Convert first index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const override;

         /**
          * @brief Convert first index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         int iS(const int i, const int j, const int k) const override;

         /**
          * @brief Convert first index (2D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         int i(const int i, const int j, const int idxI, const int idxJ) const override;

         /**
          * @brief Convert first index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         int iS(const int i, const int j) const override;

         /**
          * @brief Convert first index (1D)
          *
          * @param i    Array index off first dimension
          */
         int i(const int i) const override;

         /**
          * @brief Convert second index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const override;

         /**
          * @brief Convert second index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         int jS(const int i, const int j, const int k) const override;

         /**
          * @brief Convert second index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         int j(const int i, const int j, const int idxI, const int idxJ) const override;

         /**
          * @brief Convert second index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         int jS(const int i, const int j) const override;

         /**
          * @brief Convert second index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const override;

         /**
          * @brief Convert third index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         int kS(const int i, const int j, const int k) const override;

      private:
         /**
          * @brief Shared transform resolution on forward side
          */
         SharedCTransformResolution mspTResFwd;

         /**
          * @brief Shared transform resolution backward side
          */
         SharedCTransformResolution mspTResBwd;

         /**
          * @brief Simulation resolution
          */
         int mSimN;
   };

   inline int PMIndexConv::centralPadding(const int idx, const int k) const
   {
      // Compute central padding if necessary
      int pad = 0;
      if(idx >= this->mSimN/2 + (this->mSimN % 2))
      {
         pad = this->mspTResBwd->dim<Dimensions::Data::DATB1D>(k) - this->mSimN;
      }

      return pad;
   }

   inline int PMIndexConv::i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const
   {
      return idxK + this->centralPadding(idxK, j);
   }

   inline int PMIndexConv::iS(const int i, const int j, const int k) const
   {
      return this->i(i,j,k,i,j,k);
   }

   inline int PMIndexConv::i(const int i, const int j, const int idxI, const int idxJ) const
   {
      return j;
   }

   inline int PMIndexConv::iS(const int i, const int j) const
   {
      return this->i(i,j,i,j);
   }

   inline int PMIndexConv::i(const int i) const
   {
      return i;
   }

   inline int PMIndexConv::j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const
   {
      return i;
   }

   inline int PMIndexConv::jS(const int i, const int j, const int k) const
   {
      return this->j(i,j,k,i,j,k);
   }

   inline int PMIndexConv::j(const int i, const int j, const int idxI, const int idxJ) const
   {
      return i;
   }

   inline int PMIndexConv::jS(const int i, const int j) const
   {
      return this->j(i,j,i,j);
   }

   inline int PMIndexConv::k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const
   {
      return j;
   }

   inline int PMIndexConv::kS(const int i, const int j, const int k) const
   {
      return this->k(i,j,k,i,j,k);
   }

}
}

#endif // QUICC_PARALLEL_PMINDEXCONV_HPP
