/**
 * @file SHlIndexConv.hpp
 * @brief Implementation of the index converter for spherical harmonics with l spectral ordering
 */

#ifndef QUICC_PARALLEL_SHLINDEXCONV_HPP
#define QUICC_PARALLEL_SHLINDEXCONV_HPP

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
    * @brief Implementation of the index converter for spherical harmonics with l spectral ordering
    */
   class SHlIndexConv: public IIndexConv
   {
      public:
         /**
          * @brief Constructor
          */
         SHlIndexConv();

         /**
          * @brief Destructor
          */
         virtual ~SHlIndexConv();

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
          * @brief Convert second index (2D)
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
   };

   inline int SHlIndexConv::centralPadding(const int, const int) const
   {
      return 0;
   }

   inline int SHlIndexConv::i(const int, const int, const int, const int, const int idxJ, const int idxK) const
   {
      return idxK - idxJ;
   }

   inline int SHlIndexConv::iS(const int i, const int j, const int k) const
   {
      return this->i(i,j,k,i,j,k);
   }

   inline int SHlIndexConv::i(const int, const int, const int idxI, const int idxJ) const
   {
      return idxJ - idxI;
   }

   inline int SHlIndexConv::iS(const int i, const int j) const
   {
      return this->i(i,j,i,j);
   }

   inline int SHlIndexConv::i(const int i) const
   {
      return i;
   }

   inline int SHlIndexConv::j(const int i, const int, const int, const int, const int, const int) const
   {
      return i;
   }

   inline int SHlIndexConv::jS(const int i, const int j, const int k) const
   {
      return this->j(i,j,k,i,j,k);
   }

   inline int SHlIndexConv::j(const int i, const int, const int, const int) const
   {
      return i;
   }

   inline int SHlIndexConv::jS(const int i, const int j) const
   {
      return this->j(i,j,i,j);
   }

   inline int SHlIndexConv::k(const int, const int j, const int, const int, const int, const int) const
   {
      return j;
   }

   inline int SHlIndexConv::kS(const int i, const int j, const int k) const
   {
      return this->k(i,j,k,i,j,k);
   }

}
}

#endif // QUICC_PARALLEL_SHLINDEXCONV_HPP
