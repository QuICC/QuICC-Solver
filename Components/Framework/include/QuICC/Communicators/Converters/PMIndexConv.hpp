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
    *
    * For complex FFT, because of the ordering of positive and negative frequency, padding happens in the middle of the array. This class computes the mapping from (IJK) to (JKI) for complex Fourier modes.
    */
   class PMIndexConv: public IIndexConv
   {
      public:
         /**
          * @brief Constructor
          */
         PMIndexConv() = default;

         /**
          * @brief Destructor
          */
         virtual ~PMIndexConv() = default;

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
          * Takes data index (i,j,k) and physical mode indexes (idxI, idxJ, idxK) from input data type, to compute corresponding data index i for output data
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
          * @brief Convert first index (2D)
          *
          * Takes data index (i,j) and physical mode indexes (idxI, idxJ) from input data type, to compute corresponding data index i for output data
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         int i(const int i, const int j, const int idxI, const int idxJ) const override;

         /**
          * @brief Convert first index (1D)
          *
          * No reordering is possible in 1D
          *
          * @param i    Array index off first dimension
          */
         int i(const int i) const override;

         /**
          * @brief Convert second index (3D)
          *
          * Takes data index (i,j,k) and physical mode indexes (idxI, idxJ, idxK) from input data type, to compute corresponding data index j for output data
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
          * @brief Convert second index (3D)
          *
          * Takes data index (i,j) and physical mode indexes (idxI, idxJ) from input data type, to compute corresponding data index j for output data
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         int j(const int i, const int j, const int idxI, const int idxJ) const override;

         /**
          * @brief Convert second index (3D)
          *
          * Takes data index (i,j,k) and physical mode indexes (idxI, idxJ, idxK) from input data type, to compute corresponding data index k for output data
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const override;

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

   inline int PMIndexConv::i(const int i, const int j, const int idxI, const int idxJ) const
   {
      return j;
   }

   inline int PMIndexConv::i(const int i) const
   {
      return i;
   }

   inline int PMIndexConv::j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const
   {
      return i;
   }

   inline int PMIndexConv::j(const int i, const int j, const int idxI, const int idxJ) const
   {
      return i;
   }

   inline int PMIndexConv::k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const
   {
      return j;
   }

}
}

#endif // QUICC_PARALLEL_PMINDEXCONV_HPP
