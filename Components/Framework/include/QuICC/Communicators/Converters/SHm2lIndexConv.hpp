/**
 * @file SHm2lIndexConv.hpp
 * @brief Implementation of the index converter for spherical harmonics with m spectral ordering to l spectral ordering
 */

#ifndef QUICC_PARALLEL_SHM2LINDEXCONV_HPP
#define QUICC_PARALLEL_SHM2LINDEXCONV_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Communicators/Converters/IIndexConv.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the index converter for spherical harmonics with m spectral ordering to l spectral ordering
    *
    * The input data has ordering (NLM), ie. m ordering, and the corresponding indexes for (NML), i.e. m-ordering are computed.
    */
   class SHm2lIndexConv: public IIndexConv
   {
      public:
         /**
          * @brief Constructor
          */
         SHm2lIndexConv() = default;

         /**
          * @brief Destructor
          */
         virtual ~SHm2lIndexConv() = default;

         /**
          * @brief Initialize index conversions
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
          * No reodering is possible in 1D
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
          * @brief Convert second index (2D)
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
          * @brief Store minimal harmonic degree
          */
         int mMinL;
   };

   inline int SHm2lIndexConv::centralPadding(const int, const int) const
   {
      return 0;
   }

   inline int SHm2lIndexConv::i(const int, const int, const int, const int n, const int, const int) const
   {
      return n;
   }

   inline int SHm2lIndexConv::i(const int, const int, const int globalI, const int) const
   {
      return globalI;
   }

   inline int SHm2lIndexConv::i(const int ii) const
   {
      return ii;
   }

   inline int SHm2lIndexConv::j(const int, const int, const int ik, const int, const int, const int) const
   {
      return ik;
   }

   inline int SHm2lIndexConv::j(const int, const int ij, const int, const int) const
   {
      return ij;
   }

   inline int SHm2lIndexConv::k(const int, const int , const int , const int n, const int l, const int m) const
   {
      return l-this->mMinL;
   }

}
}

#endif // QUICC_PARALLEL_SHM2LINDEXCONV_HPP
