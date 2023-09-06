/**
 * @file PassthroughIndexConv.hpp
 * @brief Implementation of the passthrough index converter
 */

#ifndef QUICC_PARALLEL_PASSTHROUGHINDEXCONV_HPP
#define QUICC_PARALLEL_PASSTHROUGHINDEXCONV_HPP

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
    * @brief Implementation of the passthrough index converter
    *
    * Input and output data have the same data ordering. Resolution for the first index might different (dealiasing).
    */
   class PassthroughIndexConv: public IIndexConv
   {
      public:
         /**
          * @brief Constructor
          */
         PassthroughIndexConv() = default;

         /**
          * @brief Destructor
          */
         virtual ~PassthroughIndexConv() = default;

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
          * Takes data index (i,j,k) and physical mode indexes (idxI, idxJ, idxK) from input data type, to compute corresponding data index j for output data
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
          * Takes data index (i,j) and physical mode indexes (idxI, idxJ) from input data type, to compute corresponding data index j for output data
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const override;
   };

   inline int PassthroughIndexConv::centralPadding(const int, const int) const
   {
      return 0;
   }

   inline int PassthroughIndexConv::i(const int i_, const int, const int, const int, const int, const int) const
   {
      return i_;
   }

   inline int PassthroughIndexConv::i(const int i_, const int, const int, const int) const
   {
      return i_;
   }

   inline int PassthroughIndexConv::i(const int i_) const
   {
      return i_;
   }

   inline int PassthroughIndexConv::j(const int, const int j_, const int, const int, const int, const int) const
   {
      return j_;
   }

   inline int PassthroughIndexConv::j(const int, const int j_, const int, const int) const
   {
      return j_;
   }

   inline int PassthroughIndexConv::k(const int, const int , const int k_, const int, const int, const int) const
   {
      return k_;
   }

}
}

#endif // QUICC_PARALLEL_PASSTHROUGHINDEXCONV_HPP
