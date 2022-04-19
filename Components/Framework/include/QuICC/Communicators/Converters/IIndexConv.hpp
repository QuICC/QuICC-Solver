/**
 * @file IIndexConv.hpp
 * @brief Interface for generic index converter
 */

#ifndef QUICC_PARALLEL_IINDEXCONV_HPP
#define QUICC_PARALLEL_IINDEXCONV_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Interface for generic index converter
    */
   class IIndexConv
   {
      public:
         /**
          * @brief Constructor
          */
         IIndexConv();

         /**
          * @brief Destructor
          */
         virtual ~IIndexConv();

         /**
          * @brief Initialize index conversions
          */
         virtual void init(const Resolution& res, const Dimensions::Transform::Id id);

         /**
          * @brief Compute shift due to central padding
          *
          * @param idx  Reference index
          * @param k    Array index of third dimension
          */
         virtual int centralPadding(const int idx, const int k) const = 0;

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
         virtual int i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const = 0;

         /**
          * @brief Convert first index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         virtual int iS(const int i, const int j, const int k) const = 0;

         /**
          * @brief Convert first index (2D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         virtual int i(const int i, const int j, const int idxI, const int idxJ) const = 0;

         /**
          * @brief Convert first index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         virtual int iS(const int i, const int j) const = 0;

         /**
          * @brief Convert first index (1D)
          *
          * @param i    Array index off first dimension
          */
         virtual int i(const int i) const = 0;

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
         virtual int j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const = 0;

         /**
          * @brief Convert second index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         virtual int jS(const int i, const int j, const int k) const = 0;

         /**
          * @brief Convert second index (2D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         virtual int j(const int i, const int j, const int idxI, const int idxJ) const = 0;

         /**
          * @brief Convert second index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         virtual int jS(const int i, const int j) const = 0;

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
         virtual int k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK) const = 0;

         /**
          * @brief Convert third index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         virtual int kS(const int i, const int j, const int k) const = 0;
   };

}
}

#endif // QUICC_PARALLEL_IINDEXCONV_HPP
