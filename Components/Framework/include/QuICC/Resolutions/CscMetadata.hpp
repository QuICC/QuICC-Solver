/**
 * @file CscMetadata.hpp
 * @brief Compressed index metadata
 */

#ifndef QUICC_CSCMETADATA_HPP
#define QUICC_CSCMETADATA_HPP

// System includes
//
#include <vector>
#include <cstdint>

namespace QuICC {
   
   /**
    * @brief Small struct for passing compressed index metadata
    */
   struct CscMetadata
   {
      /// Global 1D size
      std::uint32_t global1D;
      /// Global 2D size
      std::uint32_t global2D;
      /// Global 3D size
      std::uint32_t global3D;

      /// 1D size of each 2D profile
      std::vector<std::uint32_t> dim1D;

      /// Pointers for 2D
      std::vector<std::uint32_t> ptr2D;

      /// Global indexes for 2D
      std::vector<std::uint32_t> idx2D;
   };

} // namespace QuICC

#endif // QUICC_CSCMETADATA_HPP
