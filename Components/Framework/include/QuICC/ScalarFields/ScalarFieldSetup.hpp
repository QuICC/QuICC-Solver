/**
 * @file ScalarFieldSetup.hpp
 * @brief Single configuration class for the different scalar fields
 */

#ifndef QUICC_DATATYPES_SCALARFIELDSETUP_HPP
#define QUICC_DATATYPES_SCALARFIELDSETUP_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "QuICC/Resolutions/CscMetadata.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Single configuration class for the different scalar fields
    */
   class ScalarFieldSetup
   {
      public:
         /**
          * @brief Constructor for 3D scalar field
          */
         ScalarFieldSetup(std::shared_ptr<CscMetadata> spMeta, std::shared_ptr<Memory::memory_resource> mem);

         /**
          * @brief Destructor
          */
         virtual ~ScalarFieldSetup() = default;

         /**
          * @brief Get number of rows in data storage
          */
         int dataRows() const;

         /**
          * @brief Get number of rows in data storage
          */
         int dataCols() const;

         /**
          * @brief Get column index corresponding to the given 2D and 3D indexes
          */
         int colIdx(const int j, const int k = 0) const;

         /**
          * @brief Get index of start of block for 3D index
          */
         int blockIdx(const int k) const;

         /**
          * @brief Get number of rows of block for 3D index
          */
         int blockRows(const int k) const;

         /**
          * @brief Get number of columns of block for 3D index
          */
         int blockCols(const int k) const;

         /**
          * @brief Get number of blocks
          */
         int nBlock() const;

         /**
          * @brief Get ptr and idx vectors for View
          */
         std::shared_ptr<CscMetadata> viewMeta() const;

         /**
          * @brief Get memory resource
          */
         std::shared_ptr<Memory::memory_resource> mem() const;

      protected:

      private:
         /**
          * @brief Array of dimensions for the first data dimension
          */
         SharedArrayI   mspDim1D;

         /**
          * @brief Array of dimensions for the first data dimension
          */
         SharedArrayI   mspDim2D;

         /**
          * @brief Array of dimensions for the first data dimension
          */
         int mDim3D;

         /**
          * @brief Number of rows of 2D storage
          */
         int mDataRows;

         /**
          * @brief Number of columns of 2D storage
          */
         int mDataCols;

         /**
          * @brief Csc Metadata
          */
         std::shared_ptr<CscMetadata> mspMeta;

         /**
          * @brief Memory resources
          */
         std::shared_ptr<Memory::memory_resource> mMem;
   };

} // namespace Datatypes
} // namespace QuICC

#endif // QUICC_DATATYPES_SCALARFIELDSETUP_HPP
