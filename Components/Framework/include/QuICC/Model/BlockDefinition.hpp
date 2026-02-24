/**
 * @file BlockDefinition.hpp
 * @brief Class for providing a model block definition
 */

#ifndef QUICC_MODEL_BLOCKDEFINITION_HPP
#define QUICC_MODEL_BLOCKDEFINITION_HPP

// System includes
//
#include <memory>
#include <vector>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Model {

namespace details {

/**
 * @brief Base class for proving options for system block builder
 */
struct BlockOptions
{
   /**
    * @brief default ctor
    */
   BlockOptions() = default;

   /**
    * @brief default dtor
    */
   virtual ~BlockOptions() = default;
};

/**
 * @brief Operator block description
 */
struct BlockDescription
{
   /// Starting row shift
   int nRowShift = 0;
   /// Starting column shift
   int nColShift = 0;
   /// Options to build block
   std::shared_ptr<BlockOptions> opts;
   /// Builder for real part
   SparseMatrix (*realOp)(const int nNr, const int nNc, const int j,
      std::shared_ptr<BlockOptions> opts,
      const NonDimensional::NdMap& nds) = nullptr;
   /// Builder for imaginary part
   SparseMatrix (*imagOp)(const int nNr, const int nNc, const int j,
      std::shared_ptr<BlockOptions> opts,
      const NonDimensional::NdMap& nds) = nullptr;
};

/**
 * @brief Class full block definition
 */
struct BlockDefinition
{
   /// Field ID of row
   SpectralFieldId rowId;
   /// Field ID of column
   SpectralFieldId colId;

   /// Block is complex?
   bool isComplex;
   /// Use galerkin description?
   bool isGalerkin;

   std::vector<BlockDescription> descr;
};

} // namespace details
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BLOCKDEFINITION_HPP
