/**
 * @file ICppModelBackend.hpp
 * @brief Interface for a model backend
 */

#ifndef QUICC_MODEL_ICPPMODELBACKEND_HPP
#define QUICC_MODEL_ICPPMODELBACKEND_HPP

// System includes
//
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Model/EquationInfo.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/NonDimensional/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Model {

namespace details {

/**
 * @brief Struct holding general information on full system size
 */
struct SystemInfo
{
   /// Size of full system
   int systemSize;
   /// Rows per block
   int blockRows;
   /// Columns per block
   int blockCols;
   /// Starting row
   int startRow;
   /// Starting column
   int startCol;

   /**
    * @brief ctor
    *
    * @param size System of full system
    * @param rows Number of rows per block
    * @param cols Number of columns per block
    * @param row  Starting row
    * @param col  Starting col
    */
   SystemInfo(const int size, const int rows, const int cols, const int row,
      const int col) :
       systemSize(size),
       blockRows(rows),
       blockCols(cols),
       startRow(row),
       startCol(col) {};
};

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
} // namespace details

/**
 * @brief Interface for a model backend
 */
class ICppModelBackend : public IModelBackend
{
public:
   /**
    * @brief Constructor
    */
   ICppModelBackend() = default;

   /**
    * @brief Destructor
    */
   virtual ~ICppModelBackend() = default;

protected:
   /**
    * @brief Operators are complex?
    *
    * @param fId  Field ID
    */
   virtual bool isComplex(const SpectralFieldId& fId) const = 0;

   /**
    * @brief Get coupled fields
    *
    * @param fId  Field ID
    */
   virtual SpectralFieldIds implicitFields(
      const SpectralFieldId& fId) const = 0;

   /**
    * @brief Apply tau line for boundary condition
    *
    * @param mat     Input/Output matrix to apply tau line to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param j       2D index
    * @param opts    Additional options
    * @param res     Resolution object
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    * @param isSplitOperator  Is second operator of split 4th order system?
    */
   virtual void applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int j,
      std::shared_ptr<details::BlockOptions> opts, const Resolution& res,
      const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator) const = 0;

   /**
    * @brief Apply galerkin stencil for boundary condition
    *
    * @param mat     Input/Output matrix to apply stencil to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param jr      Row space index
    * @param jc      Column space index
    * @param opts    Additional options
    * @param res     Resolution object
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    */
   virtual void applyGalerkinStencil(SparseMatrix& decMat,
      const SpectralFieldId& rowId, const SpectralFieldId& colId, const int jr,
      const int jc, std::shared_ptr<details::BlockOptions> opts,
      const Resolution& res, const BcMap& bcs,
      const NonDimensional::NdMap& nds) const = 0;

   /**
    * @brief Number of boundary conditions
    *
    * @fId  Field ID
    */
   virtual int nBc(const SpectralFieldId& fId) const = 0;

   /**
    * @brief Get operator block information
    *
    * @param tN      Tau radial size
    * @param gN      Galerkin radial truncation
    * @param shift   Shift in each direction due to Galerkin basis
    * @param fId     ID of the field
    * @param res     Resolution object
    * @param j       Index on which fast resolution may depend
    * @param bcs     Boundary conditions
    */
   void blockInfo(int& tN, int& gN, ArrayI& shift, int& rhs,
      const SpectralFieldId& fId, const Resolution& res, const MHDFloat j,
      const BcMap& bcs) const;

   /**
    * @brief Compute size information of full system
    *
    * @param rowId  Equation Field ID
    * @param colId  Field Id
    * @param fields List of fields Id
    * @param j0     First 2D index
    * @param maxJ   Max 2D index
    * @param res    Resolution object
    * @param bcs    Boundary conditions
    * @param isGalerkin Use Galerkin scheme?
    * @param dropRows?  Drop Tau line rows
    */
   details::SystemInfo systemInfo(const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const SpectralFieldIds& fields,
      const int j0, const int maxJ, const Resolution& res, const BcMap& bcs,
      const bool isGalerkin, const bool dropRows) const;

   /**
    * @brief Get operator information
    *
    * @param fId  Field ID
    * @param j0   First 2D index
    * @param maxJ Max 2D index
    * @param res  Resolution object
    * @param bcs  Boundary conditions
    * @param isGalerkin Use Galerkin scheme?
    */
   int blockSize(const SpectralFieldId& fId, const int j0, const int maxJ,
      const Resolution& res, const BcMap& bcs, const bool isGalerkin) const;

   /**
    * @brief Get operator block shape
    *
    * @param rowId   Equation Field ID
    * @param colId    Field Id
    * @param j0   First 2D index
    * @param maxJ Max 2D index
    * @param res  Resolution object
    * @param bcs  Boundary conditions
    * @param isGalerkin Use Galerkin scheme?
    * @param dropRows?  Drop Tau line rows
    */
   std::pair<int, int> blockShape(const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int j0, const int maxJ,
      const Resolution& res, const BcMap& bcs, const bool isGalerkin,
      const bool dropRows) const;

   /**
    * @brief Build 2D matrix block with fixed width from description
    *
    * @param decMat  Ouput matrix
    * @param cols    Number of columns
    * @param isComplexBlock block is complex?
    * @param descr   Block description
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param matIdx  Matrix ID
    * @param bcType  Type of boundary condition
    * @param res     Resolution object
    * @param j0      First 2D index
    * @param maxJ    Max 2D index
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    * @param isSplitOperator  Set operator of split system
    * @param ignoreStart  Ignore start shift
    */
   void buildFixedBlock(DecoupledZSparse& decMat, const int cols,
      const bool isComplex, const std::vector<details::BlockDescription>& descr,
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const SpectralFieldIds& fields, const int matIdx,
      const std::size_t bcType, const Resolution& res, const int j0,
      const int maxJ, const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator, const bool ignoreStart = false) const;

   /**
    * @brief Build 2D matrix block from description
    *
    * @param decMat  Ouput matrix
    * @param descr   Block description
    * @param rowId   Field ID of block matrix row
    * @param colId   Field ID of block matrix column
    * @param matIdx  Matrix ID
    * @param bcType  Type of boundary condition
    * @param res     Resolution object
    * @param j0      First 2D index
    * @param maxJ    Max 2D index
    * @param bcs     Boundary conditions for each field
    * @param nds     Nondimension parameters
    * @param isSplitOperator  Set operator of split system
    * @param ignoreStart  Ignore start shift
    */
   void buildBlock(DecoupledZSparse& decMat,
      const std::vector<details::BlockDescription>& descr,
      const SpectralFieldId& rowId, const SpectralFieldId& colId,
      const SpectralFieldIds& fields, const int matIdx,
      const std::size_t bcType, const Resolution& res, const int j0,
      const int maxJ, const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator, const bool ignoreStart = false) const;

   /**
    * @brief Add block matrix to full system matrix
    *
    * @param mat Input/Output matrix to add matrix block to
    * @param block matrix block to add
    * @param rowShift   Start row of matrix block
    * @param colShift   Start column of matrix block
    * @param coeff      Scaling coefficient of matrix block
    */
   void addBlock(SparseMatrix& mat, const SparseMatrix& block,
      const int rowShift, const int colShift, const MHDFloat coeff = 1.0) const;

private:
};

} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_ICPPMODELBACKEND_HPP
