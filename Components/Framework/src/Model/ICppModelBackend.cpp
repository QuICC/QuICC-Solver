/**
 * @file ICppModelBackend.cpp
 * @brief Source of the interface for a C++ model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Model/ICppModelBackend.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace Model {

void ICppModelBackend::blockInfo(int& tN, int& gN, ArrayI& shift, int& rhs,
   const SpectralFieldId& fId, const Resolution& res, const MHDFloat j,
   const BcMap& bcs) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, j)(0);
   tN = nN;

   int shiftI = this->nBc(fId);
   if (this->useGalerkin())
   {
      gN = (nN - shiftI);
   }
   else
   {
      shiftI = 0;
      gN = nN;
   }

   // Set galerkin shifts
   shift(0) = shiftI;
   shift(1) = 0;
   shift(2) = 0;

   rhs = 1;
}

int ICppModelBackend::blockSize(const SpectralFieldId& fId, const int j0,
   const int maxJ, const Resolution& res, const BcMap& bcs,
   const bool isGalerkin) const
{
   // Compute size
   auto s = 0;
   for (int j = j0; j <= maxJ; j++)
   {
      int tN, gN, rhs;
      ArrayI shift(3);
      this->blockInfo(tN, gN, shift, rhs, fId, res, j, bcs);
      if (isGalerkin)
      {
         s += gN;
      }
      else
      {
         s += tN;
      }
   }

   return s;
}

std::pair<int, int> ICppModelBackend::blockShape(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int j0, const int maxJ,
   const Resolution& res, const BcMap& bcs, const bool isGalerkin,
   const bool dropRows) const
{
   // Compute number of rows
   auto rows =
      this->blockSize(rowId, j0, maxJ, res, bcs, isGalerkin || dropRows);

   // Compute number of cols
   int cols = this->blockSize(colId, j0, maxJ, res, bcs, isGalerkin);

   return std::make_pair(rows, cols);
}

details::SystemInfo ICppModelBackend::systemInfo(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const SpectralFieldIds& fields, const int j0,
   const int maxJ, const Resolution& res, const BcMap& bcs,
   const bool isGalerkin, const bool dropRows) const
{
   auto shape =
      this->blockShape(rowId, colId, j0, maxJ, res, bcs, isGalerkin, dropRows);

   int sysN = 0;
   bool rowCount = true;
   bool colCount = true;
   int rowIdx = 0;
   int colIdx = 0;
   for (auto it = fields.begin(); it != fields.end(); ++it)
   {
      int s = this->blockSize(*it, j0, maxJ, res, bcs, isGalerkin);
      sysN += s;

      // Get block index of rowId
      if (rowCount && rowId != *it)
      {
         rowIdx += s;
      }
      else if (rowId == *it)
      {
         rowCount = false;
      }

      // Get block index of colId
      if (colCount && colId != *it)
      {
         colIdx += s;
      }
      else if (colId == *it)
      {
         colCount = false;
      }
   }

   details::SystemInfo info(sysN, shape.first, shape.second, rowIdx, colIdx);
   return info;
}

void ICppModelBackend::addBlock(SparseMatrix& mat, const SparseMatrix& block,
   const int rowShift, const int colShift, const MHDFloat coeff) const
{
   std::vector<Eigen::Triplet<MHDFloat>> triplets;
   triplets.reserve(block.nonZeros());
   for (int k = 0; k < block.outerSize(); ++k)
   {
      for (SparseMatrix::InnerIterator it(block, k); it; ++it)
      {
         triplets.emplace_back(Eigen::Triplet<MHDFloat>(it.row() + rowShift,
            it.col() + colShift, coeff * it.value()));
      }
   }
   SparseMatrix full(mat.rows(), mat.cols());
   full.setFromTriplets(triplets.begin(), triplets.end());
   mat += full;
}

void ICppModelBackend::buildFixedBlock(DecoupledZSparse& decMat, const int cols,
   const bool isComplexBlock,
   const std::vector<details::BlockDescription>& descr,
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const SpectralFieldIds& fields, const int matIdx, const std::size_t bcType,
   const Resolution& res, const int j0, const int maxJ, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator,
   const bool ignoreStart) const
{
   // Compute system size
   const auto sysInfo = systemInfo(rowId, colId, fields, j0, maxJ, res, bcs,
      this->useGalerkin(), false);
   const auto& sysN = sysInfo.systemSize;
   auto baseRowShift = sysInfo.startRow;
   auto baseColShift = sysInfo.startCol;
   if (ignoreStart)
   {
      baseRowShift = 0;
      baseColShift = 0;
   }

   // Resize matrices the first time
   if (decMat.real().size() == 0)
   {
      decMat.real().resize(sysN, cols);
      if (isComplexBlock)
      {
         decMat.imag().resize(sysN, cols);
      }
   }
   assert(decMat.real().rows() == sysN);
   assert(decMat.real().cols() == cols);
   if (isComplexBlock)
   {
      assert(decMat.imag().rows() == sysN);
      assert(decMat.imag().cols() == cols);
   }

   int tN, gN, rhs;
   ArrayI shift(3);

   bool needStencil = (this->useGalerkin());
   bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id() &&
                   !this->useGalerkin());

   for (auto&& d: descr)
   {
      assert(d.nRowShift == 0 || d.nColShift == 0);

      // Shift starting row
      int rowShift = baseRowShift;
      for (int s = 0; s < d.nRowShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, rowId, res, j0 + s, bcs);
         rowShift += gN;
      }

      // Shift starting col
      int colShift = baseColShift;
      for (int s = 0; s < d.nColShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, colId, res, j0 + s, bcs);
         colShift += cols;
      }

      int jShift = -d.nRowShift + d.nColShift;

      for (int j = j0 + d.nRowShift; j <= maxJ - d.nColShift; j++)
      {
         auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, j)(0);
         auto nNc = res.counter().dimensions(Dimensions::Space::SPECTRAL,
            j + jShift)(0);

         //
         // Build real part of block
         if (d.realOp)
         {
            auto bMat = d.realOp(nNr, nNc, j, d.opts, nds);
            assert(bMat.cols() == cols);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, j, j + jShift,
                  d.opts, res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, j + jShift, d.opts, res, bcs,
                  nds, isSplitOperator);
            }
            this->addBlock(decMat.real(), bMat, rowShift, colShift);
         }

         //
         // Build imaginary part of block
         if (d.imagOp)
         {
            auto bMat = d.imagOp(nNr, nNc, j, d.opts, nds);
            assert(bMat.cols() == cols);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, j, j + jShift,
                  d.opts, res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, j + jShift, d.opts, res, bcs,
                  nds, isSplitOperator);
            }
            this->addBlock(decMat.imag(), bMat, rowShift, colShift);
         }

         // Shift to next block
         this->blockInfo(tN, gN, shift, rhs, rowId, res, j, bcs);
         rowShift += gN;

         this->blockInfo(tN, gN, shift, rhs, colId, res, j + jShift, bcs);
         colShift += cols;
      }
   }
}

void ICppModelBackend::buildBlock(DecoupledZSparse& decMat,
   const std::vector<details::BlockDescription>& descr,
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const SpectralFieldIds& fields, const int matIdx, const std::size_t bcType,
   const Resolution& res, const int j0, const int maxJ, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator,
   const bool ignoreStart) const
{
   // Compute system size
   const auto sysInfo = systemInfo(rowId, colId, fields, j0, maxJ, res, bcs,
      this->useGalerkin(), false);
   const auto& sysN = sysInfo.systemSize;
   auto baseRowShift = sysInfo.startRow;
   auto baseColShift = sysInfo.startCol;
   if (ignoreStart)
   {
      baseRowShift = 0;
      baseColShift = 0;
   }

   // Resize matrices the first time
   if (decMat.real().size() == 0)
   {
      decMat.real().resize(sysN, sysN);
      if (this->isComplex(rowId))
      {
         decMat.imag().resize(sysN, sysN);
      }
   }
   assert(decMat.real().rows() == sysN);
   assert(decMat.real().cols() == sysN);
   if (this->isComplex(rowId))
   {
      assert(decMat.imag().rows() == sysN);
      assert(decMat.imag().cols() == sysN);
   }

   int tN, gN, rhs;
   ArrayI shift(3);

   bool needStencil = (this->useGalerkin());
   bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id() &&
                   !this->useGalerkin());

   for (auto&& d: descr)
   {
      assert((d.imagOp && decMat.imag().size() > 0) || !d.imagOp);
      assert(d.nRowShift == 0 || d.nColShift == 0);

      // Shift starting row
      int rowShift = baseRowShift;
      for (int s = 0; s < d.nRowShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, rowId, res, j0 + s, bcs);
         rowShift += gN;
      }

      // Shift starting col
      int colShift = baseColShift;
      for (int s = 0; s < d.nColShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, colId, res, j0 + s, bcs);
         colShift += gN;
      }

      int jShift = -d.nRowShift + d.nColShift;

      for (int j = j0 + d.nRowShift; j <= maxJ - d.nColShift; j++)
      {
         auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, j)(0);
         auto nNc = res.counter().dimensions(Dimensions::Space::SPECTRAL,
            j + jShift)(0);

         //
         // Build real part of block
         if (d.realOp)
         {
            auto bMat = d.realOp(nNr, nNc, j, d.opts, nds);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, j, j + jShift,
                  d.opts, res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, j + jShift, d.opts, res, bcs,
                  nds, isSplitOperator);
            }
            this->addBlock(decMat.real(), bMat, rowShift, colShift);
         }

         //
         // Build imaginary part of block
         if (d.imagOp)
         {
            auto bMat = d.imagOp(nNr, nNc, j, d.opts, nds);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, j, j + jShift,
                  d.opts, res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, j + jShift, d.opts, res, bcs,
                  nds, isSplitOperator);
            }
            this->addBlock(decMat.imag(), bMat, rowShift, colShift);
         }

         // Shift to next block
         this->blockInfo(tN, gN, shift, rhs, rowId, res, j, bcs);
         rowShift += gN;

         this->blockInfo(tN, gN, shift, rhs, colId, res, j + jShift, bcs);
         colShift += gN;
      }
   }
}

} // namespace Model
} // namespace QuICC
