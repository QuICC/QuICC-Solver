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
   const int nTauLines, const int nN) const
{
   tN = nN;

   int shiftI = nTauLines;
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

int ICppModelBackend::blockSize(const int nTauLines, const std::vector<int>& nNs,
   const bool isGalerkin) const
{
   // Compute size
   auto s = 0;
   for (auto&& nN: nNs)
   {
      int tN, gN, rhs;
      ArrayI shift(3);
      this->blockInfo(tN, gN, shift, rhs, nTauLines, nN);
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

std::pair<int, int> ICppModelBackend::blockShape(const int nTauLinesRow, const int nTauLinesCol, const std::vector<int>& nNs, const bool isGalerkin,
   const bool dropRows) const
{
   // Compute number of rows
   auto rows =
      this->blockSize(nTauLinesRow, nNs, isGalerkin || dropRows);

   // Compute number of cols
   int cols = this->blockSize(nTauLinesCol, nNs, isGalerkin);

   return std::make_pair(rows, cols);
}

details::SystemInfo ICppModelBackend::systemInfo(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const SpectralFieldIds& fields, const std::vector<int>& nNs,
   const bool isGalerkin, const bool dropRows) const
{
   auto nTauLinesRow = this->nBc(rowId);
   auto nTauLinesCol = this->nBc(colId);

   auto shape =
      this->blockShape(nTauLinesRow, nTauLinesCol, nNs, isGalerkin, dropRows);

   int sysN = 0;
   bool rowCount = true;
   bool colCount = true;
   int rowIdx = 0;
   int colIdx = 0;
   for (auto it = fields.begin(); it != fields.end(); ++it)
   {
      auto nTauLines = this->nBc(*it);
      int s = this->blockSize(nTauLines, nNs, isGalerkin);
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

std::tuple<int,int,int,int> ICppModelBackend::blockSystemInfo(const SpectralFieldId& rowId, const SpectralFieldId& colId, const SpectralFieldIds& fields, const std::vector<int>& nNs, const bool ignoreStart, const int fixedCols) const
{
   // Compute system size
   const auto sysInfo = systemInfo(rowId, colId, fields, nNs,
      this->useGalerkin(), false);
   const auto& sysN = sysInfo.systemSize;
   int sysCols = sysN;
   if(fixedCols > 0)
   {
      sysCols = fixedCols;
   }
   auto baseRowShift = sysInfo.startRow;
   auto baseColShift = sysInfo.startCol;
   if (ignoreStart)
   {
      baseRowShift = 0;
      baseColShift = 0;
   }

   std::tuple<int,int,int,int> info = std::make_tuple(sysN, sysCols, baseRowShift, baseColShift);
   return info;
}

void ICppModelBackend::computeBlockShift(int& blockShift, const int s0, const int nShift, const int nTauLines, const std::vector<int>& nNs, const int fixedCols) const
{
   int tN, gN, rhs;
   ArrayI shift(3);

   // Shift starting block
   if(fixedCols > 0)
   {
      for (int s = s0; s < s0 + nShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, nTauLines, nNs.at(s));
         blockShift += fixedCols;
      }
   }
   else
   {
      for (int s = s0; s < s0 + nShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, nTauLines, nNs.at(s));
         blockShift += gN;
      }
   }
}

void ICppModelBackend::buildBlock(DecoupledZSparse& decMat,
   const bool isComplexBlock,
   const std::vector<details::BlockDescription>& descr,
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const SpectralFieldIds& fields, const int matIdx, const std::size_t bcType,
   const Resolution& res, const int j0, const int maxJ, const std::vector<int>& nNs, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator, const int fixedCols,
   const bool ignoreStart) const
{

   auto nTauLinesRow = this->nBc(rowId);
   auto nTauLinesCol = this->nBc(colId);

   // Compute system size
   const auto sysInfo = this->blockSystemInfo(rowId, colId, fields, nNs, ignoreStart, fixedCols);
   const auto sysRows = std::get<0>(sysInfo);
   const auto sysCols = std::get<1>(sysInfo);
   const auto baseRowShift = std::get<2>(sysInfo);
   const auto baseColShift = std::get<3>(sysInfo);

   // Resize matrices the first time
   if (decMat.real().size() == 0)
   {
      decMat.real().resize(sysRows, sysCols);
      if (isComplexBlock)
      {
         decMat.imag().resize(sysRows, sysCols);
      }
   }
   assert(decMat.real().rows() == sysRows);
   assert(decMat.real().cols() == sysCols);
   if (isComplexBlock)
   {
      assert(decMat.imag().rows() == sysRows);
      assert(decMat.imag().cols() == sysCols);
   }

   bool needStencil = (this->useGalerkin());
   bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id() &&
                   !this->useGalerkin());

   for (auto&& d: descr)
   {
      assert(d.nRowShift == 0 || d.nColShift == 0);

      // Shift starting row
      int rowShift = baseRowShift;
      this->computeBlockShift(rowShift, 0, d.nRowShift, nTauLinesRow, nNs, -1);

      // Shift starting col
      int colShift = baseColShift;
      this->computeBlockShift(colShift, 0, d.nColShift, nTauLinesCol, nNs, fixedCols);

      int jShift = -d.nRowShift + d.nColShift;

      for (int j = j0 + d.nRowShift; j <= maxJ - d.nColShift; j++)
      {
         auto nNr = nNs.at(j - j0);
         auto nNc = nNs.at(j + jShift - j0);

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
         int s_ = j-j0;
         this->computeBlockShift(rowShift, s_, 1, nTauLinesRow, nNs, -1);

         s_ = j + jShift - j0;
         this->computeBlockShift(colShift, s_, 1, nTauLinesCol, nNs, fixedCols);
      }
   }
}

} // namespace Model
} // namespace QuICC
