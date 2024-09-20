/**
 * @file ISphericalModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Model/ISphericalModelBackend.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace Model {

void ISphericalModelBackend::blockInfo(int& tN, int& gN, ArrayI& shift,
   int& rhs, const SpectralFieldId& fId, const Resolution& res,
   const MHDFloat l, const BcMap& bcs) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
   tN = nN;

   int shiftR = this->nBc(fId);
   if (this->useGalerkin())
   {
      gN = (nN - shiftR);
   }
   else
   {
      shiftR = 0;
      gN = nN;
   }

   // Set galerkin shifts
   shift(0) = shiftR;
   shift(1) = 0;
   shift(2) = 0;

   rhs = 1;
}

void ISphericalModelBackend::buildFixedBlock(DecoupledZSparse& decMat,
   const int cols, const bool isComplexBlock,
   const std::vector<details::BlockDescription>& descr,
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const SpectralFieldIds& fields, const int matIdx, const std::size_t bcType,
   const Resolution& res, const int l0, const int maxL, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator,
   const bool ignoreStart) const
{
   // Compute system size
   const auto sysInfo = systemInfo(rowId, colId, fields, l0, maxL, res, bcs,
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
         this->blockInfo(tN, gN, shift, rhs, rowId, res, l0 + s, bcs);
         rowShift += gN;
      }

      // Shift starting col
      int colShift = baseColShift;
      for (int s = 0; s < d.nColShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, colId, res, l0 + s, bcs);
         colShift += cols;
      }

      int lShift = -d.nRowShift + d.nColShift;

      for (int l = l0 + d.nRowShift; l <= maxL - d.nColShift; l++)
      {
         auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
         auto nNc = res.counter().dimensions(Dimensions::Space::SPECTRAL,
            l + lShift)(0);

         //
         // Build real part of block
         if (d.realOp)
         {
            auto bMat = d.realOp(nNr, nNc, l, d.opts, nds);
            assert(bMat.cols() == cols);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, l, l + lShift,
                  res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, l + lShift, res, bcs, nds,
                  isSplitOperator);
            }
            this->addBlock(decMat.real(), bMat, rowShift, colShift);
         }

         //
         // Build imaginary part of block
         if (d.imagOp)
         {
            auto bMat = d.imagOp(nNr, nNc, l, d.opts, nds);
            assert(bMat.cols() == cols);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, l, l + lShift,
                  res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, l + lShift, res, bcs, nds,
                  isSplitOperator);
            }
            this->addBlock(decMat.imag(), bMat, rowShift, colShift);
         }

         // Shift to next block
         this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
         rowShift += gN;

         this->blockInfo(tN, gN, shift, rhs, colId, res, l + lShift, bcs);
         colShift += cols;
      }
   }
}

void ISphericalModelBackend::buildBlock(DecoupledZSparse& decMat,
   const std::vector<details::BlockDescription>& descr,
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const SpectralFieldIds& fields, const int matIdx, const std::size_t bcType,
   const Resolution& res, const int l0, const int maxL, const BcMap& bcs,
   const NonDimensional::NdMap& nds, const bool isSplitOperator,
   const bool ignoreStart) const
{
   // Compute system size
   const auto sysInfo = systemInfo(rowId, colId, fields, l0, maxL, res, bcs,
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
         this->blockInfo(tN, gN, shift, rhs, rowId, res, l0 + s, bcs);
         rowShift += gN;
      }

      // Shift starting col
      int colShift = baseColShift;
      for (int s = 0; s < d.nColShift; s++)
      {
         this->blockInfo(tN, gN, shift, rhs, colId, res, l0 + s, bcs);
         colShift += gN;
      }

      int lShift = -d.nRowShift + d.nColShift;

      for (int l = l0 + d.nRowShift; l <= maxL - d.nColShift; l++)
      {
         auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
         auto nNc = res.counter().dimensions(Dimensions::Space::SPECTRAL,
            l + lShift)(0);

         //
         // Build real part of block
         if (d.realOp)
         {
            auto bMat = d.realOp(nNr, nNc, l, d.opts, nds);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, l, l + lShift,
                  res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, l + lShift, res, bcs, nds,
                  isSplitOperator);
            }
            this->addBlock(decMat.real(), bMat, rowShift, colShift);
         }

         //
         // Build imaginary part of block
         if (d.imagOp)
         {
            auto bMat = d.imagOp(nNr, nNc, l, d.opts, nds);

            if (needStencil)
            {
               this->applyGalerkinStencil(bMat, rowId, colId, l, l + lShift,
                  res, bcs, nds);
            }
            else if (needTau)
            {
               this->applyTau(bMat, rowId, colId, l + lShift, res, bcs, nds,
                  isSplitOperator);
            }
            this->addBlock(decMat.imag(), bMat, rowShift, colShift);
         }

         // Shift to next block
         this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
         rowShift += gN;

         this->blockInfo(tN, gN, shift, rhs, colId, res, l + lShift, bcs);
         colShift += gN;
      }
   }
}

void ISphericalModelBackend::addBlock(SparseMatrix& mat,
   const SparseMatrix& block, const int rowShift, const int colShift,
   const MHDFloat coeff) const
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

int ISphericalModelBackend::blockSize(const SpectralFieldId& fId, const int l0,
   const int maxL, const Resolution& res, const BcMap& bcs,
   const bool isGalerkin) const
{
   // Compute size
   auto s = 0;
   for (int l = l0; l <= maxL; l++)
   {
      int tN, gN, rhs;
      ArrayI shift(3);
      this->blockInfo(tN, gN, shift, rhs, fId, res, l, bcs);
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

std::pair<int, int> ISphericalModelBackend::blockShape(
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int l0,
   const int maxL, const Resolution& res, const BcMap& bcs,
   const bool isGalerkin, const bool dropRows) const
{
   // Compute number of rows
   auto rows =
      this->blockSize(rowId, l0, maxL, res, bcs, isGalerkin || dropRows);

   // Compute number of cols
   int cols = this->blockSize(colId, l0, maxL, res, bcs, isGalerkin);

   return std::make_pair(rows, cols);
}

details::SystemInfo ISphericalModelBackend::systemInfo(
   const SpectralFieldId& rowId, const SpectralFieldId& colId,
   const SpectralFieldIds& fields, const int l0, const int maxL,
   const Resolution& res, const BcMap& bcs, const bool isGalerkin,
   const bool dropRows) const
{
   auto shape =
      this->blockShape(rowId, colId, l0, maxL, res, bcs, isGalerkin, dropRows);

   int sysN = 0;
   bool rowCount = true;
   bool colCount = true;
   int rowIdx = 0;
   int colIdx = 0;
   for (auto it = fields.begin(); it != fields.end(); ++it)
   {
      int s = this->blockSize(*it, l0, maxL, res, bcs, isGalerkin);
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

} // namespace Model
} // namespace QuICC
