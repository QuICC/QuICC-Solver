/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for test model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/TestSuite/Framework/Io/TestBackend.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

TestBackend::TestBackend()
{}

std::vector<std::string> TestBackend::fieldNames() const
{
   std::vector<std::string> names = {
      PhysicalNames::Velocity().tag(),
      PhysicalNames::Magnetic().tag(),
      PhysicalNames::Temperature().tag()
   };

   return names;
}

std::vector<std::string> TestBackend::paramNames() const
{
   std::vector<std::string> names;

   return names;
}

std::vector<bool> TestBackend::isPeriodicBox() const
{
   std::vector<bool> periodic = {false, false, false};

   return periodic;
}


std::map<std::string, MHDFloat> TestBackend::automaticParameters(
   const std::map<std::string, MHDFloat>& cfg) const
{
   std::map<std::string, MHDFloat> params;

   return params;
}

void TestBackend::equationInfo(QuICC::Model::EquationInfo& info, const SpectralFieldId& fId,
   const Resolution& res) const
{
   // Operators are real
   info.isComplex = false;

   // Splitting 4th poloidal equation into two systems
   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::POL))
   {
      info.isSplitEquation = this->useSplitEquation();
   }
   else
   {
      info.isSplitEquation = false;
   }

   // Implicit coupled fields
   info.im.clear();

   // Explicit linear terms
   info.exL.clear();

   // Explicit nonlinear terms
   info.exNL.clear();

   // Explicit nextstep terms
   info.exNS.clear();

   // Index mode
   info.indexMode =
      static_cast<int>(Equations::CouplingIndexType::SLOWEST_MULTI_RHS);
}

void TestBackend::operatorInfo(QuICC::Model::OperatorInfo& info, const SpectralFieldId& fId,
   const Resolution& res, const Equations::Tools::ICoupling& coupling,
   const BcMap& bcs) const
{
   // Loop overall matrices/eigs
   for (int idx = 0; idx < info.tauN.size(); ++idx)
   {
      auto eigs = coupling.getIndexes(res, idx);

      int tN, gN, rhs;
      ArrayI shift(3);

      this->blockInfo(tN, gN, shift, rhs, fId, res, eigs.at(0), bcs);

      info.tauN(idx) = tN;
      info.galN(idx) = gN;
      info.galShift.row(idx) = shift;
      info.rhsCols(idx) = rhs;

      // Compute system size
      int sN = 0;
      for (auto f: this->implicitFields(fId))
      {
         this->blockInfo(tN, gN, shift, rhs, f, res, eigs.at(0), bcs);
         sN += gN;
      }

      if (sN == 0)
      {
         sN = info.galN(idx);
      }

      info.sysN(idx) = sN;
   }
}

void TestBackend::modelMatrix(DecoupledZSparse& rModelMatrix,
   const std::size_t opId,
   const Equations::CouplingInformation::FieldId_range imRange,
   const int matIdx, const std::size_t bcType, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   throw std::logic_error("modelMatrix: Should not be required");
}

void TestBackend::galerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   throw std::logic_error("galerkinStencil: Should not be required");
}

void TestBackend::explicitBlock(DecoupledZSparse& mat,
   const SpectralFieldId& fId, const std::size_t opId,
   const SpectralFieldId fieldId, const int matIdx, const Resolution& res,
   const std::vector<MHDFloat>& eigs, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   throw std::logic_error("explicitBlock: Should not be required");
}

bool TestBackend::isComplex(const SpectralFieldId& fId) const
{
   return 0;
}

int TestBackend::nBc(const SpectralFieldId& fId) const
{
   return 0;
}

TestBackend::SpectralFieldIds TestBackend::implicitFields(
   const SpectralFieldId& fId) const
{
   SpectralFieldIds fields = {fId};

   return fields;
}

} // namespace Io
} // namespace Framewok
} // namespace TestSuite
} // namespace QuICC
