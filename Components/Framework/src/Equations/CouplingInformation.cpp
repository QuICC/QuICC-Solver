/**
 * @file CouplingInformation.cpp
 * @brief Source of the base implementation of a scalar equation
 */

// Debug includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/CouplingInformation.hpp"

// Project includes
//
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"

namespace QuICC {

namespace Equations {

   CouplingInformation::CouplingInformation()
      : mEquationType(TRIVIAL), mHasNonlinear(false), mHasQuasiInverse(false), mHasSource(false), mHasBoundaryValue(false), mIsComplex(true), mIsSplitEquation(false), mIsGalerkin(false), mIndexType(CouplingIndexType::SLOWEST_SINGLE_RHS), mNSystems(0), mFieldIndex(-1), mSolverIndex(-1), mFieldStart(-1)
   {
   }

   CouplingInformation::~CouplingInformation()
   {
   }

   CouplingInformation::EquationTypeId CouplingInformation::equationType() const
   {
      return this->mEquationType;
   }

   bool CouplingInformation::hasNonlinear() const
   {
      return this->mHasNonlinear;
   }

   bool CouplingInformation::hasQuasiInverse() const
   {
      return this->mHasQuasiInverse;
   }

   bool CouplingInformation::hasSource() const
   {
      return this->mHasSource;
   }

   bool CouplingInformation::hasBoundaryValue() const
   {
      return this->mHasBoundaryValue;
   }

   bool CouplingInformation::isComplex() const
   {
      return this->mIsComplex;
   }

   bool CouplingInformation::isSplitEquation() const
   {
      return this->mIsSplitEquation;
   }

   bool CouplingInformation::isGalerkin() const
   {
      return this->mIsGalerkin;
   }

   CouplingIndexType CouplingInformation::indexType() const
   {
      return this->mIndexType;
   }

   const Tools::ICoupling& CouplingInformation::couplingTools() const
   {
      return *(this->mspCouplingTools);
   }

   int CouplingInformation::nBlocks() const
   {
      return this->mImplicitFields.size();
   }

   int CouplingInformation::nSystems() const
   {
      return this->mNSystems;
   }

   int CouplingInformation::tauN(const int idx) const
   {
      return this->mTauNs(idx);
   }

   int CouplingInformation::galerkinN(const int idx) const
   {
      return this->mGalerkinNs(idx);
   }

   int CouplingInformation::galerkinShift(const int idx, const int dim) const
   {
      return this->mGalerkinShifts(idx, dim);
   }

   int CouplingInformation::systemN(const int idx) const
   {
      return this->mSystemNs(idx);
   }

   int CouplingInformation::fieldIndex() const
   {
      return this->mFieldIndex;
   }

   int CouplingInformation::solverIndex() const
   {
      return this->mSolverIndex;
   }

   int CouplingInformation::fieldStart() const
   {
      return this->mFieldStart;
   }

   int CouplingInformation::rhsCols(const int idx) const
   {
      return this->mRhsCols(idx);
   }

   void CouplingInformation::sortImplicitFields(const std::size_t fieldId, const FieldComponents::Spectral::Id compId)
   {
      // Sort the implicit fields
      std::sort(this->mImplicitFields.begin(), this->mImplicitFields.end());

      // Extract the position of the equation field
      FieldId_iterator pos = std::find(this->mImplicitFields.begin(), this->mImplicitFields.end(), std::make_pair(fieldId, compId));
      assert((this->equationType() == TRIVIAL && pos == this->mImplicitFields.end()) || (this->equationType() == WRAPPER && pos == this->mImplicitFields.end()) || pos != this->mImplicitFields.end());

      // Set initial field index
      this->mFieldIndex = pos - this->mImplicitFields.begin();

      // Report position of equation field
      DebuggerMacro_showValue("Coupling information field index: ", 1, this->mFieldIndex);
   }

   void CouplingInformation::addImplicitField(const std::size_t fieldId, const FieldComponents::Spectral::Id compId)
   {
      this->mImplicitFields.push_back(std::make_pair(fieldId,compId));
   }

   void CouplingInformation::addExplicitField(const std::size_t fieldId, const FieldComponents::Spectral::Id compId, const std::size_t opId)
   {
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         this->mExplicitLFields.push_back(std::make_pair(fieldId,compId));

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         this->mExplicitNLFields.push_back(std::make_pair(fieldId,compId));

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         this->mExplicitNSFields.push_back(std::make_pair(fieldId,compId));
      }

   }

   void CouplingInformation::setGeneral(const CouplingInformation::EquationTypeId typeId, const bool isComplex, const int fieldStart, const bool isSplitEquation)
   {
      this->mEquationType = typeId;

      this->mIsComplex = isComplex;

      this->mFieldStart = fieldStart;

      this->mIsSplitEquation = isSplitEquation;
   }

   void CouplingInformation::setNonlinear(const bool hasNonlinear, const bool hasQuasiInverse)
   {
      this->mHasNonlinear = hasNonlinear;

      this->mHasQuasiInverse = hasQuasiInverse;
   }

   void CouplingInformation::setSource(const bool hasSource)
   {
      this->mHasSource = hasSource;
   }

   void CouplingInformation::setBoundaryValue(const bool hasBoundaryValue)
   {
      this->mHasBoundaryValue = hasBoundaryValue;
   }

   void CouplingInformation::setSizes(const int nSystems, const ArrayI& tauNs, const ArrayI& galerkinNs, const MatrixI& galerkinShifts, const ArrayI& rhsCols, const ArrayI& systemNs)
   {
      this->mNSystems = nSystems;

      this->mTauNs = tauNs;

      this->mGalerkinNs = galerkinNs;

      this->mGalerkinShifts = galerkinShifts;

      this->mRhsCols = rhsCols;

      this->mIsGalerkin = (this->mGalerkinShifts.sum() > 0);

      this->mSystemNs = systemNs;
   }

   void CouplingInformation::setSolverIndex(const int idx)
   {
      this->mSolverIndex = idx;
   }

   void CouplingInformation::setIndexType(const CouplingIndexType id, Tools::SharedICoupling spCoupling)
   {
      this->mIndexType = id;

      this->mspCouplingTools = spCoupling;
   }

   CouplingInformation::FieldId_range CouplingInformation::implicitRange() const
   {
      return std::make_pair(this->mImplicitFields.begin(), this->mImplicitFields.end());
   }

   CouplingInformation::FieldId_range CouplingInformation::explicitRange(const std::size_t opId) const
   {
      CouplingInformation::FieldId_range range;
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         range = std::make_pair(this->mExplicitLFields.begin(), this->mExplicitLFields.end());

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         range = std::make_pair(this->mExplicitNLFields.begin(), this->mExplicitNLFields.end());

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         range = std::make_pair(this->mExplicitNSFields.begin(), this->mExplicitNSFields.end());
      }

      return range;
   }
}
}
