/**
 * @file EquationData.cpp
 * @brief Source of building block for the implementation of a time dependend evolution equation
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Equations/EquationData.hpp"

// Project includes
//
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"

namespace QuICC {

namespace Equations {

   EquationData::EquationData(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend)
      : mForwardPathsType(FWD_IS_CUSTOM), mspEqParams(spEqParams), mspSpatialScheme(spScheme), mspBackend(spBackend), mTime(-1.0)
   {
   }

   EquationData::~EquationData()
   {
   }

   void EquationData::setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 0);

      this->mScalars.insert(std::make_pair(name, spField));
   }

   void EquationData::setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 0);

      this->mVectors.insert(std::make_pair(name, spField));
   }

   Framework::Selector::VariantSharedScalarVariable EquationData::spScalar(std::size_t name) const
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 1);

      return this->mScalars.find(name)->second;
   }

   Framework::Selector::VariantSharedVectorVariable EquationData::spVector(std::size_t name) const
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 1);

      return this->mVectors.find(name)->second;
   }

   const EquationParameters& EquationData::eqParams() const
   {
      return *this->mspEqParams;
   }

   const SharedCEquationParameters EquationData::spEqParams() const
   {
      return this->mspEqParams;
   }

   const SimulationBoundary& EquationData::bcIds() const
   {
      return *this->mspBcIds;
   }

   const SpatialScheme::ISpatialScheme&  EquationData::ss() const
   {
      // Safety assert
      assert(this->mspSpatialScheme);

      return *this->mspSpatialScheme;
   }

   const Model::IModelBackend& EquationData::backend() const
   {
      return *this->mspBackend;
   }

   void EquationData::cleanupBackend()
   {
      this->mspBackend.reset();
   }

   const SparseMatrix& EquationData::galerkinStencil(const FieldComponents::Spectral::Id compId, const int j) const
   {
      // Safety assert
      assert(this->mGStencils.count(compId) > 0);

      return this->mGStencils.find(compId)->second.at(j);
   }

   bool EquationData::hasQID(const FieldComponents::Spectral::Id compId) const
   {
      return (this->mQIDMatrices.count(compId) > 0);
   }

   bool EquationData::hasQIZ(const FieldComponents::Spectral::Id compId) const
   {
      return (this->mQIZMatrices.count(compId) > 0);
   }

   bool EquationData::hasExplicitDTerm(const std::size_t opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::ExplicitLinear::id())
      {
         return (this->mELDMatrices.count(key) > 0);

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         return (this->mENLDMatrices.count(key) > 0);

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         return (this->mENSDMatrices.count(key) > 0);

      } else
      {
         throw std::logic_error("Requested inexistant explicit matrices");
      }
   }

   bool EquationData::hasExplicitZTerm(const std::size_t opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::ExplicitLinear::id())
      {
         return (this->mELZMatrices.count(key) > 0);

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         return (this->mENLZMatrices.count(key) > 0);

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         return (this->mENSZMatrices.count(key) > 0);

      } else
      {
         throw std::logic_error("Requested inexistant explicit matrices");
      }
   }

   template <> const SparseMatrix& EquationData::quasiInverse<SparseMatrix>(const FieldComponents::Spectral::Id compId, const int j) const
   {
      // Safety assert
      assert(this->mQIDMatrices.count(compId) > 0);

      return this->mQIDMatrices.find(compId)->second.at(j);
   }

   template <> const SparseMatrixZ& EquationData::quasiInverse<SparseMatrixZ>(const FieldComponents::Spectral::Id compId, const int j) const
   {
      // Safety assert
      assert(this->mQIZMatrices.count(compId) > 0);

      return this->mQIZMatrices.find(compId)->second.at(j);
   }

   template <> const SparseMatrix& EquationData::explicitOperator<SparseMatrix>(const std::size_t opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::ExplicitLinear::id())
      {
         // Safety assert
         assert(this->mELDMatrices.count(key) > 0);

         return this->mELDMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         // Safety assert
         assert(this->mENLDMatrices.count(key) > 0);

         return this->mENLDMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         // Safety assert
         assert(this->mENSDMatrices.count(key) > 0);

         return this->mENSDMatrices.find(key)->second.at(j);

      } else
      {
         throw std::logic_error("Requested inexistant explicit matrices");
      }
   }

   template <> const SparseMatrixZ& EquationData::explicitOperator<SparseMatrixZ>(const std::size_t opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::ExplicitLinear::id())
      {
         // Safety assert
         assert(this->mELZMatrices.count(key) > 0);

         return this->mELZMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         // Safety assert
         assert(this->mENLZMatrices.count(key) > 0);

         return this->mENLZMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         // Safety assert
         assert(this->mENSZMatrices.count(key) > 0);

         return this->mENSZMatrices.find(key)->second.at(j);

      } else
      {
         throw std::logic_error("Requested inexistant explicit matrices");
      }
   }

   std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> >& EquationData::rEDMatrices(const std::size_t opId)
   {
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         return this->mELDMatrices;

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         return this->mENLDMatrices;

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         return this->mENSDMatrices;

      } else
      {
         throw std::logic_error("Requested inexistant explicit matrices");
      }
   }

   std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> >& EquationData::rEZMatrices(const std::size_t opId)
   {
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         return this->mELZMatrices;

      } else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         return this->mENLZMatrices;

      } else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         return this->mENSZMatrices;

      } else
      {
         throw std::logic_error("Requested inexistant explicit matrices");
      }
   }

   const CouplingInformation& EquationData::couplingInfo(const FieldComponents::Spectral::Id compId) const
   {
      // Safety assert
      assert(this->mCouplingInfos.count(compId) > 0);

      return this->mCouplingInfos.find(compId)->second;
   }

   void EquationData::setForwardPathsType(const EquationData::ForwardPathsId id)
   {
      this->mForwardPathsType = id;
   }

   void EquationData::setName(std::size_t name)
   {
      this->mName = name;
   }

   std::size_t   EquationData::name() const
   {
      return this->mName;
   }

   const VariableRequirement&  EquationData::requirements() const
   {
      return this->mRequirements;
   }

   const VariableRequirement&  EquationData::imposedRequirements() const
   {
      return this->mImposedRequirements;
   }

   const FieldRequirement&  EquationData::requirements(std::size_t id) const
   {
      return this->mRequirements.field(id);
   }

   const FieldRequirement&  EquationData::imposedRequirements(std::size_t id) const
   {
      return this->mImposedRequirements.field(id);
   }

   FieldRequirement&  EquationData::updateFieldRequirements(std::size_t id)
   {
      return this->mRequirements.rField(id);
   }

   void EquationData::setSolverIndex(const FieldComponents::Spectral::Id compId, const int idx)
   {
      // Safety assert
      assert(this->mCouplingInfos.count(compId) > 0);

      this->mCouplingInfos.find(compId)->second.setSolverIndex(idx);
   }

   void EquationData::setSolveTiming(const std::size_t timeId)
   {
      this->mSolveTiming = timeId;
   }

   std::size_t  EquationData::solveTiming() const
   {
      return this->mSolveTiming;
   }

   MHDFloat  EquationData::time() const
   {
      return this->mTime;
   }

   void  EquationData::setTime(const MHDFloat time, const bool finished)
   {
      this->mTime = time;
   }

   CouplingInformation::EquationTypeId EquationData::equationType() const
   {
      auto it = this->mCouplingInfos.cbegin();
      CouplingInformation::EquationTypeId id = it->second.equationType();
      ++it;

      for(; it != this->mCouplingInfos.cend(); ++it)
      {
         if(it->second.equationType() != id)
         {
            throw std::logic_error("Different equation types for vector field components is not supported");
         }
      }

      return id;
   }

   const std::vector<std::pair<FieldComponents::Spectral::Id,std::size_t> >& EquationData::nlComponents() const
   {
      return this->mNLComponents;
   }

   void EquationData::addNLComponent(const FieldComponents::Spectral::Id compId, const std::size_t flag)
   {
      this->mNLComponents.push_back(std::make_pair(compId,flag));
   }

   MHDFloat incrementTimeAverage(const MHDComplex avg, const MHDFloat newData, const MHDFloat time, const MHDFloat timestep)
   {
      throw std::logic_error("Setup is wrong, should not have been called!");

      return std::numeric_limits<MHDFloat>::quiet_NaN();
   }

   MHDFloat noupdateTimeAverage(const MHDComplex avg, const MHDFloat newData)
   {
      throw std::logic_error("Setup is wrong, should not have been called!");

      return std::numeric_limits<MHDFloat>::quiet_NaN();
   }

}
}
