/**
 * @file EnergyInnerProduct.cpp
 * @brief Source of test functor for matrix A
 */

// System includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Timestep/Exponential/Functors/EnergyFunctor.hpp"
#include "Timestep/Exponential/Functors/EnergyInnerProduct.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "Timestep/Exponential/Functors/DoNothingFunctor.hpp"
#include "Timestep/Exponential/Functors/OutputFunctor.hpp"
#include "Timestep/Exponential/EpirkTimestepper.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

EnergyInnerProduct::EnergyInnerProduct(const std::size_t regId, const std::size_t regCol)
   : mAn(0), mspData(nullptr), mRegId(regId), mRegCol(regCol), mpIdMap(nullptr), _mem(nullptr)
{
   this->mpNFunc = std::make_shared<Functors::DoNothingFunctor>();

   this->mFieldIds.push_back(std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR));
   this->mFieldIds.push_back(std::make_pair(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::TOR));
   this->mFieldIds.push_back(std::make_pair(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::POL));
}

void EnergyInnerProduct::configure(std::shared_ptr<Functors::FunctorData> spData, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem)
{
   if(!this->mspData)
   {
      auto spSetup = spData->res().spSpectralSetup();
      this->mspField = std::make_shared<Framework::Selector::ComplexScalarField>(spSetup);

      this->mpIdMap = idMap;
      this->_mem = mem;

      // Copy minimal information for spData
      this->mspData = std::make_shared<FunctorData>();
      this->mspData->spRes = spData->spRes;
      this->mspData->cInfos = spData->cInfos;
      this->mspData->stencils = spData->stencils;
      this->mspData->solups = spData->solups;

      // Set field pointers to local storage
      for(const auto& [id,p]: spData->fields)
      {
         this->mspData->fields.emplace(id, this->mspField.get());
      }

      this->mpEFunc = std::make_shared<Functors::EnergyFunctor>(spData->spRes, spData->spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL), _mem);
   }
}

Array EnergyInnerProduct::operator()(const Matrix& u, const int i0, const int i1, const Matrix& v, const int j, const int n) const
{
   assert(n == this->mAn);

   int cols = i1 - i0 + 1;
   Array dot = Array::Zero(cols);
   auto&& oFunc = *this->mpOutFunc;
   auto&& eFunc = *this->mpEFunc;

   for(auto&& id: mFieldIds)
   {
      for(auto&& qstId: qstIds(id))
      {
         // Copy data to handle
         this->mpHandle->col(this->mRegCol) = v.block(0, j, n, 1);

         // Convert component from v(:,j) to ScalarField
         oFunc(id);

         eFunc.init(this->mspField->globalView());

         // Get energy basis for v
         eFunc.convert(eFunc.mEModsView, this->mspField->globalView(), qstId);

         for(int i = i0; i <= i1; i++)
         {
            // Copy data to handle
            this->mpHandle->col(this->mRegCol) = u.block(0, i, n, 1);

            // Convert component from u(:,i) to ScalarField
            oFunc(id);

            // Compute energy inner product
            dot(i - i0) += eFunc(this->mspField->globalView(), eFunc.mEModsView, qstId);
         }
      }
   }

#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, dot.data(), dot.size(), Environment::MpiTypes::type<MHDFloat>(), MPI_SUM, MPI_COMM_WORLD);
#endif

   // Add dot product from augmented part
   int p = u.rows() - n;
   dot += u.block(n, i0, p, cols).transpose() * v.col(j).bottomRows(p);

   return dot;
}

double EnergyInnerProduct::operator()(const Matrix& u, const int i, const Matrix& v, const int j, const int n) const
{
   Array dot = this->operator()(u, i, i, v, j,n);
   assert(dot.rows() == 1);

   return dot(0,0);
}

double EnergyInnerProduct::norm(const Matrix& u, const int i, const int n) const
{
   double norm = 0.0;

   auto&& oFunc = *this->mpOutFunc;
   auto&& eFunc = *this->mpEFunc;

   for(auto&& id: mFieldIds)
   {
      // Copy data to handle
      this->mpHandle->col(this->mRegCol) = u.block(0, i, n, 1);

      // Convert component from u(:,j) to ScalarField
      oFunc(id);

      for(auto&& qstId: qstIds(id))
      {
         // Compute energy norm
         norm += eFunc(this->mspField->globalView(), qstId);
      }
   }

#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, &norm, 1, Environment::MpiTypes::type<MHDFloat>(), MPI_SUM, MPI_COMM_WORLD);
#endif

   // Add squared 2-norm of augmented part
   int p = u.rows() - n;
   norm += u.col(i).bottomRows(p).squaredNorm();

   norm = std::sqrt(norm);
   return norm;
}

void EnergyInnerProduct::setWorkspaceHandle(Matrix& mat)
{
   this->mpHandle = &mat;
   this->mAn = mat.rows();
}

void EnergyInnerProduct::setStepper(std::shared_ptr<TsFunctor> pStepper)
{
   // Timestepper wrapper
   this->mpTsFunc = pStepper;

   // Output functors
   this->mpOviewFunc = std::make_shared<OviewFunctor>(this->mspData, this->mpTsFunc, this->mRegId, this->mRegCol);
   this->mpOutFunc = std::make_shared<Functors::OutputFunctor<OviewFunctor, Functors::DoNothingFunctor>>(this->mspData, this->mpOviewFunc, this->mpNFunc, this->mpIdMap, this->_mem);
}

std::vector<SpectralFieldId> EnergyInnerProduct::qstIds(const SpectralFieldId& id) const
{
   std::vector<SpectralFieldId> qst;
   if(id.second == FieldComponents::Spectral::SCALAR)
   {
      qst.push_back(id);
   }
   else if(id.second == FieldComponents::Spectral::TOR)
   {
      qst.emplace_back(id.first, FieldComponents::Spectral::T);
   }
   else if(id.second == FieldComponents::Spectral::POL)
   {
      qst.emplace_back(id.first, FieldComponents::Spectral::Q);
      qst.emplace_back(id.first, FieldComponents::Spectral::S);
   }
   else
   {
      throw std::logic_error("Unknown conversion to QST IDs");
   }

   return qst;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
