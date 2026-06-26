/**
 * @file AugmentedJacobianFunctor.cpp
 * @brief Source of test functor for matrix A
 */

// System includes
//

// Project includes
//
#include "QuICC/Register/Intermediate.hpp"
#include "QuICC/Register/Temporary.hpp"
#include "Timestep/Exponential/AugmentedJacobianFunctor.hpp"
#include "Timestep/Exponential/Functors/FunctorData.hpp"
#include "Timestep/Exponential/Functors/DoNothingFunctor.hpp"
#include "Timestep/Exponential/Functors/ProcessRangeFunctor.hpp"
#include "Timestep/Exponential/Functors/InputFunctor.hpp"
#include "Timestep/Exponential/Functors/OutputFunctor.hpp"
#include "Timestep/Exponential/Functors/ApplyConstraintFunctor.hpp"
#include "Timestep/Exponential/EpirkTimestepper.hpp"
#include "Timestep/Exponential/Functors/CallExplicitPrognosticFunctor.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Timestep {

namespace Exponential {

AugmentedJacobianFunctor::AugmentedJacobianFunctor(std::shared_ptr<Functors::FunctorData> spData, const MHDFloat dt, const std::size_t regId, const std::size_t regCol, const std::set<int>& fixedIts, Pseudospectral::Coordinator* pPseudo, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem)
   : mcFixedIts(fixedIts), mcEps(1e-8), mAn(0), mBn(0), mN(0), mspData(spData), mDt(dt), mRegId(regId), mRegCol(regCol), mpHandle(nullptr), mpPseudo(pPseudo), matB(0,0), mpIdMap(idMap), _mem(mem)
{
   this->mpNFunc = std::make_shared<Functors::DoNothingFunctor>();

}

void AugmentedJacobianFunctor::updateTimestep(const MHDFloat dt)
{
   this->mDt = dt;
}

void AugmentedJacobianFunctor::operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in) const
{
   if(QuICCEnv().allowsIO())
   {
      std::cerr << "\t" << "- compute Jacobian" << std::endl;
   }

   DebuggerMacro_msg("Applying augmented Jacobian", 2);

   assert(in.rows() == this->mN);
   assert(out.rows() == this->mN);
   assert(in.cols() == 1);
   assert(out.cols() == in.cols());
   assert(this->mpHandle->rows() == this->mAn);

   if(this->mAn == 0 || this->mBn == 0)
   {
      throw std::logic_error("Operators have not been initialized");
   }

   // Copy data into handle
   this->mpHandle->col(0) = (this->mcEps*this->mDt)*in.topRows(this->mAn);

   this->applyJacobian();
   
   // Copy data from handle
   out.topRows(this->mAn) = (1./this->mcEps)*this->mpHandle->col(0);

   // Add part from augmented matrix
   DebuggerMacro_msg("Applying B matrix from augmented Jacobian", 3);
   out.topRows(this->mAn) += this->matB * in.bottomRows(this->mBn);
   out.block(this->mAn, 0, this->mBn-1, out.cols()) = in.block(this->mAn + 1, 0, this->mBn-1, in.cols());
   out.bottomRows(1).array() = 0;
}

void AugmentedJacobianFunctor::updateB(const Matrix& matB)
{
   this->matB = matB;
   this->mBn = matB.cols();
   this->mN = this->mAn + this->mBn;
}

void AugmentedJacobianFunctor::applyJacobian() const
{
   int lastIt = *this->mcFixedIts.rbegin();

   DebuggerMacro_msg("Applying Jacobian to stepper handle", 3);

   // Transfer timestep output back to equations
   Functors::ProcessRangeFunctor processOut(this->mpNFunc, this->mpOutFunc, this->mpNFunc, lastIt);
   processOut(this->mspData->eqInfos);

   // Clear RHS
   this->mpHandle->setZero();

   this->mpPseudo->evolveAfterPrognostic(this->mcFixedIts, false);

   auto progFunc = std::make_shared<Functors::CallExplicitPrognosticFunctor<TsFunctor>>(this->mspData, this->mpTsFunc, this->mRegId, this->mRegCol, lastIt, this->mpIdMap, this->_mem);
   this->mpPseudo->evolveUntilPrognostic(this->mcFixedIts, false, progFunc);

   // Update the equation input to the timestepper
   Functors::ProcessRangeFunctor processIn(this->mpIbefFunc, this->mpInFunc, this->mpNFunc, lastIt);
   processIn(this->mspData->eqInfos);
}

void AugmentedJacobianFunctor::setStepper(std::shared_ptr<TsFunctor> pStepper)
{
   // Timestepper wrapper
   this->mpTsFunc = pStepper;

   // Output functors
   this->mpOviewFunc = std::make_shared<OviewFunctor>(this->mspData, this->mpTsFunc, this->mRegId, this->mRegCol);
   this->mpOcorrFunc = std::make_shared<OcorrFunctor>(this->mspData, this->mpTsFunc, this->mRegId, this->mRegCol);
   this->mpOutFunc = std::make_shared<Functors::OutputFunctor<OviewFunctor, OcorrFunctor>>(this->mspData, this->mpOviewFunc, this->mpOcorrFunc, this->mpIdMap, this->_mem);

   // Input functors
   this->mpIbefFunc = std::make_shared<Functors::ApplyConstraintFunctor>(this->mspData, SolveTiming::Before::id());
   this->mpIviewFunc = std::make_shared<IviewFunctor>(this->mspData, this->mpTsFunc, this->mRegId, this->mRegCol);
   this->mpInFunc = std::make_shared<Functors::InputFunctor<IviewFunctor>>(this->mspData, this->mpIviewFunc, this->mpIdMap, this->_mem);
}

void AugmentedJacobianFunctor::setMatrixHandle(Matrix& mat)
{
   this->mpHandle = &mat;
   this->mAn = mat.rows();
   this->mN = this->mAn + this->mBn;
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
