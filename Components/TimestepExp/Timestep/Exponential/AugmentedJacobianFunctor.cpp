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
#include "Timestep/Exponential/InterfaceFunctors.hpp"
#include "Timestep/Exponential/EpirkTimestepper.hpp"
#include "Timestep/Exponential/Functors/CallExplicitPrognosticFunctor.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Timestep {

namespace Exponential {

AugmentedJacobianFunctor::AugmentedJacobianFunctor(const MHDFloat dt, const std::size_t regId, const std::size_t regCol, const int fixedIt, Pseudospectral::Coordinator* pPseudo, std::shared_ptr<IdMap> idMap, std::shared_ptr<Memory::memory_resource> mem)
   : mcFixedIt(fixedIt), mcEps(1e-8), mAn(0), mBn(0), mN(0), mDt(dt), mRegId(regId), mRegCol(regCol), mpHandle(nullptr), mpPseudo(pPseudo), matB(0,0), mpIdMap(idMap), _mem(mem), mpScalEq(nullptr), mpVectEq(nullptr)
{
   this->mpNFunc = std::make_shared<DoNothingFunctor>();

}

void AugmentedJacobianFunctor::updateTimestep(const MHDFloat dt)
{
   this->mDt = dt;
}

void AugmentedJacobianFunctor::operator()(Eigen::Ref<Matrix> out, Eigen::Ref<Matrix> in) const
{
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
   std::set<int> itIds = {this->mcFixedIt};

   DebuggerMacro_msg("Applying Jacobian to stepper handle", 3);
   assert(this->mpScalEq);
   auto scalEq = *mpScalEq;
   assert(this->mpVectEq);
   auto vectEq = *mpVectEq;

   // Transfer timestep output back to equations
   ProcessRangeFunctor processOut(this->mpNFunc, this->mpOutFunc, this->mpNFunc, this->mcFixedIt);
   processOut(scalEq);
   processOut(vectEq);

   // Clear RHS
   this->mpHandle->setZero();

   this->mpPseudo->evolveAfterPrognostic(itIds, false);

   auto progFunc = std::make_shared<Functors::CallExplicitPrognosticFunctor<TsFunctor>>(this->mpTsFunc, this->mRegId, this->mRegCol, this->mcFixedIt, this->mpIdMap, this->_mem);
   this->mpPseudo->evolveUntilPrognostic(itIds, false, progFunc);

   // Update the equation input to the timestepper
   ProcessRangeFunctor processIn(this->mpIbefFunc, this->mpInFunc, this->mpNFunc, this->mcFixedIt);
   processIn(scalEq);
   processIn(vectEq);
}

void AugmentedJacobianFunctor::setStepper(std::shared_ptr<TsFunctor> pStepper)
{
   // Timestepper wrapper
   this->mpTsFunc = pStepper;

   // Output functors
   this->mpOviewFunc = std::make_shared<OviewFunctor>(this->mpTsFunc, this->mRegId, this->mRegCol);
   this->mpOcorrFunc = std::make_shared<OcorrFunctor>(this->mpTsFunc, this->mRegId, this->mRegCol);
   this->mpOutFunc = std::make_shared<OutputFunctor<OviewFunctor, OcorrFunctor>>(this->mpOviewFunc, this->mpOcorrFunc, this->mpIdMap, this->_mem);

   // Input functors
   this->mpIbefFunc = std::make_shared<ApplyConstraintFunctor>(SolveTiming::Before::id());
   this->mpIviewFunc = std::make_shared<IviewFunctor>(this->mpTsFunc, this->mRegId, this->mRegCol);
   this->mpInFunc = std::make_shared<InputFunctor<IviewFunctor>>(this->mpIviewFunc, this->mpIdMap, this->_mem);
}

void AugmentedJacobianFunctor::setEquations(const Timestep::Interface::ScalarEquation_range& scalEq, const Timestep::Interface::VectorEquation_range& vectEq)
{
   this->mpScalEq = &scalEq;
   this->mpVectEq = &vectEq;
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
