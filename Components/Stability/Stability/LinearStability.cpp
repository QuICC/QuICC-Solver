/**
 * @file LinearStability.cpp
 * @brief Source of the high level simulation
 */

// System includes
//
#include <algorithm>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Equations/Tools/EquationTools.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "Stability/LinearStability.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace internal {
bool sortDecreasingReal(MHDComplex a, MHDComplex b)
{
   if (std::real(a) == std::real(b))
   {
      return std::imag(a) > std::imag(b);
   }
   return std::real(a) > std::real(b);
}
} // namespace internal

LinearStability::LinearStability(const std::vector<MHDFloat>& eigs,
   SharedResolution spRes,
   const Equations::EquationParameters::NDMapType& params,
   const std::map<std::size_t, std::size_t>& bcs,
   std::shared_ptr<Model::IModelBackend> spModel) :
    mcUseMumps(true),
    mNeedInit(true),
    mEigs(eigs),
    mspRes(spRes),
    mParams(params),
    mBcs(bcs),
    mspModel(spModel)
{}

LinearStability::~LinearStability()
{
   // Free work space
   PetscCallVoid(EPSDestroy(&this->mEps));
   PetscCallVoid(MatDestroy(&this->mA));
   PetscCallVoid(MatDestroy(&this->mB));
}

void LinearStability::buildMatrices(SparseMatrixZ& matA, SparseMatrixZ& matB,
   SparseMatrixZ& matC, const std::vector<MHDFloat>& eigs,
   const Equations::EquationParameters::NDMapType& nds)
{
   // Fields
   auto fId = std::make_pair(PhysicalNames::Velocity::id(),
      FieldComponents::Spectral::TOR);

   const int matIdx = 0;
   const auto& res = *this->mspRes;
   Model::EquationInfo eqInfo;
   this->model().equationInfo(eqInfo, fId, res);

   auto imRange = std::make_pair(eqInfo.im.begin(), eqInfo.im.end());

   // Build matrix A (linear operator)
   auto opId = ModelOperator::ImplicitLinear::id();
   auto bcType = ModelOperatorBoundary::SolverNoTau::id();
   DecoupledZSparse matT;
   this->model().modelMatrix(matT, opId, imRange, matIdx, bcType, res, eigs,
      this->mBcs, nds);
   matA = matT.real().cast<MHDComplex>() + matT.imag() * Math::cI;

   // Build matrix B (mass matrix)
   opId = ModelOperator::Time::id();
   matT.setZero();
   this->model().modelMatrix(matT, opId, imRange, matIdx, bcType, res, eigs,
      this->mBcs, nds);
   matB = matT.real().cast<MHDComplex>() + matT.imag() * Math::cI;

   // Build boundary matrix if needed
   if (this->model().useGalerkin())
   {
      matC.resize(0, 0);
   }
   else
   {
      opId = ModelOperator::Boundary::id();
      bcType = ModelOperatorBoundary::SolverHasBc::id();
      matT.setZero();
      this->model().modelMatrix(matT, opId, imRange, matIdx, bcType, res, eigs,
         this->mBcs, nds);
      matC = matT.real().cast<MHDComplex>() + matT.imag() * Math::cI;
   }
}

void LinearStability::eigenvalues(std::vector<MHDComplex>& evs, const int nev,
   const MHDFloat Ra)
{
   // Update rayleigh number
   this->mParams[NonDimensional::Rayleigh::id()] =
      std::make_shared<NonDimensional::Rayleigh>(Ra);

   SparseMatrixZ matA;
   SparseMatrixZ matB;
   SparseMatrixZ matC;
   this->buildMatrices(matA, matB, matC, this->mEigs, this->mParams);

   // Add tau lines if needed
   if (matC.size() != 0)
   {
      matA += matC;
   }

   // Convert matrices for SLEPc/PETSC
   this->convertMatrices(matA, matB);

   // Solve GEVP
   this->solveGEVP(evs, nev);

   // Sort in decreasing real part order
   std::sort(evs.begin(), evs.end(), internal::sortDecreasingReal);

   this->mNeedInit = false;

#ifdef QUICC_STABILITY_VERBOSE
   // print details results
   this->printDetails();
#endif
}

MHDFloat LinearStability::operator()(const MHDFloat Ra)
{
   const int nev = 5;
   std::vector<MHDComplex> evs(nev);
   this->eigenvalues(evs, nev, Ra);

   return evs.at(0).real();
}

MHDFloat LinearStability::operator()(const MHDFloat Ra, std::vector<MHDComplex>& evs)
{
   this->eigenvalues(evs, evs.size(), Ra);

   return evs.at(0).real();
}

void LinearStability::convertMatrices(const SparseMatrixZ& matA,
   const SparseMatrixZ& matB)
{
   PetscFunctionBeginUser;

   InsertMode mode;
   // Build PETSc matrix A
   if (this->mNeedInit)
   {
      PetscInt rows = matA.rows();
      PetscInt cols = matA.cols();
      PetscInt nnz = matA.nonZeros();
      PetscCallVoid(
         MatCreateSeqAIJ(PETSC_COMM_WORLD, rows, cols, nnz, NULL, &this->mA));
      mode = INSERT_VALUES;
   }
   else
   {
      PetscCallVoid(MatZeroEntries(this->mA));
      mode = ADD_VALUES;
   }
   for (int k = 0; k < matA.outerSize(); ++k)
   {
      for (SparseMatrixZ::InnerIterator it(matA, k); it; ++it)
      {
         PetscInt i = it.row();
         PetscInt j = it.col();
         if (it.value() == 0.0)
         {
            std::cerr << "WARNING: Matrix has explicit zero!" << std::endl;
         }
         PetscCallVoid(MatSetValues(this->mA, 1, &i, 1, &j, &it.value(), mode));
      }
   }
   PetscCallVoid(MatAssemblyBegin(this->mA, MAT_FINAL_ASSEMBLY));
   PetscCallVoid(MatAssemblyEnd(this->mA, MAT_FINAL_ASSEMBLY));

   // Build PETSc matrix B
   if (this->mNeedInit)
   {
      PetscInt rows = matB.rows();
      PetscInt cols = matB.cols();
      PetscInt nnz = matB.nonZeros();
      PetscCallVoid(
         MatCreateSeqAIJ(PETSC_COMM_WORLD, rows, cols, nnz, NULL, &this->mB));
   }
   else
   {
      PetscCallVoid(MatZeroEntries(this->mB));
      mode = ADD_VALUES;
   }

   for (int k = 0; k < matB.outerSize(); ++k)
   {
      for (SparseMatrixZ::InnerIterator it(matB, k); it; ++it)
      {
         PetscInt i = it.row();
         PetscInt j = it.col();
         PetscCallVoid(MatSetValues(this->mB, 1, &i, 1, &j, &it.value(), mode));
      }
   }
   PetscCallVoid(MatAssemblyBegin(this->mB, MAT_FINAL_ASSEMBLY));
   PetscCallVoid(MatAssemblyEnd(this->mB, MAT_FINAL_ASSEMBLY));
}

void LinearStability::solveGEVP(std::vector<MHDComplex>& evs, const int nev)
{
   ST st;

   PetscFunctionBeginUser;

   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Create the eigensolver and set various options
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   /*
      Create eigensolver context
      */
   if (this->mNeedInit)
   {
      PetscCallVoid(EPSCreate(PETSC_COMM_WORLD, &this->mEps));
   }

   /*
      Set operators. In this case, it is a generalized eigenvalue problem
      */
   PetscCallVoid(EPSSetOperators(this->mEps, this->mA, this->mB));

   PetscCallVoid(EPSGetST(this->mEps, &st));
   PetscCallVoid(STSetType(st, STSINVERT));

   // Use MUMPS
   if (this->mcUseMumps)
   {
      KSP ksp;
      PC pc;
      PetscCallVoid(STGetKSP(st, &ksp));
      PetscCallVoid(KSPSetType(ksp, KSPPREONLY));
      PetscCallVoid(KSPGetPC(ksp, &pc));
      PetscCallVoid(PCSetType(pc, PCLU));
      PetscCallVoid(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
      // next line is required to force the creation of the ST operator and its
      // passing to KSP */
      PetscCallVoid(STGetOperator(st, NULL));
      PetscCallVoid(PCFactorSetUpMatSolverType(pc));
      // Example to show how to pass additional options to Mumps solver:
      // Mat K;
      // PetscCallVoid(PCFactorGetMatrix(pc,&K));
      // PetscCallVoid(MatMumpsSetIcntl(K,14,50)); // Memory increase
      // PetscCallVoid(MatMumpsSetCntl(K,3,1e-12)); // Zero pivot detection
   }

   PetscCallVoid(
      EPSSetDimensions(this->mEps, nev, PETSC_DEFAULT, PETSC_DEFAULT));
   PetscCallVoid(EPSSetTarget(this->mEps, 0.0));
   PetscCallVoid(EPSSetWhichEigenpairs(this->mEps, EPS_TARGET_REAL));

   // Solve eigensystem
   PetscCallVoid(EPSSolve(this->mEps));

   evs.resize(nev);
   MHDComplex tmp;
   for (int i = 0; i < nev; i++)
   {
      PetscCallVoid(EPSGetEigenvalue(this->mEps, i, &evs.at(i), &tmp));
   }
}

void LinearStability::printDetails()
{
   ST st;
   KSP ksp;
   PC pc;
   EPSType type;
   STType sttype;
   KSPType ksptype;
   PCType pctype;
   MatSolverType pcsolvertype;
   PetscInt its, lits, nev, maxit;
   PetscReal tol;

   PetscFunctionBeginUser;

   // Optional: Get some information from the solver and display it
   PetscCallVoid(EPSGetIterationNumber(this->mEps, &its));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD,
      " Number of iterations of the method: %" PetscInt_FMT "\n", its));
   PetscCallVoid(EPSGetST(this->mEps, &st));
   PetscCallVoid(STGetKSP(st, &ksp));
   PetscCallVoid(KSPGetTotalIterations(ksp, &lits));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD,
      " Number of linear iterations of the method: %" PetscInt_FMT "\n", lits));
   PetscCallVoid(KSPGetPC(ksp, &pc));
   PetscCallVoid(EPSGetType(this->mEps, &type));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n", type));
   PetscCallVoid(STGetType(st, &sttype));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD, " ST method: %s\n", sttype));
   PetscCallVoid(KSPGetType(ksp, &ksptype));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD, " KSP method: %s\n", ksptype));
   PetscCallVoid(PCGetType(pc, &pctype));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD, " PC method: %s\n", pctype));
   PetscCallVoid(PCFactorGetMatSolverType(pc, &pcsolvertype));
   PetscCallVoid(
      PetscPrintf(PETSC_COMM_WORLD, " PC Solver: %s\n\n", pcsolvertype));
   PetscCallVoid(EPSGetDimensions(this->mEps, &nev, NULL, NULL));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD,
      " Number of requested eigenvalues: %" PetscInt_FMT "\n", nev));
   PetscCallVoid(EPSGetTolerances(this->mEps, &tol, &maxit));
   PetscCallVoid(PetscPrintf(PETSC_COMM_WORLD,
      " Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n", (double)tol,
      maxit));

   // Show detailed info
   PetscCallVoid(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,
      PETSC_VIEWER_ASCII_INFO_DETAIL));
   PetscCallVoid(EPSConvergedReasonView(this->mEps, PETSC_VIEWER_STDOUT_WORLD));
   PetscCallVoid(
      EPSErrorView(this->mEps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD));
   PetscCallVoid(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
}

const Model::IModelBackend& LinearStability::model() const
{
   return *this->mspModel;
}

} // namespace QuICC
