/**
 * @file LinearStability.cpp
 * @brief Source of the high level simulation
 */

// System includes
//
#include <algorithm>
#include <limits>

#include "QuICC/Equations/EquationParameters.hpp"
#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
#include <unsupported/Eigen/SparseExtra>
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX

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
#include "QuICC/NonDimensional/Omega.hpp"
#include "QuICC/NonDimensional/Sort.hpp"
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

bool sortDecreasingRealIdx(std::pair<MHDComplex, int> a,
   std::pair<MHDComplex, int> b)
{
   if (std::real(a.first) == std::real(b.first))
   {
      return std::imag(a.first) > std::imag(b.first);
   }
   return std::real(a.first) > std::real(b.first);
}

} // namespace internal

LinearStability::LinearStability(const std::size_t idc,
   const std::vector<MHDFloat>& eigs, SharedResolution spRes,
   const Equations::EquationParameters::NDMapType& params,
   const std::map<std::size_t, std::size_t>& bcs,
   std::shared_ptr<Model::IModelBackend> spModel) :
    mcUseMumps(true),
    mNeedInit(true),
    mIdc(idc),
    mEigs(eigs),
    mspRes(spRes),
    mParams(params),
    mBcs(bcs),
    mspModel(spModel),
    mTarget(0.0)
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
#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   Eigen::saveMarket(matT.real(), "A_re.mtx");
   Eigen::saveMarket(matT.imag(), "A_im.mtx");
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   if (eqInfo.isComplex)
   {
      matA = matT.real().cast<MHDComplex>() + matT.imag() * Math::cI;
   }
   else
   {
      matA = matT.real().cast<MHDComplex>();
   }

   // Build matrix B (mass matrix)
   opId = ModelOperator::Time::id();
   matT.setZero();
   this->model().modelMatrix(matT, opId, imRange, matIdx, bcType, res, eigs,
      this->mBcs, nds);
#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   Eigen::saveMarket(matT.real(), "B_re.mtx");
   Eigen::saveMarket(matT.imag(), "B_im.mtx");
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
   if (eqInfo.isComplex)
   {
      matB = matT.real().cast<MHDComplex>() + matT.imag() * Math::cI;
   }
   else
   {
      matB = matT.real().cast<MHDComplex>();
   }

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
#ifdef QUICC_DEBUG_OUTPUT_MODEL_MATRIX
      Eigen::saveMarket(matT.real(), "C_re.mtx");
      Eigen::saveMarket(matT.imag(), "C_im.mtx");
#endif // QUICC_DEBUG_OUTPUT_MODEL_MATRIX
      if (eqInfo.isComplex)
      {
         matC = matT.real().cast<MHDComplex>() + matT.imag() * Math::cI;
      }
      else
      {
         matC = matT.real().cast<MHDComplex>();
      }
   }

   // Set target
   if (nds.count(NonDimensional::Omega::id()) > 0)
   {
      this->mTarget =
         MHDComplex(0, nds.at(NonDimensional::Omega::id())->value());
   }
}

std::pair<int, int> LinearStability::setupGEVP(const MHDFloat vc)
{
   // Update critical parameter
   this->mParams[this->mIdc] = std::make_shared<NonDimensional::INumber>(vc,
      this->mParams[this->mIdc]->tag());

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

   auto dims = std::make_pair(matA.rows(), matA.cols());

   return dims;
}

void LinearStability::eigenpairs(std::vector<MHDComplex>& evs,
   std::vector<std::vector<MHDComplex>>& efs, const int nev, const MHDFloat vc)
{
   auto dims = this->setupGEVP(vc);

   std::vector<Vec> petscEfs;
   if (efs.size() == static_cast<std::size_t>(nev))
   {
      for (int i = 0; i < nev; i++)
      {
         efs.at(i).resize(dims.first);
         petscEfs.push_back(Vec());
         PetscCallVoid(MatCreateVecs(this->mA, &petscEfs.back(), nullptr));
      }
   }
   else if (efs.size() > 0 && efs.size() != static_cast<std::size_t>(nev))
   {
      throw std::logic_error(
         "Eigenvector storage initialized with wrong size: " +
         std::to_string(efs.size()) + " vs " + std::to_string(nev));
   }

   // Solve GEVP
   this->solveGEVP(evs, petscEfs, nev);

   // Setup sorting data
   std::vector<std::pair<MHDComplex, int>> evs_idx;
   for (std::size_t i = 0; i < evs.size(); i++)
   {
      evs_idx.push_back(std::make_pair(evs.at(i), i));
   }

   // Sort eigenvalues
   switch (static_cast<int>(this->mParams[NonDimensional::Sort::id()]->value()))
   {
   case 1: {
      // Sort in decreasing real part order
      std::sort(evs_idx.begin(), evs_idx.end(),
         internal::sortDecreasingRealIdx);
   }
   }

   // Extract sorted eigenpairs
   for (std::size_t i = 0; i < evs_idx.size(); i++)
   {
      const int& i_ = evs_idx.at(i).second;
      evs.at(i) = evs_idx.at(i).first;
      if (efs.size() > 0)
      {
         const PetscScalar* val;
         PetscCallVoid(VecGetArrayRead(petscEfs.at(i_), &val));
         for (std::size_t j = 0; j < efs.at(i).size(); j++)
         {
            efs.at(i).at(j) = val[j];
         }
         PetscCallVoid(VecRestoreArrayRead(petscEfs.at(i_), &val));
      }
   }

   this->mNeedInit = false;

#ifdef QUICC_STABILITY_VERBOSE
   // print details results
   this->printDetails();
#endif

   // Destroy PETSc Vec
   for (auto& ef: petscEfs)
   {
      PetscCallVoid(VecDestroy(&ef));
   }
}

MHDFloat LinearStability::operator()(const MHDFloat vc)
{
   const int nev = 5;
   std::vector<MHDComplex> evs(nev);
   std::vector<std::vector<MHDComplex>> efs;
   this->eigenpairs(evs, efs, nev, vc);

   return evs.at(0).real();
}

MHDFloat LinearStability::operator()(const MHDFloat vc,
   std::vector<MHDComplex>& evs)
{
   std::vector<std::vector<MHDComplex>> efs;
   this->eigenpairs(evs, efs, evs.size(), vc);

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

void LinearStability::solveGEVP(std::vector<MHDComplex>& evs,
   std::vector<Vec>& efs, const int nev)
{
   ST st;
   PetscInt rows, cols;
   PetscCallVoid(MatGetSize(this->mA, &rows, &cols));
   bool useShift = (nev < rows);

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
   PetscCallVoid(EPSSetProblemType(this->mEps, EPS_GNHEP));

   if (useShift)
   {
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
         // next line is required to force the creation of the ST operator and
         // its passing to KSP */
         PetscCallVoid(STGetOperator(st, NULL));
         PetscCallVoid(PCFactorSetUpMatSolverType(pc));
         // Example to show how to pass additional options to Mumps solver:
         Mat K;
         PetscCallVoid(PCFactorGetMatrix(pc, &K));
         PetscCallVoid(MatMumpsSetIcntl(K, 14, 50)); // Memory increase
         // PetscCallVoid(MatMumpsSetCntl(K,3,1e-12)); // Zero pivot detection
      }
   }

   PetscCallVoid(
      EPSSetDimensions(this->mEps, nev, PETSC_DEFAULT, PETSC_DEFAULT));
   if (useShift)
   {
      PetscCallVoid(EPSSetTarget(this->mEps, this->mTarget));
      if (this->mTarget.imag() == 0)
      {
         PetscCallVoid(EPSSetWhichEigenpairs(this->mEps, EPS_TARGET_REAL));
      }
      else if (this->mTarget.real() == 0)
      {
         PetscCallVoid(EPSSetWhichEigenpairs(this->mEps, EPS_TARGET_IMAGINARY));
      }
      else
      {
         PetscCallVoid(EPSSetWhichEigenpairs(this->mEps, EPS_TARGET_MAGNITUDE));
      }
   }

   // Solve eigensystem
   PetscCallVoid(EPSSolve(this->mEps));
   PetscInt t;
   PetscCallVoid(EPSGetConverged(this->mEps, &t));
   int nconv = static_cast<int>(t);

   evs.resize(nev);
   if (efs.size() == 0)
   {
      for (int i = 0; i < std::min(nconv, nev); i++)
      {
         PetscCallVoid(EPSGetEigenvalue(this->mEps, i, &evs.at(i), nullptr));
      }
   }
   else
   {
      for (int i = 0; i < std::min(nconv, nev); i++)
      {
         PetscCallVoid(EPSGetEigenpair(this->mEps, i, &evs.at(i), nullptr,
            efs.at(i), nullptr));
      }
   }

   // Mark unconverged eigenvales
   for (int i = nconv; i < nev; i++)
   {
      evs.at(i) = std::numeric_limits<MHDComplex>::max();
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
