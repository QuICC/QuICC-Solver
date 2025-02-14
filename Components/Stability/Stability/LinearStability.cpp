/**
 * @file LinearStability.cpp
 * @brief Source of the high level simulation
 */

// System includes
//
#include <algorithm>
#include <limits>
#include <random>

#include "QuICC/Equations/EquationParameters.hpp"
#include <unsupported/Eigen/SparseExtra>

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
#include "QuICC/NonDimensional/GrowthRate.hpp"
#include "QuICC/NonDimensional/Sort.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "Stability/LinearStability.hpp"
#include "Stability/Options.hpp"
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

LinearStability::LinearStability(const std::vector<MHDFloat>& eigs, SharedResolution spRes,
   const Equations::EquationParameters::NDMapType& params,
   const std::map<std::size_t, std::size_t>& bcs,
   std::shared_ptr<Model::IModelBackend> spModel, std::shared_ptr<const Stability::Options> opt) :
    mcUseMumps(true),
    mNeedInit(true),
    mIdc(0),
    mEigs(eigs),
    mspRes(spRes),
    mParams(params),
    mBcs(bcs),
    mspModel(spModel),
    mTarget(0.0),
    mOptions(opt)
{}

LinearStability::~LinearStability()
{
   // Free work space
   PetscCallVoid(EPSDestroy(&this->mEps));
   PetscCallVoid(MatDestroy(&this->mA));
   PetscCallVoid(MatDestroy(&this->mB));
}

const Stability::Options& LinearStability::options() const
{
   return *this->mOptions;
}

void LinearStability::setCriticalId(const std::size_t idc)
{
   this->mIdc = idc;
}

void LinearStability::buildMatrices(DecoupledZSparse& matA, DecoupledZSparse& matB,
   const std::vector<MHDFloat>& eigs,
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
   this->model().modelMatrix(matA, opId, imRange, matIdx, bcType, res, eigs,
      this->mBcs, nds);
   if(this->options().writeMtx)
   {
      Eigen::saveMarket(matA.real(), "A_re.mtx");
      Eigen::saveMarket(matA.imag(), "A_im.mtx");
   }

   // Build matrix B (mass matrix)
   opId = ModelOperator::Time::id();
   this->model().modelMatrix(matB, opId, imRange, matIdx, bcType, res, eigs,
      this->mBcs, nds);
   if(this->options().writeMtx)
   {
      Eigen::saveMarket(matB.real(), "B_re.mtx");
      Eigen::saveMarket(matB.imag(), "B_im.mtx");
   }

   // Build boundary matrix if needed
   if (!this->model().useGalerkin())
   {
      DecoupledZSparse matC;
      opId = ModelOperator::Boundary::id();
      bcType = ModelOperatorBoundary::SolverHasBc::id();
      this->model().modelMatrix(matC, opId, imRange, matIdx, bcType, res, eigs,
         this->mBcs, nds);
      if(this->options().writeMtx)
      {
         Eigen::saveMarket(matC.real(), "C_re.mtx");
         Eigen::saveMarket(matC.imag(), "C_im.mtx");
      }

      // Add BC to matA}
      matA.real() += matC.real();
      matA.imag() += matC.imag();
   }

   // Set target
   if (nds.count(NonDimensional::Omega::id()) > 0 || nds.count(NonDimensional::GrowthRate::id()) > 0)
   {
      MHDFloat re = 0.0;
      if(nds.count(NonDimensional::GrowthRate::id()) > 0)
      {
         re = nds.at(NonDimensional::GrowthRate::id())->value();
      }
      MHDFloat im = 0.0;
      if(nds.count(NonDimensional::Omega::id()) > 0)
      {
         im = nds.at(NonDimensional::Omega::id())->value();
      }

      this->mTarget = MHDComplex(re, im);
   }
}

void LinearStability::castMatrices(SparseMatrixZ& matA, SparseMatrixZ& matB, const DecoupledZSparse& decA, const DecoupledZSparse& decB)
{
   if (decA.imag().size() > 0)
   {
      matA = decA.real().cast<MHDComplex>() + decA.imag() * Math::cI;
   }
   else
   {
      matA = decA.real().cast<MHDComplex>();
   }

   if (decB.imag().size() > 0)
   {
      matB = decB.real().cast<MHDComplex>() + decB.imag() * Math::cI;
   }
   else
   {
      matB = decB.real().cast<MHDComplex>();
   }
}

std::pair<int, int> LinearStability::setupGEVP(const MHDFloat vc)
{
   if(this->mIdc != 0)
   {
      // Update critical parameter
      this->mParams[this->mIdc] = std::make_shared<NonDimensional::INumber>(vc,
         this->mParams[this->mIdc]->tag());
   }

   DecoupledZSparse decA;
   DecoupledZSparse decB;
   this->buildMatrices(decA, decB, this->mEigs, this->mParams);
   SparseMatrixZ matA;
   SparseMatrixZ matB;
   this->castMatrices(matA, matB, decA, decB);

   // Convert matrices for SLEPc/PETSC
   this->convertMatrices(matA, matB);

   auto dims = std::make_pair(matA.rows(), matA.cols());

   std::cerr << "Finshed setting up matrices. Starting solver..." << std::endl;

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

   if(this->options().verboseDiagnostics)
   {
      // print details results
      this->printDetails();
   }

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

   // Allocate PETSc matrix
   auto allocatePetscMat = [](auto& petscMat, const auto& eigenMat)
   {
      PetscInt rows = eigenMat.rows();
      PetscInt cols = eigenMat.cols();
      PetscInt tnz = 0;
      std::vector<PetscInt> nnz(rows, 0);
      for (int k = 0; k < eigenMat.outerSize(); ++k)
      {
         for (SparseMatrixZ::InnerIterator it(eigenMat, k); it; ++it)
         {
            ++nnz.at(it.row());
            ++tnz;
         }
      }
      if(tnz != eigenMat.nonZeros())
      {
         throw std::logic_error("Counting NNZ per row failed");
      }

      PetscCallVoid(
         MatCreateSeqAIJ(PETSC_COMM_WORLD, rows, cols, tnz, nnz.data(), &petscMat));
   };

   // Set PETSc matrix values
   auto setPetscMat = [](auto& petscMat, const auto& eigenMat, const auto& mode)
   {
      for (int k = 0; k < eigenMat.outerSize(); ++k)
      {
         for (SparseMatrixZ::InnerIterator it(eigenMat, k); it; ++it)
         {
            PetscInt i = it.row();
            PetscInt j = it.col();
            if (it.value() == 0.0)
            {
               std::cerr << "WARNING: Matrix has explicit zero!" << std::endl;
            }
            PetscCallVoid(MatSetValues(petscMat, 1, &i, 1, &j, &it.value(), mode));
         }
      }
      PetscCallVoid(MatAssemblyBegin(petscMat, MAT_FINAL_ASSEMBLY));
      PetscCallVoid(MatAssemblyEnd(petscMat, MAT_FINAL_ASSEMBLY));
   };

   InsertMode mode;
   // Build PETSc matrix A
   if (this->mNeedInit)
   {
      allocatePetscMat(this->mA, matA);
      mode = INSERT_VALUES;
   }
   else
   {
      PetscCallVoid(MatZeroEntries(this->mA));
      mode = ADD_VALUES;
   }
   setPetscMat(this->mA, matA, mode);

   // Build PETSc matrix B
   if (this->mNeedInit)
   {
      allocatePetscMat(this->mB, matB);
      mode = INSERT_VALUES;
   }
   else
   {
      PetscCallVoid(MatZeroEntries(this->mB));
      mode = ADD_VALUES;
   }
   setPetscMat(this->mB, matB, mode);
}

void LinearStability::setCustomGuess()
{
   if(this->options().guessType == 0)
   {
      std::vector<int> parity = {0, 1, 1, 0};
      this->setParityGuess(parity);
   }
   else if(this->options().guessType == 1)
   {
      std::vector<int> parity = {1, 0, 0, 1};
      this->setParityGuess(parity);
   }
   else
   {
      throw std::logic_error("Unknown initial guess type");
   }
}

void LinearStability::setParityGuess(const std::vector<int>& parity)
{
   const auto& res = *this->mspRes;
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);

   Vec guess;
   PetscCallVoid(MatCreateVecs(this->mA, &guess, nullptr));

   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<> dis(-1.0, 1.0);
   std::size_t idx = 0;
   int m = this->mspRes->sim().dim(Dimensions::Simulation::SIM3D,
         Dimensions::Space::SPECTRAL) -
      1;
   int p;
   PetscScalar val;
   // Loop over fields
   for(int c = 0; c < parity.size(); c++)
   {
      p = parity.at(c);
      for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
      {
         int k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
         if (k_ == m)
         {
            for (int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k);
                  j++)
            {
               int j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
               if(j_ != 0 && j % 2 == p)
               {
                  val = 1.0;
               }
               else
               {
                  val = 0.0;
               }
               for (int i = 0;
                     i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
               {
                  PetscCallVoid(VecSetValue(guess,idx,val*dis(gen),INSERT_VALUES));
                  idx++;
               }
            }
         }
      }
   }
   PetscCallVoid(EPSSetInitialSpace(this->mEps,1,&guess));
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
   PetscCallVoid(EPSSetTolerances(this->mEps, this->options().tolerance, this->options().maxIteration));
   PetscCallVoid(EPSSetBalance(this->mEps,  EPS_BALANCE_TWOSIDE, PETSC_DETERMINE, PETSC_DETERMINE));

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

   // Use custom initial guess
   if(this->options().useCustomGuess)
   {
      this->setCustomGuess();
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
