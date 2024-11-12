/**
 * @file LinearStability.hpp
 * @brief Functor to compute linear stability
 */

#ifndef QUICC_LINEARSTABILITY_HPP
#define QUICC_LINEARSTABILITY_HPP

// System includes
//
#include <memory>
#include <slepceps.h>

// Project includes
//
#include "QuICC/Equations/EquationParameters.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

/**
 * @brief Functor to compute linear stability
 */
class LinearStability
{
public:
   /**
    * @brief Constructor
    *
    * @param eigs    Indexes of matrix to solve (eg. m for rotating spherical setup)
    * @param spRes   Resolution
    * @param params  Nondimensional parameters
    * @param bcMap   Boundary conditions
    * @param spModel Model backend
    */
   LinearStability(const std::vector<MHDFloat>& eigs, SharedResolution spRes,
      const Equations::EquationParameters::NDMapType& params,
      const std::map<std::size_t, std::size_t>& bcMap,
      std::shared_ptr<Model::IModelBackend> spModel);

   /**
    * @brief Simple empty destructor
    */
   virtual ~LinearStability();

   /**
    * @brief Compute growth rate for given Rayleigh number
    */
   MHDFloat operator()(const MHDFloat Ra);

   /**
    * @brief Compute growth rate for given Rayleigh number
    *
    * @param evs  Storage for returning eigenvalues
    */
   MHDFloat operator()(const MHDFloat Ra, std::vector<MHDComplex>& evs);

   /**
    * @brief Compute eigenvalues
    *
    * @param evs  Eigenvalues
    * @param nev  Number of eigenvalues
    */
   void eigenvalues(std::vector<MHDComplex>& evs, const int nev,
      const MHDFloat Ra);

protected:
   /**
    * @brief Get resolution
    */
   const Resolution& res() const;

   /**
    * @brief Get model backend
    */
   const Model::IModelBackend& model() const;

private:
   /**
    * @brief Solve Generalized eigenvalue problem (GEVP)
    *
    * @param matA Linear operator
    * @param matB Mass matrix
    * @param matC Boundary condition matrix
    * @param eigs Indexes for matrix to solve
    * @param nds  Nondimensional parameters
    */
   void buildMatrices(SparseMatrixZ& matA, SparseMatrixZ& matB,
      SparseMatrixZ& matC, const std::vector<MHDFloat>& eigs,
      const Equations::EquationParameters::NDMapType& nds);

   /**
    * @brief Convert matrices to used with SLEPc/PETSc
    *
    * @param matA Linear operator
    * @param matB Mass matrix
    */
   void convertMatrices(const SparseMatrixZ& matA, const SparseMatrixZ& matB);

   /**
    * @brief Solve Generalized eigenvalue problem (GEVP)
    *
    * @param evs  Eigenvalues
    * @param nev  Number of eigenvalues
    */
   void solveGEVP(std::vector<MHDComplex>& evs, const int nev);

   /**
    * @brief Solve Generalized eigenvalue problem (GEVP)
    *
    * @param matA Linear operator
    * @param matB Mass matrix
    */
   void buildMatrix(const SparseMatrixZ& matA, const SparseMatrixZ& matB);

   /**
    * @brief Print solver details
    */
   void printDetails();

   /**
    * @brief Use MUMPS solver
    */
   const bool mcUseMumps;

   /**
    * @brief Need initialization?
    */
   bool mNeedInit;

   /**
    * @brief Indexes for independent dimension(s)
    *
    * example: m for spherical rtc
    */
   std::vector<MHDFloat> mEigs;

   /**
    * @brief Shared resolution
    */
   SharedResolution mspRes;

   /**
    * @brief Shared Equation parameters
    */
   Equations::EquationParameters::NDMapType mParams;

   /**
    * @brief Boundary condition map
    */
   std::map<std::size_t, std::size_t> mBcs;

   /**
    * @brief model backend
    */
   std::shared_ptr<Model::IModelBackend> mspModel;

   /**
    * @brief SLEPc/PETSc matrix A
    */
   Mat mA;

   /**
    * @brief SLEPc/PETSc matrix A
    */
   Mat mB;

   /**
    * @brief SLEPc/PETSc EPS object
    */
   EPS mEps;
};

/// Typedef for a shared pointer of a Simulation
typedef std::shared_ptr<LinearStability> SharedLinearStability;
} // namespace QuICC

#endif // QUICC_LINEARSTABILITY_HPP
