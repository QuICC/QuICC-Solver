/**
 * @file TestHelper.cpp
 * @brief Setup helpers for the Timestep tests
 */

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>
#include <limits>

// Project includes
//
#include <iostream>

#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D2.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I2Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl2.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "TestSuite/Timestep/PredictorCorrector/TestHelper.hpp"
namespace QuICC {

namespace TestSuite {

namespace Timestep {

namespace PredictorCorrector {

void initQuasiInverse(Test& test, SparseMatrix& qi)
{
   const auto& nN = test.nN;
   if (test.basisId == Test::BasisId::WORLAND)
   {
      if (test.equationId == Test::EquationId::DIFFUSION)
      {
         const auto& w_a = Polynomial::Worland::worland_chebyshev_t::ALPHA;
         const auto& w_db = Polynomial::Worland::worland_chebyshev_t::DBETA;

         const auto& l = test.l;
         SparseSM::Worland::I2 i2(nN, nN, w_a, w_db, l);
         SparseSM::Worland::Id qid(nN, nN, w_a, w_db, l, nN - 1);
         qi = i2.mat() - i2.mat() * qid.mat();
      }
      else if (test.equationId == Test::EquationId::BIDIFFUSION)
      {
         const auto& w_a = Polynomial::Worland::worland_chebyshev_t::ALPHA;
         const auto& w_db = Polynomial::Worland::worland_chebyshev_t::DBETA;

         const auto& l = test.l;
         SparseSM::Worland::I4 i4(nN, nN, w_a, w_db, l);
         SparseSM::Worland::Id qid(nN, nN, w_a, w_db, l, nN - 2);
         qi = i4.mat() - i4.mat() * qid.mat();
      }
   }
   else if (test.basisId == Test::BasisId::CHEBYSHEV)
   {
      if (test.equationId == Test::EquationId::DIFFUSION)
      {
         const MHDFloat lb = -1.0;
         const MHDFloat ub = 1.0;
         SparseSM::Chebyshev::LinearMap::I2 i2(nN, nN, lb, ub);
         qi = i2.mat();
      }
      else if (test.equationId == Test::EquationId::BIDIFFUSION)
      {
         const MHDFloat lb = -1.0;
         const MHDFloat ub = 1.0;
         SparseSM::Chebyshev::LinearMap::I4 i4(nN, nN, lb, ub);
         qi = i4.mat();
      }
   }
}

void setForcing(Test& test, Matrix& forcing, const MHDFloat t)
{
   auto ref = test.getForcing();
   const auto& a = test.timeCoeffs(0);
   const auto& b = test.timeCoeffs(1);

   for (int i = 0; i < ref.rows(); i++)
   {
      forcing.row(i).array() = ref(i, 0) * std::cos(a * t);
      forcing.row(i).array() += ref(i, 1) * std::sin(a * t);
      forcing.row(i).array() += ref(i, 2) * std::cos(b * t);
      forcing.row(i).array() += ref(i, 3) * std::sin(b * t);
   }
}

void checkSolution(Test& test, const Matrix& data, const int tIdx)
{
   auto ref = test.getReference();

   auto err = (data - ref.col(tIdx)).array().abs().maxCoeff();
   std::cerr << test.tEnd(tIdx) << "  " << err << std::endl;
}

void createOperators(Test& test, std::map<std::size_t, DecoupledZSparse>& ops)
{
   const auto& nN = test.nN;

   if (test.basisId == Test::BasisId::WORLAND)
   {
      if (test.equationId == Test::EquationId::DIFFUSION)
      {
         const auto& l = test.l;
         const auto& w_a = Polynomial::Worland::worland_chebyshev_t::ALPHA;
         const auto& w_db = Polynomial::Worland::worland_chebyshev_t::DBETA;
         // Compute model's linear operator (without Tau lines)
         ops.insert(std::make_pair(ModelOperator::ImplicitLinear::id(),
            DecoupledZSparse()));
         SparseSM::Worland::I2Lapl opL(nN, nN, w_a, w_db, l);
         ops.find(ModelOperator::ImplicitLinear::id())->second.real() =
            opL.mat();
         // Compute model's time operator (Mass matrix)
         ops.insert(
            std::make_pair(ModelOperator::Time::id(), DecoupledZSparse()));
         SparseSM::Worland::I2 opT(nN, nN, w_a, w_db, l);
         SparseSM::Worland::Id qid(nN, nN, w_a, w_db, l, nN - 1);
         ops.find(ModelOperator::Time::id())->second.real() =
            opT.mat() - opT.mat() * qid.mat();
         // Compute model's boundary operator
         ops.insert(
            std::make_pair(ModelOperator::Boundary::id(), DecoupledZSparse()));
         SparseMatrix bc(nN, nN);
         std::vector<Eigen::Triplet<MHDFloat>> coeffs;
         SparseSM::Worland::Boundary::Value bcValue(w_a, w_db, l);
         auto bcVal = bcValue.compute(nN - 1).cast<MHDFloat>();
         for (int i = 0; i < nN; i++)
         {
            coeffs.push_back(Eigen::Triplet<MHDFloat>(0, i, bcVal(i)));
         }
         bc.setFromTriplets(coeffs.begin(), coeffs.end());
         ops.find(ModelOperator::Boundary::id())->second.real() = bc;
      }
      else if (test.equationId == Test::EquationId::BIDIFFUSION)
      {
         const auto& l = test.l;
         const auto& w_a = Polynomial::Worland::worland_chebyshev_t::ALPHA;
         const auto& w_db = Polynomial::Worland::worland_chebyshev_t::DBETA;
         // Compute model's linear operator (without Tau lines)
         ops.insert(std::make_pair(ModelOperator::ImplicitLinear::id(),
            DecoupledZSparse()));
         SparseSM::Worland::I4Lapl2 opL(nN, nN, w_a, w_db, l);
         ops.find(ModelOperator::ImplicitLinear::id())->second.real() =
            opL.mat();
         // Compute model's time operator (Mass matrix)
         ops.insert(
            std::make_pair(ModelOperator::Time::id(), DecoupledZSparse()));
         SparseSM::Worland::I4Lapl opT(nN, nN, w_a, w_db, l);
         ops.find(ModelOperator::Time::id())->second.real() = opT.mat();
         // Compute model's boundary operator
         ops.insert(
            std::make_pair(ModelOperator::Boundary::id(), DecoupledZSparse()));
         SparseMatrix bc(nN, nN);
         std::vector<Eigen::Triplet<MHDFloat>> coeffs;
         SparseSM::Worland::Boundary::Value bcValue(w_a, w_db, l);
         auto bcVal = bcValue.compute(nN - 1).cast<MHDFloat>();
         SparseSM::Worland::Boundary::D2 bcDiff2(w_a, w_db, l);
         auto bcD2 = bcDiff2.compute(nN - 1).cast<MHDFloat>();
         for (int i = 0; i < nN; i++)
         {
            coeffs.push_back(Eigen::Triplet<MHDFloat>(0, i, bcVal(i)));
            coeffs.push_back(Eigen::Triplet<MHDFloat>(1, i, bcD2(i)));
         }
         bc.setFromTriplets(coeffs.begin(), coeffs.end());
         ops.find(ModelOperator::Boundary::id())->second.real() = bc;
      }
   }
   else if (test.basisId == Test::BasisId::CHEBYSHEV)
   {
      const MHDFloat lb = -1.0;
      const MHDFloat ub = 1.0;

      // Compute model's linear operator (without Tau lines)
      ops.insert(std::make_pair(ModelOperator::ImplicitLinear::id(),
         DecoupledZSparse()));
      SparseMatrix opL(nN, nN);
      std::vector<Eigen::Triplet<MHDFloat>> coeffs;
      for (int i = 2; i < nN; i++)
      {
         coeffs.push_back(Eigen::Triplet<MHDFloat>(i, i, 1.0));
      }
      opL.setFromTriplets(coeffs.begin(), coeffs.end());
      ops.find(ModelOperator::ImplicitLinear::id())->second.real() = opL;
      // Compute model's time operator (Mass matrix)
      ops.insert(std::make_pair(ModelOperator::Time::id(), DecoupledZSparse()));
      SparseSM::Chebyshev::LinearMap::I2 opT(nN, nN, lb, ub);
      ops.find(ModelOperator::Time::id())->second.real() = opT.mat();
      // Compute model's boundary operator
      ops.insert(
         std::make_pair(ModelOperator::Boundary::id(), DecoupledZSparse()));
      SparseMatrix bc(nN, nN);
      coeffs.clear();
      coeffs.push_back(Eigen::Triplet<MHDFloat>(0, 0, 1.0));
      coeffs.push_back(Eigen::Triplet<MHDFloat>(1, 0, 1.0));
      for (int i = 1; i < nN; i++)
      {
         coeffs.push_back(Eigen::Triplet<MHDFloat>(0, i, 2.0));
         coeffs.push_back(
            Eigen::Triplet<MHDFloat>(1, i, std::pow(-1.0, i) * 2.0));
      }
      bc.setFromTriplets(coeffs.begin(), coeffs.end());
      ops.find(ModelOperator::Boundary::id())->second.real() = bc;
   }
}

void initSolution(Test& test, DecoupledZMatrix& sol)
{
   auto ref = test.getInitial();
   sol.real() = ref;
   sol.imag() = ref;
}

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace TestSuite
} // namespace QuICC
