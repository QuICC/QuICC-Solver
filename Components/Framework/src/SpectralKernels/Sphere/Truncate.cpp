/**
 * @file Truncate.cpp
 * @brief Source of Truncate Worland expansion with minimal energy loss
 */

// System includes
//

// Class include
//
#include "QuICC/SpectralKernels/Sphere/Truncate.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Polynomial/Worland/Utils.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/NoPenetration.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD2.hpp"
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphere.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolverTools.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

namespace Sphere {

   Truncate::Truncate(const bool isComplex)
      : ISpectralKernel(isComplex)
   {
   }

   void Truncate::init(const FieldComponents::Spectral::Id comp, const std::string& wType, const std::size_t bcId, const int outN, const int outL, const int outM)
   {
      this->mComp = comp;
      this->mWType = wType;
      this->mBcId = bcId;
      this->mOutN = outN;
      this->mOutL = outL;
      this->mOutM = outM;
   }

   MHDVariant Truncate::compute(const int i, const int j, const int k) const
   {
      if(this->mIsComplex)
      {
         return MHDComplex(0,0);
      }
      else
      {
         return 0.0;
      }
   }

   void Truncate::apply(const std::size_t timeId)
   {
      if(timeId == SolveTiming::After::id())
      {
         const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         int nR = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL);

         Internal::MHDFloat inAlpha, inDBeta;
         Internal::Array inGrid(nR), inWeights(nR);
         Polynomial::Worland::Utils::selectJacobi(this->mWType, inAlpha, inDBeta, inGrid, inWeights);

         Internal::MHDFloat outAlpha, outDBeta;
         Internal::Array outGrid(nR), outWeights(nR);
         std::string outWType;
         if(this->mComp == FieldComponents::Spectral::TOR)
         {
            outWType = "SphEnergy";
         }
         else if(this->mComp == FieldComponents::Spectral::POL)
         {
            outWType = "SphEnergy";
         }
         else if(this->mComp == FieldComponents::Spectral::SCALAR)
         {
            outWType = "SphEnergy";
         }
         else
         {
            throw std::logic_error("Unknown component");
         }
         Polynomial::Worland::Utils::selectJacobi(outWType, outAlpha, outDBeta, outGrid, outWeights);

         const auto hasMOrdering= this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

         Polynomial::Worland::Wnl inWnl(inAlpha, inDBeta);
         Polynomial::Worland::Wnl outWnl(outAlpha, outDBeta);

         auto applyImpl = [&](auto&& comp)
         {
            if(hasMOrdering)
            {
               typedef typename std::remove_cv<decltype(comp.point(0,0,0))>::type  DataType;
               typedef Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>  MatType;

               // Loop over harmonic order m
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
               {
                  int m_ = tRes.idx<Dimensions::Data::DAT3D>(k);

                  // Loop over harmonic degree l
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);

                     if(m_ <= this->mOutM && l_ <= this->mOutL)
                     {
                        int inNPoly = tRes.dim<Dimensions::Data::DATB1D>(j,k);
                        int outNPoly = this->mOutN + 1; // Deal with triangular truncation

                        if(inNPoly > 0)
                        {
                           Internal::Matrix iproj(outGrid.size(), inNPoly);
                           Internal::Matrix iintg(outGrid.size(), inNPoly);

                           // Convert to energy norm basis
                           inWnl.compute<Internal::MHDFloat>(iproj, inNPoly, l_, outGrid, Internal::Array(), Polynomial::Worland::Evaluator::Set());
                           Matrix proj = iproj.cast<MHDFloat>();
                           outWnl.compute<Internal::MHDFloat>(iintg, inNPoly, l_, outGrid, outWeights, Polynomial::Worland::Evaluator::Set());
                           Matrix intg = iintg.cast<MHDFloat>();

                           MatType tv = comp.profile(j, k);
                           MatType iv(inNPoly, tv.cols());
                           for(int i = 0; i < inNPoly; i++)
                           {
                              iv(i,0) = tv(i,0);
                           }
                           iv = intg.transpose() * (proj * iv);

                           // Truncate Galerkin expansion
                           SparseMatrix matS = stencil(inNPoly, outAlpha, outDBeta, l_, true);
                           if(this->mBcId == 0 || matS.rows() == 0)
                           {
                              tv = iv.topRows(outNPoly);
                              iv = tv;
                           }
                           else
                           {
                              Framework::Selector::SparseSolver<SparseMatrix> solver;
                              solver.compute(matS);
                              tv = iv.topRows(matS.rows());
                              iv = tv;
                              iv.setZero();
                              Solver::details::solveWrapper(iv, solver, tv);

                              int nPoly = outNPoly - (inNPoly - matS.rows());
                              matS = stencil(outNPoly, outAlpha, outDBeta, l_, false);
                              tv = iv.topRows(nPoly);
                              iv = matS * tv;
                           }

                           // Project back to original basis
                           iintg.resize(inGrid.size(), outNPoly);
                           iproj.resize(inGrid.size(), outNPoly);

                           inWnl.compute<Internal::MHDFloat>(iintg, outNPoly, l_, inGrid, inWeights, Polynomial::Worland::Evaluator::Set());
                           intg = iintg.cast<MHDFloat>();
                           outWnl.compute<Internal::MHDFloat>(iproj, outNPoly, l_, inGrid, Internal::Array(), Polynomial::Worland::Evaluator::Set());
                           proj = iproj.cast<MHDFloat>();

                           iv = intg.transpose() * (proj * iv);

                           tv = comp.profile(j, k);
                           tv.setZero();
                           tv.topRows(outNPoly) = iv;
                           comp.setProfile(tv, j, k);
                        }
                     }
                     else
                     {
                        MatType tv = comp.profile(j, k);
                        tv.setZero();
                        comp.setProfile(tv, j, k);
                     }
                  }
               }
            }
            else
            {
               throw std::logic_error("Not yet implemented for L ordering");

               typedef typename std::remove_cv<decltype(comp.point(0,0,0))>::type  DataType;
               typedef Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>  MatType;

               // Loop over harmonic degree l
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
               {
                  int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);

                  // Loop over harmonic order m
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     int m_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);

                     if(m_ <= this->mOutM && l_ <= this->mOutL)
                     {
                        int inNPoly = tRes.dim<Dimensions::Data::DATB1D>(j,k);
                        int outNPoly = this->mOutN + 1; // Deal with triangular truncation

                        if(inNPoly > 0)
                        {
                           Internal::Matrix iproj(outGrid.size(), inNPoly);
                           Internal::Matrix iintg(outGrid.size(), outNPoly);

                           inWnl.compute<Internal::MHDFloat>(iproj, inNPoly, l_, outGrid, Internal::Array(), Polynomial::Worland::Evaluator::Set());
                           Matrix proj = iproj.cast<MHDFloat>();
                           outWnl.compute<Internal::MHDFloat>(iintg, outNPoly, l_, outGrid, outWeights, Polynomial::Worland::Evaluator::Set());
                           Matrix intg = iintg.cast<MHDFloat>();

                           MatType tv = comp.profile(j, k);
                           MatType iv(inNPoly, tv.cols());
                           for(int i = 0; i < inNPoly; i++)
                           {
                              iv(i,0) = tv(i,0);
                           }
                           tv.setZero();
                           iv = intg.transpose() * (proj * iv);
                           tv.topRows(outNPoly) = iv;
                           comp.setProfile(tv, j, k);
                        }
                     }
                     else
                     {
                        MatType tv = comp.profile(j, k);
                        tv.setZero();
                        comp.setProfile(tv, j, k);
                     }
                  }
               }
            }
         };

         if(this->mScalars.size() == 1)
         {
            auto& s = this->mScalars.begin()->second;
            std::visit(
                  [&](auto&& field)
                  {
                     auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                     applyImpl(comp);
                  }, s);
         }
         else if(this->mVectors.size() == 1)
         {
            auto& v = this->mVectors.begin()->second;
            std::visit(
                  [&](auto&& field)
                  {
                     auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                     applyImpl(comp);
                  }, v);
         }
         else
         {
            throw std::logic_error("Missing field");
         }
      }
   }

   SparseMatrix Truncate::stencil(const int tN, const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta, const int l, const bool isSquare)
   {
      const std::size_t& bcId = this->mBcId;
      const FieldComponents::Spectral::Id& comp = this->mComp;

      SparseMatrix matS;
      int bc = 0;
      if(
         bcId == Bc::Name::FixedTemperature::id() ||
         (comp == FieldComponents::Spectral::TOR && bcId == Bc::Name::Insulating::id()) ||
         (comp == FieldComponents::Spectral::POL && bcId == Bc::Name::NoPenetration::id())
         )
      {
         bc = 1;
         SparseSM::Worland::Stencil::Value S(tN, tN - bc, alpha, dBeta, l);
         matS = S.mat();
      }
      else if(bcId == Bc::Name::FixedFlux::id())
      {
         bc = 1;
         SparseSM::Worland::Stencil::D1 S(tN, tN - bc, alpha, dBeta, l);
         matS = S.mat();
      }
      else if(
         (comp == FieldComponents::Spectral::POL && bcId == Bc::Name::Insulating::id())
         )
      {
         bc = 1;
         SparseSM::Worland::Stencil::InsulatingSphere S(tN, tN - bc, alpha, dBeta, l);
         matS = S.mat();
      }
      else if(bcId == Bc::Name::NoSlip::id())
      {
         bc = 2;
         SparseSM::Worland::Stencil::ValueD1 S(tN, tN - bc, alpha, dBeta, l);
         matS = S.mat();
      }
      else if(bcId == Bc::Name::StressFree::id())
      {
         bc = 2;
         SparseSM::Worland::Stencil::ValueD2 S(tN, tN - bc, alpha, dBeta, l);
         matS = S.mat();
      }
      else if(
         (comp == FieldComponents::Spectral::TOR && bcId == Bc::Name::NoPenetration::id())
         )
      {
         matS.resize(0,0);
      }
      else
      {
         throw std::logic_error("Unknown boundary condition");
      }
      matS.makeCompressed();

      if (isSquare && matS.rows() > 0)
      {
         SparseSM::Worland::Id qId(tN - bc, tN, alpha, dBeta, l);
         matS = qId.mat() * matS;
      }

      return matS;
   }

} // Sphere
} // Kernel
} // Spectral
} // QuICC
