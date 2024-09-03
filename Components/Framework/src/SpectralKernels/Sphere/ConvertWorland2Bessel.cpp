/**
 * @file ConvertWorland2Bessel.cpp
 * @brief Source of converter from Worland basis to another Worland type basis
 */

// System includes
//

// Class include
//
#include "QuICC/SpectralKernels/Sphere/ConvertWorland2Bessel.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Bessel/Value.hpp"
#include "QuICC/Polynomial/Bessel/Insulating.hpp"
#include "QuICC/Polynomial/Bessel/NoSlip.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/SparseSM/Id.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/Value.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/D1.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/D2.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/R1D1DivR1.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

namespace Sphere {

   ConvertWorland2Bessel::ConvertWorland2Bessel(const bool isComplex)
      : ISpectralKernel(isComplex)
   {
   }

   void ConvertWorland2Bessel::init(const FieldComponents::Spectral::Id comp, const WorlandKind inWType, const BesselKind outBType, const MHDFloat scale)
   {
      this->mComp = comp;
      this->mInWType = inWType;
      this->mOutBType = outBType;
      this->mScale = scale;
   }

   MHDVariant ConvertWorland2Bessel::compute(const int i, const int j, const int k) const
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

   void ConvertWorland2Bessel::apply(const std::size_t timeId)
   {
      if(timeId == SolveTiming::After::id())
      {
         auto selectJacobi = [](const WorlandKind t, Internal::MHDFloat& a, Internal::MHDFloat& db, Internal::Array& igrid, Internal::Array& iweights)
         {
            auto defineWorland = [](const auto& wt, Internal::MHDFloat& a, Internal::MHDFloat& db, Internal::Array& igrid, Internal::Array& iweights)
            {
               // Set alpha and dBeta
               a = wt.ALPHA;
               db = wt.DBETA;

               // Compute quadrature
               int nR = igrid.size();
               if(nR > 0)
               {
                  typename std::remove_reference<decltype(wt)>::type::Rule wquad;
                  wquad.computeQuadrature(igrid, iweights, nR);
               }
            };

            if(t == WorlandKind::Chebyshev)
            {
               Polynomial::Worland::worland_chebyshev_t wt;
               defineWorland(wt, a, db, igrid, iweights);
            }
            else if(t == WorlandKind::Legendre)
            {
               Polynomial::Worland::worland_legendre_t wt;
               defineWorland(wt, a, db, igrid, iweights);
            }
            else if(t == WorlandKind::SphEnergy)
            {
               Polynomial::Worland::worland_sphenergy_t wt;
               defineWorland(wt, a, db, igrid, iweights);
            }
            else if(t == WorlandKind::CylEnergy)
            {
               Polynomial::Worland::worland_cylenergy_t wt;
               defineWorland(wt, a, db, igrid, iweights);
            }
            else
            {
               throw std::logic_error("Unknown Worland type");
            }
         };

         Internal::MHDFloat inAlpha, inDBeta;
         Internal::Array igrid(0), iweights(0);
         selectJacobi(this->mInWType, inAlpha, inDBeta, igrid, iweights);

         // Use SphEnergy grid for conversion
         const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         Internal::MHDFloat tmp;
         int nR = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL);
         igrid.resize(nR);
         iweights.resize(nR);
         selectJacobi(WorlandKind::SphEnergy, tmp, tmp, igrid, iweights);

         const auto hasMOrdering= this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

         Polynomial::Worland::Wnl inWnl(inAlpha, inDBeta);

         auto applyImpl = [&](auto& comp, auto& outJnl)
         {
            typedef typename std::remove_cv<decltype(comp.point(0,0,0))>::type  DataType;
            typedef Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>  MatType;
            if(hasMOrdering)
            {
               // Loop over harmonic order m
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
               {
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     int l_ = tRes.idx<Dimensions::Data::DAT2D>(j, k);
                     int nPoly = tRes.dim<Dimensions::Data::DATB1D>(j,k);

                     if(nPoly > 0)
                     {
                        Internal::Matrix tmp(igrid.size(), nPoly);
                        Matrix proj(igrid.size(), nPoly);
                        Matrix intg(igrid.size(), nPoly);

                        inWnl.compute<Internal::MHDFloat>(tmp, nPoly, l_, igrid, Internal::Array(), Polynomial::Worland::Evaluator::Set());
                        proj = tmp.cast<MHDFloat>();
                        tmp.setZero();
                        outJnl.template compute<Internal::MHDFloat>(tmp, nPoly, l_, igrid, iweights);
                        intg = tmp.cast<MHDFloat>();

                        MatType tv = comp.profile(j, k);
                        MatType iv(nPoly, tv.cols());
                        for(int i = 0; i < nPoly; i++)
                        {
                           iv(i,0) = tv(i,0);
                        }
                        tv.setZero();
                        iv = intg.transpose() * (proj * iv);
                        tv.topRows(nPoly) = iv;
                        tv *= this->mScale;
                        comp.setProfile(tv, j, k);
                     }
                  }
               }
            }
            else
            {
               // Loop over harmonic order l
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); ++k)
               {
                  int l_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     int nPoly = tRes.dim<Dimensions::Data::DATB1D>(j,k);
                     SparseSM::Bessel::Boundary::Operator bcOp(nPoly, nPoly, SparseSM::Bessel::BesselKind::NOSLIP, l_, true);
                     bcOp.addRow<SparseSM::Bessel::Boundary::Value>();

                     if(nPoly > 0)
                     {
                        Internal::Matrix tmp(igrid.size(), nPoly);
                        Matrix proj(igrid.size(), nPoly);
                        Matrix intg(igrid.size(), nPoly);
                        Matrix projJ(igrid.size(), nPoly);

                        inWnl.compute<Internal::MHDFloat>(tmp, nPoly, l_, igrid, Internal::Array(), Polynomial::Worland::Evaluator::Set());
                        proj = tmp.cast<MHDFloat>();
                        tmp.setZero();
                        outJnl.template compute<Internal::MHDFloat>(tmp, nPoly, l_, igrid, iweights);
                        intg = tmp.cast<MHDFloat>();
                        tmp.setZero();
                        outJnl.template compute<Internal::MHDFloat>(tmp, nPoly, l_, igrid, Internal::Array());
                        projJ = tmp.cast<MHDFloat>();

                        MatType tv = comp.profile(j, k);
                        MatType iv(nPoly, tv.cols());
                        for(int i = 0; i < nPoly; i++)
                        {
                           iv(i,0) = tv(i,0);
                        }
                        tv.setZero();
                        iv = intg.transpose() * (proj * iv);
                        tv.topRows(nPoly) = iv;
                        if(this->mComp == FieldComponents::Spectral::POL)
                        {
                           tv.topRows(1) += -(bcOp.mat() * tv).topRows(1)/bcOp.mat().coeffRef(0,0);
                        }
                        tv *= this->mScale;
                        comp.setProfile(tv, j, k);
                     }
                  }
               }
            }
         };

         if(this->mScalars.size() == 1)
         {
            auto& s = this->mScalars.begin()->second;
            if(this->mOutBType == BesselKind::Value)
            {
               Polynomial::Bessel::Value<Polynomial::Bessel::SphJnl> outJnl;
               std::visit(
                     [&](auto&& field)
                     {
                        auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                        applyImpl(comp, outJnl);
                     }, s);
            }
            else if(this->mOutBType == BesselKind::Insulating)
            {
               Polynomial::Bessel::Insulating<Polynomial::Bessel::SphJnl> outJnl;
               std::visit(
                     [&](auto&& field)
                     {
                        auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                        applyImpl(comp, outJnl);
                     }, s);
            }
            else if(this->mOutBType == BesselKind::NoSlip)
            {
               Polynomial::Bessel::NoSlip<Polynomial::Bessel::SphJnl> outJnl;
               std::visit(
                     [&](auto&& field)
                     {
                        auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                        applyImpl(comp, outJnl);
                     }, s);
            }
         }
         else if(this->mVectors.size() == 1)
         {
            auto& v = this->mVectors.begin()->second;
            if(this->mOutBType == BesselKind::Value)
            {
               Polynomial::Bessel::Value<Polynomial::Bessel::SphJnl> outJnl;
               std::visit(
                     [&](auto&& field)
                     {
                        auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                        applyImpl(comp, outJnl);
                     }, v);
            }
            else if(this->mOutBType == BesselKind::Insulating)
            {
               Polynomial::Bessel::Insulating<Polynomial::Bessel::SphJnl> outJnl;
               std::visit(
                     [&](auto&& field)
                     {
                        auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                        applyImpl(comp, outJnl);
                     }, v);
            }
            else if(this->mOutBType == BesselKind::NoSlip)
            {
               Polynomial::Bessel::NoSlip<Polynomial::Bessel::SphJnl> outJnl;
               std::visit(
                     [&](auto&& field)
                     {
                        auto&& comp = field->rDom(0).rPerturbation().rComp(this->mComp);
                        applyImpl(comp, outJnl);
                     }, v);
            }
         }
         else
         {
            throw std::logic_error("Missing field");
         }
      }
   }

} // Sphere
} // Kernel
} // Spectral
} // QuICC
