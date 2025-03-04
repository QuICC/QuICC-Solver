/**
 * @file MarginalCurve.cpp
 * @brief Source of the high level simulation
 */

// System includes
//
#include <algorithm>
#include <boost/math/tools/roots.hpp>
#include <iostream>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/NonDimensional/Nev.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/StabilityMode.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/Feature.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Variables/RequirementTools.hpp"
#include "Stability/LinearStability.hpp"
#include "Stability/MarginalCurve.hpp"
#include "Types/Math.hpp"

namespace QuICC {

/**
 * @brief Functor to wrap calls for stability calculations
 */
template <unsigned int TLevel> class func
{
public:
   /**
    * @brief Constructor
    *
    * @param name   Name to use in log
    * @param spStab Linear stability solver
    * @param spLog  Logging stream
    */
   func(const std::string& name, std::shared_ptr<LinearStability> spStab,
      std::shared_ptr<std::ofstream> spLog) :
       mcName(name), mspStab(spStab), mspLog(spLog) {};

   /**
    * @brief Destructor
    */
   ~func() = default;

   /**
    * @brief Interface to call solver
    */
   MHDFloat operator()(const MHDFloat vc)
   {
      if constexpr (TLevel == 0)
      {
         return (*this->mspStab)(vc);
      }
      else if constexpr (TLevel == 1)
      {
         const unsigned int nev = 5;
         std::vector<MHDComplex> evs(nev);
         auto ret = (*this->mspStab)(vc, evs);

         this->print(vc, evs);

         return ret;
      }
      else
      {
         throw std::logic_error("Unknown output level");
      }
   };

   /**
    * @brief Logging computed eigenvalues
    */
   void print(const MHDFloat vc, const std::vector<MHDComplex>& evs) const
   {
      std::ofstream& log = *this->mspLog;
      log << this->mcName << " = " << vc << std::endl;
      for (unsigned int i = 0; i < evs.size(); i++)
      {
         log << " \t growth: " << evs.at(i) << std::endl;
      }
   }

protected:
   /**
    * @brief Name for logging
    */
   const std::string mcName;

   /**
    * @brief Linear stability calculator
    */
   std::shared_ptr<LinearStability> mspStab;

   /**
    * @brief Logging stream
    */
   std::shared_ptr<std::ofstream> mspLog;
};

MarginalCurve::MarginalCurve() : StabilityBase()
{
   // Initialize SLEPc/PETSc
   PetscCallVoid(SlepcInitializeNoArguments());
}

MarginalCurve::~MarginalCurve()
{
   // Finalize SLEPc/PETSc
   PetscCallVoid(SlepcFinalize());
}

void MarginalCurve::saveEigenfunction(const int m, const MHDComplex ev,
   const std::vector<MHDComplex>& ef)
{
   auto&& ss = this->mspRes->sim().ss();
   if (!this->mpH5File)
   {
      this->mpH5File = std::make_shared<Io::Variable::StateFileWriter>(ss.tag(),
         ss.has(SpatialScheme::Feature::RegularSpectrum));

      VariableRequirement varInfo;
      for (auto fId: this->mspBackend->fieldIds())
      {
         if (fId == PhysicalNames::Velocity::id() ||
             fId == PhysicalNames::Magnetic::id())
         {
            auto& req = varInfo.addField(fId,
               FieldRequirement(false, ss.spectral(), ss.physical()));
            req.enableSpectral();
         }
         else if (fId == PhysicalNames::Temperature::id())
         {
            auto& req = varInfo.addField(fId,
               FieldRequirement(true, ss.spectral(), ss.physical()));
            req.enableSpectral();
         }
         else
         {
            throw std::logic_error("Unknown field for saving eigenfuction");
         }
         this->mpH5File->expect(fId);
      }

      // Initialize variables
      RequirementTools::initVariables(this->mScalars, this->mVectors, varInfo,
         this->mspRes);

      // Add scalars
      for (auto&& s: this->mScalars)
      {
         this->mpH5File->addScalar(s);
      }

      // Add vectors
      for (auto&& v: this->mVectors)
      {
         this->mpH5File->addVector(v);
      }
   }

   // Fields
   auto pId = this->mspBackend->fieldIds().at(0);
   FieldComponents::Spectral::Id comp;
   if (this->mVectors.count(pId) > 0)
   {
      comp = FieldComponents::Spectral::TOR;
   }
   else
   {
      comp = FieldComponents::Spectral::SCALAR;
   }
   auto fId = std::make_pair(pId, comp);

   const auto& res = *this->mspRes;
   Model::EquationInfo eqInfo;
   this->backend().equationInfo(eqInfo, fId, res);

   auto imRange = std::make_pair(eqInfo.im.begin(), eqInfo.im.end());

   std::size_t idx = 0;
   const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
   for (auto it = imRange.first; it != imRange.second; ++it)
   {
      if (this->mScalars.count(it->first) > 0)
      {
         std::visit(
            [&](auto& p)
            {
               // Clear previous values
               p->rDom(0).rPerturbation().setZeros();

               for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  int k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  if (k_ == m)
                  {
                     for (int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k);
                          j++)
                     {
                        for (int i = 0;
                             i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
                        {
                           p->rDom(0).rPerturbation().setPoint(ef.at(idx), i, j,
                              k);
                           idx++;
                        }
                     }
                  }
               }
            },
            this->mScalars.at(it->first));
      }

      if (this->mVectors.count(it->first) > 0)
      {
         std::visit(
            [&](auto& p)
            {
               // Clear previous values
               p->rDom(0).rPerturbation().rComp(it->second).setZeros();

               for (int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  int k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  if (k_ == m)
                  {
                     for (int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k);
                          j++)
                     {
                        for (int i = 0;
                             i < tRes.dim<Dimensions::Data::DATF1D>(j, k); i++)
                        {
                           p->rDom(0)
                              .rPerturbation()
                              .rComp(it->second)
                              .setPoint(ef.at(idx), i, j, k);
                           idx++;
                        }
                     }
                  }
               }
            },
            this->mVectors.at(it->first));
      }
   }

   if (idx != ef.size())
   {
      throw std::logic_error("Wrong sizes");
   }

   this->mpH5File->init();
   this->mpH5File->write();
   this->mpH5File->finalize();
}

void MarginalCurve::mainRun()
{
   std::vector<MHDFloat> eigs;
   std::size_t max_nev;

   if (this->mspRes->sim().ss().has(SpatialScheme::Feature::SpectralMatrix1D))
   {
      eigs.clear();
      eigs = {this->mspRes->sim().boxScale(Dimensions::Simulation::SIM2D),
         this->mspRes->sim().boxScale(Dimensions::Simulation::SIM3D)};

      int nF = this->mspBackend->fieldIds().size();
      int nN = this->mspRes->sim().dim(Dimensions::Simulation::SIM1D,
         Dimensions::Space::SPECTRAL);
      max_nev = nF * nN;
   }
   else if (this->mspRes->sim().ss().has(
               SpatialScheme::Feature::SpectralMatrix2D))
   {
      int m = this->mspRes->sim().dim(Dimensions::Simulation::SIM3D,
                 Dimensions::Space::SPECTRAL) -
              1;
      MHDFloat m_ = static_cast<MHDFloat>(m);
      eigs.clear();
      eigs = {m_};

      int nF = this->mspBackend->fieldIds().size();
      int nN = this->mspRes->sim().dim(Dimensions::Simulation::SIM1D,
         Dimensions::Space::SPECTRAL);
      int nL = (this->mspRes->sim().dim(Dimensions::Simulation::SIM2D,
                   Dimensions::Space::SPECTRAL) -
                m);
      max_nev = nF * nN * nL;
   }
   else
   {
      throw std::logic_error("3D spectral matrix not setup");
   }

   auto idc = NonDimensional::Rayleigh::id();
   auto name = NonDimensional::Coordinator::tag(idc);

   auto spLinStab = std::make_shared<LinearStability>(idc, eigs, this->mspRes,
      this->mspEqParams->map(), this->createBoundary()->map(),
      this->mspBackend);

   auto nev_ = this->mspEqParams->nd(NonDimensional::Nev::id());
   unsigned int nev = 0;
   if (nev_ <= 0)
   {
      nev = max_nev;
   }
   else
   {
      nev = static_cast<unsigned int>(nev_);
   }

   int solver_mode = this->mspEqParams->nd(NonDimensional::StabilityMode::id());

   // Single mode calculation
   if (solver_mode == 0)
   {
      std::vector<MHDComplex> evs(nev);
      std::vector<std::vector<MHDComplex>> efs(nev);
      auto vc = this->mspEqParams->nd(idc);
      spLinStab->eigenpairs(evs, efs, nev, vc);

      // Print eigenvalues
      for (auto&& e: evs)
      {
         std::cerr << e << std::endl;
      }

      for (std::size_t i = 0; i < efs.size(); i++)
      {
         if (this->mspRes->sim().ss().has(
                SpatialScheme::Feature::SpectralMatrix2D))
         {
            int k_ = static_cast<int>(eigs.at(0));
            saveEigenfunction(k_, evs.at(i), efs.at(i));
         }
         else
         {
            // Writing other eigenfunctions is not yet implemented
         }
      }
   }
   // Find critical parameter
   else
   {
      MHDFloat guess = this->mspEqParams->nd(idc);
      MHDFloat factor = 1.1; // To multiply
      int digits =
         std::numeric_limits<MHDFloat>::digits; // Maximum possible binary
                                                // digits accuracy for type T.
                                                // digits used to control how
                                                // accurate to try to make the
                                                // result.
      int get_digits =
         (digits * 3) / 4; // Near maximum (3/4) possible accuracy.

      const boost::uintmax_t maxit = 20;
      boost::uintmax_t it =
         maxit; // Initally our chosen max iterations, but updated with actual.
                // We could also have used a maximum iterations provided by any
                // policy: boost::uintmax_t max_it =
                // policies::get_max_root_iterations<Policy>();
      bool is_rising =
         true; // So if result if guess^3 is too low, try increasing guess.
      boost::math::tools::eps_tolerance<double> tol(get_digits);

      auto spLog = std::make_shared<std::ofstream>("marginal.log");

      std::pair<MHDFloat, MHDFloat> r =
         boost::math::tools::bracket_and_solve_root(
            func<1>(name, spLinStab, spLog), guess, factor, is_rising, tol, it);

      std::ofstream& logger = *spLog;

      unsigned int prec = (std::numeric_limits<MHDFloat>::digits10 * 3) / 4 + 1;
      logger << std::string(50, '-') << std::endl
             << "Critical " + name + " number converged to the bracket: "
             << std::endl
             << std::setprecision(prec) << r.first << " < " << name << " < "
             << r.second << std::endl;
   }
}

void MarginalCurve::preRun() {}

void MarginalCurve::postRun() {}

} // namespace QuICC
