/** 
 * @file WLFlBuilder.cpp
 * @brief Source of the sphere Worland + Spherical Harmonics (Associated Legendre + Fourrier) scheme implementation with spectral l ordering
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFlBuilder.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/Transform/Fft/Fourier/Mixed/Setup.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"
#include "QuICC/Transform/Poly/ALegendre/Setup.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Transform/Setup/Triangular.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SpatialScheme/Tools/TriangularSH.hpp"

namespace QuICC {

namespace SpatialScheme {

   void WLFlBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      auto  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      auto  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      auto  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::SharedTransformSetup WLFlBuilder::spSetup1D(SharedResolution spRes) const
   {
      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA1D);

      // Get physical size of polynomial transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      
      Transform::SharedTransformSetup spSetup;

      const auto& opt = this->mOptions.at(0);

      // Poly algorithm setup
      if(std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end())
      {
         auto spPolySetup = std::make_shared<Transform::Poly::Setup>(size, specSize, this->purpose());
         spSetup = spPolySetup;
      }
      // FFT algorithm setup
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end())
      {
         auto spFftSetup = std::make_shared<Transform::Fft::Worland::Setup>(size, specSize, this->purpose());
         spSetup = spFftSetup;
      }

      // Get number of transforms and list of indexes
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         auto l = tRes.idx<Dimensions::Data::DAT3D>(i);
         auto nN = spRes->counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
         spSetup->addIndex(l, tRes.dim<Dimensions::Data::DAT2D>(i), nN);
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::SharedTransformSetup WLFlBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Check options consistency
      const auto& opt = this->mOptions.at(1);
      assert(opt.size() == 1 && std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end());

      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA2D);

      // Get physical size of polynomial transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Poly::ALegendre::Setup>(size, specSize, this->purpose());

      // Get number of transforms and list of indexes
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         spSetup->addIndex(tRes.idx<Dimensions::Data::DAT3D>(i), tRes.dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::SharedTransformSetup WLFlBuilder::spSetup3D(SharedResolution spRes) const
   {
      // Check options consistency
      const auto& opt = this->mOptions.at(2);
      assert(opt.size() == 1 && std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end());

      const auto& tRes = *spRes->cpu()->dim(Dimensions::Transform::TRA3D);

      // Get size of FFT transform
      int size = tRes.dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = 0;
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); i++)
      {
         blockSize += tRes.dim<Dimensions::Data::DAT2D>(i);
      }

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Mixed::Setup>(size, blockSize, specSize, this->purpose());
      spSetup->setBoxScale(1.0);
      spSetup->lock();

      return spSetup;
   }

   WLFlBuilder::WLFlBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options)
      : IRegularSHlBuilder(dim, purpose, options)
   {
   }

   void WLFlBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(3);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mL + 1;
      traSize(2) = this->mM + 1;
      this->setTransformSpace(traSize);

      const auto& opt = this->mOptions.at(0);

      int nN_ = 0;
      int nR_ = 0;
      // Poly algorithm setup with triangular truncation
      if(
            std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end() && 
            std::find(opt.begin(), opt.end(), Transform::Setup::Triangular::id()) != opt.end())
      {
         Tools::TriangularSH t;

         // Check if resolution is compatible with triangular truncation
         if(!t.isOptimal(this->mI+1,this->mL))
         {
            throw std::logic_error("Triangular truncation requires L+1 = 2(N-1)");
         }

         // radial spectral resolution
         nN_ = this->mI+1;

         // radial grid resolution
         nR_ = (2*(t.truncationBwd(nN_+7,0)-1) + 1)/2;
         // ... check against max L requirements
         nR_ = std::max(nR_, (2*(t.truncationBwd(nN_+7,this->mL)-1) + this->mL + 1)/2);
      }
      // Poly algorithm setup
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end())
      {
         // radial spectral resolution
         nN_ = this->mI+8;

         // radial grid resolution
         nR_ = this->mI+this->mL/2 + 8;
      }
      // FFT algorithm setup
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end())
      {
         // radial spectral resolution
         nN_ = this->mI+1;

         // radial grid resolution
         nR_ = this->mI+this->mL/2 + 4;
      }

      //
      // Compute sizes
      //

      // Get dealiased Worland transform size
      int nR = Transform::Poly::Tools::dealias(nR_);

      // Get dealiased associated Legendre transform size
      int nTh = Transform::Poly::Tools::dealias(this->mL+1);

      // Get standard dealiased FFT size
      int nPh = Transform::Fft::Tools::dealiasMixedFft(this->mM+1);
      // Check for optimised FFT sizes
      nPh = Transform::Fft::Tools::optimizeFft(nPh);

      // Modify grid size for visualiation
      if(this->purpose() == GridPurpose::VISUALIZATION)
      {
         // Make space for theta = 0 and  theta = pi
         nTh += 2;
         // Make space for r = 0, r = 1
         nR += 2;
      }

      //
      // Initialise spectral space
      //

      // Initialise forward dimension of first transform
      this->setDimension(traSize(0), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(traSize(0), Dimensions::Transform::SPECTRAL, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(2), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::SPECTRAL, Dimensions::Data::DAT3D);

      //
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(nN_, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

      // Initialise second dimension of first transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Initialise third dimension of first transform
      this->setDimension(traSize(1), Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);

      //
      // Initialise second transform
      //

      // Initialise forward dimension of second transform
      this->setDimension(nTh, Dimensions::Transform::TRA2D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of second transform
      this->setDimension(this->mL + 1, Dimensions::Transform::TRA2D, Dimensions::Data::DATB1D);

      // Initialise second dimension of second transform
      this->setDimension(nR, Dimensions::Transform::TRA2D, Dimensions::Data::DAT2D);

      // Initialise third dimension of second transform
      this->setDimension(traSize(2), Dimensions::Transform::TRA2D, Dimensions::Data::DAT3D);

      //
      // Initialise third transform
      //

      // Initialise forward dimension of third transform
      this->setDimension(nPh, Dimensions::Transform::TRA3D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of third transform
      this->setDimension(nPh/2 + 1, Dimensions::Transform::TRA3D, Dimensions::Data::DATB1D);

      // Initialise second dimension of third transform
      this->setDimension(nTh, Dimensions::Transform::TRA3D, Dimensions::Data::DAT2D);

      // Initialise third dimension of third transform
      this->setDimension(nR, Dimensions::Transform::TRA3D, Dimensions::Data::DAT3D);
   }

   void WLFlBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void WLFlBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void WLFlBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

} // SpatialScheme
} // QuICC
