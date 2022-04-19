/** 
 * @file WLFlfftBuilder.cpp
 * @brief Source of the sphere Worland(poly) + Spherical Harmonics (Associated Legendre(poly) + Fourrier) scheme implementation with spectral l ordering
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/3D/WLFlfftBuilder.hpp"

// Project includes
//
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   void WLFlfftBuilder::addTransformSetups(SharedResolution spRes) const
   {
      // Add setup for first transform
      Transform::Fft::Worland::SharedSetup  spS1D = this->spSetup1D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA1D, spS1D);

      // Add setup for second transform
      Transform::Poly::ALegendre::SharedSetup  spS2D = this->spSetup2D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA2D, spS2D);

      // Add setup for third transform
      Transform::Fft::Fourier::Mixed::SharedSetup  spS3D = this->spSetup3D(spRes);
      spRes->addTransformSetup(Dimensions::Transform::TRA3D, spS3D);
   }

   Transform::Fft::Worland::SharedSetup WLFlfftBuilder::spSetup1D(SharedResolution spRes) const
   {
      // Get physical size of polynomial transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Fft::Worland::Setup>(size, specSize, this->purpose());

      // Get number of transforms and list of indexes
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         spSetup->addIndex(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i), spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::Poly::ALegendre::SharedSetup WLFlfftBuilder::spSetup2D(SharedResolution spRes) const
   {
      // Get physical size of polynomial transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the polynomial transform
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      auto spSetup = std::make_shared<Transform::Poly::ALegendre::Setup>(size, specSize, this->purpose());

      // Get number of transforms and list of indexes
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         spSetup->addIndex(spRes->cpu()->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DAT3D>(i), spRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(i));
      }

      spSetup->lock();

      return spSetup;
   }

   Transform::Fft::Fourier::Mixed::SharedSetup WLFlfftBuilder::spSetup3D(SharedResolution spRes) const
   {
      // Get size of FFT transform
      int size = spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>();

      // Get spectral size of the FFT
      int specSize = spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get number of transforms
      int blockSize = 0;
      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         blockSize += spRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(i);
      }

      auto spSetup = std::make_shared<Transform::Fft::Fourier::Mixed::Setup>(size, blockSize, specSize, this->purpose());
      spSetup->setBoxScale(1.0);
      spSetup->lock();

      return spSetup;
   }

   WLFlfftBuilder::WLFlfftBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IRegularSHlBuilder(dim, purpose)
   {
   }

   WLFlfftBuilder::~WLFlfftBuilder()
   {
   }

   void WLFlfftBuilder::setDimensions()
   {
      //
      // Set transform space sizes
      //
      ArrayI traSize(3);
      traSize(0) = this->mI + 1;
      traSize(1) = this->mL + 1;
      traSize(2) = this->mM + 1;
      this->setTransformSpace(traSize);

      //
      // Compute sizes
      //

      // Get dealiased Worland transform size
      int nR = Transform::Poly::Tools::dealias(this->mI+this->mL/2+4);

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
      // Initialise first transform
      //

      // Initialise forward dimension of first transform
      this->setDimension(nR, Dimensions::Transform::TRA1D, Dimensions::Data::DATF1D);

      // Initialise backward dimension of first transform
      this->setDimension(this->mI+1, Dimensions::Transform::TRA1D, Dimensions::Data::DATB1D);

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

   void WLFlfftBuilder::setCosts()
   {
      // Set first transform cost
      this->setCost(1.0, Dimensions::Transform::TRA1D);

      // Set second transform cost
      this->setCost(1.0, Dimensions::Transform::TRA2D);

      // Set third transform cost
      this->setCost(1.0, Dimensions::Transform::TRA3D);
   }

   void WLFlfftBuilder::setScalings()
   {
      // Set first transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA1D);

      // Set second transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA2D);

      // Set third transform scaling
      this->setScaling(1.0, Dimensions::Transform::TRA3D);
   }

   void WLFlfftBuilder::setMemoryScore()
   {
      // Set first transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA1D);

      // Set second transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA2D);

      // Set third transform memory footprint
      this->setMemory(1.0, Dimensions::Transform::TRA3D);
   }

   bool WLFlfftBuilder::applicable() const
   {
      return true;
   }
}
}
