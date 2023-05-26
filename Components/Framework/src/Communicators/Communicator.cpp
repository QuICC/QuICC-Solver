/**
 * @file Communicator.cpp
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Communicators/Communicator.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   Communicator::Communicator()
   {
   }

   Communicator::~Communicator()
   {
   }

   void Communicator::init(const Dimensions::Transform::Id id, SharedFwdSetupType spSetupFwd, SharedBwdSetupType spSetupBwd)
   {
      // Initialise storage for id
      this->addStorage(id);
      this->storage(id).init(spSetupFwd, spSetupBwd);

      #ifdef QUICC_STORAGEPROFILE
         MHDFloat mem = this->storage(id).requiredStorage();
         switch(id)
         {
            case Dimensions::Transform::TRA1D:
               StorageProfilerMacro_update(Debug::StorageProfiler::TEMPTRA1D, mem);
               break;
            case Dimensions::Transform::TRA2D:
               StorageProfilerMacro_update(Debug::StorageProfiler::TEMPTRA2D, mem);
               break;
            case Dimensions::Transform::TRA3D:
               StorageProfilerMacro_update(Debug::StorageProfiler::TEMPTRA3D, mem);
               break;
            case Dimensions::Transform::SPECTRAL:
               StorageProfilerMacro_update(Debug::StorageProfiler::TEMPSPECTRAL, mem);
               break;
         }
         StorageProfilerMacro_update(Debug::StorageProfiler::TEMPORARIES, mem);
      #endif // QUICC_STORAGEPROFILE
   }

   void Communicator::initConverter(SharedResolution spRes, const std::vector<ArrayI>& packs, Splitting::Locations::Id split)
   {
      /////////////////////////////////////////////////////////////////////////////////
      // Initialise spectral/1D converter
      //
      #ifdef QUICC_MPI
         if(split == Splitting::Locations::FIRST || split == Splitting::Locations::BOTH)
         {
            this->createMpiConverter<Dimensions::Transform::TRA1D>(spRes, packs.at(0), packs.at(1));
         }
         else
         {
            this->createSerialConverter<Dimensions::Transform::TRA1D>(spRes);
         }
      #else
         this->createSerialConverter<Dimensions::Transform::TRA1D>(spRes);
      #endif
      this->converter(Dimensions::Transform::TRA1D).setup();


      /////////////////////////////////////////////////////////////////////////////////
      // Initialise 1D/2D converter
      //

      //
      // Serial implementation
      //
      if(split == QuICC::Splitting::Locations::NONE)
      {
         this->createSerialConverter<Dimensions::Transform::TRA2D>(spRes);
      }

      //
      // MPI implementation
      //
      #ifdef QUICC_MPI
         // Load splitting has been done on first dimension
         if(split == Splitting::Locations::FIRST || split == Splitting::Locations::COUPLED2D)
         {
            assert(packs.size() > 1);
            this->createMpiConverter<Dimensions::Transform::TRA2D>(spRes, packs.at(0), packs.at(1));

         // Load splitting has been done on second dimension
         } else if(split == Splitting::Locations::SECOND)
         {
            this->createSerialConverter<Dimensions::Transform::TRA2D>(spRes);

         // Load splitting has been done on two dimensions
         } // else if(split == Splitting::Locations::BOTH)
         // {
         //    Initialisation of this part is done at a higher level to optimise storage/buffer use
         // }
      #endif // QUICC_MPI

      // If only one dimension is split. Else the setup will be done at a later stage to optimise memory/buffer usage.
      if(split != Splitting::Locations::BOTH)
      {
         // Setup converter
         this->converter(Dimensions::Transform::TRA2D).setup();
      }

      /////////////////////////////////////////////////////////////////////////////////
      // Initialise 2D/3D converter
      //

      //
      // Serial implementation
      //
      if(split == QuICC::Splitting::Locations::NONE)
      {
         // Initialise serial 3D converter
         this->createSerialConverter<Dimensions::Transform::TRA3D>(spRes);
      }

      //
      // MPI implementation
      //
#ifdef QUICC_MPI
      // Load splitting has been done on first dimension
      if(split == Splitting::Locations::FIRST || split == Splitting::Locations::COUPLED2D)
      {
         this->createSerialConverter<Dimensions::Transform::TRA3D>(spRes);

      // Load splitting has been done on second dimension
      } else if(split == Splitting::Locations::SECOND)
      {
         assert(packs.size() == 4);
         this->createMpiConverter<Dimensions::Transform::TRA3D>(spRes, packs.at(2), packs.at(3));

         // Load splitting has been done on two dimensions
      } else if(split == Splitting::Locations::BOTH)
      {
         assert(packs.size() == 4);
         this->createMpiConverter<Dimensions::Transform::TRA2D>(spRes, packs.at(0), packs.at(1));
         this->createMpiConverter<Dimensions::Transform::TRA3D>(spRes, packs.at(2), packs.at(3));
         // A better implementation using only 3 buffers is possible but requires some cross transform synchronization
      }
#endif // QUICC_MPI

      // Setup converter
      this->converter(Dimensions::Transform::TRA3D).setup();

      // If both dimensions are split. In the other cases setup() has already been called.
      if(split == Splitting::Locations::BOTH)
      {
         this->converter(Dimensions::Transform::TRA2D).setup();
      }

#ifdef QUICC_STORAGEPROFILE
      // Do (MPI) storage profiling on 3D converter
      this->converter(Dimensions::Transform::TRA2D).profileStorage();

      // Do (MPI) storage profiling on 3D converter
      this->converter(Dimensions::Transform::TRA3D).profileStorage();
#endif // QUICC_STORAGEPROFILE
   }

}
}
