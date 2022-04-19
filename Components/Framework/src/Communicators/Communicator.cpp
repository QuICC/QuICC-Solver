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
         }
         StorageProfilerMacro_update(Debug::StorageProfiler::TEMPORARIES, mem);
      #endif // QUICC_STORAGEPROFILE
   }

   void Communicator::initConverter(SharedResolution spRes, const std::vector<ArrayI>& packs, Splitting::Locations::Id split)
   {
      /////////////////////////////////////////////////////////////////////////////////
      // Initialise 1D/2D converter
      //

      //
      // MPI code needs to come first
      // Separating out MPI and serial code makes sure Serial version can be compiled without MPI installed
      //
      #ifdef QUICC_MPI
         // Load splitting has been done on first dimension
         if(split == Splitting::Locations::FIRST || split == Splitting::Locations::COUPLED2D)
         {
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
      //
      // Serial implementation
      //
      #else
         this->createSerialConverter<Dimensions::Transform::TRA2D>(spRes);
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
      // MPI code needs to come first
      // Separating out MPI and serial code makes sure Serial version can be compiled without MPI installed
      //
#ifdef QUICC_MPI
      // Load splitting has been done on first dimension
      if(split == Splitting::Locations::FIRST || split == Splitting::Locations::COUPLED2D)
      {
         this->createSerialConverter<Dimensions::Transform::TRA3D>(spRes);

      // Load splitting has been done on second dimension
      } else if(split == Splitting::Locations::SECOND)
      {
         this->createMpiConverter<Dimensions::Transform::TRA3D>(spRes, packs.at(2), packs.at(3));

         // Load splitting has been done on two dimensions
      } else if(split == Splitting::Locations::BOTH)
      {
         this->createMpiConverter<Dimensions::Transform::TRA2D>(spRes, packs.at(0), packs.at(1));
         this->createMpiConverter<Dimensions::Transform::TRA3D>(spRes, packs.at(2), packs.at(3));
         // A better implementation using only 3 buffers is possible but requires some cross transform synchronization
      }
      //
      // Serial implementation
      //
#else
      // Initialise serial 3D converter
      this->createSerialConverter<Dimensions::Transform::TRA3D>(spRes);
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
