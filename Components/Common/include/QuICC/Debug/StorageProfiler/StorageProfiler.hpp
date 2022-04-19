/**
 * @file StorageProfiler.hpp
 * @brief Implementation of a storage profiler
 *
 *  \mhdTodo Cleanup the StorageProfiler class similarly to the Profiler class
 */

#ifndef QUICC_DEBUG_STORAGEPROFILER_HPP
#define QUICC_DEBUG_STORAGEPROFILER_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Debug {

   /**
    * @brief Implementation of a storage profiler
    */
   class StorageProfiler
   {
      public:
         // Index for level 0 measure point
         static const unsigned int LVL0 = 1000;
         // Index for level 1 measure point
         static const unsigned int LVL1 = LVL0/10;
         // Index for level 2 measure point
         static const unsigned int LVL2 = LVL1/10;
         // Index for level 3 measure point
         static const unsigned int LVL3 = LVL2/10;

         /**
          * @brief List of storage locations
          */
         enum StoragePoint {
            VARIABLES = 1*LVL0,
               VARIABLESSPEC = VARIABLES + 1*LVL1,
               VARIABLESPHYS = VARIABLES + 2*LVL1,
            TRANSFORMS = 2*LVL0,
               TRASHELLCHEBYSHEV = TRANSFORMS + 1*LVL1,
               TRAANNULUSCHEBYSHEV = TRANSFORMS + 2*LVL1,
               TRACARTESIANCHEBYSHEV = TRANSFORMS + 3*LVL1,
               TRASPHEREWORLAND = TRANSFORMS + 4*LVL1,
               TRACYLINDERWORLAND = TRANSFORMS + 5*LVL1,
               TRAALEGENDRE = TRANSFORMS + 6*LVL1,
               TRAMIXEDFFT = TRANSFORMS + 7*LVL1,
               TRACOMPLEXFFT = TRANSFORMS + 8*LVL1,
            TEMPORARIES = 3*LVL0,
               TEMPTRA1D = TEMPORARIES + 1*LVL1,
               TEMPTRA2D = TEMPORARIES + 2*LVL1,
               TEMPTRA3D = TEMPORARIES + 3*LVL1,
            EQUATIONS = 4*LVL0,
               TRIVIALEQUATION = EQUATIONS + 1*LVL1,
               DIAGNOSTICEQUATION = EQUATIONS + 2*LVL1,
               PROGNOSTICEQUATION = EQUATIONS + 3*LVL1,
            MPI = 5*LVL0,
               MPIBUFFERS = MPI + 1*LVL1,
               MPITYPES = MPI + 2*LVL1,
                  MPIFTYPES = MPITYPES + 1*LVL2,
                  MPIBTYPES = MPITYPES + 2*LVL2,
               MPICOMM = MPI + 3*LVL1,
            IO = 6*LVL0,
         };

         /**
          * @brief Initialise the profiler base
          */
         static void init();

         /**
          * @brief Get memory requirement
          *
          * @param point measure point
          */
         static MHDFloat requirement(StoragePoint point);

         /**
          * @brief Reset the profiling timings
          */
         static void reset();

         /**
          * @brief Check if measure point is active
          *
          * @param point measure point
          */
         static bool isActive(StoragePoint point);

         /**
          * @brief Get map of used points to name
          */
         static std::map<StoragePoint,std::string> namePoints();

         /**
          * @brief Create enum ID to name string map
          *
          * @param map  Storage for the map to build
          */
         static std::string pointName(StoragePoint point);

         /**
          * @brief Update required storage
          *
          * @param point   Storage location
          * @param memory  Required storage
          */
         static void update(StoragePoint point, MHDFloat memory);

         /**
          * @brief Analyze the measured storage requirements
          *
          * @param min  Array of MPI minimal values
          * @param max  Array of MPI maximal values
          */
         static void analyze(Array& reqs, Array& min, Array& max);

         /**
          * @brief Get number of measures
          */
         static size_t size();

         /**
          * @brief Get the requirements
          *
          * @param reqs Requirements for the measure points
          */
         static void getRequirements(Array& reqs);

      protected:
         /**
          * @brief Constructor
          */
         StorageProfiler();

         /**
          * @brief Destructor
          */
         ~StorageProfiler();


      private:
         /**
          * @brief Profiling timings
          */
         static const unsigned int mActiveLvl;

         /**
          * @brief Storage requirements
          */
         static std::map<StoragePoint, MHDFloat> mRequirements;
   };

}
}

#endif // QUICC_DEBUG_STORAGEPROFILER_HPP
