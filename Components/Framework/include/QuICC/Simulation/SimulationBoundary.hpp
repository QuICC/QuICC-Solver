/**
 * @file SimulationBoundary.hpp
 * @brief Implementation of a simple simulation wide boundary condition interface
 */

#ifndef QUICC_SIMULATIONBOUNDARY_HPP
#define QUICC_SIMULATIONBOUNDARY_HPP

// Configuration includes
//

// System includes
//
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/IPhysicalName.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a simple simulation wide boundary condition interface
    */
   class SimulationBoundary
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationBoundary(const std::map<std::string,std::size_t>& bcIds);

         /**
          * @brief Destructor
          */
         ~SimulationBoundary();

         /**
          * @brief Get tag map
          */
         const std::map<std::size_t,std::size_t>& map() const;

         /**
          * @brief Get tag map
          */
         std::size_t bcId(const std::size_t id) const;

      protected:

      private:
         /**
          * @brief Convert tag map to ID map
          *
          * @param bcIds   Tag map
          */
         void convert(const std::map<std::string,std::size_t>& bcIds);

         /**
          * @brief Storage for the boundary conditions
          */
         std::map<std::size_t,std::size_t> mBcs;
   };

   /// Typedef for a shared pointer to a SimulationBoundary object
   typedef std::shared_ptr<SimulationBoundary>   SharedSimulationBoundary;
}

#endif // QUICC_SIMULATIONBOUNDARY_HPP
