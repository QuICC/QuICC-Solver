/**
 * @file SparseSolverBase.hpp
 * @brief Implementation of the base for the solver structures
 */

#ifndef QUICC_SOLVER_SPARSESOLVERBASE_HPP
#define QUICC_SOLVER_SPARSESOLVERBASE_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of the base for the solver structures
    */
   class SparseSolverBase
   {
      public:
         /// Typedef to simplify notation for the field data
         typedef std::vector<SpectralFieldId> FieldIdVector;

         /// Typedef for an iterator for the field data
         typedef FieldIdVector::const_iterator  FieldId_iterator;

         /// Typedef for a range iterator for the field coupling data
         typedef std::pair<FieldId_iterator,FieldId_iterator>  FieldId_range;

         /**
          * @brief Constructor
          *
          * @param start   Starting index (for example without m=0)
          * @param timeId  Solver timing with respect to timestepping
          */
         SparseSolverBase(const int start, const std::size_t timeId);

         /**
          * @brief Destructor
          */
         virtual ~SparseSolverBase();

         /**
          * @brief Add storage information
          *
          * @param id         Field ID
          * @param idx        Field index in solver matrix
          * @param startRow   Sizes required to compute start row sizes
          */
         void addInformation(const SpectralFieldId& id, const int idx, const ArrayI& startRow);

         /**
          * @brief Initialise the startRow based on added information
          */
         void initStartRow();

         /**
          * @brief Get start row
          */
         int startRow(const SpectralFieldId& id, const int i) const;

         /**
          * @brief Range of stored fields
          */
         FieldId_range fieldRange() const;

         /**
          * @brief Is operator initialized?
          */
         bool isInitialized() const;

         /**
          * @brief Set operator to initialized
          */
         void setInitialized();

         /**
          * @brief Solve timing
          */
         std::size_t solveTiming() const;

      protected:
         /**
          * @brief Starting index
          */
         int mZeroIdx;

         /**
          * @brief Solver timing
          */
         std::size_t   mSolveTiming;

         /**
          * @brief Storage for the field Ids
          */
         std::vector<SpectralFieldId>   mFieldIds;

         /**
          * @brief Storage for the storage information
          */
         std::map<SpectralFieldId, std::pair<int, ArrayI> > mInformation;

         /**
          * @brief Flag for operator initialization
          */
         bool mIsInitialized;

      private:
   };

   /// Typedef for a shared pointer of a SparseSolverBase
   typedef std::shared_ptr<SparseSolverBase>  SharedSparseSolverBase;
}
}

#endif // QUICC_SOLVER_SPARSESOLVERBASE_HPP
