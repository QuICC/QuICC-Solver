/** 
 * @file IScheme.hpp
 * @brief Interface for a generic implicit/explicit Runge-Kutta scheme (Cavaglieri & Bewley, 2015)
 */

#ifndef QUICC_TIMESTEP_ISCHEME_HPP
#define QUICC_TIMESTEP_ISCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Generic interface for timestepping scheme
    */
   class IScheme
   {
      public:
         /**
          * @brief Constructor
          */
         IScheme();

         /**
          * @brief Destructor
          */
         virtual ~IScheme();

         /**
          * @brief Number of substeps for final step (this is +1 compared to theoretical value due to implementation)
          */
         virtual int steps() const = 0; 

         /**
          * @brief Order of the scheme
          */
         virtual int order() const = 0;

         /**
          * @brief Name of the scheme
          */
         virtual std::string name() const = 0;

         /**
          * @brief Initialize Butcher's tableau
          */
         virtual void init() = 0;

         /**
          * @brief Scheme has embedded lower order scheme?
          */
         virtual bool hasEmbedded() const = 0;

         /**
          * @brief Use scheme's embedded lower order scheme?
          */
         virtual bool useEmbedded() const;

         /**
          * @brief Enable embedded scheme
          */
         void enableEmbedded();
         
      protected:

      private:
         /**
          * @brief Use embedded scheme?
          */
         bool mUseEmbedded;

   };

   /// Typedef for a shared pointer Scheme
   typedef std::shared_ptr<IScheme> SharedIScheme;

}
}

#endif // QUICC_TIMESTEP_ISCHEME_HPP
