/**
 * @file EquationParameters.hpp
 * @brief Definition of the non-dimensional parameters
 */

#ifndef QUICC_EQUATIONS_EQUATIONPARAMETERS_HPP
#define QUICC_EQUATIONS_EQUATIONPARAMETERS_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/NonDimensional/INumber.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief This class provides the constant coefficients apparearing in the equations
    */
   class EquationParameters
   {
      public:
         // Typedef for the nondimensional number map
         typedef std::map<std::size_t,NonDimensional::SharedINumber> NDMapType;

         /**
          * @brief Constructor
          */
         EquationParameters();

         /**
          * @brief Destructor
          */
         virtual ~EquationParameters();

         /**
          * @brief Get the IDs of the equation parameters
          */
         std::vector<std::size_t>  ids() const;

         /**
          * @brief Get the names of the equation parameters
          */
         std::vector<std::string>  names() const;

         /**
          * @brief Initialise the values from given parameters
          *
          * @param parameters Parameter values read from configuration
          */
         void init(const std::map<std::string, MHDFloat>& parameters);

         /**
          * @brief Get a nondimensional parameter
          *
          * @param name Name of the parameter
          */
         MHDFloat nd(std::size_t id) const;

         /**
          * @brief Get full map of parameters
          */
         const NDMapType& map() const;

         /**
          * @brief Update full map of parameters
          */
         NDMapType& rMap();

      protected:
         /**
          * @brief Storage for the nondimensional parameters
          */
         NDMapType   mND;

      private:
   };

   /// Typedef for a shared pointer to an EquationParameters object
   typedef std::shared_ptr<EquationParameters>   SharedEquationParameters;
   typedef std::shared_ptr<const EquationParameters>   SharedCEquationParameters;
}
}

#endif // QUICC_EQUATIONS_EQUATIONPARAMETERS_HPP
