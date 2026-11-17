/**
 * @file ExpEuler.hpp
 * @brief Implementation of an exponential Euler scheme of
 * order 2
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_EXPEULER_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_EXPEULER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Timestep/Exponential/IExpScheme.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

/**
 * @brief Implementation of an implicit/explicit Euler scheme of
 * order 1
 */
class ExpEuler: public IExpScheme
{
public:
   /**
    * @brief Constructor
    */
   ExpEuler();

   /**
    * @brief Destructor
    */
   virtual ~ExpEuler() = default;

   /**
    * @brief Number of substeps for final step (this is +1 compared to
    * theoretical value due to implementation)
    */
   int steps() const final;

   /**
    * @brief Order of the scheme
    */
   int order() const final;

   /**
    * @brief Scheme has embedded lower order scheme?
    */
   bool hasEmbedded() const final;

   /**
    * @brief Name of the scheme
    */
   std::string name() const final;

protected:

private:
};

/// Typedef for a shared pointer ExpEuler
typedef std::shared_ptr<ExpEuler> SharedExpEuler;

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_EXPEULER_HPP
