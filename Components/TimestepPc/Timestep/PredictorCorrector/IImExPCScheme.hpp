/**
 * @file IImExPCScheme.hpp
 * @brief Interface for a generic implicit/explicit predictor-corrector scheme
 */

#ifndef QUICC_TIMESTEP_PREDICTORCORRECTOR_IIMEXPCSCHEME_HPP
#define QUICC_TIMESTEP_PREDICTORCORRECTOR_IIMEXPCSCHEME_HPP

// System includes
//

// Project includes
//
#include "QuICC/Timestep/IScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

namespace PredictorCorrector {

/**
 * @brief Interface of generic of an implicit/explicit predictor-corrector
 * scheme
 */
class IImExPCScheme : public IScheme
{
public:
   /**
    * @brief Constructor
    */
   IImExPCScheme() = default;

   /**
    * @brief Destructor
    */
   virtual ~IImExPCScheme() = default;

   /**
    * @brief Implicit coefficient
    */
   virtual MHDFloat aIm(const int i) const = 0;

   /**
    * @brief Step fraction
    */
   virtual MHDFloat cEx(const int i) const = 0;

protected:
private:
};

/// Typedef for a shared pointer IImExPCScheme
typedef std::shared_ptr<IImExPCScheme> SharedIImExPCScheme;

} // namespace PredictorCorrector
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_PREDICTORCORRECTOR_IIMEXPCSCHEME_HPP
