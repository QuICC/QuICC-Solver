/**
 * @file BuildTimestepMatrixWrapper.cpp
 * @brief Implementation of functors used in Interface
 */

// System includes
//

// Project includes
//
#include "Timestep/Exponential/BuildTimestepMatrix.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/Dispatchers.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/ModelOperator/QuasiInverse.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoBc.hpp"
#include "QuICC/Tag/Operator/Qi.hpp"
#include "QuICC/Tag/Operator/Lhs.hpp"
#ifdef QUICC_DEBUG
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#endif
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Timestep {

namespace Exponential {

void buildTimestepMatrix(std::map<std::size_t, DecoupledZSparse>& ops, FieldComponents::Spectral::Id comp, const int idx, const Resolution& res, const Model::IModelBackend& backend, const Equations::CouplingInformation cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams)
{
   auto buildOp = [&](const std::size_t opId, const std::size_t tId, const std::size_t bcId)
   {
      auto ret = ops.insert(std::make_pair(tId, DecoupledZSparse()));
      auto& mat = ret.first->second;
      Equations::dispatchModelMatrix(mat, opId, comp, idx, bcId, res, backend, cinfo, bcIds, eqParams);
   };

   using namespace ModelOperator;
   using namespace ModelOperatorBoundary;
   // Compute model's mass operator (with boundary conditions)
   buildOp(Time::id(), Tag::Operator::Lhs::id(), SolverNoTau::id());

   // Compute model's quasi-inverse operator (without boundary conditions)
   buildOp(QuasiInverse::id(), Tag::Operator::Qi::id(), SolverNoBc::id());
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
