/**
 * @file BuildTimestepMatrixWrapper.hpp
 * @brief Implementation of functors used in Interface
 */

#pragma once

// System includes
//
#include <memory>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
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

namespace QuICC {

namespace Timestep {

namespace Exponential {

/**
 * @brief Wrapper to build timestepping matrices
 */
void buildTimestepMatrixWrapper(std::map<std::size_t, DecoupledZSparse>& ops, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

inline void buildTimestepMatrixWrapper(std::map<std::size_t, DecoupledZSparse>& ops, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp,
   const int idx)
{
   auto buildOp = [&](const std::size_t opId, const std::size_t tId, const std::size_t bcId)
   {
      auto ret = ops.insert(std::make_pair(tId, DecoupledZSparse()));
      spEq->buildModelMatrix(ret.first->second, opId, comp, idx, bcId);
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
