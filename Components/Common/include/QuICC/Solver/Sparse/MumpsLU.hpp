/** 
 * @file MumpsLU.hpp
 * @brief Simple typedef to use MumpsLU solver.
 */

#ifndef QUICC_SOLVER_SPARSE_MUMPSLU_HPP
#define QUICC_SOLVER_SPARSE_MUMPSLU_HPP

#include "../External/Interfaces/MumpsLU.hpp"

namespace QuICC {

namespace Solver {

namespace Sparse {

      typedef <class TMatrix>
         using MumpsLU = Eigen::MumpsLU<TMatrix> Type;
}
}
}

#endif // QUICC_SOLVER_SPARSE_MUMPSLU_HPP
