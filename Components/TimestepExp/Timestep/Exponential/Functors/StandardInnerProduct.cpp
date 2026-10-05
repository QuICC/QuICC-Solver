/**
 * @file StandardInnerProduct.cpp
 * @brief Source of test functor for matrix A
 */

// System includes
//

// Project includes
//
#include "Timestep/Exponential/Functors/StandardInnerProduct.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

Array StandardInnerProduct::operator()(const Matrix& u, const int i0, const int i1, const Matrix& v, const int j, const int n) const
{
   // Dot product of field values
   int cols = i1 - i0 + 1;
   Array dot = u.block(0, i0, n, cols).transpose()*v.col(j).topRows(n);

#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, dot.data(), dot.size(), Environment::MpiTypes::type<MHDFloat>(), MPI_SUM, MPI_COMM_WORLD);
#endif

   // Add dot product from augmented part
   int p = u.rows() - n;
   dot += u.block(n, i0, p, cols).transpose() * v.col(j).bottomRows(p);

   return dot;
}

MHDFloat StandardInnerProduct::operator()(const Matrix& u, const int i, const Matrix& v, const int j, const int n) const
{
   Array dot = this->operator()(u, i, i, v, j,n);
   assert(dot.rows() == 1);

   return dot(0,0);
}

MHDFloat StandardInnerProduct::norm(const Matrix& u, const int i, const int n) const
{
   MHDFloat norm = u.col(i).topRows(n).squaredNorm();

#ifdef QUICC_MPI
   MPI_Allreduce(MPI_IN_PLACE, &norm, 1, Environment::MpiTypes::type<MHDFloat>(), MPI_SUM, MPI_COMM_WORLD);
#endif

   // Add squared 2-norm of augmented part
   int p = u.rows() - n;
   norm += u.col(i).bottomRows(p).squaredNorm();

   norm = std::sqrt(norm);
   return norm;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
