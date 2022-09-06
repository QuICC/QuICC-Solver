/** 
 * @file MpiConverterImpl.hpp
 * @brief Selector for configured MPI converter implementation
 */

#ifndef QUICC_FRAMEWORK_SELECTOR_MPICONVERTERIMPL_HPP
#define QUICC_FRAMEWORK_SELECTOR_MPICONVERTERIMPL_HPP

#if defined QUICC_MPICOMM_SENDRECV
   #include "QuICC/Communicators/Converters/MpiConverterSendRecv.hpp"
#elif defined QUICC_MPICOMM_ALLTOALL
   #include "QuICC/Communicators/Converters/MpiConverterAllToAll.hpp"
#endif //defined QUICC_MPICOMM_SENDRECV

namespace QuICC {

namespace Framework {

namespace Selector {

#if defined QUICC_MPICOMM_SENDRECV
   typedef Parallel::MpiConverterSendRecv MpiConverterImpl;
#elif defined QUICC_MPICOMM_ALLTOALL
   typedef Parallel::MpiConverterAllToAll MpiConverterImpl;
#endif //defined QUICC_MPICOMM_SENDRECV

}
}
}

#endif // QUICC_FRAMEWORK_SELECTOR_MPICONVERTERIMPL_HPP
