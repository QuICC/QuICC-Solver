/** 
 * @file TimeSchemeSelector.hpp
 * @brief Definition of some useful typedefs for the timestep scheme
 */

#ifndef QUICC_TIMESTEP_TIMESCHEMESELECTOR_HPP
#define QUICC_TIMESTEP_TIMESCHEMESELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

// Configure code to use ImExRKCB2
#ifdef QUICC_TIMESTEPPER_IMEXRKCB2

   #include "QuICC/Timestep/ImExRKCB2.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB2 TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB2

// Configure code to use ImExRKCB3a
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3A

   #include "QuICC/Timestep/ImExRKCB3a.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB3a TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB3A

// Configure code to use ImExRKCB3b
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3B

   #include "QuICC/Timestep/ImExRKCB3b.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB3b TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB3B

// Configure code to use ImExRKCB3c
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3C

   #include "QuICC/Timestep/ImExRKCB3c.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB3c TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB3C

// Configure code to use ImExRKCB3d
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3D

   #include "QuICC/Timestep/ImExRKCB3d.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB3d TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB3D

// Configure code to use ImExRKCB3e
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3E

   #include "QuICC/Timestep/ImExRKCB3e.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB3e TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB3E

// Configure code to use ImExRKCB3f
#ifdef QUICC_TIMESTEPPER_IMEXRKCB3F

   #include "QuICC/Timestep/ImExRKCB3f.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB3f TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB3F

// Configure code to use ImExRKCB4
#ifdef QUICC_TIMESTEPPER_IMEXRKCB4

   #include "QuICC/Timestep/ImExRKCB4.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRKCB4 TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRKCB4

// Configure code to use ImExRK3
#ifdef QUICC_TIMESTEPPER_IMEXRK3

   #include "QuICC/Timestep/ImExRK3.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExRK3 TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXRK3

// Configure code to use ImExSBDF2
#ifdef QUICC_TIMESTEPPER_IMEXSBDF2

   #include "QuICC/Timestep/ImExSBDF2.hpp"

   namespace QuICC {

      namespace Timestep {

         typedef ImExSBDF2 TimeSchemeSelector;

      }
   }
#endif //QUICC_TIMESTEPPER_IMEXSBDF2

#endif // QUICC_TIMESTEP_TIMESCHEMESELECTOR_HPP
