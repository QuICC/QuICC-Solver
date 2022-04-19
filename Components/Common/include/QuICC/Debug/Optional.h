#ifndef QUICC_DEBUG_OPTIONAL_H
#define QUICC_DEBUG_OPTIONAL_H

#ifdef QUICC_OPTIONAL_ENABLED
   #define QUICC_OPTIONAL(x) x
#else
   #define QUICC_OPTIONAL(x)
#endif

#undef QUICC_OPTIONAL_ENABLED

// Level 1
#ifdef QUICC_OPTIONAL_LVL1_ENABLED
   #define QUICC_OPTIONAL_LVL1(x) x
#else
   #define QUICC_OPTIONAL_LVL1(x)
#endif

#undef QUICC_OPTIONAL_LVL1_ENABLED

// Level 2
#ifdef QUICC_OPTIONAL_LVL2_ENABLED
   #define QUICC_OPTIONAL_LVL2(x) x
#else
   #define QUICC_OPTIONAL_LVL2(x)
#endif

#undef QUICC_OPTIONAL_LVL2_ENABLED

// Level 3
#ifdef QUICC_OPTIONAL_LVL3_ENABLED
   #define QUICC_OPTIONAL_LVL3(x) x
#else
   #define QUICC_OPTIONAL_LVL3(x)
#endif

#undef QUICC_OPTIONAL_LVL3_ENABLED

// Level 4
#ifdef QUICC_OPTIONAL_LVL4_ENABLED
   #define QUICC_OPTIONAL_LVL4(x) x
#else
   #define QUICC_OPTIONAL_LVL4(x)
#endif

#undef QUICC_OPTIONAL_LVL4_ENABLED

#endif // QUICC_DEBUG_OPTIONAL_H
