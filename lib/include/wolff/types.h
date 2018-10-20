
#ifndef WOLFF_TYPES_H
#define WOLFF_TYPES_H

#include <inttypes.h>

typedef uint_fast32_t v_t;     // vertex and edge indices
typedef uint_fast8_t q_t;      // enumerated states
typedef uint_fast8_t D_t;      // dimension
typedef uint_fast16_t L_t;     // linear size
typedef uint_fast64_t N_t;     // cycle iterator

#define MAX_v UINT_FAST32_MAX
#define MAX_q UINT_FAST8_MAX
#define MAX_D UINT_FAST8_MAX
#define MAX_L UINT_FAST16_MAX
#define MAX_N UINT_FAST64_MAX

#define PRIv PRIuFAST32
#define PRIq PRIuFAST8
#define PRID PRIuFAST8
#define PRIL PRIuFAST16
#define PRIN PRIuFAST64

#define SCNv SCNuFAST32
#define SCNq SCNuFAST8
#define SCND SCNuFAST8
#define SCNL SCNuFAST16
#define SCNN SCNuFAST64

#endif

