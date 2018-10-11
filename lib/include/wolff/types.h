
#include <inttypes.h>

typedef uint_fast32_t v_t;
typedef uint_fast8_t q_t;
typedef uint_fast16_t R_t;
typedef uint_fast8_t D_t;
typedef uint_fast16_t L_t;
typedef uint_fast64_t count_t;
typedef int_fast64_t h_t;

#define MAX_v UINT_FAST32_MAX
#define MAX_Q UINT_FAST8_MAX
#define MAX_R UINT_FAST16_MAX
#define MAX_D UINT_FAST8_MAX
#define MAX_L UINT_FAST16_MAX
#define MAX_COUNT UINT_FAST64_MAX
#define MAX_H INT_FAST64_MAX
#define MIN_H INT_FAST64_MIN

#define PRIv PRIuFAST32
#define PRIq PRIuFAST8
#define PRIR PRIuFAST16
#define PRID PRIuFAST8
#define PRIL PRIuFAST16
#define PRIcount PRIuFAST64
#define PRIh PRIdFAST64

#define SCNv SCNuFAST32
#define SCNq SCNuFAST8
#define SCNR SCNuFAST16
#define SCND SCNuFAST8
#define SCNL SCNuFAST16
#define SCNcount SCNuFAST64
#define SCNh SCNdFAST64

