
#include <convex.h>
#include <wstp.h>

extern int WSMain(int, char **);

void convexminorant(double *, int);

:Begin:
:Function:      convexminorant
:Pattern:       GetConvexMinorant[ list:{___Real} ]
:Arguments:     { list }
:ArgumentTypes: { Real64List }
:ReturnType:    Manual
:End:

:Evaluate:      GetConvexMinorant[ sequence___Float]:= GetConvexMinorant[ {sequence} ]

void convexminorant(double * Gammas, int len) {
  int i;
  for (i = 0; i < len; i++) {
    if (Gammas[i] <= 0) {
       break;
    }
  }
  double *m = get_convex_minorant(i, Gammas);
  WSPutReal64List(stdlink, m, i);
  free(m);
}

int main(int argc, char **argv) {
  return WSMain(argc, argv);
}

