#pragma once

#include "types.h"

double hue_to_R(double theta) {
  if (((M_PI / 3 <= theta) && (theta < 2 * M_PI / 3)) || ((4 * M_PI / 3 <= theta) && (theta < 5 * M_PI / 3))) {
    return 1.0 - fabs(fmod(theta / (2 * M_PI / 6), 2) - 1.0);
  } else if (((0 <= theta) && (theta < M_PI / 3)) || ((5 * M_PI / 3 <= theta) && (theta <= 2 * M_PI))) {
    return 1.0;
  } else {
    return 0.0;
  }
}

double hue_to_G(double theta) {
  if (((0 <= theta) && (theta < M_PI / 3)) || ((M_PI <= theta) && (theta < 4 * M_PI / 3))) {
    return 1.0 - fabs(fmod(theta / (2 * M_PI / 6), 2) - 1.0);
  } else if (((M_PI / 3 <= theta) && (theta < 2 * M_PI / 3)) || ((2 * M_PI / 3 <= theta) && (theta < M_PI))) {
    return 1.0;
  } else {
    return 0.0;
  }
}

double hue_to_B(double theta) {
  if (((2 * M_PI / 3 <= theta) && (theta < M_PI)) || ((5 * M_PI / 3 <= theta) && (theta <= 2 * M_PI))) {
    return 1.0 - fabs(fmod(theta / (2 * M_PI / 6), 2) - 1.0);
  } else if (((M_PI <= theta) && (theta < 4 * M_PI / 3)) || ((4 * M_PI / 3 <= theta) && (theta < 5 * M_PI / 3))) {
    return 1.0;
  } else {
    return 0.0;
  }
}

