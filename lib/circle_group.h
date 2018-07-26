#pragma once

#include "angle.h"

class circle_group_t {
  public:
    bool is_reflection;
    double x;

    circle_group_t() : is_reflection(false), x(0) {}
    circle_group_t(bool x, double y) : is_reflection(x), x(y) {}

    angle_t act(const angle_t& theta) const {
      if (is_reflection) {
        return angle_t(fmod(2 * M_PI + x - theta.x, 2 * M_PI));
      } else {
        return angle_t(fmod(x + theta.x, 2 * M_PI));
      }
    }

    circle_group_t act(const circle_group_t& r) const {
      if (is_reflection) {
        return circle_group_t(!r.is_reflection, fmod(2 * M_PI + x - r.x, 2 * M_PI));
      } else {
        return circle_group_t(r.is_reflection, fmod(x + r.x, 2 * M_PI));
      }
    }

    angle_t act_inverse(const angle_t& theta) const {
      if (is_reflection) {
        return act(theta);
      } else {
        return angle_t(fmod(2 * M_PI + theta.x - x, 2 * M_PI));
      }
    }

    circle_group_t act_inverse(const circle_group_t& r) const {
      if (is_reflection) {
        return act(r);
      } else {
        return circle_group_t(r.is_reflection, fmod(2 * M_PI + r.x - x, 2 * M_PI));
      }
    }
};


