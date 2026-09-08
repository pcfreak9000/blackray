#pragma once

#include <cmath>

static constexpr Real a1 = 1.0 / 4.0;
static constexpr Real b1 = 3.0 / 32.0;
static constexpr Real b2 = 9.0 / 32.0;
static constexpr Real c1 = 1932.0 / 2197.0;
static constexpr Real c2 = -7200.0 / 2197.0;
static constexpr Real c3 = 7296.0 / 2197.0;
static constexpr Real d1 = 439.0 / 216.0;
static constexpr Real d2 = -8.0;
static constexpr Real d3 = 3680.0 / 513.0;
static constexpr Real d4 = -845.0 / 4104.0;
static constexpr Real e1 = -8.0 / 27.0;
static constexpr Real e2 = 2.0;
static constexpr Real e3 = -3544.0 / 2565.0;
static constexpr Real e4 = 1859.0 / 4104.0;
static constexpr Real e5 = -11.0 / 40.0;
static constexpr Real f1 = 25.0 / 216.0;
static constexpr Real f2 = 0.0;
static constexpr Real f3 = 1408.0 / 2565.0;
static constexpr Real f4 = 2197.0 / 4104.0;
static constexpr Real f5 = -1.0 / 5.0;
static constexpr Real g1 = 16.0 / 135.0;
static constexpr Real g2 = 0.0;
static constexpr Real g3 = 6656.0 / 12825.0;
static constexpr Real g4 = 28561.0 / 56430.0;
static constexpr Real g5 = -9.0 / 50.0;
static constexpr Real g6 = 2.0 / 55.0;

static constexpr Real errmin = 1.0e-8;
static constexpr Real errmax = 1.0e-6;

template<size_t N>
void adaptiveRK45(const Real &b, Real *const vars, Real &h,
    const bool &freeze_h, void (*diffeqs)(const Real&, const Real *const, Real[])) {
  //evolve the system one step such that the error stays below certain limits. For this, adaptively increase or decrease step size
  Real diffs[N], vars_4th[N], vars_temp[N], k1[N], k2[N], k3[N], k4[N], k5[N]; //, k6[5];
  do {
    int check = 0;

    /* ----- compute RK1 ----- */
    //maybe this diffeqs only needs computation once and not everytime to adjust h?
    diffeqs(b, vars, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k1[i] = h * diffs[i];
      vars_temp[i] = vars[i] + a1 * k1[i];
    }

    /* ----- compute RK2 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k2[i] = h * diffs[i];
      vars_temp[i] = vars[i] + b1 * k1[i] + b2 * k2[i];
    }

    /* ----- compute RK3 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k3[i] = h * diffs[i];
      vars_temp[i] = vars[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i];
    }

    /* ----- compute RK4 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k4[i] = h * diffs[i];
      vars_temp[i] = vars[i] + d1 * k1[i] + d2 * k2[i] + d3 * k3[i]
          + d4 * k4[i];
    }

    /* ----- compute RK5 ----- */

    diffeqs(b, vars_temp, diffs);
#pragma omp simd
    for (size_t i = 0; i < N; i++) {
      k5[i] = h * diffs[i];
      vars_temp[i] = vars[i] + e1 * k1[i] + e2 * k2[i] + e3 * k3[i] + e4 * k4[i]
          + e5 * k5[i];
    }

    /* ----- compute RK6 ----- */

    diffeqs(b, vars_temp, diffs);
//this for loop is integrated into the next for-loop because it is redundant
//    for (int i = 0; i <= 4; i++)
//      k6[i] = h * diffs[i];

    /* ----- local error ----- */

    for (size_t i = 0; i < N; i++) {
      vars_4th[i] = vars[i] + f1 * k1[i] + f2 * k2[i] + f3 * k3[i] + f4 * k4[i]
          + f5 * k5[i];
      //vars_5th[i] =
      Real varfith = vars[i] + g1 * k1[i] + g2 * k2[i] + g3 * k3[i] + g4 * k4[i]
          + g5 * k5[i] + g6 * diffs[i] * h; //k6[i];

      Real err = std::fabs(
          (vars_4th[i] - varfith) / std::max(vars_4th[i], vars[i]));

      if (err > errmax && !freeze_h) check = 1;
      else if (err < errmin && check != 1 && !freeze_h) check = -1;
    }

    if (check == 1) {
      h /= 2.0;
    } else {
      if (check == -1) h *= 2.0;
      //apply the new step to the variables
      for (size_t i = 0; i < N; i++) {
        vars[i] = vars_4th[i];
      }
      break;
    }

  } while (true);
}
