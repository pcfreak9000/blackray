#define DEBUG_DIV 1000.0
#define MAX_ITER 3000
#define NO_INTERSECT -1
#define INTERSECT 0
#define TOO_MANY_INTERSECT -2

#define SQR(x) ((x)*(x))

long double interpolate(long double a, long double b, long double f) {
  return (1.0 - f) * a + f * b;
}

long double checkIntersect(long double x1, long double y1, long double x2,
    long double y2, long double x3, long double y3, long double x4,
    long double y4) {
  long double num = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
  long double denum = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
  long double t = num / denum;
  if (t >= 0 && t <= 1) {
    return t;
  }
  return NO_INTERSECT;
}

//finds any intersection. it is not guranteed that it is the closest one, i.e. the first hit. Sufficiently small stepsize should circumvent the problem,
//as well as retaking the step with a smaller stepsize.
int get_interpolated_sp(const long double x1, const long double y1,
    const long double x2, const long double y2, const SurfacePoint *diskdata,
    const size_t ddsize, SurfacePoint &out, bool mirrory) {
  for (size_t i = 0; i < ddsize - 1; i++) {
    SurfacePoint p0 = diskdata[i];
    SurfacePoint p1 = diskdata[i + 1];
    long double mul = mirrory ? -1 : 1;
    long double pt = checkIntersect(p0.x, mul * p0.y / DEBUG_DIV, p1.x,
        mul * p1.y / DEBUG_DIV, x1, y1, x2, y2);
    int intersects = 0;
    if (pt != NO_INTERSECT) {
      long double xi = p0.x + pt * (p1.x - p0.x);
      long double yi = mul * (p0.y + pt * (p1.y - p0.y));
      out.x = xi;
      out.y = yi / DEBUG_DIV;
      out.u0 = interpolate(p0.u0, p1.u0, pt);
      out.u1 = interpolate(p0.u1, p1.u1, pt);
      out.u2 = interpolate(p0.u2, p1.u2, pt);
      out.u3 = interpolate(p0.u3, p1.u3, pt);
      return INTERSECT;
    }
  }
  return NO_INTERSECT;
}

void raytrace(long double xobs, long double yobs, long double iobs,
    long double rin, long double disk_length_combined, RayHit &hit,
    int &stop_integration, const SurfacePoint *diskdata, const size_t ddsize) {
  long double dobs;
  long double xobs2, yobs2;
  long double atol, rtol;
  long double hstart;
  long double t0, r0, th0, phi0;
  long double kt0, kr0, kth0, kphi0;
  long double r02, s0, s02;
  long double fact1, fact2, fact3;
  long double t, r, th, phi;
  long double kt, kr, kth, kphi;
  long double rau, thau, phiau, krau, kthau;
  long double kyem;
  long double const0, const1;
  long double v1, v2;
  long double h, hnext;
  long double Delta;
  long double spin2 = spin * spin;

  long double carter, cosem, c02;
  long double b;

  long double met[4][4];
  long double diffs[5], vars[5], vars_temp[5], vars_4th[5], vars_5th[5], k1[5],
      k2[5], k3[5], k4[5], k5[5], k6[5];
  long double xem[4];
  long double gfactor;
  long double err, errmin, errmax;

  int check, check2 = 0;
  int i;

  int div;

  /* ----- Set computational parameters ----- */
  dobs = 1.0e+8; /* distance of the observer */
  errmin = 1.0e-8;
  errmax = 1.0e-6;
  atol = 1.0e-10;
  //rtol = 1.0e-10;
  rtol = 1.0e-8;
  long double thtol = 1.0e-8;
  int count, iter;

  hstart = -1.0;

  /* ----- compute photon initial conditions ----- */
  xobs2 = xobs * xobs;
  yobs2 = yobs * yobs;

  fact1 = yobs * sin(iobs) + dobs * cos(iobs);
  fact2 = dobs * sin(iobs) - yobs * cos(iobs);

  r02 = xobs2 + yobs2 + dobs * dobs;

  r0 = sqrt(r02);
  th0 = acos(fact1 / r0);
  phi0 = atan2(xobs, fact2);

  s0 = sin(th0);
  s02 = s0 * s0;

  kr0 = dobs / r0;
  kth0 = -(cos(iobs) - dobs * fact1 / r02) / sqrt(r02 - fact1 * fact1);
  kphi0 = -xobs * sin(iobs) / (xobs2 + fact2 * fact2);

  metric(r0, th0, met);

  fact3 = sqrt(
      met[0][3] * met[0][3] * kphi0 * kphi0
          - met[0][0]
              * (met[1][1] * kr0 * kr0 + met[2][2] * kth0 * kth0
                  + met[3][3] * kphi0 * kphi0));

  kt0 = -(met[0][3] * kphi0 + fact3) / met[0][0];

  b = -(met[3][3] * kphi0 + met[0][3] * kt0)
      / (met[0][0] * kt0 + met[0][3] * kphi0);

  kr0 /= fact3;
  kth0 /= fact3;

  /* ----- carter constant ----- */

  c02 = 1. - s02;

  carter = yobs2 - spin2 * c02 + xobs2 * c02;
  carter = sqrt(carter);

  /* ----- solve geodesic equations ----- */

  r = r0;
  th = th0;
  phi = phi0;

  kr = kr0;
  kth = kth0;

  const0 = kt0;
  const1 = r02 * s02 * kphi0 / kt0;

  stop_integration = 0;

  h = hstart;
  count = 0;
  iter = 0;

  long double a1 = 1.0 / 4.0;
  long double b1 = 3.0 / 32.0;
  long double b2 = 9.0 / 32.0;
  long double c1 = 1932.0 / 2197.0;
  long double c2 = -7200.0 / 2197.0;
  long double c3 = 7296.0 / 2197.0;
  long double d1 = 439.0 / 216.0;
  long double d2 = -8.0;
  long double d3 = 3680.0 / 513.0;
  long double d4 = -845.0 / 4104.0;
  long double e1 = -8.0 / 27.0;
  long double e2 = 2.0;
  long double e3 = -3544.0 / 2565.0;
  long double e4 = 1859.0 / 4104.0;
  long double e5 = -11.0 / 40.0;
  long double f1 = 25.0 / 216.0;
  long double f2 = 0.0;
  long double f3 = 1408.0 / 2565.0;
  long double f4 = 2197.0 / 4104.0;
  long double f5 = -1.0 / 5.0;
  long double g1 = 16.0 / 135.0;
  long double g2 = 0.0;
  long double g3 = 6656.0 / 12825.0;
  long double g4 = 28561.0 / 56430.0;
  long double g5 = -9.0 / 50.0;
  long double g6 = 2.0 / 55.0;
  SurfacePoint spi;
  long double prevh=-1.0;
  do {
    iter++;
    vars[0] = r;
    vars[1] = th;
    vars[2] = phi;
    vars[3] = kr;
    vars[4] = kth;
//evolve the system one step such that the error stays below certain limits. For this, adaptively increase or decrease step size
    do {
      check = 0;

      /* ----- compute RK1 ----- */

      diffeqs(b, vars, diffs);
      for (i = 0; i <= 4; i++) {
        k1[i] = h * diffs[i];
        vars_temp[i] = vars[i] + a1 * k1[i];
      }

      /* ----- compute RK2 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k2[i] = h * diffs[i];
        vars_temp[i] = vars[i] + b1 * k1[i] + b2 * k2[i];
      }

      /* ----- compute RK3 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k3[i] = h * diffs[i];
        vars_temp[i] = vars[i] + c1 * k1[i] + c2 * k2[i] + c3 * k3[i];
      }

      /* ----- compute RK4 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k4[i] = h * diffs[i];
        vars_temp[i] = vars[i] + d1 * k1[i] + d2 * k2[i] + d3 * k3[i]
            + d4 * k4[i];
      }

      /* ----- compute RK5 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++) {
        k5[i] = h * diffs[i];
        vars_temp[i] = vars[i] + e1 * k1[i] + e2 * k2[i] + e3 * k3[i]
            + e4 * k4[i] + e5 * k5[i];
      }

      /* ----- compute RK6 ----- */

      diffeqs(b, vars_temp, diffs);
      for (i = 0; i <= 4; i++)
        k6[i] = h * diffs[i];

      /* ----- local error ----- */

      for (i = 0; i <= 4; i++) {
        vars_4th[i] = vars[i] + f1 * k1[i] + f2 * k2[i] + f3 * k3[i]
            + f4 * k4[i] + f5 * k5[i];
        vars_5th[i] = vars[i] + g1 * k1[i] + g2 * k2[i] + g3 * k3[i]
            + g4 * k4[i] + g5 * k5[i] + g6 * k6[i];

        err = fabs((vars_4th[i] - vars_5th[i]) / max(vars_4th[i], vars[i]));

        if (err > errmax && check2 == 0)
          check = 1;
        else if (err < errmin && check != 1 && check2 == 0)
          check = -1;
      }

      if (check == 1)
        h /= 2.0;
      else if (check == -1)
        h *= 2.0;

    } while (check == 1);
    //std::cout << h << std::endl;

    /* ----- solutions to the fourth-order RKN method ----- */
//apply the new step to the variables
    rau = r;
    thau = th;
    phiau = phi;
    krau = kr;
    kthau = kth;

    r = vars_4th[0];
    th = vars_4th[1];
    phi = vars_4th[2];
    kr = vars_4th[3];
    kth = vars_4th[4];

    Delta = r * r - 2.0 * r + spin2;
//check if the new position ends the integration
    if (Delta < 1.e-3) {
      stop_integration = 4; // printf("photon crosses the horizon\n"); /* the
                            // photon crosses the horizon */
      break;
    }
    if (r < 1.0) {
      stop_integration = 5; // printf("photon crosses the horizon\n"); /* the
                            // photon crosses the horizon */
      break;
    }

    if (r != r) {
      stop_integration = 6; // printf("numerical problem\n");          /*
                            // numerical problems! */
      break;
    }

    if (r > 1.05 * dobs) {
      stop_integration = 7; // printf("photon escaped to infinity\n");   /* the
                            // photon escapes to infinity */
      break;
    }
    if (iter > MAX_ITER) {
      stop_integration = 255;
      break;
    }
//check if the new position intersects the accretion disk
    //convert coordinates of current and previous position via BL-cartesian conversion
    long double xcoord = std::sqrt(r * r + spin2) * sin(th);
    long double ycoord = r * cos(th);
    long double xcoordprev = std::sqrt(rau * rau + spin2) * sin(thau);
    long double ycoordprev = rau * cos(thau);
    SurfacePoint spia, spib;
    int resa = get_interpolated_sp(xcoordprev, ycoordprev, xcoord, ycoord,
        diskdata, ddsize, spia, false);
    int resb = get_interpolated_sp(xcoordprev, ycoordprev, xcoord, ycoord,
        diskdata, ddsize, spib, true);
    int res=NO_INTERSECT;
    int index=0;
    resb = NO_INTERSECT;
    if (resa == INTERSECT && resb == INTERSECT) {
      long double dista = std::sqrt(
      SQR(spia.x-xcoordprev) + SQR(spia.y - ycoordprev));
      long double distb = std::sqrt(
      SQR(spib.x-xcoordprev) + SQR(spib.y - ycoordprev));
      if (dista < distb) {
        res = resa;
        spi = spia;
        index = 1;
      } else {
        res = resb;
        spi = spib;
        index = 128;
      }
    } else if (resa == INTERSECT) {
      res = resa;
      spi = spia;
      index = 1;
    } else if (resb == INTERSECT) {
      res = resb;
      spi = spib;
      index = 128;
    } else {
      index = 0;
      res = NO_INTERSECT;
    }
    //deal with (possible) intersection
//    if((cos(th))<0.0){
//      if(r>6.0&&r<50.0){
//        stop_integration = 1;
//      }
//    }
#ifndef xxx
    if (res == INTERSECT) {
     // std::cout << "int" << std::endl;
      if(check2!=1){
        prevh = h;
      }
      check2 = 1;//don't adapt stepsize anymore, this is now done manually to reach certain tolerances
      if (fabs(th - thau) <= thtol && fabs(r - rau) <= rtol) {
        count++;
      }
      if (count > 0) {
        //if(xcoord>=6.0&&xcoord<=50.0){
        if (spi.y > 0.0 && spi.x > 0.0) {
          // next step is redshift calculation with data from the intersection
          // point. we also need the interpolated 4-vel etc
          stop_integration = index;
          //break; /* the photon hits the disk */
        } else {
          check2 = 0;
          count = 0;
          h = prevh;
          //std::cout << "reset" << std::endl;
          // this is a simplification, we can only be sure about the final
          // whereabouts of the photon if it ends up beyond the event horizon ore
          // is ejected to infinity
          //stop_integration = 2; /* the photon misses the disk */
        }

      } else {
        r = rau;
        th = thau;
        phi = phiau;
        kr = krau;
        kth = kthau;
        h /= 2.0;
      }
    }
#endif
  } while (stop_integration == 0);

  if (stop_integration == 1 || stop_integration == 128) {
    xem[1] = r;
    //we also need the density at the point of the hit... for what????

    //to calculate the redshift, we need the photon momentum k (which is present with kr and kth, kt=-E=kt0, kphi=L=kphi0) the observer 4-vel,
    //which is (1,0,0,0), and the interpolated 4-vel of the disk. With this, we can calculate the gfactor.
    //do we need the metric or is it contained in the kvector already???
    //redshift(xem[1], const1, gfactor);
    gfactor = kt0
        / (kt0 * spi.u0 + kr * spi.u1 + kth * spi.u2 + kphi0 * spi.u3);

    /*Non Kerr PRD 90, 064002 (2014) Eq. 34*/
    cosem = carter * gfactor / sqrt(xem[1] * xem[1] + epsi3 / xem[1]);
  } else {
    xem[1] = 0.0;
    gfactor = 0.0;
    cosem = 0.0;
  }
  hit.cosem = cosem;
  hit.r = xem[1];
  hit.gfactor = gfactor;
  hit.hc = 0.0;
//  traces[0] = xem[1];
//  traces[1] = cosem;
//  // traces[2] = xem[3];
//  traces[3] = gfactor;
}
