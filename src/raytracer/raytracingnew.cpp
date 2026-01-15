
long double interpolate(long double a, long double b, long double f) {
  return (1.0 - f) * a + f * b;
}

Real checkIntersect(long double x1, long double y1, long double x2,
    long double y2, long double x3, long double y3, long double x4,
    long double y4) {
  long double axl = std::min(x1, x2);
  long double axh = std::max(x1, x2);
  long double ayl = std::min(y1, y2);
  long double ayh = std::max(y1, y2);
  long double bxl = std::min(x3, x4);
  long double bxh = std::max(x3, x4);
  long double byl = std::min(y3, y4);
  long double byh = std::max(y3, y4);
  if (!(axl < bxh && axh > bxl && ayl < byh && ayh > byl)) {
    return NO_INTERSECT;
  }
  long double num = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
  long double denum = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
  long double t = num / denum;
  if (t >= 0 && t <= 1) {
    return t;
  }
  return NO_INTERSECT;
}


//finds any intersection. it is not guranteed (but very likely) that it is the closest one, i.e. the first hit. Sufficiently small stepsize should circuvent the problem,
//as well as retaking the step with a smaller stepsize.
//this might interpolate between a mesh cell right beyond the horizon, which might contain invalid data, and the first cell outside the horizon,
//for disk shapes which come very close to the horizon
//at the moment, this does not seem to be a problem, but keep this in mind, especially if and when doing a clean rewrite of this
int get_interpolated_sp(const long double x1, const long double y1,
    const long double x2, const long double y2, QuadTree* quadtree, SurfacePoint &out, int& index) {
  SurfaceElement* elem;
  Real result = quadtree->check_intersect(x1,y1,x2,y2,&elem);
  if(result != NO_INTERSECT) {
    index = elem->index;
    long double xi = (elem->sp0->x) + result * ((elem->sp1->x) - (elem->sp0->x));
    long double yi = (elem->sp0->y) + result * ((elem->sp1->y) - (elem->sp0->y));
    out.x = xi;
    out.y = yi;
    out.density = interpolate(elem->sp0->density, elem->sp1->density, result);
    //should be zero if either p is zero because linear interpolation of these velocities at that place is probably not physically
    out.u0 = interpolate(elem->sp0->u0, elem->sp1->u0, result);
    out.u1 = interpolate(elem->sp0->u1, elem->sp1->u1, result);
    out.u2 = interpolate(elem->sp0->u2, elem->sp1->u2, result);
    out.u3 = interpolate(elem->sp0->u3, elem->sp1->u3, result);
    return INTERSECT;
  }
  index = 0;
  return NO_INTERSECT;
}


void scalarProduct(Real met[4][4], Real* fvec0, Real* fvec1, Real& scal) {
  scal = 0.0;
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      scal += met[i][j]*fvec0[i]*fvec1[j];
    }
  }
}

void correct4VelNorm(Real met[4][4], Real norm, Real* fvel) {
  Real dif = norm+1;
  Real deltaut = dif/met[0][0];
  fvel[0] = sqrt(fvel[0]*fvel[0]-deltaut);
}

void raytrace(long double xobs, long double yobs, long double iobs,
    long double rin, long double disk_length_combined, RayHit &hit,
    int &stop_integration, SurfacePoint **diskdata, const size_t ddsize, QuadTree* tree, Real checkr) {
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
  rtol = 1.0e3;
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

  Real obsuarray[4] = {1.0, 0.0, 0.0, 0.0};
  Real obskarray[4] = {kt0, kr0, kth0, kphi0};
  Real obsenergy;
  scalarProduct(met,obsuarray,obskarray,obsenergy);

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
  long double prevh = -1.0;



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

      if (check == 1) {
        h /= 2.0;
#ifdef ITER_WARN
        if (iter > MAX_ITER - 10) {
          std::cout << "descale0" << std::endl;
        }
#endif
      } else if (check == -1)
        h *= 2.0;

    } while (check == 1);

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
    if (Delta < 1.0e-3) {
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
#ifdef ITER_WARN
    if (iter > MAX_ITER - 10) {
      std::cout << "Reaching max iter with..." << std::endl;
      std::cout << h << " " << iter << std::endl;
      std::cout << r << " " << th << " " << phi << std::endl;
      std::cout << rau << " " << thau << " " << phiau << std::endl;
    }
#endif
    if (iter > MAX_ITER) {
      stop_integration = 255;
      break;
    }

    //not at all close to disk so we don't need to perform the checks below
    if(r > checkr) continue;

    //check if the new position intersects the accretion disk
    //convert coordinates of current and previous position via a BL-cartesian conversion
    long double xcoord = std::sqrt(r * r + spin2) * sin(th);
    long double ycoord = r * cos(th);
    long double xcoordprev = std::sqrt(rau * rau + spin2) * sin(thau);
    long double ycoordprev = rau * cos(thau);
    int res = NO_INTERSECT;
    int index = 0;

    res = get_interpolated_sp(xcoordprev,ycoordprev,xcoord,ycoord,tree,spi,index);
    //deal with (possible) intersection
    if (res == INTERSECT) {
#ifdef ITER_WARN
      if (iter > MAX_ITER - 10) {
        std::cout << "int" << std::endl;
      }
#endif
      if (check2 != 1) {
        prevh = h;
      }
      check2 = 1; //don't adapt stepsize anymore, this is now done manually to reach certain tolerances
      if (fabs(th - thau) <= thtol) {
        count++;
      }
      if (count > 0) {
        if (std::abs(spi.y) > 0.0 && spi.x > 0.0) {
          // next step is redshift calculation with data from the intersection
          // point. we also need the interpolated 4-vel etc
          stop_integration = index;
          //break; /* the photon hits the disk */
        } else {
          check2 = 0;
          count = 0;
          h = prevh;
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
#ifdef ITER_WARN
        if (iter > MAX_ITER - 10) {
          std::cout << "descale1" << std::endl;
        }
#endif
      }
    }
  } while (stop_integration == 0);

  if (stop_integration >= 128 && stop_integration <= 131) {
    xem[1] = r;

    //to calculate the redshift, we need the photon momentum k (which is present with kr and kth, kt=-E=kt0, kphi=L=kphi0) the observer 4-vel,
    //which is (1,0,0,0), and the interpolated 4-vel of the disk. With this, we can calculate the gfactor.
    metric(r, th, met);


    //Real x = std::sqrt(r);
    //Real p_ut = (0.0 + CUBE(x))/std::sqrt(CUBE(x)*(2*0.0+CUBE(x)-3*x));
    //Real p_uph = 1/std::sqrt(CUBE(x)*(2*0.0+CUBE(x)-3*x));
    //long double uarray[4] = {p_ut,0.0,0.0,p_uph};
    //long double uarray[4] = {1,0,0,0};
    Real uarray[4] = {spi.u0, spi.u1, spi.u2, spi.u3};
    //to fix any inconsistencies introduced by linear interpolation or the change of coordinate chart (KS->BL)
    //or code differences between Athena++ and Blackray or simply numerical issues in the entire pipeline
    //the fix is done by recalculating (only) the time component of the 4-velocity so that the normalization is correct, i.e. far closer to -1.
    //the highest delta |spi.u0-fixedu0| is approximately 0.05 for an average disk.
    Real norm;
    scalarProduct(met, uarray, uarray, norm);
    correct4VelNorm(met, norm, uarray);
    Real newnorm;
    scalarProduct(met, uarray, uarray, newnorm);
#ifdef DEBUG_FVEL_NORM
    if(norm > -0.97 || norm < -1.03) {
      std::cout << "4-Vel norm deviates significantly, ignoring ray" << std::endl;
      std::cout << "Old norm: " << norm << std::endl;
      std::cout << "Fixed norm: " << newnorm << std::endl;
      std::cout << "Delta components: " << spi.u0-uarray[0] << " " << spi.u1-uarray[1] << " " << spi.u2-uarray[2] << " " << spi.u3-uarray[3] << " " << std::endl;
      stop_integration = 6;
    }
#endif
    //if we can't fix this mess, invalidate ray
    if(newnorm > -0.95 || newnorm < -1.05) {
      std::cout << "even fixed 4-vel norm deviates significantly, ignoring ray" << std::endl;
      stop_integration = 6;
      cosem = 0.0;
      gfactor = 1.0;
      xem[1] = r;
    }else{


      Real g_tt, g_pp, g_tp;
      g_tt = met[0][0];
      g_pp = met[3][3];
      g_tp = met[0][3];
      Real denom = (g_tt * g_pp - g_tp * g_tp);
      Real ktcalc = -(g_pp + b * g_tp) / denom;
      Real kphicalc = (g_tp + b * g_tt) / denom;
      Real karray[4] = {ktcalc, kr, kth, kphicalc};

  #ifdef DEBUG_FMOM_NORM
      Real knorm;
      scalarProduct(met, karray, karray, knorm);
      if(knorm > 0.03 || knorm < -0.03) {
        std::cout << "4-Momentum norm deviates significantly" << std::endl;
        std::cout << "Norm: " << knorm << std::endl;
        stop_integration = 6;
      }
  #endif

      Real emenergy;
      scalarProduct(met, uarray, karray, emenergy);
      gfactor = obsenergy/emenergy;

      //cosem stays artifical
      Real gfactorforcosem;
      redshift(xem[1],const1,gfactorforcosem);
      /*Non Kerr PRD 90, 064002 (2014) Eq. 34*/
      cosem = carter * gfactorforcosem / sqrt(xem[1] * xem[1] + epsi3 / xem[1]);
      //Workaround for redshift function giving nan...
      if(std::isnan(cosem)) {
        cosem = carter * gfactor / sqrt(xem[1] * xem[1] + epsi3 / xem[1]);
        if(cosem > 1.05){
          std::cout << "Cosem was nan, then fixed cosem was > 1.05, ignoring ray: " << cosem << std::endl;
          stop_integration = 6;
          xem[1] = r;
          gfactor = 1.0;
          cosem = 0.0;
        } else if(cosem > 1.0) {
          cosem = 1.0;
        }
      }
    }
  } else {
    xem[1] = r;
    gfactor = 1.0;
    cosem = 0.0;
  }
  if(gfactor < 0.0){
    std::cout << "gfactor is < 0.0, ignoring ray" << std::endl;
    stop_integration = 6;
    xem[1] = r;
    gfactor = 1.0;
    cosem = 0.0;
  }
#ifdef DEBUG_COSEM
  if(std::isnan(cosem)) {
    std::cout << "Cosem is nan, ignoring ray" << std::endl;
    stop_integration = 6;
    xem[1] = r;
    gfactor = 1.0;
    cosem = 0.0;
  }
  if(cosem > 1.05) {
    std::cout << "Cosem > 1.05 detected, ignoring ray: " << cosem << std::endl;
    stop_integration = 6;
    xem[1] = r;
    gfactor = 1.0;
    cosem = 0.0;
  }
  if(cosem > 1.0) {
    std::cout << "1.05 >= Cosem > 1.0 detected, clamping to 1.0: " << cosem << std::endl;
    cosem = 1.0;
  }
#endif
  hit.cosem = cosem;
  hit.r = xem[1];
  hit.gfactor = gfactor;
  hit.hc = 0.0;

}
