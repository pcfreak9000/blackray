#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>

#include "problem.hpp"

#include "raytracingnew.hpp"
#include "metric.hpp"
#include "find_isco.hpp"
#include "quadtree.hpp"
#include "environment.hpp"

//oof
Real chi,eps3,Mdot,eta;
Real johannsen_isco;
#include "johannsen.h"


#define RayNum 100000
#define imax 400
static constexpr Real dobs = 2000; /* distance of the observer */

struct Record {
	Real emisDelta;
	Real incDelta;
	Real rCylindCoord;
	Real rSurf;
	Real thetaSurf;
	Real rBeforeHit;
	Real thetaBeforeHit;
	Real rAfterHit;
};

class del_ic: public InitialCondition {
public:
	Real delta;
	Real stepDelta;
};

static std::vector<Record> records;
static size_t arr_index = 0;
static Real height=5.0;

static Real rDisk[imax], rDiskNew[imax], RayCount[imax];
static Real rHit_prev;
static int j;

void setupProblem(int argc, char *argv[], Env *env, size_t &ray_count_total) {
	/* ----- SETUP ENVIRONMENT ----- */
	Real maxx;
	Real maxy;
	const Real spin2 = SQR(spin);
	chi = spin;
	eps3 = epsi3;
	Mdot=0.0;
//	std::unique_ptr<QuadTree> tree = readFileToTree(diskdatafile, maxx, maxy);
//  if (!tree) return 1;
	Real maxr_xdir = std::sqrt(SQR(maxx+10.0) - spin2) + 10.0;
	Real maxr_ydir = std::sqrt(SQR(maxy + 10.0)) + 10.0;
	Real checkr = maxr_ydir;
	if (maxr_xdir > maxr_ydir)
		checkr = maxr_xdir;
	Real isco;

	find_isco(15.0, isco); /* Depends upon the properties of BH */
//  GRMHDDisk disk(tree.get(), checkr);
	//env.addEntity(&disk);
//  env->addEntity(std::make_unique<GRMHDDisk>(std::move(tree), checkr));
	env->addEntity(std::make_unique<ThinDisk>(isco, 1050));
	//env->addEntity(std::make_unique<PlungingRegion>(isco));
	johannsen_isco = kerr_rms(spin);
    eta = 1.0 - specific_energy(johannsen_isco);
	for (int i = 0; i < imax - 1; i++) {
		rDisk[i] = std::pow((double) i / (imax - 2), 3) * (1000. - johannsen_isco)
				+ johannsen_isco;
		RayCount[i] = 0;
	}
	rDisk[imax - 1] = 1050;

	ray_count_total = RayNum;
	records.reserve(RayNum);
}

std::unique_ptr<InitialCondition> initialcondition(const size_t &ray_index) {
	std::unique_ptr<del_ic> cond = std::make_unique<del_ic>();
	Real met[4][4];

	Real Emis_delta[2];
	Emis_delta[1] = M_PI - 1e-8;
	Emis_delta[0] = 1E-8;
	Real delta = Emis_delta[0];
	Real stepDelta = (Emis_delta[1] - Emis_delta[0]) / ((double) RayNum - 1.);
	//delta -= stepDelta;
	delta = Emis_delta[0] + stepDelta * (ray_index - 1);
	const Real r0 = height;  //height
	const Real th0 = 1E-8;
	const Real phi0 = 0.0;
	metric(r0, th0, met);
	Real g_tt = met[0][0];
	Real g_rr = met[1][1];
	Real g_thth = met[2][2];
	Real g_pp = met[3][3];
	Real g_tp = met[0][3];

	const Real bkr0 = std::cos(delta) / std::sqrt(g_rr);
	const Real bkth0 = std::sin(delta) / std::sqrt(g_thth);
	const Real kphi0 = 0.0;
	const Real E = sqrt(
			g_tp * g_tp * kphi0 * kphi0
					- g_tt
							* (g_rr * SQR(bkr0) + g_thth * SQR(bkth0)
									+ g_pp * kphi0 * kphi0));
	const Real kr0 = bkr0 / E;
	const Real kth0 = bkth0 / E;

	const Real kt0 = -(g_tp * kphi0 + E) / g_tt;
	const Real b = -(g_pp * kphi0 + g_tp * kt0) / (g_tt * kt0 + g_tp * kphi0);

	cond->pvec0.r = r0;
	cond->pvec0.th = th0;
	cond->pvec0.phi = phi0;

	cond->pvec0.kr = kr0;
	cond->pvec0.kth = kth0;
	cond->b = b;
	cond->stepDelta = stepDelta;
	cond->delta = delta;
	cond->dobs = dobs;
	return cond;
}

void notifyDone(const int &stopping_condition, const InitialCondition *const ic,
		const RayHit &hit, const size_t &ray_index) {
	Record rec;
	if (hit.entity != nullptr) {
		//std::cout << stopping_condition << std::endl;
		const del_ic *const dic = (del_ic*) ic;
		Real rHit = (hit.pvec.r + hit.pvecau.r) * 0.5;
		rec.emisDelta = dic->delta;
		rec.incDelta = cal_delta_inc(dic->delta, hit.pvec.th, height, rHit,hit.pvec.vars,hit.pvecau.vars);
		rec.rCylindCoord = std::fabs(rHit * std::sin(hit.pvec.th));
		rec.rSurf = std::fabs(rHit);
		rec.thetaSurf = hit.pvec.th;
		// phflag[arr_index] = 1

		rec.rBeforeHit = std::fabs(hit.pvecau.r);
		rec.thetaBeforeHit = hit.pvecau.th;
		rec.rAfterHit = std::fabs(hit.pvec.r);

#pragma omp ordered
		{
			arr_index++;
			//std::cout << rHit << std::endl;
			//std::cout << rDisk[j] << " " << rDisk[j+1] << std::endl;
			while(rHit >= rDisk[j+1]) j++;
			//std::cout << j << std::endl;
			if (rDisk[j] <= rHit && rHit < rDisk[j + 1]) {
				RayCount[j] = std::sin(dic->delta) * dic->stepDelta
						/ std::fabs(rHit - rHit_prev);
				rDiskNew[j] = rHit;
				j = j + 1;
			}

			rHit_prev = rHit;

			records.push_back(rec);
		}
	}
}

void finishProblem(const char *tempdir, const char *outtxt) {
	std::cout << "Writing data..." << std::endl;
	/* ----- file stuff ----- */
	Real III[imax-1], ED[imax-1], HD[imax-1];
	for(int i=0, j=0; i<imax-1; i++) {

	      bool nonono = false;
	      j=0;
	        while(!(records[j].rCylindCoord < rDisk[i] && rDisk[i] < records[j+1].rCylindCoord)){
	          j++;
	          if(j>RayNum-1){
	            nonono = true;
	            j--;
	            break;
	          }
	        }
	        Real A;
	        Real LorentzF;
	        if(!nonono){
	        //assert(rCylindCoord[j] < rDisk[i] && rDisk[i] < rCylindCoord[j+1]);

	        Real frac = (records[j+1].rCylindCoord - rDisk[i]) / (records[j+1].rCylindCoord - records[j].rCylindCoord);
	        Real thetaDisk = 0.5 * (1-frac) * (records[j].thetaSurf + records[j].thetaBeforeHit) + 0.5 * frac * (records[j+1].thetaSurf + records[j+1].thetaBeforeHit);

	        Real thetaAfterHitj = (1-frac) * records[j].thetaSurf + frac * records[j+1].thetaSurf;
	        Real thetaBeforeHitj = (1-frac)*records[j].thetaBeforeHit + frac*records[j+1].thetaBeforeHit;
	        Real rAfterHitj = (1-frac)*records[j].rAfterHit + frac*records[j+1].rAfterHit;
	        Real rBeforeHitj = (1-frac)*records[j].rBeforeHit + frac*records[j+1].rBeforeHit;

	        Real denom2 = rAfterHitj * std::sin(thetaAfterHitj) - rBeforeHitj * std::sin(thetaBeforeHitj);
	        Real dthdrhau;
	        if (denom2 != 0)
	            dthdrhau = (thetaAfterHitj - thetaBeforeHitj) / denom2;


	        ED[i] = (1-frac) * records[j].emisDelta + frac * records[j+1].emisDelta;
	        HD[i] = (1-frac) * records[j].incDelta + frac * records[j+1].incDelta;

	        if (ED[i] > M_PI/2.0) HD[i] *= (-1.0);
	        Real g[4][4];
	        Real dg_dr[4][4];
	        metric(rDiskNew[i], thetaDisk, g);
	        metric_rderivatives(rDiskNew[i], thetaDisk, dg_dr);
	        if (Mdot == 0) A = 2 * M_PI * std::sqrt(g[1][1] * g[3][3]);
	        else A = 2 * M_PI * sqrt(g[1][1] + g[2][2] * dthdrhau * dthdrhau * sin(thetaDisk) * sin(thetaDisk)) * sqrt(g[3][3] / sin(thetaDisk) / sin(thetaDisk));

	        //A *= (rDisk[i+1] - rDisk[i]); A => dA/dr
	        Real g_src[4][4];
	        metric(height, 0.0, g_src);
	        Real g_tt_source = g_src[0][0];
	        Real Omega = (-dg_dr[0][3] + std::sqrt(dg_dr[0][3]*dg_dr[0][3] - dg_dr[0][0]*dg_dr[3][3])) / dg_dr[3][3];
	        LorentzF = std::pow(std::pow(Omega + g[0][3] / g[3][3], 2) * g[3][3] * g[3][3] / (g[0][0]*g[3][3] - g[0][3]*g[0][3]) + 1.0, -0.5);
	        Real redshift = std::sqrt(g_tt_source / (g[0][0] + 2 * g[0][3] * Omega + g[3][3] * Omega * Omega));
	        }
	        III[i] =nonono?0.0: RayCount[i] / (A * LorentzF);
	        // III[i] *= pow(redshift, gamma);

	    }

		FILE *out = fopen("data_lp.dat","w");
	    for (int i = 0; i<imax-1; i++) {
	        std::cout << III[i] << std::endl;
	        std::cout << rDiskNew[i] << std::endl;

	        if (std::isnan(III[i])) III[i] = 0.;
	        // fprintf(out, "%f %f %f %e %f %f\n", a, height, rDisk[i], III[i], ED[i], HD[i]);
	        //fprintf(out, "%f %e %f %f\n", rDisk[i], III[i], ED[i], HD[i]);
	        fprintf(out, "%f %e\n", (double)rDiskNew[i], (double)III[i]);

	    }
	    fprintf(out, "\n");
	    fclose(out);
//  std::cout << "Integrated " << (records.size()) << " rays of which "
//      << photon_index << " hit the disk" << std::endl;
    std::cout << "Finishing..." << std::endl;
}

