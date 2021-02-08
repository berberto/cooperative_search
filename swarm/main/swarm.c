#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl_sf_bessel.h>
#include "random.h"

/* Parameters defined in the MATLAB script generating the control */
#include "parameters.h"

double gammat = gamma;
double lambda;

const double dx = 1.0*(xmax-xmin)/(Lx-1);
const double dy = 1.0*(ymax-ymin)/(Ly-1);

/* Parameters of the simulation */
#define N 1000     /* # particles (trajectories to simulate) */
#define Ntp 1000       /* # trajectories to print */
#define ndt 10      /* # time steps for the Langevin equation every dt */
#define T 1.         /* simulation time */


double x[Lx], y[Ly];            /* lattice positions */
double ux[Lx*Ly], uy[Lx*Ly];    /* drift field */
double px[N], py[N];            /* positions of particles */


/* Find lattice position of i-th particle */
int findindex (double X, double Y){
    int ii, jj;

    ii = (int)rint((X-xmin)/dx);
    jj = (int)rint((Y-ymin)/dy);

    return jj*Lx + ii;
}

/* Initial density */
double rho0 (double x, double y){

    /* (R^2-(x-x0)^2)*Theta[R^2-(x-x0)^2] */
    /* double rsq = 1. - ((x-x0)*(x-x0) + (y-y0)*(y-y0))/(s0*s0);    
    if (rsq > 0.) return rsq;
    else return 0.; */
    
    /* exp[-((x-x0)/s0)^4] */
    return exp(-(pow((x - x0),2.)+pow((y - y0),2.))/pow(s0,4.))/(2.78416*s0*s0);
}


int main (int argc, char *argv[]){

    if (argc < 5){
        printf("Enter\n\
            1) value of alpha\n\
            2) drift file \"..._con.dat\"\n\
            3) output directory (within \"output/\")\n\
            4) free (0) or controlled (1)\n");
		exit(EXIT_FAILURE);
    } else if (atoi(argv[4])!=0 && atoi(argv[4])!=1){
        printf("Fourth argument must be 0 or 1\n");
    }

    double alpha = atof(argv[1]);
    int control = atoi(argv[4]); /* Controlled dynamics or not */;

    int i, idx[N], idxn, absorbed[N];
    int ts;
    
    double ddt = dt/ndt;
    double r[3], dB[2*N], tempx, tempy;
    double cx, cy;

    char filename[64], buff[64], cmd[64];
    FILE *traj[Ntp];

    FILE *drift;
    drift   = fopen(argv[2], "r");

    if (alpha!=0){
        control=2;
        gammat = gamma/(1.-2*D*alpha*gamma);
    }
    lambda = sqrt(q/(2*gammat))/D;

    printf("Results in output/%s/trajs_%d\n\n", argv[3], control);
    printf("chem coeff: 2Dgt/g = %lf\ninv length: lambda = %lf\n", 2*D*gammat/gamma, lambda);

    
    sprintf(cmd, "mkdir -p output/%s/trajs_%d", argv[3], control);
    system(cmd);

    for (i=0; i<Ntp; i++){
        sprintf(filename, "output/%s/trajs_%d/traj_%03d.dat", argv[3], control, i);
        traj[i] = fopen(filename, "w");
    }

    /* Initialization of the randomizer */
    int seed;
    seed = rlxd_seed();
	rlxd_init(1,seed);

    /* Initialize particles and costs */
    i = 0;
    while(i<N){
        /* Initialize particles' positions */

        /*  by rejection with distr rho0  */
        ranlxd(r,3);
        r[0] = (2.*r[0] - 1.)*4*s0 + x0;
        r[1] = (2.*r[1] - 1.)*4*s0 + y0;
        if (r[2] < rho0(r[0],r[1])){
            px[i] = r[0];
            py[i] = r[1];
            absorbed[i] = 0;
            idx[i] = findindex(px[i], py[i]);
            i++;
        }

        /* from one point */
        /** px[i] = x0;
        py[i] = y0;
        timec[i] = 0.;
        contc[i] = 0.;
        absorbed[i] = 0;
        idx[i] = findindex(px[i], py[i]);
        i++;
        **/
    }

    /* Save lattice points */
    for (i=0; i<Lx; i++) x[i] = xmin + dx*i;
    for (i=0; i<Ly; i++) y[i] = ymin + dy*i;

    /* Simulating trajectories */
    ts = 0;
    while (ts*ddt < T){

        /* If time step is integer multiple of ndt */
        /* read the drift file and fill ux and uy vectors */
        if (ts % ndt == 0){
            for (i=0; i<Lx*Ly; i++){
                fgets(buff, 32, drift);
                sscanf(buff, "%lf %lf", &ux[i], &uy[i]);
                if (!control){
                    if (ux[i]!=1000. || uy[i]!=1000.){
                        ux[i] = 0.;
                        uy[i] = 0.;
                    }
                }
            }
            fgets(buff, 32, drift);
        }
        
        /* evolve trajectories */
        for (i=0; i<N; i++){

            /* print trajectories on separate files */
            if (i<Ntp) fprintf(traj[i],"%lf\t%lf\t%lf\n", ts*ddt, px[i], py[i]);
            
            /* if particle has been absorbed, do nothing to it and go to the next */
            if (absorbed[i]) continue;

            /* Move the particle */
            /* ..using lattice data */
            cx = ux[idx[i]];
            cy = uy[idx[i]];
            /* ..using exact solution */
            /* cx = uex(px[i],py[i]);
            cy = uey(px[i],py[i]); */

            gauss_dble(dB, 2);
            tempx = px[i] + cx*ddt + sqrt(2*D*ddt)*dB[0];
            tempy = py[i] + cy*ddt + sqrt(2*D*ddt)*dB[1];

            /* if the particle has reached the target, mark it as absorbed */
            if (tempx*tempx + tempy*tempy <= R*R){
                px[i] = tempx;
                py[i] = tempy;
                absorbed[i] = 1;
                continue;
            }

            /* index of the new position */
            idxn = findindex(tempx,tempy);

            /* if particle would end up outside the domain, try until it stays -- reflecting boundaries */
            /* ...only for control from lattice data */
            while((ux[idxn]==1000. && uy[idxn]==1000.) || tempx > xmax || tempx < xmin || tempy > ymax || tempy < ymin){
                gauss_dble(r,2);
                tempx = px[i] + ux[idx[i]]*ddt + sqrt(2*D*ddt)*r[0];
                tempy = py[i] + uy[idx[i]]*ddt + sqrt(2*D*ddt)*r[1];
                idxn = findindex(tempx,tempy);
            }

            px[i] = tempx;
            py[i] = tempy;
            idx[i] = idxn;

        }
        ++ts;
    }
    for (i=0; i<Ntp; i++) fclose(traj[i]);

}
