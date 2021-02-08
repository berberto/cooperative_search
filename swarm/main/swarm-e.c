#define MAIN_PROGRAM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl_sf_bessel.h>
#include "random.h"
#include <mpi.h>

/* Parameters defined in the MATLAB script generating the control */
#include "parameters.h"

/* Parameters of the simulation */
#define N 1000      /* # particles (trajectories to simulate) */
#define Ntp 1000       /* # trajectories to print */
#define ndt 10      /* # time steps for the Langevin equation every dt */
#define T 4.         /* simulation time */

double alpha, gammat, lambda;

double uex (double X, double Y){
    double r = sqrt(X*X + Y*Y);
    /* return -2.*D*X/r; */
    return -2.*X/r*lambda*D*gammat/gamma*gsl_sf_bessel_K1(lambda*r)/gsl_sf_bessel_K0(lambda*r) ;
}

double uey (double X, double Y){
    double r = sqrt(X*X + Y*Y);
    /* return -2.*D*Y/r; */
    return -2.*Y/r*lambda*D*gammat/gamma*gsl_sf_bessel_K1(lambda*r)/gsl_sf_bessel_K0(lambda*r);
}

int main (int argc, char *argv[]){

  	/* Initialize the MPI environment. The two arguments to MPI Init are not
	 * currently used by MPI implementations, but are there in case future
	 * implementations might need the arguments. */
	MPI_Init(NULL, NULL);

	/* Get the number of processes */
	int Nruns;
	MPI_Comm_size(MPI_COMM_WORLD, &Nruns);

	/* Get the rank of the process */
	int run;
	MPI_Comm_rank(MPI_COMM_WORLD, &run);


    if (argc < 3){
        if (run == 0){
            printf("Enter\n\
            1) value of alpha\n\
            2) output directory (within \"output/\")\n");
        }
        MPI_Finalize();
		exit(EXIT_FAILURE);
    }

    int i, idx[N], idxn, absorbed[N];
    int ts;
    
    double ddt = dt/ndt;
    double r[3], dB[2], tempx, tempy;
    double px[N], py[N];
    double cost_exp[N], cost[N], contc[N], timec[N];
    double cx, cy;

    char filename[64], buff[64], cmd[64], outdir[64];
    FILE *outcost, *outcont, *outtime;
    FILE *traj[Ntp];
    
    alpha = atof(argv[1]);
    gammat = gamma/(1.-2*D*alpha*gamma);
    lambda = sqrt(q/(2*gammat))/D;

    sprintf(outdir,"%s", argv[2]);

    if (run == 0){
        printf("Results in output/%s/trajs_%.2lf\n\n", outdir, alpha);
        printf("chem coeff: 2Dgt/g = %lf\ninv length: lambda = %lf\n", 2*D*gammat/gamma, lambda);
        sprintf(cmd, "mkdir -p output/%s/trajs_%.2lf", outdir, alpha);
        system(cmd);
    }
    
    sprintf(filename, "output/%s/trajs_%.2lf/cost_tot_%03d.dat", outdir, alpha, run);
    outcost = fopen(filename, "w");

    sprintf(filename, "output/%s/trajs_%.2lf/cost_cont_%03d.dat", outdir, alpha, run);
    outcont = fopen(filename, "w");

    sprintf(filename, "output/%s/trajs_%.2lf/cost_time_%03d.dat", outdir, alpha, run);
    outtime = fopen(filename, "w");


    if (run == 0){
        for (i=0; i<Ntp; i++){
            sprintf(filename, "output/%s/trajs_%.2lf/traj_%03d.dat", outdir, alpha, i);
            traj[i] = fopen(filename, "w");
        }
    }

    /* Initialization of the randomizer */
    int seed;
    seed = rlxd_seed();
	printf("%d --> start.\t%d\n", run, seed);
	rlxd_init(1,seed);

    /* Initialize particles and costs */
    i = 0;
    for(i=0;i<N;i++){
        px[i] = x0;
        py[i] = y0;
        timec[i] = 0.;
        contc[i] = 0.;
        absorbed[i] = 0;
    }

    /* Simulating trajectories */
    ts = 0;
    while (ts*ddt < T){

        /* evolve trajectories */
        for (i=0; i<N; i++){

            /* print trajectories on separate files */
            if (run==0 && i<Ntp) fprintf(traj[i],"%lf\t%lf\t%lf\n", ts*ddt, px[i], py[i]);
            
            /* if particle has been absorbed, do nothing to it and go to the next */
            if (absorbed[i]) continue;

            /* Move the particle */
            cx = uex(px[i],py[i]);
            cy = uey(px[i],py[i]);

            gauss_dble(dB, 2);
            px[i] = px[i] + cx*ddt + sqrt(2*D*ddt)*dB[0];
            py[i] = py[i] + cy*ddt + sqrt(2*D*ddt)*dB[1];

            /* Add running costs */
            timec[i] += q*ddt;
            contc[i] += gamma*.5*(cx*cx + cy*cy)*ddt;

            /* if the particle has reached the target, mark it as absorbed */
            if (px[i]*px[i] + py[i]*py[i] <= R*R) absorbed[i] = 1;
        }
        ++ts;
    }
    if (run == 0) for (i=0; i<Ntp; i++) fclose(traj[i]);

    for (i=0; i<N; i++){
        cost[i] = contc[i] + timec[i];
        fprintf(outcost, "%lf\n", cost[i]);
        fprintf(outcont, "%lf\n", contc[i]);
        fprintf(outtime, "%lf\n", timec[i]);
    }
    fclose(outcost);
    fclose(outcont);
    fclose(outtime);

    printf("%d --> end\n",run ); fflush(stdout);
    
	/* Finalize the MPI environment. No more MPI calls can be made after this */
	MPI_Finalize();
	
	exit(EXIT_SUCCESS);
}
