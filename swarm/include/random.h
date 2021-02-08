
#ifndef RANDOM_H
#define RANDOM_H


#ifndef RANLXS_C
extern void ranlxs(float r[],int n);
extern void rlxs_init(int level,int seed);
extern int rlxs_size(void);
extern void rlxs_get(int state[]);
extern void rlxs_reset(int state[]);
#endif

#ifndef RANLXD_C
extern void ranlxd(double r[],int n);
extern int rlxd_seed();
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);
#endif

#ifndef GAUSS_C
extern void gauss(float r[],int n);
extern void gauss_dble(double r[],int n);
extern double gaussdistr(double x);
#endif

#ifndef GAUSSALT_C
extern void wgauss_dble(double r[], int n, double sigma);
extern double wgaussdistr(double x, double sigma);
extern void stdgauss_dble(double r[],int n);
extern double stdgaussdistr(double x);
#endif

#ifndef EXPO_C
extern void expo(float r[], int n);
extern void expo_dble(double r[], int n);
extern double expodistr(double x);
extern double wexpo(float r[], int n, float alpha);
extern double wexpo_dble(double r[], int n, double alpha);
extern double wexpodistr(double x, double alpha);
#endif

#ifndef SYMEXPO_C
extern void symexpo(float r[], int n);
extern void symexpo_dble(double r[], int n);
extern double symexpodistr(double x);
extern double symwexpo(float r[], int n, float alpha);
extern double symwexpo_dble(double r[], int n, double alpha);
extern double symwexpodistr(double x, double alpha);
#endif

#ifndef ROOTEXPO_C
extern void rootexpo(float r[], int n);
extern void rootexpo_dble(double r[], int n);
extern double rootexpodistr(double x);
#endif

#ifndef SYMROOTEXP_C
extern void symrootexpo(float r[], int n);
extern void symrootexpo_dble(double r[], int n);
extern double symrootexpodistr(double x);
#endif

#endif

