#include "fargo3d.h"
#include <stdlib.h>  // Para abs
#include <math.h>    // Para fabs

void InitPlanetarySystemFromPar(PlanetarySystem *sys) {

  // Parámetros del planeta
    real x0 = 60.0;  // Coordenada x inicial en AU
    real y0 = 10.0;    // Coordenada y inicial en AU
    real v0 = 1.0*sqrt(2.0 * G * MSTAR / sqrt(x0 * x0 + y0 * y0));  // Velocidad de escape
    real theta = M_PI ;  // Ángulo de lanzamiento en radianes

    real vx0 = v0 * cos(theta);
    real vy0 = v0 * sin(theta);

  

  // Inicialización de la estructura del sistema planetario
  sys->nb = 1; // Número de planetas
  sys->mass[0] = 0.01 * MSTAR; // Masa del planeta en masas solares
  sys->x[0] = x0; // Posición inicial en x
  sys->y[0] = y0; // Posición inicial en y
  sys->z[0] = 0.0; // Posición inicial en z
  sys->vx[0] = vx0; // vx0; // Velocidad inicial en x
  sys->vy[0] = vy0; // vy0; // Velocidad inicial en y
  sys->vz[0] = 0.0; // Velocidad inicial en z

  printf("Planet initialized at x=%g, y=%g, with velocity: vy=%g, vx=%g\n", sys->x[0], sys->y[0], sys->vy[0], sys->vx[0]);
}


void Init() {
  
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  int i,j,k;
  real r, omega;
  real soundspeed;
  
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
  real *rho  = Density->field_cpu;
  
#ifdef ADIABATIC
  real *e   = Energy->field_cpu;
#endif
#ifdef ISOTHERMAL
  real *cs   = Energy->field_cpu;
#endif

  i = j = k = 0;
  
  for (j=0; j<Ny+2*NGHY; j++) {
    for (i=0; i<Nx+2*NGHX; i++) {
      
      r = Ymed(j);
      omega = sqrt(G*MSTAR/r/r/r);
      
      rho[l] = SIGMA0*pow(r/R0,-SIGMASLOPE)*(1.0+NOISE*(drand48()-.5));
      soundspeed  = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;

#ifdef ISOTHERMAL
      cs[l] = soundspeed;
#endif
#ifdef ADIABATIC
      e[l] = pow(soundspeed,2)*rho[l]/(GAMMA-1.0);
#endif
      
      vphi[l] = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
			     (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
      vphi[l] -= OMEGAFRAME*r;
      vphi[l] *= (1.+ASPECTRATIO*NOISE*(drand48()-.5));
      
      vr[l]    = soundspeed*NOISE*(drand48()-.5);
    }
  } 


  //Sys = InitPlanetarySystem(PLANETCONFIG); // Leer el archivo de configuración del planeta
  InitPlanetarySystemFromPar(Sys);         // Inicializar el sistema planetario con parámetros del .par
  ListPlanets();

}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   Init();
}
