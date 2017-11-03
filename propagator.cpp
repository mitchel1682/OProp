#include <iostream>
#include <math.h>
#include "nrlmsise/nrlmsise-00.h"
#include <armadillo>

using namespace std;
using namespace arma;

//important globals:
double msc; //spacecraft mass, in Kg
double a_max; //maximum drag area, m^2.

//Atmospheric Drag Globals:
struct nrlmsise_output output[1];
struct nrlmsise_input input[1];
struct nrlmsise_flags flags;

void update_time(double t){
	
	double t_alt = t;
	double d_years = floor(t_alt/(86400.0*365.0));
	t_alt = t_alt - (d_years*86400.0*365.0);
	
	double d_days = floor(t_alt/86400.0);
	t_alt = t_alt - d_days*86400.0;
	
	//update time vector for atmospheric model:
	input[0].doy= d_days + input[0].doy; 					 /* day of year */
	input[0].year= d_years + input[0].year; 				 /* year, currently ignored */
	input[0].sec= t_alt + input[0].sec; 				 /* seconds in day (UT) */
	//input[0].alt= 400;  				 /* altitude in kilometers */
	//input[0].g_lat= 0;  				 /* geodetic latitude */
	//input[0].g_long= 0;				     /* geodetic longitude */
	input[0].lst=input[0].sec/3600 + input[0].g_long/15;	/* local apparent solar time (hours), see note below */
	input[0].f107A=70;
	input[0].f107=70;
	input[0].ap=4;
}

mat dydt(double t, mat y){
	
	//establish accelerations
	mat r_ECI(3,1); r_ECI<< y(0,0)<< endr << y(1,0) << endr << y(2,0) << endr;
	mat v_ECI(3,1); v_ECI<< y(3,0)<< endr << y(4,0) << endr << y(5,0) << endr;
	
	//keplerian gravity:
	double mu_e = 3.98574405096E+14; //earth's standard gravitational parameter, SI Units.
	double r_scalar = norm(r_ECI,2);
	double a_scalar = (-1.0*mu_e)/(pow(r_scalar,3));
	mat a_grav(3,1); a_grav = a_scalar*(r_ECI);
	
	//Equatorial Bulge (J2) Effects:
	double re = 6378137; //radius of the earth - meters
	double J2 = 1.0826E-3;
	double eq_x = ((mu_e*r_ECI(0,0))/(pow(r_scalar,3)))*J2*(3.0/2.0)*(pow((re/r_scalar),2))*(5.0*(pow(r_ECI(2,0),2))/(pow(r_scalar,2))-1.0);
	double eq_y = ((mu_e*r_ECI(1,0))/(pow(r_scalar,3)))*J2*(3.0/2.0)*(pow((re/r_scalar),2))*(5.0*(pow(r_ECI(2,0),2))/(pow(r_scalar,2))-1.0);
	double eq_z = ((-mu_e*r_ECI(2,0))/(pow(r_scalar,3)))*J2*(3.0/2.0)*(pow((re/r_scalar),3))*(3.0 - 5.0*(pow(r_ECI(2,0),2))/(pow(r_scalar,2)));
	mat a_J2(3,1); a_J2 << eq_x << endr << eq_y << endr << eq_z << endr;
	
	//Atmospheric Drag Effects:
	input[0].alt= (norm(r_ECI,2)- 6378137)/1000;				 /* altitude in kilometers */
	input[0].g_lat= 0;  				 /* geodetic latitude */
	input[0].g_long= 0;				     /* geodetic longitude */

	gtd7(&input[0], &flags, &output[0]);	
	double rho = output[0].d[5]*1000.0;
	
	
	mat w_e(3,1); w_e << 0.0 << endr << 0.0 << endr << 7.292115E-5 << endr;
	mat skew_w(3,3); skew_w << 0.0 << w_e(2,0) << -1.0*w_e(1,0) << endr << -1.0*w_e(2,0) << 0.0 << w_e(0,0) << endr << w_e(1,0) << -1.0*w_e(0,0) << 0.0 << endr;
	
	mat r_dot_gas(3,1); r_dot_gas = v_ECI + skew_w*r_ECI;
	double r_dot_sc = norm(r_dot_gas,2);
	double Cd = 2.00; //from Wertz.
	
	mat a_drag(3,1); a_drag = (-1.0/2.0)*Cd*rho*(a_max/msc)*(r_dot_sc)*(r_dot_gas);
	
	//update acceleration
	mat a(3,1); 
	a = a_grav + a_J2 + a_drag;  

	//update velocity
	mat dydt(6,1);
	dydt <<y(3,0) << endr
	<< y(4,0) << endr
	<< y(5,0) << endr
	<< a(0,0) << endr
	<< a(1,0) << endr
	<< a(2,0) << endr;
	
	return dydt;
}

mat RK4(double t, double h, mat y){
	int ksize = y.n_rows;
	mat k1(ksize,1); mat k2(ksize,1); mat k3(ksize,1); mat k4(ksize,1);
	
	k1 = dydt(t, y);
	k2 = dydt((t+(h/2)),(y+((h/2)*k1)));
	k3 = dydt((t+(h/2)),(y+((h/2)*k2)));
	k4 = dydt((t+h),(y+(h*k3)));
	
	mat y_2(ksize,1);
	y_2 = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
	return y_2;
}
// main () is where program execution begins
int main()
{
	mat ro(3,1);
	mat vo(3,1);
	double year_i;
	double doy_i;
	double second_i;
	
	double dt;
	double t_end;
	double dt_out;
	
	FILE * IC;
	printf("Maximum Drag Orbital Propagator - M. McDonald 10/31/17 \n");
	IC = fopen("IC.txt","r");
	fscanf(IC,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&a_max,&msc,&ro(0,0),&ro(1,0),&ro(2,0),&vo(0,0),&vo(1,0),&vo(2,0),&year_i,&doy_i,&second_i,&dt,&t_end,&dt_out);
	fclose(IC);
	printf("Sim Parameters: \n");
	printf("total time: %f days \n",(t_end/86400.0));
	printf("timestep %f seconds \n",dt);
	printf("Output Decimation: %f seconds \n",dt_out);
	printf("Initial Conditions: \n");
	printf("S/C Mass: %f Kg \n ",msc);
	printf("Max Area: %f m^2 \n",a_max);
	printf("Deployment Date: %f %f %f\n",year_i,doy_i,second_i);
	printf("r_ECI_i: %f %f %f \n",ro(0,0),ro(1,0),ro(2,0));
	printf("v_ECI_i: %f %f %f \n",vo(0,0),vo(1,0),vo(2,0));
	
	mat y(6,1); y << ro(0,0) << endr << ro(1,0) << endr << ro(2,0) << endr << vo(0,0) << endr << vo(1,0) << endr << vo(2,0) << endr;
  
	double t = 0; //start time, seconds
	long int tot_i = t_end/dt; //total iterations (integer value) (limited run)
  
	//File Setup:
	FILE * statevector;
	statevector = fopen("statevector.txt","w");
	fprintf(statevector,"%f %f %f %f %f %f %f \n",t,y(0,0),y(1,0),y(2,0),y(3,0),y(4,0), y(5,0));
	fclose(statevector);
	double t_out = 0; //starting output counter

  	//update time vector for atmospheric model:
	input[0].doy= doy_i; 					 /* day of year */
	input[0].year= year_i; 				 /* year, currently ignored */
	input[0].sec= second_i; 				 /* seconds in day (UT) */
	//input[0].alt= 400;  				 /* altitude in kilometers */
	//input[0].g_lat= 0;  				 /* geodetic latitude */
	//input[0].g_long= 0;				     /* geodetic longitude */
	input[0].lst=input[0].sec/3600 + input[0].g_long/15;	/* local apparent solar time (hours), see note below */
	input[0].f107A=70;
	input[0].f107=70;
	input[0].ap=4;
  for (int i = 0; i<tot_i; i++) {
	y = RK4(t,dt,y);
	t = t + dt;
	update_time(t);
	//recording data
	if (t_out<dt_out){
		t_out = t_out + dt;
	}
	else{
		t_out = 0;
		statevector = fopen("statevector.txt","a");
		fprintf(statevector,"%f %f %f %f %f %f %f \n",t,y(0,0),y(1,0),y(2,0),y(3,0),y(4,0), y(5,0));
		fclose(statevector);
		
		double alt = (sqrt(y(0,0)*y(0,0) + y(1,0)*y(1,0) + y(2,0)*y(2,0)) - 6378137)/1000;
		int days = floor(t/86400);
		printf("days completed: %i, altitude: %f Km \n",days,alt);
		if (alt <120){
			break;
		}
	}
  
  }

  return 0;
}