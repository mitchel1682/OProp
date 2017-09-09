#include <iostream>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

mat dydt(double t, double dt, mat y){
	double msc = 1.0; //Mass of s/c/, Kilograms
	
	//establish forces
	//test 1.
	mat Ft(3,1); Ft << 5.0 << endr << 0 << endr << 0 << endr;
	
	//sum forces:
	mat F(3,1); F = Ft;
	
	//extract acceleration
	mat a(3,1);
	a = F/msc;
	
	//update velocity
	
	mat dydt(6,1);
	dydt << (a(0,0)*dt) + y(3,0) << endr
	<< (a(1,0)*dt) + y(4,0) << endr
	<< (a(2,0)*dt) + y(5,0) << endr
	<< a(0,0) << endr
	<< a(1,0) << endr
	<< a(2,0) << endr;
	
	return dydt;
}

mat RK4(double t, double h, mat y){
	int ksize = y.n_rows;
	mat k1(ksize,1); mat k2(ksize,1); mat k3(ksize,1); mat k4(ksize,1);
	
	k1 = dydt(t, h, y);
	k2 = dydt((t+(h/2)),h,(y+((h/2)*k1)));
	k3 = dydt((t+(h/2)),h,(y+((h/2)*k2)));
	k4 = dydt((t+h),h,(y+(h*k3)));
	
	mat y_2(ksize,1);
	y_2 = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
	return y_2;
}
// main () is where program execution begins
int main()
{
  //double re = 6.371*(1e+6);
  //double alt = 400*(1e+3);
  
  mat ro(3,1);	ro <<   0.0 << endr << 0.0 << endr << 0.0 << endr;
  mat vo(3,1); 	vo <<	0.0	<< endr << 0.0 << endr << 0.0 << endr;

  mat y(6,1); y << ro(0,0) << endr << ro(1,0) << endr << ro(2,0) << endr << vo(0,0) << endr << vo(1,0) << endr << vo(2,0) << endr;
  
  double t = 0; //start time, seconds
  double dt = .001; //time stepsize, seconds
  
  double t_end = 20; //timespan of simulation (whole number only) (limited run)
  int tot_i = t_end/dt; //total iterations (integer value) (limited run)
  
  //File Setup:
  FILE * statevector;
  statevector = fopen("statevector.txt","w");
  fprintf(statevector,"%f %f %f %f %f %f %f \n",t,y(0,0),y(1,0),y(2,0),y(3,0),y(4,0), y(5,0));
  fclose(statevector);
  double t_out = 0; //starting output counter
  double dt_out = 1; //output value every set of seconds. (whole number only)

  //mat y2 = dydt(t, h, y);
  //y2.print("y2:");
  for (int i = 0; i<tot_i; i++) {
	y = RK4(t,dt,y);
	t = t + dt;
	
	//recording data
	if (t_out<dt_out){
		t_out = t_out + dt;
	}
	else{
		t_out = 0;
		statevector = fopen("statevector.txt","a");
		fprintf(statevector,"%f %f %f %f %f %f %f \n",t,y(0,0),y(1,0),y(2,0),y(3,0),y(4,0), y(5,0));
		fclose(statevector);
	}
  
  }
  
  
  //y.print("y:");
  //ro.print("ro:");
  //cout << "Hello World " << re <<endl;      // prints Hello World

  return 0;
}