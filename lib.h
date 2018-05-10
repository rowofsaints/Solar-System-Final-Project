// This is our header file which will have the classes and methods we will use

#include<iostream>
#include<vector>

using namespace std;


class Planet{
	
	public:
		double x,y,z;
		double vx,vy,vz;
		double fx,fy,fz;
		string name;
		double mass;
		double f; // What is f in threebody.py? its the total force. 
//  void set(double x,double y,double z,double vx,double vy,double vz,double mass);		
	// Planet Methods
        Planet();
		Planet(double x1,double y1,double z1, \
		       double vx1,double vy1,double vz1,\
		       double mass);

/*
		//changing data
		void setPosition(double x,double y,double z);
		void setVelocity(double vx, double vy, double vz);
		void setForces(double fx, double fy, double fz);
		void setName(string name);
		void setMass(double m);
		
		//accessing data
		string getName();
		double getMass();
		double * getPosition();
		double * getVelocity();
		double * getForces();


		//update methods
		void updateForces();
		void updateVelocities();
		void updatePosition();
*/
	//More functions to come...

}; //end Class Planet


double MoonTheta(double time);
double MoonPhi(double theta);
bool firstBurn(Planet mod, double time, double *a, double *dt);
bool secondBurn(Planet mod, double time, double a, double*dt);
void burnBaby(Planet *mod,double time);
void burnBabyBurn(Planet *mod, double a,double time);

void save_planets(vector<Planet> solarsystem, double time, ofstream out_file);

void adaptDT(double x1, double y1, double x2, double y2, double a, double *dt);
//What other classes do we need?
