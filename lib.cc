//This is the library folder, containing the definitions of our functions and classes from the header file

#include "lib.h"
#include<iostream>
#include<vector>
#include<fstream>
#include "math.h"
using namespace std;

void save_planets(vector<Planet> solarsystem, double time,ofstream out_file){
  out_file << time << endl;
  for(int i=0;i<solarsystem.size();i++){
    out_file << solarsystem[i].x << " " << solarsystem[i].y << " "<< solarsystem[i].z <<endl;
  }

}




Planet ::Planet(){x=0;y=0;z=0;vx=0;vy=0;vz=0;mass=0;};

Planet ::Planet(double x1,double y1,double z1, \
		double vx1,double vy1,double vz1,\
		double mass)
	{
	x = x1;
	y = y1;
	z = z1;
	vx = vx1;
	vy = vy1;
	vz = vz1;
//	 setPosition(x,y,z);

//	 setVelocity(vx,vy,vz);

//	 setForces(fx,fy,fz);

//	 setName(name);

//	 setMass(mass); 	

	}; //Planet initialization
/*
 

void Planet :: setPosition(double X,double Y,double Z)
		{x[0] = X; x[1] = Y; x[2] = Z;};

void Planet :: setVelocity(double vx, double vy, double vz)
		{v[0] = vx; v[1] = vy; v[2] = vz;};

void Planet :: setForces(double fx, double fy, double fz)
		{f[0]=fx; f[1] = fy; f[2] = fz;};

void Planet :: setName(string name)
		{nm = name;};

void Planet :: setMass(double m)
		{mass = m;};


string Planet :: getName()
		{return nm;};

double Planet :: getMass()
		{return mass;}

double * Planet :: getPosition()
		{return x;};

double * Planet :: getVelocity()
		{return v;};

double * Planet :: getForces()
		{return f;};


void Planet :: updateForces()
		{ };
void Planet :: updateVelocities()
		{ };
void Planet :: updatePosition()
		{ };
*/

double MoonTheta(double time)
	{ return time*2*3.1415926535/(2.3718*pow(10,6));};

double MoonPhi(double theta)
	{ return (5.15*3.1415926535/180.0)*cos(theta);};

bool firstBurn(Planet mod,double time, double *a, double *dt)
	{
	double theta;
	double T;
	double Mtheta; double phi; double Mx, My; double theta2;
	double dist = 	(382848220+837800);
	double ellip_x,ellip_y; 
	double arrival;
	theta = atan((mod.y/mod.x)); //angle of the Module around earth
	
	if (mod.x<0)
	{
	ellip_x = (dist*cos(theta));//compare this to the moon position at future time T
	ellip_y = (dist*sin(theta));
	}

	else        
	{
	theta = theta + 3.1415926535;
	ellip_x = (dist*cos(theta));//compare this to the moon position at future time T
	ellip_y = (dist*sin(theta));
	}
	
	*a = 0.5*sqrt(pow(mod.x,2) + pow(mod.y,2)) + .5*sqrt(pow(ellip_x,2) + pow(ellip_y,2));
	T = 3.1415926535*(pow(*a,1.5))/sqrt(3.986*pow(10,14)); //half period of transfer orbit
	
	
	Mtheta = MoonTheta(time+ .93*T);	//future angle and positions of the moon
	phi = MoonPhi(Mtheta);
	Mx = cos(phi)*cos(Mtheta)*384400000;
	My = cos(phi)*sin(Mtheta)*384400000;

	theta2 = abs(atan(My/Mx));
	
	adaptDT(ellip_x,ellip_y, Mx+cos(Mtheta)*837800, My+837800*sin(Mtheta), .1, dt);
	
	arrival = sqrt(pow(ellip_x-(Mx+cos(Mtheta)*837800),2) + pow(ellip_y-(My+837800*sin(Mtheta)),2));
	if(arrival <= 1.6*pow(10,6) /*&&abs(abs(tan(theta2))-abs(tan(theta))) <=.1(*/) {
		//cout << "ellip_x: " << ellip_x << " ellip_y = " << ellip_y << " Mx = " << Mx << " My = " << My <<endl;
				  *dt = 1; 	
				  cout<< "Time at which they'll arrive "<< T + time<< endl;
				  cout << "Ellip_x: "<<ellip_x << " Ellip_y: " << ellip_y<< endl;
				  return true;}

	else		        return false;

	
	}; //end firstBurn


	void adaptDT(double x1, double y1, double x2, double y2, double a, double *dt)
	{
	double distance =  sqrt(pow(x1-x2,2)+pow(y1-y2,2));
	
	if(distance > pow(10,8) && *dt != a) {*dt = a;} 
	if(distance <= 1.6*pow(10,6)) {*dt = 1;  cout << "less than 5.1*10^5. " <<distance << endl;}
	else if(distance <=5*pow(10,7)) {*dt = .001;}//cout << "less than 10^6. dist: "<<distance << endl;}
	else if(distance <= pow(10,8)) {*dt = .01; }//cout << "less than 10^7. dist: "<<distance << endl;}
	else {};

    //else if(distance <= pow(10,7)) {*dt = .001; count++;cout << "Less than 10^7. count : " <<count << endl;}

   

	
	};


	void burnBaby(Planet *mod, double time)
	{double theta;
	double a;
	double T;
	double Mtheta; double Mx, My; double phi;
	double dist = 	(382848220+837800);
	double ellip_x,ellip_y; 
	double arrival;
	theta = atan((mod->y/mod->x)); //angle of the Module around earth
	
	if (mod->x<0)
	{
	ellip_x = (dist*cos(theta));//compare this to the moon position at future time T
	ellip_y = (dist*sin(theta));
	}

	else        
	{
	theta = theta + 3.1415926535;
	ellip_x = (dist*cos(theta));//compare this to the moon position at future time T
	ellip_y = (dist*sin(theta));
	}
	
	a = 0.5*sqrt(pow(mod->x,2) + pow(mod->y,2)) + .5*sqrt(pow(ellip_x,2) + pow(ellip_y,2));

	T = 3.1415926535*(pow(a,1.5))/sqrt(3.986*pow(10,14)); //half period of transfer orbit
	
	
	Mtheta = MoonTheta(time+T);	//future angle and positions of the moon
	phi = MoonPhi(Mtheta);
	Mx = cos(phi)*cos(Mtheta)*384400000;
	My = cos(phi)*sin(Mtheta)*384400000;
	//Mx = 382848220*cos(Mtheta);
	//My = 382848220*sin(Mtheta);
		


	double va1, va2, rloga,relogc;

	rloga = sqrt(pow(mod->x,2) + pow(mod->y,2));
	relogc = 2*a - rloga;

	cout << "rloga: " << rloga << " relogc: " << relogc << endl;

	va1 = sqrt(pow(mod->vx,2)+ pow(mod->vy,2));
	cout << "VA1 = " << va1 << endl;
	va2 = sqrt(2*3.986*pow(10,14))*sqrt((rloga * relogc)/(rloga + relogc))/rloga;
	cout << "VA2 = " << va2 << endl;
	theta = theta - 3.1415926535/2;
	cout << "DeltaV = " << va2-va1 << endl;

	cout << "vx = " << mod->vx << endl;
	cout << "vy = " << mod->vy << endl;
	mod->vx += (va2-va1)*cos(theta);
	mod->vy += (va2-va1)*sin(theta);
	cout << "vx = " << mod->vx << endl;
	cout << "vy = " << mod->vy << endl;

	
	};


	bool secondBurn(Planet mod, double time, double a, double *dt)

	{double Mtheta; double Mx, My, Mz; double phi;double rloga; double Vmag; double dot;
	vector<double> newPlane(3,0);
	Mtheta = MoonTheta(time);	//future angle and positions of the moon
	phi = MoonPhi(Mtheta);
	Mx = cos(phi)*cos(Mtheta)*384400000;
	My = cos(phi)*sin(Mtheta)*384400000;
	Mz= sin(phi)*384400000;

	//Mx = 382848220*cos(Mtheta);
	//My = 382848220*sin(Mtheta);
	//adaptDT(mod.x, mod.y, Mx,My, 1, dt);


	//Trying new test condition
	rloga = sqrt(pow(mod.x - Mx,2) + pow(mod.y - My,2)+ pow(mod.z - Mz,2));
	newPlane[0] = (Mx - mod.x)/rloga;newPlane[1] = (My - mod.y)/rloga;newPlane[2] = (Mz - mod.z)/rloga;
	Vmag = sqrt(pow(mod.vx,2) + pow(mod.vy,2) +pow(mod.vz,2));
	dot = newPlane[0]*mod.vx/Vmag + newPlane[1]*mod.vy/Vmag + newPlane[2]*mod.vz/Vmag;





	double arrival = sqrt(pow(mod.x-Mx,2) + pow(mod.y-My,2) + pow(mod.z-Mz,2));
	if(arrival <= 3*pow(10,7) && abs(dot)<=.05) { *dt = 1;return true;}

	else {return false;}
	};// Second Burn


	void burnBabyBurn(Planet *mod, double a, double time) 
	{double theta;
	double T;
	double Mtheta; double Mx, My,Mz; double phi;
	double dist = 	(382848220+837800);
	double ellip_x,ellip_y; 
	double arrival;
	double h;
	vector<double> newPlane(3,0);
	double Vmag;

	Mtheta = MoonTheta(time);	//future angle and positions of the moon
	phi = MoonPhi(Mtheta);
	Mx = cos(phi)*cos(Mtheta)*384400000;
	My = cos(phi)*sin(Mtheta)*384400000;
	Mz = sin(phi)*384400000;


	

	theta = atan((mod->y)/(mod->x)); //angle of the Module around earth

	
	double va1, va2, rloga,relogc,DeltaV,MoonV,burnPhi, e,rp;

	rloga = sqrt(pow(mod->x - Mx,2) + pow(mod->y - My,2)+ pow(mod->z - Mz,2));
	e = 0.8;
	rp = rloga*(1-e)/(1+e);
	newPlane[0] = (Mx - mod->x)/rloga;newPlane[1] = (My - mod->y)/rloga;newPlane[2] = (Mz - mod->z)/rloga;
	

	cout << "distance from Module and Moon " << rloga << endl;

	Vmag = sqrt(pow(mod->vx,2) + pow(mod->vy,2) +pow(mod->vz,2));
	cout << "Magnitude of Module's Velocity: " << Vmag << endl;
	

	cout << "dot Product: " << newPlane[0]*mod->vx/Vmag + newPlane[1]*mod->vy/Vmag + newPlane[2]*mod->vz/Vmag << endl;

	h = sqrt(rp*4.905*pow(10,12)*1.8);
	MoonV = sqrt(3.986*pow(10,14)/384400000);
	cout <<"Magnitude of Moon's Velocity: "<< MoonV << endl;

	DeltaV = h/rloga - (MoonV - Vmag);
	burnPhi = asin((Mz - mod->z)/rloga);

	cout <<"Relative velocity between Moon and Module: " << MoonV - Vmag << endl;

	cout << " DeltaV = "<< DeltaV <<endl;

	/*
	va1 = sqrt(pow(mod->vx,2)+ pow(mod->vy,2));
	cout << "VA1 = " << va1 << endl;
	va2 = sqrt(rloga*4.905*pow(10,12))/sqrt(pow(mod->x - Mx,2) + pow(mod->y - My,2));
	
	cout << "VA2 = " << va2 << endl;
	theta = theta - 3.1415926535/2;
	cout << "DeltaV = " << va2-va1 << endl;

	*/
	cout << "vx = " << mod->vx << endl;
	cout << "vy = " << mod->vy << endl;
	cout << "vz = " << mod->vz << endl;
	mod->vx += (mod->vx/Vmag)*DeltaV - MoonV*sin(Mtheta);
	mod->vy += (mod->vy/Vmag)*DeltaV + MoonV*cos(Mtheta);
	mod->vz += (mod->vz/Vmag)*DeltaV;
	cout << "vx = " << mod->vx << endl;
	cout << "vy = " << mod->vy << endl;
	cout << "vz = " << mod->vz << endl;
	DeltaV = sqrt(pow((mod->vx/Vmag)*DeltaV - MoonV*sin(Mtheta),2) + pow((mod->vy/Vmag)*DeltaV + MoonV*cos(Mtheta),2) + pow((mod->vz/Vmag)*DeltaV,2)); 
	cout <<"Actual total Delta V forSecond Burn: " << DeltaV<< endl;
}; //Second burn
