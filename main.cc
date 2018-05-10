//This is going to be the main folder


#include "lib.h"
#include<iostream>
#include<vector>
#include<math.h>
#include<sstream>
#include<fstream>
#include<list>

using namespace std;


int main()
{

/* The planets are in the following order:
0 - Earth
1 - Moon
2 - Logistics module 
3 - DSG
*/

//Now to read in the data from the csv files
//3/22/caleb

ifstream initial;
std::string initial_line;
//initial.open("earthmoon.txt");
initial.open("planetinitial.txt");
ofstream out_file;
out_file.open("simulation.txt");
vector<Planet> solarsystem;
if(!initial.is_open())
	std::cout << "Error opening planetinitial.\n";

int count = 0;
   for(int l=0;l<11;l++){
	double x,y,z,vx,vy,vz,mass,extra;
	string name;
        std::getline(initial,initial_line);
	std::istringstream(initial_line) >> x >> y >> z >> vx >> vy >> vz >> mass >> name;
	Planet a;
	a.x = x;
	a.y = y;
	a.z = z;
	a.vx = vx;
	a.vy = vy;
	a.vz = vz;
	a.mass = mass;
	solarsystem.push_back(a);
	count++;
	cout << "x = " << z << endl;
	cout << "vy = " << vz << endl;
	cout << "mass = " << mass << endl;
}
initial.close();
//simulate here
//solarsystem[i].x or .vx

double a;
double G =6.67384* pow(10,-11);
double rij,fij,fxij,fyij,fzij;
double sum_fx = 0;
double sum_fy = 0;
double sum_fz = 0;
double dt = 1E4;
double time = 0;



while(time < 3.6E9){

for(int i=1;i<11;i++){ //change back to i <
	sum_fx = 0;
	sum_fy = 0;
	sum_fz = 0;
    for(int j=0;j< 11;j++) {
	if(i!=j){
	   rij = sqrt( pow(solarsystem[i].x - solarsystem[j].x,2)+pow(solarsystem[i].y - solarsystem[j].y,2) + pow(solarsystem[i].z - solarsystem[j].z,2) );
	   fij = (G * solarsystem[i].mass * solarsystem[j].mass) / pow(rij,2);	

           fxij = -fij * (solarsystem[i].x - solarsystem[j].x)/rij;
           fyij = -fij * (solarsystem[i].y - solarsystem[j].y)/rij;
           fzij = -fij * (solarsystem[i].z - solarsystem[j].z)/rij;
       	}//end if
	else{
	fxij=0;
	fyij=0;
	fzij=0;
	} // end else
	sum_fx = sum_fx + fxij;
	sum_fy = sum_fy + fyij;
	sum_fz = sum_fz + fzij;
    }//end j

	solarsystem[i].vx = solarsystem[i].vx + dt*sum_fx/solarsystem[i].mass;
	solarsystem[i].vy = solarsystem[i].vy + dt*sum_fy/solarsystem[i].mass;
	solarsystem[i].vz = solarsystem[i].vz + dt*sum_fz/solarsystem[i].mass;
}// end i

  out_file << time << endl;

for(int i=0;i<solarsystem.size();i++){  
	time = time+dt;
	solarsystem[i].x = solarsystem[i].x + dt*solarsystem[i].vx;
	solarsystem[i].y = solarsystem[i].y + dt*solarsystem[i].vy;
	solarsystem[i].z = solarsystem[i].z + dt*solarsystem[i].vz;
	//Writing Results
    out_file << solarsystem[i].x << " " << solarsystem[i].y << " "<< solarsystem[i].z << " "<<solarsystem[i].vx<< " " << solarsystem[i].vy << " " << solarsystem[i].vz << endl;
}// Writing Results Loops

}//While



return 0;
}//main 

