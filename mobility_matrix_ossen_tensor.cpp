#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <dirent.h>
# include "structure_definitions.h"
# include "defs.h"
using namespace std;


// random numbers using random_device option with normal distribution
/*void createInitialPosition_N_particles(std::string fileName, int N) {
      std::random_device rd, rd1;
      std::mt19937 genx(rd()),geny(rd1());
      std::ofstream outFile(fileName);
      std::normal_distribution<> d(0,1);
      for(int i=0;i<N;i++) {
      outFile<<d(genx)<<'\t'<<d(geny)<<std::endl;
      }
      outFile.close();
}
*/

// random numbers using rand function
void createInitialPosition_N_particles(std::string fileName,std::string fileName2, int N, double Lx, double Ly, double Lz) {
	double x,y,z, vx, vy, vz;
 	srand (time(NULL)); // initialize random seed
 	std::ofstream outFile(fileName);
 	std::ofstream outFile2(fileName2);
 	for(int i=0;i<N;i++) {
 		x=((double) rand() / (RAND_MAX/Lx))-Lx/2;  // create particle position from -Lx/2 to Lx/2
		y=((double) rand() / (RAND_MAX/Ly))-Ly/2;
		z=((double) rand() / (RAND_MAX/Lz))-Lz/2;
		vx= ((double) rand()/(RAND_MAX)-0.5);
		vy= ((double) rand()/(RAND_MAX)-0.5);
		vz= ((double) rand()/(RAND_MAX)-0.5);
		outFile<<x<<'\t'<<y<<'\t'<<z<<std::endl;
		outFile2<<vx<<'\t'<<vy<<'\t'<<vz<<std::endl;
 	}
 	outFile.close();
}

std::vector<int> radialDistFunc(double XYZ[][3], double Lx,double Ly, double Lz, double dr, int N) {
    std::vector<int> rdf((int) floor(sqrt(pow(Lx/2,2)+pow(Ly/2,2)+pow(Lz/2,2)))/dr,0);
	double r;
	for(int j=0;j<N;j++) {
		r=sqrt(pow(XYZ[j][0],2)+pow(XYZ[j][1],2)+pow(XYZ[j][2],2));
	    rdf[(int) floor(r/dr)]+=1;                        // put each particle in a bin according to its position from origin ie. (0,0)
	}
	return rdf;
}

// forceUpdate fucntion included as force.h header file
void verlet( vector<SubData>& particle ) {
	
	for(int i=0;i<NrParticles;i++) 
	{
		particle[i].vel+=particle[i].frc*(0.5*dt*inv_mass);
		particle[i].pos+=particle[i].vel*dt;
		particle[i].pos.PBC(box,rbox);

	}
}

void verletB(vector<SubData>& particle, double vel_scale) {
	if(0) 
		{
		for(int i=0;i<NrParticles;i++) 
			{
				particle[i].vel+=particle[i].frc*(0.5*dt*inv_mass);
				particle[i].vel=(particle[i].vel)*vel_scale;
			}
       	} 
	else 
		{
		for(int i=0;i<NrParticles;i++) 
			{
				particle[i].vel+=particle[i].frc*(0.5*dt*inv_mass);
			}
       	}
}

int main() {
// current date/time based on current system
   time_t now = time(0);
   struct tm *ltm = localtime(&now);
   cout << "start time"<< '\t'<< ltm->tm_hour << ":";
   cout << ltm->tm_min << ":";
   cout << ltm->tm_sec << endl;
         
int if_create_particles = 0, ifrestart=1;
         
double kb=1 , T0=0.3, tauT=0.1;
double Temp=0;
double shear_rate = 0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="../1",dataFileName_new="../1" ;
int NrParticles=4;
double simu_time=dt;
int step=0, nSteps=10000, frame=10;
double vel_scale;
int if_Periodic =1;
vector<vector<vector<vector<int>>>>  grid(ceil ( box.comp[x] / (r_cut) ),vector<vector<vector<int>>>(ceil ( box.comp[y] / (r_cut) ),vector<vector<int>>(ceil ( box.comp[z] / (r_cut) ),vector<int>(50,0))));
std::cout<<cellx<<'\t'<<celly<<'\t'<<cellz<<std::endl;
double  T_Energy, K_Energy, P_Energy, p_energy=0;
vctr3D dR, dr2;
double R, r2;
double dr=0.05; // step size for RDF calculation
// std::vector<int> RDF((int)  floor(sqrt((Lx/2)*(Lx/2)+(Ly/2)*(Ly/2)+(Lz/2)*(Lz/2)))/dr,0), RDF1((int)  floor(sqrt(Lx*Lx+Ly*Ly))/dr,0);
double KE_rot=0;

vector<SubData>  particle(NrParticles);
vector<vector<mtrx3D>>  Mobility_Tnsr_tt(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr_tr(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr_rt(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr_rr(NrParticles,vector<mtrx3D>(NrParticles));
vector<vector<mtrx3D>>  Mobility_Tnsr(NrParticles,vector<mtrx3D>(NrParticles));

// variables for mobility tensor calculation
double eta_0=0.01;
vctr3D Rij;
double Rij2, Rij2_inv, temp, temp1, temp2, temp3, tau;
//

if(ifrestart)	{
std::string fileName=dataFileName+"/End_positions.dat";

//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
if(!dataFile.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);    
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
		currentLine >> particle[i].radius;
		cout<<particle[i].radius<<endl;
    }
}	
	fileName=dataFileName+"/Velocities.dat";
	std::ifstream dataFile1(fileName);

	if(!dataFile1.good()) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading Vx, Vy, Vz Velocities"<<std::endl;
    std::string line;
     
    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile1,line);
    	std::istringstream currentLine(line);
        currentLine >> particle[i].vel.comp[0];
        currentLine >> particle[i].vel.comp[1];
        currentLine >> particle[i].vel.comp[2];

    }
}	
} else {

	std::string fileName="../XYZ.dat";
	std::string fileName2="../Initial_Velocities.dat";

if (if_create_particles) {
createInitialPosition_N_particles(fileName,fileName2,NrParticles,Lx,Ly,Lz);
}
//read x,y positions from XY.dat
std::ifstream dataFile(fileName);
std::ifstream dataFile2(fileName2);

if(!dataFile.good() && !dataFile2.good() ) {
	std::cerr<<"Given file is corrupt /n"<<std::endl;
}
else {
	std::cout<<"Reading X, Y, Z co-ordinates"<<std::endl;
    std::string line;
    std::string line1;


    for (int i=0;i<NrParticles;i++) {
		std::getline(dataFile,line);
    	std::istringstream currentLine(line);
        currentLine >> particle[i].pos.comp[0];
        currentLine >> particle[i].pos.comp[1];
        currentLine >> particle[i].pos.comp[2];
        currentLine >> particle[i].radius;
		std::getline(dataFile2,line1);
    	std::istringstream currentLine1(line1);
        currentLine1 >> particle[i].vel.comp[0];
        currentLine1 >> particle[i].vel.comp[1];
        currentLine1 >> particle[i].vel.comp[2];  
        Temp+=0.5*m*(particle[i].vel.comp[0]*particle[i].vel.comp[0]
				   + particle[i].vel.comp[1]*particle[i].vel.comp[1]
				   + particle[i].vel.comp[2]*particle[i].vel.comp[2]);
              
    }
    	Temp=(Temp)/(1.5*NrParticles*kb);
		vel_scale = sqrt(T0/Temp);
		std::cout<<Temp<<'\t'<<vel_scale<<std::endl;

}

}		
//delete all files before writing data

// following snippet taken from stakcflow link  http://stackoverflow.com/questions/11007494/how-to-delete-all-files-in-a-folder-but-not-delete-the-folder-c-linux
if (ifrestart) {
dataFileName=dataFileName_new;
}
const char *dataFileNamePointer = dataFileName.c_str();  // covnert the datafilename to a char pointer ans pass it to the snippet below which delete all files in that folder before running the simulation
if (!ifrestart) {
struct dirent *next_file;
DIR *theFolder;
char filepath[256];
theFolder = opendir(dataFileNamePointer);
while (( next_file = readdir(theFolder)) )
	{
    // build the full path for each file in the folder
    sprintf(filepath, "%s/%s",dataFileNamePointer, next_file->d_name);
    remove(filepath);
    }
//
}

/* sort particles into cells */

// 	Garcia de la Torre-Bloomfield Tensor calculation
// 	based on the paper 
//	Improved Calculation of Rotational Diffusion and Intrinsic Viscosity of Bead Models for
//	Macromolecules and Nanoparticles, J. Phys. Chem. B 2007, 111, 955-961 .


/*
 * snippet to call the HYDRO++ routine 
 * 
std::ofstream outFile7("square.dat");
outFile7<<1<<",    !Unit of length for coordinates and radii, cm (10 A)"<<endl;
outFile7<<4<<",        !Number of beads"<<endl;

for (int i=0; i<NrParticles; i++)
	{
		outFile7<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<'\t'<<particle[i].radius<<std::endl;
	}
system("../diffusion_tensor/hydro++10-lnx.exe < ../diffusion_tensor/input.txt");
* 
*/

for (int i=0; i<NrParticles; i++)
	{
		for (int j=0; j<NrParticles; j++)
			{
				cout<<(i==j)<<endl;
				Rij=particle[i].pos-particle[j].pos;
				vctr3D col1(0, -Rij.comp[2], Rij.comp[1]);
				vctr3D col2(Rij.comp[2],0,-Rij.comp[0]);
				vctr3D col3(-Rij.comp[1],Rij.comp[0],0);
				mtrx3D epsilon_rij(col1 , col2, col3);
				Rij2=Rij.norm2();
				Rij2_inv=1/Rij2;
				temp1=1.0/(8.0*M_PI*eta_0);
				tau = 1.0/(6.0*M_PI*eta_0*particle[i].radius);
				temp=temp1/(sqrt(Rij2));
				temp2=temp1/(particle[i].radius*particle[i].radius*particle[i].radius);
				temp3=temp/(2.0*Rij2);
				if(i==j) {
				Mobility_Tnsr_tt[i][j]	=	 	Unit_diag * tau;
											
				Mobility_Tnsr_rr[i][j]	=		Unit_diag*temp2 ;
				Mobility_Tnsr_rt[i][j]	= 		null33D ;
				Mobility_Tnsr_tr[i][j]	= 		null33D ;

			 } else {
				Mobility_Tnsr_tt[i][j]	=	(	Unit_diag
											+	(Rij^Rij)*Rij2_inv
											+	(Unit_diag*(1.0/3.0)-(Rij^Rij)*Rij2_inv)*(particle[i].radius*particle[i].radius+particle[j].radius*particle[j].radius)*Rij2_inv
											)	*	temp;
				
				Mobility_Tnsr_rr[i][j]	=		(Unit_diag*(-1.0) + (Rij^Rij)*Rij2_inv*3.0)*temp3 ;
				
				Mobility_Tnsr_rt[i][j]	=	  	epsilon_rij*(-2.0)*temp3;
				Mobility_Tnsr_tr[i][j]	= 	Mobility_Tnsr_rt[i][j];
												
			 } 
			 	Mobility_Tnsr[i][j]		=	 Mobility_Tnsr_tt[i][j] -  Mobility_Tnsr_tr[i][j]*Mobility_Tnsr_rr[i][j]*Mobility_Tnsr_rt[i][j] ;

				
			}	
	}	
				
	
				
std::ofstream outFile10(dataFileName+"/molibity_tensor_tt.dat");
std::ofstream outFile11(dataFileName+"/molibity_tensor_tr.dat");
std::ofstream outFile12(dataFileName+"/molibity_tensor_rr.dat");
std::ofstream outFile13(dataFileName+"/molibity_tensor.dat");

		// save position, Kinetic energy, Potential energy, Forces every 'frame' steps
for (int i=0; i<NrParticles; i++)
	{
		for (int j=0; j<NrParticles; j++)
			{
				outFile10 << setw(10) << Mobility_Tnsr_tt[i][j].comp[0][0] << "  " << setw(10) << Mobility_Tnsr_tt[i][j].comp[0][1] << "  " << setw(10) << Mobility_Tnsr_tt[i][j].comp[0][2] << endl;
				outFile10 << setw(10) << Mobility_Tnsr_tt[i][j].comp[1][0] << "  " << setw(10) << Mobility_Tnsr_tt[i][j].comp[1][1] << "  " << setw(10) << Mobility_Tnsr_tt[i][j].comp[1][2] << endl;
				outFile10 << setw(10) << Mobility_Tnsr_tt[i][j].comp[2][0] << "  " << setw(10) << Mobility_Tnsr_tt[i][j].comp[2][1] << "  " << setw(10) << Mobility_Tnsr_tt[i][j].comp[2][2] << endl;
				
				outFile11 << setw(10) << Mobility_Tnsr_tr[i][j].comp[0][0] << "  " << setw(10) << Mobility_Tnsr_tr[i][j].comp[0][1] << "  " << setw(10) << Mobility_Tnsr_tr[i][j].comp[0][2] << endl;
				outFile11 << setw(10) << Mobility_Tnsr_tr[i][j].comp[1][0] << "  " << setw(10) << Mobility_Tnsr_tr[i][j].comp[1][1] << "  " << setw(10) << Mobility_Tnsr_tr[i][j].comp[1][2] << endl;
				outFile11 << setw(10) << Mobility_Tnsr_tr[i][j].comp[2][0] << "  " << setw(10) << Mobility_Tnsr_tr[i][j].comp[2][1] << "  " << setw(10) << Mobility_Tnsr_tr[i][j].comp[2][2] << endl;
				
				outFile12 << setw(10) << Mobility_Tnsr_rr[i][j].comp[0][0] << "  " << setw(10) << Mobility_Tnsr_rr[i][j].comp[0][1] << "  " << setw(10) << Mobility_Tnsr_rr[i][j].comp[0][2] << endl;
				outFile12 << setw(10) << Mobility_Tnsr_rr[i][j].comp[1][0] << "  " << setw(10) << Mobility_Tnsr_rr[i][j].comp[1][1] << "  " << setw(10) << Mobility_Tnsr_rr[i][j].comp[1][2] << endl;
				outFile12 << setw(10) << Mobility_Tnsr_rr[i][j].comp[2][0] << "  " << setw(10) << Mobility_Tnsr_rr[i][j].comp[2][1] << "  " << setw(10) << Mobility_Tnsr_rr[i][j].comp[2][2] << endl;

				outFile13 << setw(10) << Mobility_Tnsr[i][j].comp[0][0] << " \t " << setw(10) << Mobility_Tnsr[i][j].comp[0][1] << " \t " << setw(10) << Mobility_Tnsr[i][j].comp[0][2] << endl;
				outFile13 << setw(10) << Mobility_Tnsr[i][j].comp[1][0] << " \t " << setw(10) << Mobility_Tnsr[i][j].comp[1][1] << " \t " << setw(10) << Mobility_Tnsr[i][j].comp[1][2] << endl;
				outFile13 << setw(10) << Mobility_Tnsr[i][j].comp[2][0] << " \t " << setw(10) << Mobility_Tnsr[i][j].comp[2][1] << " \t " << setw(10) << Mobility_Tnsr[i][j].comp[2][2] << endl;

			}
	}					

outFile10.close();

std::ofstream outFile8(dataFileName+"/logfile");
	outFile8<<"NrParticles"<<'\t'<<NrParticles<<std::endl;
	outFile8<<"mass"<<'\t'<<m<<std::endl;
	outFile8<<"kb"<<'\t'<<kb<<std::endl;
	outFile8<<"T0"<<'\t'<<T0<<std::endl;
	outFile8<<"box"<<'\t'<<box.comp[0]<<'\t'<<box.comp[1]<<'\t'<<box.comp[2]<<std::endl;
	outFile8<<"shear rate"<<'\t'<<shear_rate<<std::endl;
	outFile8<<"R_cut"<<'\t'<<r_cut<<std::endl;
	outFile8<<"rs"<<'\t'<<rs<<std::endl;
	outFile8<<"epsilon"<<'\t'<<epsilon<<std::endl;
	outFile8<<"sigma"<<'\t'<<sigma<<std::endl;
outFile8.close();

     // get time now
	now = time(0);
	ltm = localtime(&now);
	cout << "end time"<< '\t'<< ltm->tm_hour << ":";
	cout << ltm->tm_min << ":";
	cout << ltm->tm_sec << endl;
return 0;

}
