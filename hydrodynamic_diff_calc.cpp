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
int main() {
// current date/time based on current system
	time_t now = time(0);
	struct tm *ltm = localtime(&now);
	cout << "start time"<< '\t'<< ltm->tm_hour << ":";
	cout << ltm->tm_min << ":";
	cout << ltm->tm_sec << endl;
   
	int NrParticles = 4; 
   
	int mu_6N[6*NrParticles][6*NrParticles];
   
	mtrx3D Mobility_Tnsr_tt;
	mtrx3D Mobility_Tnsr_tr;
	mtrx3D Mobility_Tnsr_rt;
	mtrx3D Mobility_Tnsr_rr;
   
   
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
				Mobility_Tnsr_tt	=	 	Unit_diag * tau;
											
				Mobility_Tnsr_rr	=		Unit_diag*temp2 ;
				Mobility_Tnsr_rt	= 		null33D ;
				Mobility_Tnsr_tr	= 		null33D ;

			 } else {
				Mobility_Tnsr_tt	=	(	Unit_diag
											+	(Rij^Rij)*Rij2_inv
											+	(Unit_diag*(1.0/3.0)-(Rij^Rij)*Rij2_inv)*(particle[i].radius*particle[i].radius+particle[j].radius*particle[j].radius)*Rij2_inv
											)	*	temp;
				
				Mobility_Tnsr_rr	=		(Unit_diag*(-1.0) + (Rij^Rij)*Rij2_inv*3.0)*temp3 ;
				
				Mobility_Tnsr_rt	=	  	epsilon_rij*(-2.0)*temp3;
				Mobility_Tnsr_tr	= 	Mobility_Tnsr_rt[i][j];
												
			 } 
			
			for (int k=0; k<3; k++)
				{
					for (int l=0; l<3; l++)
						{
							mu_6N[k+3*i][l+3*i] 								=	 Mobility_Tnsr_tt.comp[k][l];
							mu_6N[k+3*i][l+3*i+3*NrParticles] 					=	 Mobility_Tnsr_tr.comp[k][l];
							mu_6N[k+3*i+3*NrParticles][l+3*i] 					=	 Mobility_Tnsr_rt.comp[k][l];
							mu_6N[k+3*i+3*NrParticles][l+3*i+3*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];
						}
				}	
				
			}	
	}	
