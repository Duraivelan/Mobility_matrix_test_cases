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


extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
}


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
         
int if_create_particles = 0, ifrestart=1;
         
double kb=1 , T0=0.3, tauT=0.1;
double Temp=0;
double shear_rate = 0; //shear rate
int ifshear = 0;// set equal to 1 for shear
std::string dataFileName="1",dataFileName_new="1" ;
int NrParticles=2;
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

// variables for mobility tensor calculation
double eta_0=0.1;
vctr3D Rij , rij_unit ;
double Rij2, Rij2_inv, temp, temp1, temp2, temp3, tau ;
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
      
	mtrx3D Mobility_Tnsr_tt;
	mtrx3D Mobility_Tnsr_tr;
	mtrx3D Mobility_Tnsr_rt;
	mtrx3D Mobility_Tnsr_rr;
   
   double mu_6N[36*NrParticles*NrParticles] ;  		// grand mobility matrix
   double zeta_6N[36*NrParticles*NrParticles] ;  	// grand resistance matrix
   double xi_6x6[6*6] ; 							// generalized friction matrix
   
   for (int i=0; i<36; i++)
		{
			xi_6x6[i] = 0.0; 
		}        
		
std::ofstream outFile1(dataFileName+"/data.dat");
		
for (int i=0; i<NrParticles; i++)
	{
		for (int j=0; j<NrParticles; j++)
			{
				cout<<(i==j)<<endl;
				Rij=particle[i].pos-particle[j].pos;
			//	rij_unit = 	Rij*(1.0/sqrt(Rij2)); 
				mtrx3D Pij(Rij, Rij) ;
				Rij2=Rij.norm2();
				vctr3D col1(0.0, -Rij.comp[2], Rij.comp[1]);
				vctr3D col2(Rij.comp[2],0.0,-Rij.comp[0]);
				vctr3D col3(-Rij.comp[1],Rij.comp[0],0.0);
				mtrx3D epsilon_rij(col1 , col2 , col3);
				Rij2_inv=1.0/Rij2;
				temp1=1.0/(8.0*M_PI*eta_0);
				tau = 1.0/(6.0*M_PI*eta_0*particle[i].radius);
				temp=temp1/(sqrt(Rij2));
				temp2=temp1/(particle[i].radius*particle[i].radius*particle[i].radius);
				temp3=temp/(2.0*(Rij2));
				if(i==j) {
				Mobility_Tnsr_tt	=	 	Unit_diag * tau ;
											
				Mobility_Tnsr_rr	=		Unit_diag * temp2 ;
				Mobility_Tnsr_rt	= 		null33D ;
				Mobility_Tnsr_tr	= 		null33D ;

			 } else {
				Mobility_Tnsr_tt	=	(	Unit_diag
											+	(Pij)*Rij2_inv
											+	(Unit_diag*(1.0/3.0)-(Pij)*Rij2_inv)*(particle[i].radius*particle[i].radius+particle[j].radius*particle[j].radius)*Rij2_inv
											)	*	temp;
				
				Mobility_Tnsr_rr	=		(Unit_diag*(-1.0) + (Pij)*Rij2_inv*3.0)*temp3 ;
				
				Mobility_Tnsr_tr	=	  	epsilon_rij*(-2.0)*temp3;
				Mobility_Tnsr_rt	= 	    Mobility_Tnsr_tr*(-1.0);												
			 } 
			
			for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{
					/*		mu_6N[l+3*k+9*j+i*9*NrParticles] 								=	 Mobility_Tnsr_tt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+ 9*NrParticles*NrParticles] 	=	 Mobility_Tnsr_tr.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+18*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+27*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];		*/					
							// column major format
							zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j] 					=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles] 	=	 Mobility_Tnsr_tr.comp[k][l];
							zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j+3*NrParticles] 	=	 Mobility_Tnsr_rt.comp[k][l];
							zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles+3*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];
						}
				}	
				
			}	

		}
		
	mtrx3D Resistance_Tnsr_tt;
	mtrx3D Resistance_Tnsr_tr;
	mtrx3D Resistance_Tnsr_rt;
	mtrx3D Resistance_Tnsr_rr;	
	
	mtrx3D Friction_Tnsr_tt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_tr(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rr(0.0,0.0,0.0);
   			
	inverse ( zeta_6N , 6*NrParticles )	 ; 							
					
	for (int i=0; i<NrParticles; i++)
		{
			vctr3D col1 (  0.0 						,	particle[i].pos.comp[2] , -particle[i].pos.comp[1]	);
			vctr3D col2 ( -particle[i].pos.comp[2] 	,   0.0						,  particle[i].pos.comp[0]	);
			vctr3D col3 (  particle[i].pos.comp[1] 	,  -particle[i].pos.comp[0] , 0.0						);
			mtrx3D	Ai(	col1 , col2 , col3 ) ;   
					
			for (int j=0; j<NrParticles; j++)
				{
					
					vctr3D col4 (  0.0 						,	particle[j].pos.comp[2] , -particle[j].pos.comp[1]	);
					vctr3D col5 ( -particle[j].pos.comp[2] 	,   0.0						,  particle[j].pos.comp[0]	);
					vctr3D col6 (  particle[j].pos.comp[1] 	,  -particle[j].pos.comp[0] , 0.0						);
					mtrx3D	Aj(	col4 , col5 , col6 ) ; 
					
					for (int l=0; l<3; l++)
						{
							for (int k=0; k<3; k++)
								{
							
									Resistance_Tnsr_tt.comp[k][l] = zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j] ;
									Resistance_Tnsr_tr.comp[k][l] = zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles]	;
									Resistance_Tnsr_rt.comp[k][l] = zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j+3*NrParticles] ;
									Resistance_Tnsr_rr.comp[k][l] = zeta_6N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles+3*NrParticles] ;
								}
						}
					
					Friction_Tnsr_tt += Resistance_Tnsr_tt ;  
					Friction_Tnsr_tr += ( Resistance_Tnsr_tr    -  ( Resistance_Tnsr_tt*Aj )  ) ;
					Friction_Tnsr_rt += ( Ai*Resistance_Tnsr_tt +    Resistance_Tnsr_rt    )    ; 
					Friction_Tnsr_rr += ( Resistance_Tnsr_rr    -  ( Resistance_Tnsr_rt*Aj ) + Ai*Resistance_Tnsr_tr - Ai*Resistance_Tnsr_tt*Aj ) ; 
					
			/*		xi_6x6[0] += zeta_6N[	9*j + i*9*NrParticles] ;  
					xi_6x6[1] += zeta_6N[1+ 9*j + i*9*NrParticles] ;  
					xi_6x6[2] += zeta_6N[2+ 9*j + i*9*NrParticles] ;  
					xi_6x6[3] += zeta_6N[3+ 9*j + i*9*NrParticles] ;  
					xi_6x6[4] += zeta_6N[4+ 9*j + i*9*NrParticles] ;  
					xi_6x6[5] += zeta_6N[5+ 9*j + i*9*NrParticles] ;  
					xi_6x6[6] += zeta_6N[6+ 9*j + i*9*NrParticles] ;  
					xi_6x6[7] += zeta_6N[7+ 9*j + i*9*NrParticles] ;  
					xi_6x6[8] += zeta_6N[8+ 9*j + i*9*NrParticles] ;  					
					
					xi_6x6[9]  += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[0] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;  
					xi_6x6[10] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[1] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;  
					xi_6x6[11] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[2] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;  
					xi_6x6[12] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[3] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;  
					xi_6x6[13] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[4] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;    
					xi_6x6[14] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[5] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;  
					xi_6x6[15] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[6] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;  
					xi_6x6[16] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[7] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;    
					xi_6x6[17] += ( -zeta_6N[	9*j + i*9*NrParticles]*Aj[8] + zeta_6N[	 9*j + i*9*NrParticles + 9*NrParticles*NrParticles]) ;  
					
					xi_6x6[18] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[0] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;  
					xi_6x6[19] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[1] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;  
					xi_6x6[20] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[2] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;  
					xi_6x6[21] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[3] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;  
					xi_6x6[22] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[4] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;    
					xi_6x6[23] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[5] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;  
					xi_6x6[24] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[6] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;  
					xi_6x6[25] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[7] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;    
					xi_6x6[26] += ( zeta_6N[	9*j + i*9*NrParticles]*Ai[8] + zeta_6N[	 9*j + i*9*NrParticles + 18*NrParticles*NrParticles]) ;   

					xi_6x6[27] += zeta_6N[	 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[28] += zeta_6N[1+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[29] += zeta_6N[2+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[30] += zeta_6N[3+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[31] += zeta_6N[4+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[32] += zeta_6N[5+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[33] += zeta_6N[6+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[34] += zeta_6N[7+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  
					xi_6x6[35] += zeta_6N[8+ 9*j + i*9*NrParticles +27*NrParticles*NrParticles] ;  						*/				 
				}
		}
					//Friction_Tnsr_rt = ~Friction_Tnsr_tr;

									// column major format

					xi_6x6[0] = Friction_Tnsr_tt.comp[0][0] ;  
					xi_6x6[1] = Friction_Tnsr_tt.comp[0][1] ;  
					xi_6x6[2] = Friction_Tnsr_tt.comp[0][2] ; 
					xi_6x6[6] = Friction_Tnsr_tt.comp[1][0] ; 
					xi_6x6[7] = Friction_Tnsr_tt.comp[1][1] ;  
					xi_6x6[8] = Friction_Tnsr_tt.comp[1][2] ;  
					xi_6x6[12] = Friction_Tnsr_tt.comp[2][0] ;   
					xi_6x6[13] = Friction_Tnsr_tt.comp[2][1] ; 
					xi_6x6[14] = Friction_Tnsr_tt.comp[2][2] ; 				

					xi_6x6[18] = Friction_Tnsr_rt.comp[0][0] ;  
					xi_6x6[19] = Friction_Tnsr_rt.comp[0][1] ;  
					xi_6x6[20] = Friction_Tnsr_rt.comp[0][2] ; 
					xi_6x6[24] = Friction_Tnsr_rt.comp[1][0] ; 
					xi_6x6[25] = Friction_Tnsr_rt.comp[1][1] ;  
					xi_6x6[26] = Friction_Tnsr_rt.comp[1][2] ;  
					xi_6x6[30] = Friction_Tnsr_rt.comp[2][0] ;   
					xi_6x6[31] = Friction_Tnsr_rt.comp[2][1] ; 
					xi_6x6[32] = Friction_Tnsr_rt.comp[2][2] ; 				
										
					xi_6x6[3] = Friction_Tnsr_tr.comp[0][0] ;  
					xi_6x6[4] = Friction_Tnsr_tr.comp[0][1] ;  
					xi_6x6[5] = Friction_Tnsr_tr.comp[0][2] ; 
					xi_6x6[9] = Friction_Tnsr_tr.comp[1][0] ; 
					xi_6x6[10] = Friction_Tnsr_tr.comp[1][1] ;  
					xi_6x6[11] = Friction_Tnsr_tr.comp[1][2] ;  
					xi_6x6[15] = Friction_Tnsr_tr.comp[2][0] ;   
					xi_6x6[16] = Friction_Tnsr_tr.comp[2][1] ; 
					xi_6x6[17] = Friction_Tnsr_tr.comp[2][2] ; 				
										
					xi_6x6[21] = Friction_Tnsr_rr.comp[0][0] ;  
					xi_6x6[22] = Friction_Tnsr_rr.comp[0][1] ;  
					xi_6x6[23] = Friction_Tnsr_rr.comp[0][2] ; 
					xi_6x6[27] = Friction_Tnsr_rr.comp[1][0] ; 
					xi_6x6[28] = Friction_Tnsr_rr.comp[1][1] ;  
					xi_6x6[29] = Friction_Tnsr_rr.comp[1][2] ;  
					xi_6x6[33] = Friction_Tnsr_rr.comp[2][0] ;   
					xi_6x6[34] = Friction_Tnsr_rr.comp[2][1] ; 
					xi_6x6[35] = Friction_Tnsr_rr.comp[2][2] ; 				
		
	inverse ( xi_6x6 , 6 )	 ; 			
	for (int i=0; i<36; i++)
		{
			xi_6x6[i]*=4.0315e-14;	// multiply by kbT in erg K-1
		} 	
			
	for (int i=0; i<NrParticles; i++)
		{
			outFile1<<"particle position"<<'\t'<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<'\t'<<particle[i].radius<<std::endl ;
		} 	
		outFile1<<std::endl ;
		outFile1<<xi_6x6[0]<<'\t'<<xi_6x6[6]<<'\t'<<xi_6x6[12]<<'\t'<<xi_6x6[18]<<'\t'<<xi_6x6[24]<<'\t'<<xi_6x6[30]<<std::endl ;
		outFile1<<xi_6x6[1]<<'\t'<<xi_6x6[7]<<'\t'<<xi_6x6[13]<<'\t'<<xi_6x6[19]<<'\t'<<xi_6x6[25]<<'\t'<<xi_6x6[31]<<std::endl ;
		outFile1<<xi_6x6[2]<<'\t'<<xi_6x6[8]<<'\t'<<xi_6x6[14]<<'\t'<<xi_6x6[20]<<'\t'<<xi_6x6[26]<<'\t'<<xi_6x6[32]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_6x6[3]<<'\t'<<xi_6x6[9]<<'\t'<<xi_6x6[15]<<'\t'<<xi_6x6[21]<<'\t'<<xi_6x6[27]<<'\t'<<xi_6x6[33]<<std::endl ;
		outFile1<<xi_6x6[4]<<'\t'<<xi_6x6[10]<<'\t'<<xi_6x6[16]<<'\t'<<xi_6x6[22]<<'\t'<<xi_6x6[28]<<'\t'<<xi_6x6[34]<<std::endl ;
		outFile1<<xi_6x6[5]<<'\t'<<xi_6x6[11]<<'\t'<<xi_6x6[17]<<'\t'<<xi_6x6[23]<<'\t'<<xi_6x6[29]<<'\t'<<xi_6x6[35]<<std::endl ;


}
		
