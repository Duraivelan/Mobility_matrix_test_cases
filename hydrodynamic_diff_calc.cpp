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
int NrParticles=3;
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
vctr3D e_ab , e_ab_unit ;
double e_ab2, e_ab2_inv, temp, temp1, temp2, temp3, tau ;
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
	mtrx35D Mobility_Tnsr_td;
	mtrx35D Mobility_Tnsr_rd;
	mtrx53D Mobility_Tnsr_dt;
	mtrx53D Mobility_Tnsr_dr;
	mtrx55D Mobility_Tnsr_dd;
	
// Kronecker delta

	double kron_del[3][3] = {	
								{1.0,0.0,0.0},
								{0.0,1.0,0.0},
								{0.0,0.0,1.0}
							};
							

// Levi-Civita 

	double Levi_Civi[3][3][3] = {
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,-1.0,0.0}},
							{{0.0,0.0,-1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,1.0,0.0},{-1.0,0.0,0.0},{0.0,0.0,0.0}}
							};
	
// three and four index mobility matrices	
	double g_ijk[3][3][3];
	double h_ijk[3][3][3];
	double m_ijkl[3][3][3][3];
	
// the five base matrices for strain tensor // option 5 :  equation 419 wouter's tex version clusterdyn_110816_1556

	double e[5][3][3]= {
							{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
							{{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,-1.0}},
							{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
						};
   
   double mu_11N[121*NrParticles*NrParticles] ;  		// grand mobility matrix
   double zeta_11N[121*NrParticles*NrParticles] ;  	// grand resistance matrix
   double xi_11x11[11*11] ; 							// generalized friction matrix
   
   for (int i=0; i<121; i++)
		{
			xi_11x11[i] = 0.0; 
		}        
		
std::ofstream outFile1(dataFileName+"/data.dat");
		
for (int a=0; a<NrParticles; a++)
	{
		for (int b=a; b<NrParticles; b++)
			{
				cout<<(a==b)<<endl;
				e_ab=particle[a].pos-particle[b].pos;
				mtrx3D Pab(e_ab, e_ab) ;
				e_ab2=e_ab.norm2();
				vctr3D col1(0.0, -e_ab.comp[2], e_ab.comp[1]);
				vctr3D col2(e_ab.comp[2],0.0,-e_ab.comp[0]);
				vctr3D col3(-e_ab.comp[1],e_ab.comp[0],0.0);
				mtrx3D epsilon_e_ab(col1 , col2 , col3);
				e_ab2_inv=1.0/e_ab2;
				
				e_ab_unit = e_ab*sqrt(e_ab2_inv); 
				if(a==b) {	e_ab_unit = e_ab*0.0; }	
				temp1=1.0/(8.0*M_PI*eta_0);
				tau = 1.0/(6.0*M_PI*eta_0*particle[a].radius);
				temp=temp1/(sqrt(e_ab2));
				temp2=temp1/(particle[a].radius*particle[a].radius*particle[a].radius);
				temp3=temp/(2.0*(e_ab2));
			    double r 	= sqrt(e_ab2);
			    double r_1 	= 1.0/(r);
			    double r_2 	= 1.0/(r*r);			    
			    double r_3 	= 1.0/(r*r*r);
			    cout << r << '\t'<<r_1<< '\t'<<r_2 << '\t' << r_3 << 'a'<< a<< 'b'<< b<< endl;
				cout << e_ab_unit.comp[0] << "x-comp"<<endl; 

				// mobility scalar values - as defined in Page 46, Appendix A - Durlofsky, Louis, John F. Brady, and Georges Bossis. 
				// "Dynamic simulation of hydrodynamically interacting particles." Journal of fluid mechanics 180 (1987): 21-49.

				double x_a[2][2] = {{	1.0		,	3.0*r_1/2.0		-	1.0*r_3*particle[a].radius*particle[a].radius		},{	3.0*r_1/2.0		-	1.0*r_3*particle[a].radius*particle[a].radius		,	1.0		} }; 
				double y_a[2][2] = {{	1.0		,	(3.0*r_1/4.0)		+	(1.0*r_3*particle[a].radius*particle[a].radius/2.0)	},{	(3.0*r_1/4.0)		+	(1.0*r_3*particle[a].radius*particle[a].radius/2.0)	,	1.0		} }; 
				double y_b[2][2] = {{	0.0		,  -3.0*r_2/4.0						},{	3.0*r_2/4.0						,	0.0		} }; 
				double x_c[2][2] = {{	3.0/4.0	,  	3.0*r_3/4.0						},{	3.0*r_3/4.0						,  	3.0/4.0	} }; 
				double y_c[2][2] = {{	3.0/4.0	,  -3.0*r_3/8.0						},{-3.0*r_3/8.0						,  	3.0/4.0	} }; 
								cout << x_a[0][1] << "x-comp"<<endl; 

				if(a==b) {
				double	a_norm = 1.0/(6.0*M_PI*eta_0*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r
				double	b_norm = 1.0/(6.0*M_PI*eta_0*particle[a].radius*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r2
				double	c_norm = 1.0/(6.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r3
				Mobility_Tnsr_tr	= 		null33D ;
					
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							Mobility_Tnsr_tt.comp[i][j]		=	a_norm*(x_a[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);

							double ep_ijk_e_k = 0.0;
							
							for (int k =0 ; k < 3; k++ ) {ep_ijk_e_k+=Levi_Civi[i][j][k]*e_ab_unit.comp[k];}  
							
							Mobility_Tnsr_rt.comp[i][j]		=	b_norm*(													y_b[1][1]*ep_ijk_e_k													);
						
							Mobility_Tnsr_rr.comp[i][j]		=	c_norm*(x_c[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);

		//						cout << "expresssion"	<< endl;				
		//						cout << Mobility_Tnsr_rr.comp[i][j] << endl;		
						}
					}
	
	/*				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{									
							cout << "matrix"	<< endl;				
							cout << Mobility_Tnsr_rr.comp[i][j] << endl;
						}
					}		*/												
			
				Mobility_Tnsr_tt	=	 	Unit_diag * tau ;
											
				Mobility_Tnsr_rr	=		Unit_diag * temp2 ;
				Mobility_Tnsr_rt	= 		null33D ;
				Mobility_Tnsr_tr	= 		null33D ;
				Mobility_Tnsr_td	= 		null35D ;
				Mobility_Tnsr_dt	= 		null53D ;
							
								
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
						for (int k=0; k<3; k++)
							{
								
							}	// k
						}	//j
					}	//i				

			 } else {
				double	a_norm = 1.0/(6.0*M_PI*eta_0);						// mobility matrix a non-dimensionalized by 6*pi*mu*r
				double	b_norm = 1.0/(6.0*M_PI*eta_0);						// mobility matrix a non-dimensionalized by 6*pi*mu*r2
				double	c_norm = 1.0/(6.0*M_PI*eta_0);						// mobility matrix a non-dimensionalized by 6*pi*mu*r3		
				
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							Mobility_Tnsr_tt.comp[i][j]		=	a_norm*(x_a[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
									


							double ep_ijk_e_k = 0.0;
							
							for (int k =0 ; k < 3; k++ ) {ep_ijk_e_k+=Levi_Civi[i][j][k]*e_ab_unit.comp[k];}  
							
							Mobility_Tnsr_rt.comp[i][j]		=	b_norm*(													y_b[1][0]*ep_ijk_e_k													);
							
							Mobility_Tnsr_tr	= 	    Mobility_Tnsr_rt*(1.0);		
							
							Mobility_Tnsr_rr.comp[i][j]		=	c_norm*(x_c[0][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[0][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);


								cout << "expresssion"	<< endl;				
								cout << Mobility_Tnsr_tt.comp[i][j] << endl;	

						}
					}
									
	/*			Mobility_Tnsr_tt	=	(	Unit_diag
											+	(Pab)*e_ab2_inv
											+	(Unit_diag*(1.0/3.0)-(Pab)*e_ab2_inv)*(particle[a].radius*particle[a].radius+particle[b].radius*particle[b].radius)*e_ab2_inv
											)	*	temp;

									
				Mobility_Tnsr_rr	=		(Unit_diag*(-1.0) + (Pab)*e_ab2_inv*3.0)*temp3 ;
				
				Mobility_Tnsr_rt	=	  	epsilon_e_ab*(2.0)*temp3;
				Mobility_Tnsr_tr	= 	    Mobility_Tnsr_tr*(1.0);	
				
				
					for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{									
							cout << "matrix"	<< endl;				
							cout << Mobility_Tnsr_tt.comp[i][j] << endl;
						}
					}		
	*/		}
			for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{
					/*		mu_6N[l+3*k+9*j+i*9*NrParticles] 								=	 Mobility_Tnsr_tt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+ 9*NrParticles*NrParticles] 	=	 Mobility_Tnsr_tr.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+18*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+27*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];		*/					
							// column major format
							zeta_11N[k+6*NrParticles*l+3*a+18*NrParticles*b] 					=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k+6*NrParticles*l+3*a+18*NrParticles*b+18*NrParticles*NrParticles] 	=	 Mobility_Tnsr_tr.comp[k][l];
							zeta_11N[k+6*NrParticles*l+3*a+18*NrParticles*b+3*NrParticles] 	=	 Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k+6*NrParticles*l+3*a+18*NrParticles*b+18*NrParticles*NrParticles+3*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];							
							zeta_11N[k+6*NrParticles*l+3*b+18*NrParticles*a] 					=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k+6*NrParticles*l+3*b+18*NrParticles*a+18*NrParticles*NrParticles] 	=	 -Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k+6*NrParticles*l+3*b+18*NrParticles*a+3*NrParticles] 	=	 -Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k+6*NrParticles*l+3*b+18*NrParticles*a+18*NrParticles*NrParticles+3*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];
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
   			
	inverse ( zeta_11N , 6*NrParticles )	 ; 							
					
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
							
									Resistance_Tnsr_tt.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j] ;
									Resistance_Tnsr_tr.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles]	;
									Resistance_Tnsr_rt.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+3*NrParticles] ;
									Resistance_Tnsr_rr.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles+3*NrParticles] ;
									
								}
						}
					
					Friction_Tnsr_tt += Resistance_Tnsr_tt ;  
					Friction_Tnsr_tr += ( Resistance_Tnsr_tr    -  ( Resistance_Tnsr_tt*Aj )  ) ;
					//Friction_Tnsr_rt += ( Ai*Resistance_Tnsr_tt +    Resistance_Tnsr_rt    )    ; 
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
					Friction_Tnsr_rt = ~Friction_Tnsr_tr;

									// column major format

					xi_11x11[0] = Friction_Tnsr_tt.comp[0][0] ;  
					xi_11x11[1] = Friction_Tnsr_tt.comp[0][1] ;  
					xi_11x11[2] = Friction_Tnsr_tt.comp[0][2] ; 
					xi_11x11[6] = Friction_Tnsr_tt.comp[1][0] ; 
					xi_11x11[7] = Friction_Tnsr_tt.comp[1][1] ;  
					xi_11x11[8] = Friction_Tnsr_tt.comp[1][2] ;  
					xi_11x11[12] = Friction_Tnsr_tt.comp[2][0] ;   
					xi_11x11[13] = Friction_Tnsr_tt.comp[2][1] ; 
					xi_11x11[14] = Friction_Tnsr_tt.comp[2][2] ; 				

					xi_11x11[18] = Friction_Tnsr_rt.comp[0][0] ;  
					xi_11x11[19] = Friction_Tnsr_rt.comp[0][1] ;  
					xi_11x11[20] = Friction_Tnsr_rt.comp[0][2] ; 
					xi_11x11[24] = Friction_Tnsr_rt.comp[1][0] ; 
					xi_11x11[25] = Friction_Tnsr_rt.comp[1][1] ;  
					xi_11x11[26] = Friction_Tnsr_rt.comp[1][2] ;  
					xi_11x11[30] = Friction_Tnsr_rt.comp[2][0] ;   
					xi_11x11[31] = Friction_Tnsr_rt.comp[2][1] ; 
					xi_11x11[32] = Friction_Tnsr_rt.comp[2][2] ; 				
										
					xi_11x11[3] = Friction_Tnsr_tr.comp[0][0] ;  
					xi_11x11[4] = Friction_Tnsr_tr.comp[0][1] ;  
					xi_11x11[5] = Friction_Tnsr_tr.comp[0][2] ; 
					xi_11x11[9] = Friction_Tnsr_tr.comp[1][0] ; 
					xi_11x11[10] = Friction_Tnsr_tr.comp[1][1] ;  
					xi_11x11[11] = Friction_Tnsr_tr.comp[1][2] ;  
					xi_11x11[15] = Friction_Tnsr_tr.comp[2][0] ;   
					xi_11x11[16] = Friction_Tnsr_tr.comp[2][1] ; 
					xi_11x11[17] = Friction_Tnsr_tr.comp[2][2] ; 				
										
					xi_11x11[21] = Friction_Tnsr_rr.comp[0][0] ;  
					xi_11x11[22] = Friction_Tnsr_rr.comp[0][1] ;  
					xi_11x11[23] = Friction_Tnsr_rr.comp[0][2] ; 
					xi_11x11[27] = Friction_Tnsr_rr.comp[1][0] ; 
					xi_11x11[28] = Friction_Tnsr_rr.comp[1][1] ;  
					xi_11x11[29] = Friction_Tnsr_rr.comp[1][2] ;  
					xi_11x11[33] = Friction_Tnsr_rr.comp[2][0] ;   
					xi_11x11[34] = Friction_Tnsr_rr.comp[2][1] ; 
					xi_11x11[35] = Friction_Tnsr_rr.comp[2][2] ; 				
		
	inverse ( xi_11x11 , 6 )	 ; 			
	for (int i=0; i<36; i++)
		{
			xi_11x11[i]*=4.1419e-14;	// multiply by kbT in erg K-1
		} 	

// center of diffusion calculation based on "Hydrodynamic properties of rigid particles: comparison of different modeling and computational procedures." Biophysical journal 76.6 (1999): 3044-3057.			
// Page 3046 equation 13.
vctr3D ctr_diff ; 
double temp_mat[3*3];
temp_mat[0] =  xi_11x11[28] + xi_11x11[35]	;
temp_mat[1] = -xi_11x11[27] 				; 
temp_mat[2] = -xi_11x11[33]				; 
temp_mat[3] = -xi_11x11[27]				;
temp_mat[4] =  xi_11x11[21] + xi_11x11[35]	;
temp_mat[5] = -xi_11x11[34]				;
temp_mat[6] = -xi_11x11[33]				;
temp_mat[7] = -xi_11x11[34]				;
temp_mat[8] =  xi_11x11[28] + xi_11x11[21]	;

inverse ( temp_mat , 3 )	 ; 	

ctr_diff.comp[0] = temp_mat[0]*(xi_11x11[16] - xi_11x11[11]) + temp_mat[3]*(xi_11x11[5] - xi_11x11[15]) + temp_mat[6]*(xi_11x11[9] - xi_11x11[4]) ;
ctr_diff.comp[1] = temp_mat[1]*(xi_11x11[16] - xi_11x11[11]) + temp_mat[4]*(xi_11x11[5] - xi_11x11[15]) + temp_mat[7]*(xi_11x11[9] - xi_11x11[4]) ;
ctr_diff.comp[2] = temp_mat[2]*(xi_11x11[16] - xi_11x11[11]) + temp_mat[5]*(xi_11x11[5] - xi_11x11[15]) + temp_mat[8]*(xi_11x11[9] - xi_11x11[4]) ;

mtrx3D D_tt, D_tr, D_rt , D_rr, U_OD;

D_tt.comp[0][0] = xi_11x11[0];
D_tt.comp[1][0] = xi_11x11[1];
D_tt.comp[2][0] = xi_11x11[2];
D_tt.comp[0][1] = xi_11x11[6];
D_tt.comp[1][1] = xi_11x11[7];
D_tt.comp[2][1] = xi_11x11[8];
D_tt.comp[0][2] = xi_11x11[12];
D_tt.comp[1][2] = xi_11x11[13];
D_tt.comp[2][2] = xi_11x11[14];

D_tr.comp[0][0] = xi_11x11[3];
D_tr.comp[1][0] = xi_11x11[4];
D_tr.comp[2][0] = xi_11x11[5];
D_tr.comp[0][1] = xi_11x11[9];
D_tr.comp[1][1] = xi_11x11[10];
D_tr.comp[2][1] = xi_11x11[11];
D_tr.comp[0][2] = xi_11x11[15];
D_tr.comp[1][2] = xi_11x11[16];
D_tr.comp[2][2] = xi_11x11[17];

D_rt.comp[0][0] = xi_11x11[18];
D_rt.comp[1][0] = xi_11x11[19];
D_rt.comp[2][0] = xi_11x11[20];
D_rt.comp[0][1] = xi_11x11[24];
D_rt.comp[1][1] = xi_11x11[25];
D_rt.comp[2][1] = xi_11x11[26];
D_rt.comp[0][2] = xi_11x11[30];
D_rt.comp[1][2] = xi_11x11[31];
D_rt.comp[2][2] = xi_11x11[32];

D_rr.comp[0][0] = xi_11x11[21];
D_rr.comp[1][0] = xi_11x11[22];
D_rr.comp[2][0] = xi_11x11[23];
D_rr.comp[0][1] = xi_11x11[27];
D_rr.comp[1][1] = xi_11x11[28];
D_rr.comp[2][1] = xi_11x11[29];
D_rr.comp[0][2] = xi_11x11[33];
D_rr.comp[1][2] = xi_11x11[34];
D_rr.comp[2][2] = xi_11x11[35];

U_OD.comp[0][0] =  0.0;
U_OD.comp[1][0] =  ctr_diff.comp[2];
U_OD.comp[2][0] = -ctr_diff.comp[1];
U_OD.comp[0][1] = -ctr_diff.comp[2];
U_OD.comp[1][1] =  0.0;
U_OD.comp[2][1] =  ctr_diff.comp[0];
U_OD.comp[0][2] =  ctr_diff.comp[1];
U_OD.comp[1][2] = -ctr_diff.comp[0];
U_OD.comp[2][2] =  0.0;

mtrx3D D_tt_CoD = D_tt -  U_OD*D_rr*U_OD + D_rt*U_OD - U_OD*D_tr ; 
mtrx3D D_tr_CoD = D_tr +  D_rr*U_OD ;  // based on equations 42 from Wouter's notes "clusterdyn"
mtrx3D D_rt_CoD = D_rt -  U_OD*D_rr ;  // based on equations 43 from Wouter's notes "clusterdyn"

xi_11x11[0]  = D_tt_CoD.comp[0][0];
xi_11x11[1]  = D_tt_CoD.comp[1][0];
xi_11x11[2]  = D_tt_CoD.comp[2][0];
xi_11x11[6]  = D_tt_CoD.comp[0][1];
xi_11x11[7]  = D_tt_CoD.comp[1][1];
xi_11x11[8]  = D_tt_CoD.comp[2][1];
xi_11x11[12] = D_tt_CoD.comp[0][2];
xi_11x11[13] = D_tt_CoD.comp[1][2];
xi_11x11[14] = D_tt_CoD.comp[2][2];

xi_11x11[3]  = D_tr_CoD.comp[0][0];
xi_11x11[4]  = D_tr_CoD.comp[1][0];
xi_11x11[5]  = D_tr_CoD.comp[2][0];
xi_11x11[9]  = D_tr_CoD.comp[0][1];
xi_11x11[10] = D_tr_CoD.comp[1][1];
xi_11x11[11] = D_tr_CoD.comp[2][1];
xi_11x11[15] = D_tr_CoD.comp[0][2];
xi_11x11[16] = D_tr_CoD.comp[1][2];
xi_11x11[17] = D_tr_CoD.comp[2][2];

xi_11x11[18]  = D_rt_CoD.comp[0][0];
xi_11x11[19]  = D_rt_CoD.comp[1][0];
xi_11x11[20]  = D_rt_CoD.comp[2][0];
xi_11x11[24]  = D_rt_CoD.comp[0][1];
xi_11x11[25]  = D_rt_CoD.comp[1][1];
xi_11x11[26]  = D_rt_CoD.comp[2][1];
xi_11x11[30]  = D_rt_CoD.comp[0][2];
xi_11x11[31]  = D_rt_CoD.comp[1][2];
xi_11x11[32]  = D_rt_CoD.comp[2][2];

	for (int i=0; i<NrParticles; i++)
		{
			outFile1<<"particle position"<<'\t'<<particle[i].pos.comp[0]<<'\t'<<particle[i].pos.comp[1]<<'\t'<<particle[i].pos.comp[2]<<'\t'<<particle[i].radius<<std::endl ;
		} 	
		outFile1<<std::endl ;
		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[6]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[18]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[30]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[7]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[19]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[31]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[8]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[20]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[32]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[9]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[21]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[33]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[10]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[28]<<'\t'<<xi_11x11[34]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[17]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[29]<<'\t'<<xi_11x11[35]<<std::endl ;
		
		outFile1<<ctr_diff.comp[0]<<'\t'<<ctr_diff.comp[1]<<'\t'<<ctr_diff.comp[2]<<std::endl ;
		
		D_tt_CoD.echo(); 


}
		
