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

	mtrx3D Resistance_Tnsr_tt;
	mtrx3D Resistance_Tnsr_tr;
	mtrx3D Resistance_Tnsr_rt;
	mtrx3D Resistance_Tnsr_rr;	
	mtrx35D Resistance_Tnsr_rd;
	mtrx35D Resistance_Tnsr_td;
	mtrx53D Resistance_Tnsr_dt;
	mtrx53D Resistance_Tnsr_dr;
	mtrx55D Resistance_Tnsr_dd;
	
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

// three and four index resistance matrices	
	double G_IJK[3][3][3];
	double H_IJK[3][3][3];
	double M_IJKL[3][3][3][3];
	double H_clst_ijk[3][3][3] = {{{0}}};
		
// the five base matrices for strain tensor // option 5 :  equation 419 wouter's tex version clusterdyn_110816_1556

	double e[5][3][3]= {
							{{1.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,-1.0}},
							{{0.0,1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,1.0},{0.0,0.0,0.0},{1.0,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,1.0},{0.0,1.0,0.0}},
							{{0.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,-1.0}}
						};


	double e_l[5][3][3]= {
							{{ 2.0/3.0,0.0,0.0},{0.0,-1.0/3.0,0.0},{0.0,0.0,-1.0/3.0}},
							{{0.0,0.5,0.0},{0.5,0.0,0.0},{0.0,0.0,0.0}},
							{{0.0,0.0,0.5},{0.0,0.0,0.0},{0.5,0.0,0.0}},
							{{0.0,0.0,0.0},{0.0,0.0,0.5},{0.0,0.5,0.0}},
							{{-1.0/3.0,0.0,0.0},{0.0, 2.0/3.0,0.0},{0.0,0.0,-1.0/3.0}}
						};
   
   double mu_11N[121*NrParticles*NrParticles] ;  		// grand mobility matrix
   double zeta_11N[121*NrParticles*NrParticles] ={0.0} ;  	// grand resistance matrix
   double rho_11N[121*NrParticles*NrParticles] ;  	// grand resistance matrix
   double xi_11x11[11*11] ; 							// generalized friction matrix
   double xi_6x6[6*6] ; 							// generalized friction matrix
   
   for (int i=0; i<121; i++)
		{
			xi_11x11[i] = 0.0; 
		}        
		
std::ofstream outFile1(dataFileName+"/data.dat");

// important all lengths have been normalized by particle radius as metioned in Page 46, Appendix A - Durlofsky, Louis, John F. Brady, and Georges Bossis. 
				// "Dynamic simulation of hydrodynamically interacting particles." Journal of fluid mechanics 180 (1987): 21-49.
				// for ease of programming. 
		
for (int a=0; a<NrParticles; a++)
	{
		for (int b=a; b<NrParticles; b++)
			{
				cout<<a<<b<<endl;
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
			    double r 	= sqrt(e_ab2)/particle[a].radius;			// distance between particle vector 'r' magnitude |r| normalized by particle radius 'a' ;
			    double r_1 	= 1.0/(r);
			    double r_2 	= 1.0/(r*r);			    
			    double r_3 	= 1.0/(r*r*r);
			    double r_4 	= 1.0/(r*r*r*r);
			    double r_5 	= 1.0/(r*r*r*r*r);
			    // cout << r << '\t'<<r_1<< '\t'<<r_2 << '\t' << r_3 << 'a'<< a<< 'b'<< b<< endl;
				// cout << e_ab_unit.comp[0] << "x-comp"<<endl; 

				// mobility scalar values - as defined in Page 46, Appendix A - Durlofsky, Louis, John F. Brady, and Georges Bossis. 
				// "Dynamic simulation of hydrodynamically interacting particles." Journal of fluid mechanics 180 (1987): 21-49.

				double x_a[2][2] = {{	1.0		,	3.0*r_1/2.0		-	1.0*r_3			},{	3.0*r_1/2.0		-	1.0*r_3			,	1.0		} }; 
				double y_a[2][2] = {{	1.0		,	(3.0*r_1/4.0)	+	(1.0*r_3/2.0)	},{	3.0*r_1/4.0		+	1.0*r_3/2.0		,	1.0		} }; 
				double y_b[2][2] = {{	0.0		,  -3.0*r_2/4.0							},{	3.0*r_2/4.0							,	0.0		} }; 
				double x_c[2][2] = {{	3.0/4.0	,  	3.0*r_3/4.0							},{	3.0*r_3/4.0							,  	3.0/4.0	} }; 
				double y_c[2][2] = {{	3.0/4.0	,  -3.0*r_3/8.0							},{-3.0*r_3/8.0							,  	3.0/4.0	} }; 

				double x_g[2][2] = {{	0.0		,	9.0*r_2/4.0		-	18.0*r_4/5.0	},{-9.0*r_2/4.0		+	18.0*r_4/5.0	,	0.0		} };
				double y_g[2][2] = {{	0.0		,	6.0*r_4/5.0							},{-6.0*r_4/5.0							,	0.0		} };
				double y_h[2][2] = {{	0.0		,  -9.0*r_3/8.0							},{-9.0*r_3/8.0							,	0.0		} };
				double x_m[2][2] = {{	9.0/10.0,  -9.0*r_3/2.0		+ 	54.0*r_5/5.0	},{-9.0*r_3/2.0		+ 	54.0*r_5/5.0	,	9.0/10.0} };
				double y_m[2][2] = {{	9.0/10.0,   9.0*r_3/4.0		- 	36.0*r_5/5.0	},{ 9.0*r_3/4.0		- 	36.0*r_5/5.0	,	9.0/10.0} };
				double z_m[2][2] = {{	9.0/10.0,  					 	 9.0*r_5/5.0	},{ 				 	 9.0*r_5/5.0	,	9.0/10.0} };

				// resistance scalar values - as defined in Page 280, book of Kim and Karrila

				double X_A[2][2] = {{	1.0		,  -3.0*r_1/2.0		-	19.0*r_3/8.0	},{-3.0*r_1/2.0		-	19.0*r_3/8.0	,	1.0		} }; 	// normalised by 6*pi*r
				double Y_A[2][2] = {{	1.0		,  -3.0*r_1/4.0		-	59.0*r_3/16.0	},{-3.0*r_1/4.0		-	59.0*r_3/16.0	,	1.0		} }; 
				double Y_B[2][2] = {{	0.0		,  -r_2									},{-r_2									,	0.0		} }; 
				double X_C[2][2] = {{	4.0/3.0	,  -4.0*r_3/3.0							},{-4.0*r_3/3.0							,  	4.0/3.0	} }; 
				double Y_C[2][2] = {{	4.0/3.0	,   2.0*r_3/3.0							},{ 2.0*r_3/3.0							,  	4.0/3.0	} }; 


				double X_G[2][2] = {{	0.0		,  -15.0*r_2/4.0	/*-	39.0*r_4/16.0*/	},{0.0/*-15.0*r_2/4.0	/*-	39.0*r_4/16.0*/	,	0.0		} };
				double Y_G[2][2] = {{	0.0		,  -2.0*r_4								},{-2.0*r_4								,	1.0		} };
				double Y_H[2][2] = {{	0.0		,   5.0*r_3/3.0							},{10.0*r_3/3.0							,	0.0		} };
				double X_M[2][2] = {{	1.0		,   5.0*r_3			/*- 	51.0*r_5/8.0*/	},{ 5.0*r_3			/*- 	51.0*r_5/8.0*/	,	1.0		} };
				double Y_M[2][2] = {{	1.0		,  -5.0*r_3/2.0		/*+ 	8.0*r_5*/			},{-5.0*r_3/2.0		/*+	8.0*r_5*/			,	1.0		} };
				double Z_M[2][2] = {{	1.0		,  					0.0/*-	2.0*r_5*/			},{ 				0.0/*-	 2.0*r_5*/		,	1.0		} };
				
				//				cout << x_a[0][1] << "x-comp"<<endl; 
				
				double	a_norm = 1.0; // /(6.0*M_PI*eta_0*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r
				double	b_norm = 1.0; // /(6.0*M_PI*eta_0*particle[a].radius*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r2
				double	c_norm = 1.0; // /(6.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r3
				double	g_norm = 1.0; // /(6.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius);						//		assuming correction factor of 6*pi*mu*r3	
				double	h_norm = 1.0; // /(6.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius);						//		assuming correction factor of 6*pi*mu*r3				
			//	double	g_norm = 20.0/18.0;											//		assuming correction factor of 6*pi*mu*r3/(20*pi*mu*r3/3)	
			//	double	h_norm = 20.0/18.0;											//		assuming correction factor of 6*pi*mu*r3/(20*pi*mu*r3/3)	
				double	m_norm = 1.0; // /(6.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius);						//		assuming correction factor of 6*pi*mu*r3	

				double	A_Norm = 1.0; // (6.0*M_PI*eta_0*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r
				double	B_Norm = 1.0; // (6.0*M_PI*eta_0*particle[a].radius*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r2
				double	C_Norm = 1.0; // (6.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius);						// mobility matrix a non-dimensionalized by 6*pi*mu*r3
				double	G_Norm = 1.0; // (4.0*M_PI*eta_0*particle[a].radius*particle[a].radius);						//		assuming correction factor of 4*pi*mu*r2	
				double	H_Norm = 1.0; // (6.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius);						//		assuming correction factor of 6*pi*mu*r3				
			//	double	g_norm = 20.0/18.0;											//		assuming correction factor of 6*pi*mu*r3/(20*pi*mu*r3/3)	
			//	double	h_norm = 20.0/18.0;											//		assuming correction factor of 6*pi*mu*r3/(20*pi*mu*r3/3)	
				double	M_Norm = 1.0 ; // (20.0*M_PI*eta_0*particle[a].radius*particle[a].radius*particle[a].radius/3.0);						//		assuming correction factor of 20*pi*mu*r3/3	
				
				if(a==b) {
				cout<<a<<b<<endl;

				Mobility_Tnsr_tr	= 		null33D ;
					
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							double ep_ijk_e_k = 0.0;
							
						for (int k=0; k<3; k++)
							{	

								double ep_jkl_e_l = 0.0;
								double ep_ikl_e_l = 0.0;
								
								for (int l=0; l<3; l++)
									{
										ep_jkl_e_l	+=	Levi_Civi[j][k][l]*e_ab_unit.comp[l];
										ep_ikl_e_l	+=	Levi_Civi[i][k][l]*e_ab_unit.comp[l];
										
										m_ijkl[i][j][k][l]	=	m_norm*((3.0/2.0)*x_m[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 					-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																-(1.0/3.0)*kron_del[k][l])
																+(1.0/2.0)*y_m[1][1]*(e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	+	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																					+ e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	+ 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																					- 4.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	)
																			
																+(1.0/2.0)*z_m[1][1]*(kron_del[i][k]*kron_del[j][l]		+ 	kron_del[j][k]*kron_del[i][l]	- 	kron_del[i][j]*kron_del[k][l] 
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*kron_del[k][l]	+	kron_del[i][j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]
																- e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	- 	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																- e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	- 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																));
																																							
									}	// l
									
								ep_ijk_e_k					+=	Levi_Civi[i][j][k]*e_ab_unit.comp[k];
								
								g_ijk[i][j][k]				=	g_norm*(x_g[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		y_g[1][1]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								h_ijk[i][j][k]				= 	h_norm*(y_h[1][1]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Mobility_Tnsr_tt.comp[i][j]		=	a_norm*(x_a[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Mobility_Tnsr_rt.comp[i][j]		=	b_norm*(													y_b[1][1]*ep_ijk_e_k													);
						
							Mobility_Tnsr_rr.comp[i][j]		=	c_norm*(x_c[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);

						//		cout << "m_ijkl"	<< endl;				
						//		cout << m_ijkl[0][0][0][0] << endl;				
						}	// j
					}	// i
		/*
					for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{									
							cout << "matrix"	<< endl;				
							cout << Mobility_Tnsr_rt.comp[i][j] << endl;
						}
					}														
			
				Mobility_Tnsr_tt	=	 	Unit_diag * tau ;
											
				Mobility_Tnsr_rr	=		Unit_diag * temp2 ;
				Mobility_Tnsr_rt	= 		null33D ;
				Mobility_Tnsr_tr	= 		null33D ;
		*/

			 } else {
				 					
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							double ep_ijk_e_k = 0.0;
							
						for (int k=0; k<3; k++)
							{	

								double ep_jkl_e_l = 0.0;
								double ep_ikl_e_l = 0.0;
								
								for (int l=0; l<3; l++)
									{
										ep_jkl_e_l	+=	Levi_Civi[j][k][l]*e_ab_unit.comp[l];
										ep_ikl_e_l	+=	Levi_Civi[i][k][l]*e_ab_unit.comp[l];
										
										m_ijkl[i][j][k][l]	=	 m_norm*((3.0/2.0)*x_m[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 			-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																-(1.0/3.0)*kron_del[k][l])
																+(1.0/2.0)*y_m[1][0]*(e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	+	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																					+ e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	+ 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																					- 4.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	)
																			
																+(1.0/2.0)*z_m[1][0]*(kron_del[i][k]*kron_del[j][l]		+ 	kron_del[j][k]*kron_del[i][l]	- 	kron_del[i][j]*kron_del[k][l]
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*kron_del[k][l]	+	kron_del[i][j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]
																- e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	- 	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																- e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	- 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																));
																																							
									}	// l
									
								ep_ijk_e_k					+=	Levi_Civi[i][j][k]*e_ab_unit.comp[k];
								
								g_ijk[i][j][k]				=	g_norm*(x_g[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		y_g[1][0]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								h_ijk[i][j][k]				= 	h_norm*(y_h[1][0]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Mobility_Tnsr_tt.comp[i][j]		=	a_norm*(x_a[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_a[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Mobility_Tnsr_rt.comp[i][j]		=	b_norm*(													y_b[1][0]*ep_ijk_e_k													);
						
							Mobility_Tnsr_rr.comp[i][j]		=	c_norm*(x_c[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	y_c[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
									
							Mobility_Tnsr_tr	= 	    Mobility_Tnsr_rt*(1.0);		
		//						cout << "expresssion"	<< endl;				
		//						cout << Mobility_Tnsr_rt.comp[i][j] << endl;		
						}	// j
					}	// i	
				}				

// Resistance matrix calculation

				if(a==b) {

				Resistance_Tnsr_tr	= 		null33D ;
					
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							double ep_ijk_e_k = 0.0;
							
						for (int k=0; k<3; k++)
							{	

								double ep_jkl_e_l = 0.0;
								double ep_ikl_e_l = 0.0;
								
								for (int l=0; l<3; l++)
									{
										ep_jkl_e_l	+=	Levi_Civi[j][k][l]*e_ab_unit.comp[l];
										ep_ikl_e_l	+=	Levi_Civi[i][k][l]*e_ab_unit.comp[l];
										
										M_IJKL[i][j][k][l]	=	M_Norm*((3.0/2.0)*X_M[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 					-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																-(1.0/3.0)*kron_del[k][l])
																+(1.0/2.0)*Y_M[1][1]*(e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	+	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																					+ e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	+ 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																					- 4.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	)
																			
																+(1.0/2.0)*Z_M[1][1]*(kron_del[i][k]*kron_del[j][l]		+ 	kron_del[j][k]*kron_del[i][l]	- 	kron_del[i][j]*kron_del[k][l]
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*kron_del[k][l]	+	kron_del[i][j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]
																- e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	- 	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																- e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	- 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																));

							//	outFile1 << M_IJKL[i][j][k][l] <<	'\t'	<< m_ijkl[i][j][k][l] << '\t'	<< "m_ijkl"  << i <<j <<k<< l<< endl;																																									
																																																
									}	// l
									
								ep_ijk_e_k					+=	Levi_Civi[i][j][k]*e_ab_unit.comp[k];
								
								G_IJK[i][j][k]				=	G_Norm*(X_G[1][1]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		Y_G[1][1]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								H_IJK[i][j][k]				= 	H_Norm*(Y_H[1][1]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Resistance_Tnsr_tt.comp[i][j]		=	A_Norm*(X_A[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	Y_A[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
						//	Resistance_Tnsr_rt.comp[i][j]		=	B_Norm*(													Y_B[1][1]*ep_ijk_e_k													);
						
							Resistance_Tnsr_rr.comp[i][j]		=	C_Norm*(X_C[1][1]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	Y_C[1][1]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Resistance_Tnsr_rt	= 		null33D ;

						//		cout << "m_ijkl"	<< endl;				
						//		cout << M_IJKL[0][0][0][0] << endl;	
						}	// j
					}	// i
		/*
					for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{									
							cout << "matrix"	<< endl;				
							cout << Mobility_Tnsr_rt.comp[i][j] << endl;
						}
					}														
			
				Mobility_Tnsr_tt	=	 	Unit_diag * tau ;
											
				Mobility_Tnsr_rr	=		Unit_diag * temp2 ;
				Mobility_Tnsr_rt	= 		null33D ;
				Mobility_Tnsr_tr	= 		null33D ;
		*/

			 } else {
				 					
				for (int i=0; i<3; i++)
					{
					for (int j=0; j<3; j++)
						{
							double ep_ijk_e_k = 0.0;
							
						for (int k=0; k<3; k++)
							{	

								double ep_jkl_e_l = 0.0;
								double ep_ikl_e_l = 0.0;
								
								for (int l=0; l<3; l++)
									{
										ep_jkl_e_l	+=	Levi_Civi[j][k][l]*e_ab_unit.comp[l];
										ep_ikl_e_l	+=	Levi_Civi[i][k][l]*e_ab_unit.comp[l];
										
										M_IJKL[i][j][k][l]	=	 M_Norm*((3.0/2.0)*X_M[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 			-	(1.0/3.0)*kron_del[i][j])*(e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																-(1.0/3.0)*kron_del[k][l])
																+(1.0/2.0)*Y_M[1][0]*(e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	+	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																					+ e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	+ 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																					- 4.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	)
																			
																+(1.0/2.0)*Z_M[1][0]*(kron_del[i][k]*kron_del[j][l]		+ 	kron_del[j][k]*kron_del[i][l]	- 	kron_del[i][j]*kron_del[k][l]
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*kron_del[k][l]	+	kron_del[i][j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]	
																+ e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]*e_ab_unit.comp[l]
																- e_ab_unit.comp[i]*kron_del[j][l]*e_ab_unit.comp[k]	- 	e_ab_unit.comp[j]*kron_del[i][l]*e_ab_unit.comp[k]
																- e_ab_unit.comp[i]*kron_del[j][k]*e_ab_unit.comp[l]	- 	e_ab_unit.comp[j]*kron_del[i][k]*e_ab_unit.comp[l]
																));
																																							
									}	// l
									
								ep_ijk_e_k					+=	Levi_Civi[i][j][k]*e_ab_unit.comp[k];
								
								G_IJK[i][j][k]				=	G_Norm*(X_G[1][0]*(e_ab_unit.comp[i]*e_ab_unit.comp[j] 	-	(1.0/3.0)*kron_del[i][j])*e_ab_unit.comp[k]
																+ 		Y_G[1][0]*(e_ab_unit.comp[i]*kron_del[j][k]		+ 	e_ab_unit.comp[j]*kron_del[i][k]	-	2.0*e_ab_unit.comp[i]*e_ab_unit.comp[j]*e_ab_unit.comp[k]	)	);
										
								H_IJK[i][j][k]				= 	H_Norm*(Y_H[1][0]*(e_ab_unit.comp[i]*ep_jkl_e_l			+	e_ab_unit.comp[j]*ep_ikl_e_l										)	);
							}	// k		

							Resistance_Tnsr_tt.comp[i][j]		=	A_Norm*(X_A[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	Y_A[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
							
							Resistance_Tnsr_rt.comp[i][j]		=	B_Norm*(													Y_B[1][0]*ep_ijk_e_k													);
						
							Resistance_Tnsr_rr.comp[i][j]		=	C_Norm*(X_C[1][0]*e_ab_unit.comp[i]*e_ab_unit.comp[j]	+ 	Y_C[1][0]*(kron_del[i][j]	- e_ab_unit.comp[i]*e_ab_unit.comp[j]	)	);
									
							Resistance_Tnsr_tr	= 	    Resistance_Tnsr_rt*(1.0);		
		//						cout << "expresssion"	<< endl;				
		//						cout << Mobility_Tnsr_rt.comp[i][j] << endl;		
						}	// j
					}	// i	

Resistance_Tnsr_tt.echo();
									
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
							cout << Mobility_Tnsr_tr.comp[i][j] << endl;
						}
					}		
	*/		}
	
//	extract the reduced index mobility tesnors from g_ijk, h_ijk and m_ijkl

// based on the equation mu^dt_{\p\g} = (e^p)_{\a\b}*mu^dt_{\a\b\g}	, etc in wouter notes equation no : 422 , version : 110816_1556
	
	Mobility_Tnsr_dt	= 		null53D ;	
	Mobility_Tnsr_dr	= 		null53D ;	
	Mobility_Tnsr_dd	= 		null55D ;
		
	Resistance_Tnsr_dt	= 		null53D ;	
	Resistance_Tnsr_dr	= 		null53D ;	
	Resistance_Tnsr_dd	= 		null55D ;
	
	Mobility_Tnsr_td	= 		null35D ;	
	Mobility_Tnsr_rd	= 		null35D ;	
		
	Resistance_Tnsr_td	= 		null35D ;	
	Resistance_Tnsr_rd	= 		null35D ;	

outFile1<<"g_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		outFile1 << setw(10) << G_IJK[s][0][0] << "  " << setw(10) << G_IJK[s][0][1] << "  " << setw(10) << G_IJK[s][0][2] << endl;
		outFile1 << setw(10) << G_IJK[s][1][0] << "  " << setw(10) << G_IJK[s][1][1] << "  " << setw(10) << G_IJK[s][1][2] << endl;
		outFile1 << setw(10) << G_IJK[s][2][0] << "  " << setw(10) << G_IJK[s][2][1] << "  " << setw(10) << G_IJK[s][2][2] << endl;
    }

outFile1<<"h_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		outFile1 << setw(10) << H_IJK[s][0][0] << "  " << setw(10) << H_IJK[s][0][1] << "  " << setw(10) << H_IJK[s][0][2] << endl;
		outFile1 << setw(10) << H_IJK[s][1][0] << "  " << setw(10) << H_IJK[s][1][1] << "  " << setw(10) << H_IJK[s][1][2] << endl;
		outFile1 << setw(10) << H_IJK[s][2][0] << "  " << setw(10) << H_IJK[s][2][1] << "  " << setw(10) << H_IJK[s][2][2] << endl;
    }
	
	for (int p=0; p<5; p++)
		{
		for (int g=0; g<3; g++)
			{
			for (int a=0; a<3; a++)
				{
				for (int b=0; b<3; b++)
					{		
						Mobility_Tnsr_dt.comp[p][g]		+=		e[p][a][b]*g_ijk[a][b][g];	
						Mobility_Tnsr_dr.comp[p][g]		+=		e[p][a][b]*h_ijk[a][b][g];		
						Mobility_Tnsr_td.comp[g][p]		+=		g_ijk[g][a][b]*e[p][a][b];	
						Mobility_Tnsr_rd.comp[g][p]		+=		h_ijk[g][a][b]*e[p][a][b];
								
						Resistance_Tnsr_dt.comp[p][g]		+=		e_l[p][a][b]*G_IJK[a][b][g];	
						Resistance_Tnsr_dr.comp[p][g]		+=		e_l[p][a][b]*H_IJK[a][b][g];		
						Resistance_Tnsr_td.comp[g][p]		+=		e_l[p][a][b]*G_IJK[a][b][g];
						Resistance_Tnsr_rd.comp[g][p]		+=		e_l[p][a][b]*H_IJK[a][b][g];
					}
				}				
			}
		for (int s=0; s<5; s++)
			{
			for (int a=0; a<3; a++)
				{
				for (int b=0; b<3; b++)
					{
					for (int g=0; g<3; g++)
						{
						for (int d=0; d<3; d++)
							{							
								Mobility_Tnsr_dd.comp[p][s]		+=		e[p][a][b]*m_ijkl[a][b][g][d]*e[s][g][d];			
								Resistance_Tnsr_dd.comp[p][s]		+=		e_l[p][a][b]*M_IJKL[a][b][g][d]*e_l[s][g][d];			
							}
						}													
					}
				}
			}
		}		
	
				for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{
					/*		mu_6N[l+3*k+9*j+i*9*NrParticles] 								=	 Mobility_Tnsr_tt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+ 9*NrParticles*NrParticles] 	=	 Mobility_Tnsr_tr.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+18*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+27*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];		*/					
							// 11N column major format
							rho_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b														] 	=	 Resistance_Tnsr_tt.comp[k][l];
							rho_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles						] 	=	 Resistance_Tnsr_tr.comp[k][l];
							rho_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	3*NrParticles									] 	=	 Resistance_Tnsr_rt.comp[k][l];
							rho_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Resistance_Tnsr_rr.comp[k][l];							
							rho_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a														] 	=	 Resistance_Tnsr_tt.comp[k][l];
							rho_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles						] 	=	-Resistance_Tnsr_rt.comp[k][l];
							rho_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	3*NrParticles									] 	=	-Resistance_Tnsr_rt.comp[k][l];
							rho_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Resistance_Tnsr_rr.comp[k][l];
							// 6N column major format
/*							rho_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b														] 	=	 Resistance_Tnsr_tt.comp[k][l];
							rho_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	18*NrParticles*NrParticles						] 	=	 Resistance_Tnsr_tr.comp[k][l];
							rho_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	3*NrParticles									] 	=	 Resistance_Tnsr_rt.comp[k][l];
							rho_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	18*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Resistance_Tnsr_rr.comp[k][l];							
							rho_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a														] 	=	 Resistance_Tnsr_tt.comp[k][l];
							rho_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	18*NrParticles*NrParticles						] 	=	 -Resistance_Tnsr_rt.comp[k][l];
							rho_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	3*NrParticles									] 	=	 -Resistance_Tnsr_rt.comp[k][l];
							rho_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	18*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Resistance_Tnsr_rr.comp[k][l];
*/						}
				}

				
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<3; k++)
					{				
						// column major format
						rho_11N[k	+	11*NrParticles*l	+	3*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles						] 	=	- Resistance_Tnsr_td.comp[k][l];
						rho_11N[k	+	11*NrParticles*l	+	3*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Resistance_Tnsr_rd.comp[k][l];
						rho_11N[k	+	11*NrParticles*l	+	3*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles						] 	=	 Resistance_Tnsr_td.comp[k][l];
						rho_11N[k	+	11*NrParticles*l	+	3*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Resistance_Tnsr_rd.comp[k][l];
						
					}
				}					
			for (int l=0; l<3; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						rho_11N[k	+	11*NrParticles*l	+	5*a	+	33*NrParticles*b	+	6*NrParticles									] 	=	 Resistance_Tnsr_td.comp[l][k];
						rho_11N[k	+	11*NrParticles*l	+	5*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Resistance_Tnsr_rd.comp[l][k];
						rho_11N[k	+	11*NrParticles*l	+	5*b	+	33*NrParticles*a	+	6*NrParticles									] 	=	- Resistance_Tnsr_td.comp[l][k];
						rho_11N[k	+	11*NrParticles*l	+	5*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Resistance_Tnsr_rd.comp[l][k];						
					}
				}
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						rho_11N[k	+	11*NrParticles*l	+	5*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Resistance_Tnsr_dd.comp[k][l];
						rho_11N[k	+	11*NrParticles*l	+	5*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Resistance_Tnsr_dd.comp[k][l];
					}
				}


				
			for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{
					/*		mu_6N[l+3*k+9*j+i*9*NrParticles] 								=	 Mobility_Tnsr_tt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+ 9*NrParticles*NrParticles] 	=	 Mobility_Tnsr_tr.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+18*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rt.comp[k][l];
							mu_6N[l+3*k+9*j+i*9*NrParticles+27*NrParticles*NrParticles] 	=	 Mobility_Tnsr_rr.comp[k][l];		*/					
							// 11N column major format
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles						] 	=	 Mobility_Tnsr_tr.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	3*NrParticles									] 	=	 Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];							
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles						] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	3*NrParticles									] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	11*NrParticles*l	+	3*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];
							// 6N column major format
/*							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	18*NrParticles*NrParticles						] 	=	 Mobility_Tnsr_tr.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	3*NrParticles									] 	=	 Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*a	+	18*NrParticles*b	+	18*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];							
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a														] 	=	 Mobility_Tnsr_tt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	18*NrParticles*NrParticles						] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	3*NrParticles									] 	=	-Mobility_Tnsr_rt.comp[k][l];
							zeta_11N[k	+	6*NrParticles*l	+	3*b	+	18*NrParticles*a	+	18*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rr.comp[k][l];
*/					}
				}
				
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<3; k++)
					{				
						// column major format
						zeta_11N[k	+	11*NrParticles*l	+	3*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles						] 	=	 Mobility_Tnsr_td.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	3*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	3*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles						] 	=	 Mobility_Tnsr_td.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	3*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles	+	3*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[k][l];
						
					}
				}					
			for (int l=0; l<3; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						zeta_11N[k	+	11*NrParticles*l	+	5*a	+	33*NrParticles*b	+	6*NrParticles									] 	=	 Mobility_Tnsr_td.comp[l][k];
						zeta_11N[k	+	11*NrParticles*l	+	5*a	+	33*NrParticles*b	+	33*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[l][k];
						zeta_11N[k	+	11*NrParticles*l	+	5*b	+	33*NrParticles*a	+	6*NrParticles									] 	=	 Mobility_Tnsr_td.comp[l][k];
						zeta_11N[k	+	11*NrParticles*l	+	5*b	+	33*NrParticles*a	+	33*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_rd.comp[l][k];						
					}
				}
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						zeta_11N[k	+	11*NrParticles*l	+	5*a	+	55*NrParticles*b	+	66*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_dd.comp[k][l];
						zeta_11N[k	+	11*NrParticles*l	+	5*b	+	55*NrParticles*a	+	66*NrParticles*NrParticles	+	6*NrParticles	] 	=	 Mobility_Tnsr_dd.comp[k][l];
					}
				}		
			}	

		}
		
		
	//	cout<<zeta_11N[72]<<'\t'<<zeta_11N[84]<<'\t'<<zeta_11N[96]<<'\t'<<zeta_11N[109]<<'\t'<<zeta_11N[120]<<std::endl ;
	//	cout<<m_ijkl[0][0][0][0]<<'\t'<<m_ijkl[1][1][1][1]<<'\t'<<m_ijkl[2][2][2][2]<<'\t'<<zeta_11N[109]<<'\t'<<zeta_11N[120]<<std::endl ;

		
	mtrx3D Friction_Tnsr_tt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_tr(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rt(0.0,0.0,0.0);
	mtrx3D Friction_Tnsr_rr(0.0,0.0,0.0);
 	mtrx53D Friction_Tnsr_dt 	=	null53D;
	mtrx53D Friction_Tnsr_dr	= 	null53D;
	mtrx35D Friction_Tnsr_td	=	null35D;
	mtrx35D Friction_Tnsr_rd	=	null35D;
	mtrx55D Friction_Tnsr_dd	=	null55D;

	      for (int i=0;i<22;i++)
       { 
               for(int j=0;j<22;j++)
                       {
               outFile1<< std::setprecision(5) <<zeta_11N[i+11*NrParticles*j]<<'\t' ;
                       }
               outFile1<<std::endl; 
       }	
               outFile1<<std::endl; 
			  			
	inverse ( zeta_11N ,11*NrParticles )	 ; 	
               cout<<std::endl; 


	for (int i=0; i<11*NrParticles*11*NrParticles; i++)
		{
			zeta_11N[i] = rho_11N[i] ;
		}	
	
	      for (int i=0;i<22;i++)
       { 
               for(int j=0;j<22;j++)
                       {
               outFile1<< std::setprecision(5) <<zeta_11N[i+11*NrParticles*j]<<'\t' ;
                       }
               outFile1<<std::endl; 
       }	
               outFile1<<std::endl; 
	
		
	for (int l=0; l<5; l++)
		{
		for (int k=0; k<5; k++)
			{				
				// column major format
				Resistance_Tnsr_dd.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	5	+	55*NrParticles	+	66*NrParticles*NrParticles	+	6*NrParticles	];
			}
		}						
					
/*
	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
				for (int d=0; d<3; d++)
					{
					for (int p=0; p<5; p++)
						{
						for (int s=0; s<5; s++)
							{							
								m_ijkl[a][b][g][d]		+=		e_l[p][a][b]*Resistance_Tnsr_dd.comp[p][s]*e_l[s][g][d];		
							}
						}													
					}
				}
			}
		}	
			for (int a=0; a<3; a++)
				{
				for (int b=0; b<3; b++)
					{
					for (int g=0; g<3; g++)
						{
						for (int d=0; d<3; d++)
							{							
								outFile1 << M_IJKL[a][b][g][d] <<	'\t'	<< m_ijkl[a][b][g][d] << '\t'	<< "m_ijkl"  << a << b << g << d << endl;		
							}
						}													
					}
				}
*/	
	/*
		outFile1<<zeta_11N[0]<<'\t'<<zeta_11N[6]<<'\t'<<zeta_11N[12]<<'\t'<<zeta_11N[18]<<'\t'<<zeta_11N[24]<<'\t'<<zeta_11N[30]<<std::endl ;
		outFile1<<zeta_11N[1]<<'\t'<<zeta_11N[7]<<'\t'<<zeta_11N[13]<<'\t'<<zeta_11N[19]<<'\t'<<zeta_11N[25]<<'\t'<<zeta_11N[31]<<std::endl ;
		outFile1<<zeta_11N[2]<<'\t'<<zeta_11N[8]<<'\t'<<zeta_11N[14]<<'\t'<<zeta_11N[20]<<'\t'<<zeta_11N[26]<<'\t'<<zeta_11N[32]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<zeta_11N[3]<<'\t'<<zeta_11N[9]<<'\t'<<zeta_11N[15]<<'\t'<<zeta_11N[21]<<'\t'<<zeta_11N[27]<<'\t'<<zeta_11N[33]<<std::endl ;
		outFile1<<zeta_11N[4]<<'\t'<<zeta_11N[10]<<'\t'<<zeta_11N[16]<<'\t'<<zeta_11N[22]<<'\t'<<zeta_11N[28]<<'\t'<<zeta_11N[34]<<std::endl ;
		outFile1<<zeta_11N[5]<<'\t'<<zeta_11N[11]<<'\t'<<zeta_11N[17]<<'\t'<<zeta_11N[23]<<'\t'<<zeta_11N[29]<<'\t'<<zeta_11N[35]<<std::endl ;
		outFile1<<"asdgfersh"<<endl;
		outFile1<<rho_11N[0]<<'\t'<<rho_11N[6]<<'\t'<<rho_11N[12]<<'\t'<<rho_11N[18]<<'\t'<<rho_11N[24]<<'\t'<<rho_11N[30]<<std::endl ;
		outFile1<<rho_11N[1]<<'\t'<<rho_11N[7]<<'\t'<<rho_11N[13]<<'\t'<<rho_11N[19]<<'\t'<<rho_11N[25]<<'\t'<<rho_11N[31]<<std::endl ;
		outFile1<<rho_11N[2]<<'\t'<<rho_11N[8]<<'\t'<<rho_11N[14]<<'\t'<<rho_11N[20]<<'\t'<<rho_11N[26]<<'\t'<<rho_11N[32]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<rho_11N[3]<<'\t'<<rho_11N[9]<<'\t'<<rho_11N[15]<<'\t'<<rho_11N[21]<<'\t'<<rho_11N[27]<<'\t'<<rho_11N[33]<<std::endl ;
		outFile1<<rho_11N[4]<<'\t'<<rho_11N[10]<<'\t'<<rho_11N[16]<<'\t'<<rho_11N[22]<<'\t'<<rho_11N[28]<<'\t'<<rho_11N[34]<<std::endl ;
		outFile1<<rho_11N[5]<<'\t'<<rho_11N[11]<<'\t'<<rho_11N[17]<<'\t'<<rho_11N[23]<<'\t'<<rho_11N[29]<<'\t'<<rho_11N[35]<<std::endl ;
	*/		
	double h_clst_ijk[3][3][3] = {{{0}}};

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

// if	(del_j)_{\alpha,p}	*	E^{inf\tilde}_{p}	= 	E^inf_{\alpha\beta} *	r_\beta 

// then	(del_j)_{\alpha,p}	*	E^{inf\tilde}_{p}	= 	(e_p)_{\alpha\beta}	*	E^{inf\tilde}_{p}	*	r_\beta	// sicne E^{inf}_{\alpha\beta}	= (e_p)_{\alpha\beta}*E^{inf\tilde}_{\p}

//hence	(del_j)_{\alpha,p}							= 	(e_p)_{\alpha\beta}	*	r_\beta	

					mtrx35D	Delj;
					mtrx53D	Deli;
					mtrx53D	Unit_tnsr_Redc;
					
						for (int k=0; k<5; k++)
						{							
							for (int a=0; a<3; a++)
							{
								Delj.comp[a][k]	=	0.0	;
								Deli.comp[k][a]	=	0.0	;
								for (int b=0; b<3; b++)
								{
									Delj.comp[a][k]	+=	e_l[k][b][a]	*	particle[j].pos.comp[b]	;
									Deli.comp[k][a]	+=	e_l[k][a][b]	*	particle[i].pos.comp[b]	;
									for (int c=0; c<3; c++)
									{
										Unit_tnsr_Redc.comp[k][a]	+=	e_l[k][b][c]		* (b==c)	*	particle[i].pos.comp[a]	;	
									}
								}
							}		
						}
	
					
					for (int l=0; l<3; l++)
						{
							for (int k=0; k<3; k++)
								{
							

// 									11N format 
									Resistance_Tnsr_tt.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j														];
									Resistance_Tnsr_tr.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j	+	33*NrParticles*NrParticles						];
									Resistance_Tnsr_rt.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j	+	3*NrParticles									];
									Resistance_Tnsr_rr.comp[k][l] = zeta_11N[k	+	11*NrParticles*l	+	3*i	+	33*NrParticles*j	+	33*NrParticles*NrParticles	+	3*NrParticles	];


/*
// 									6N format
									Resistance_Tnsr_tt.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j] ;
									Resistance_Tnsr_tr.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles]	;
									Resistance_Tnsr_rt.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+3*NrParticles] ;
									Resistance_Tnsr_rr.comp[k][l] = zeta_11N[k+6*NrParticles*l+3*i+18*NrParticles*j+18*NrParticles*NrParticles+3*NrParticles] ;
								
*/								}
						}


					for (int l=0; l<5; l++)
						{
							for (int k=0; k<3; k++)
								{				
									// column major format
									Resistance_Tnsr_td.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	55*NrParticles*j	+	66*NrParticles*NrParticles						];
									Resistance_Tnsr_rd.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	3*i	+	55*NrParticles*j	+	66*NrParticles*NrParticles	+	3*NrParticles	];						
								}
						}					
					
					for (int l=0; l<3; l++)
						{
							for (int k=0; k<5; k++)
								{				
									// column major format
									Resistance_Tnsr_dt.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	5*i	+	33*NrParticles*j	+	6*NrParticles									];
									Resistance_Tnsr_dr.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	5*i	+	33*NrParticles*j	+	33*NrParticles*NrParticles	+	6*NrParticles	];				
								}
						}
					
					for (int l=0; l<5; l++)
						{
							for (int k=0; k<5; k++)
								{				
									// column major format
									Resistance_Tnsr_dd.comp[k][l] =	zeta_11N[k	+	11*NrParticles*l	+	5*i	+	55*NrParticles*j	+	66*NrParticles*NrParticles	+	6*NrParticles	];
								}
						}						



	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					g_ijk[a][b][g] = 0.0;
				for (int p=0; p<5; p++)
					{						
							g_ijk[a][b][g]	+=		e[p][a][b]*Resistance_Tnsr_td.comp[g][p];		
					}													
				}
			}
		}

outFile1<<i<<j<<endl;
outFile1<<"g_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		outFile1 << setw(10) << g_ijk[s][0][0] << "  " << setw(10) << g_ijk[s][0][1] << "  " << setw(10) << g_ijk[s][0][2] << endl;
		outFile1 << setw(10) << g_ijk[s][1][0] << "  " << setw(10) << g_ijk[s][1][1] << "  " << setw(10) << g_ijk[s][1][2] << endl;
		outFile1 << setw(10) << g_ijk[s][2][0] << "  " << setw(10) << g_ijk[s][2][1] << "  " << setw(10) << g_ijk[s][2][2] << endl;
    }
    
	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					h_ijk[a][b][g] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							h_ijk[a][b][g]	+=		e[p][a][b]*Resistance_Tnsr_rd.comp[g][p];		
					}													
				}
			}
		}

outFile1<<"h_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		outFile1 << setw(10) << h_ijk[s][0][0] << "  " << setw(10) << h_ijk[s][0][1] << "  " << setw(10) << h_ijk[s][0][2] << endl;
		outFile1 << setw(10) << h_ijk[s][1][0] << "  " << setw(10) << h_ijk[s][1][1] << "  " << setw(10) << h_ijk[s][1][2] << endl;
		outFile1 << setw(10) << h_ijk[s][2][0] << "  " << setw(10) << h_ijk[s][2][1] << "  " << setw(10) << h_ijk[s][2][2] << endl;
    }
    
	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					h_clst_ijk[a][b][g] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							h_clst_ijk[a][b][g]	+=		e[p][a][b]*Friction_Tnsr_rd.comp[g][p];		

					}													
				}
			}
		}

outFile1<<"h_clst_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		outFile1 << setw(10) << h_clst_ijk[s][0][0] << "  " << setw(10) << h_clst_ijk[s][0][1] << "  " << setw(10) << h_clst_ijk[s][0][2] << endl;
		outFile1 << setw(10) << h_clst_ijk[s][1][0] << "  " << setw(10) << h_clst_ijk[s][1][1] << "  " << setw(10) << h_clst_ijk[s][1][2] << endl;
		outFile1 << setw(10) << h_clst_ijk[s][2][0] << "  " << setw(10) << h_clst_ijk[s][2][1] << "  " << setw(10) << h_clst_ijk[s][2][2] << endl;
    }
					
					Friction_Tnsr_tt	+=		Resistance_Tnsr_tt ;  
					Friction_Tnsr_tr 	+= 	( 	Resistance_Tnsr_tr		- 	(	Resistance_Tnsr_tt*Aj	)	)	;
				//	Friction_Tnsr_rt 	+= 	(	Ai*Resistance_Tnsr_tt	+		Resistance_Tnsr_rt    	)    	; 
					Friction_Tnsr_rr 	+= 	( 	Resistance_Tnsr_rr    	-	(	Resistance_Tnsr_rt*Aj 	) 	+ 			Ai*Resistance_Tnsr_tr	-	Ai*Resistance_Tnsr_tt*Aj	)	;
				
					for (int l=0; l<5; l++)
						{
						for (int k=0; k<3; k++)
							{
								for (int m=0; m<3; m++)
									{
										Friction_Tnsr_td.comp[k][l]	+=	Resistance_Tnsr_tt.comp[k][m]*Delj.comp[m][l];
										Friction_Tnsr_rd.comp[k][l]	+= 	Resistance_Tnsr_rt.comp[k][m]*Delj.comp[m][l];
										//		Friction_Tnsr_dr.comp[l][k]	-= 	Resistance_Tnsr_dt.comp[l][m]*Aj.comp[m][k];
										for (int n=0; n<3; n++)
											{
												Friction_Tnsr_rd.comp[k][l]	+=	Ai.comp[k][n]	*	( 	Resistance_Tnsr_tt.comp[n][m]	*	Delj.comp[m][l]	)	;
											}
										Friction_Tnsr_rd.comp[k][l]	+=	Ai.comp[k][m]	*	( Resistance_Tnsr_td.comp[m][l]	)	;									
									}
								Friction_Tnsr_td.comp[k][l]	+=	Resistance_Tnsr_td.comp[k][l]	;
								Friction_Tnsr_rd.comp[k][l]	+=	Resistance_Tnsr_rd.comp[k][l]	;								
								Friction_Tnsr_dt.comp[k][l]	+=	Resistance_Tnsr_dt.comp[l][k]	;								
								Friction_Tnsr_dr.comp[l][k]	+=	Resistance_Tnsr_dr.comp[l][k]	;								
							}
						for (int k=0; k<5; k++)
						{								
							 for (int m=0; m<3; m++)
							{
								Friction_Tnsr_dd.comp[l][k]	+= 	Resistance_Tnsr_dt.comp[l][m]*Delj.comp[m][k];
										Friction_Tnsr_dd.comp[l][k]	+=	(	0.5	*	(
																				Deli.comp[l][m]	*	(	Resistance_Tnsr_td.comp[m][k]	)
																				+
																				Deli.comp[k][m]	*	(	Resistance_Tnsr_td.comp[m][l]	)
																				)
																		-	(1.0/3.0)	*	(
																				Unit_tnsr_Redc.comp[l][m]	*	(	Resistance_Tnsr_td.comp[m][k]	)
																				)
																		) ;		
									for (int n=0; n<3; n++)
									{
										Friction_Tnsr_dd.comp[l][k]	+=	(	0.5	*	(
																				Deli.comp[l][m]	*	(	Resistance_Tnsr_tt.comp[m][n]	*	Delj.comp[n][k]	)
																				+
																				Deli.comp[k][m]	*	(	Resistance_Tnsr_tt.comp[m][n]	*	Delj.comp[n][l]	)
																				)
																		-	(1.0/3.0)	*	(
																				Unit_tnsr_Redc.comp[l][m]	*	(	Resistance_Tnsr_tt.comp[m][n]	*	Delj.comp[n][k]	)
																				)
																		) ;		
									}	
							}	
						
							Friction_Tnsr_dd.comp[l][k]	+=	Resistance_Tnsr_dd.comp[l][k]	;		
						} 						

						} 


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


			for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{				
							// column major format
							xi_11x11[k	+	11*l					] 	=	 Friction_Tnsr_tt.comp[k][l];
							xi_11x11[k	+	11*l	+	33			] 	=	 Friction_Tnsr_tr.comp[k][l];
							xi_11x11[k	+	11*l	+	3			] 	=	 Friction_Tnsr_rt.comp[k][l];
							xi_11x11[k	+	11*l	+	33	+	3	] 	=	 Friction_Tnsr_rr.comp[k][l];							
						}
				}
				
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<3; k++)
					{				
						// column major format
						xi_11x11[k	+	11*l	+	66			] 	=	 -Friction_Tnsr_td.comp[k][l];
						xi_11x11[k	+	11*l	+	66	+	3	] 	=	 Friction_Tnsr_rd.comp[k][l];						
					}
				}					
			for (int l=0; l<3; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						xi_11x11[k	+	11*l	+	6			] 	=	 Friction_Tnsr_td.comp[l][k];
						xi_11x11[k	+	11*l	+	33	+	6	] 	=	 Friction_Tnsr_rd.comp[l][k];					
					}
				}
			for (int l=0; l<5; l++)
				{
				for (int k=0; k<5; k++)
					{				
						// column major format
						xi_11x11[k	+	11*l	+	66	+	6	] 	=	 Friction_Tnsr_dd.comp[k][l];
					}
				}

/*	
	 			// 6x6 format					
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
	*/
	
/*		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[33]<<'\t'<<xi_11x11[44]<<'\t'<<xi_11x11[55]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[34]<<'\t'<<xi_11x11[45]<<'\t'<<xi_11x11[56]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[35]<<'\t'<<xi_11x11[46]<<'\t'<<xi_11x11[57]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[36]<<'\t'<<xi_11x11[47]<<'\t'<<xi_11x11[58]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[37]<<'\t'<<xi_11x11[48]<<'\t'<<xi_11x11[59]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[38]<<'\t'<<xi_11x11[49]<<'\t'<<xi_11x11[60]<<std::endl ;	
		outFile1<<std::endl ;
		outFile1<<xi_11x11[6]<<'\t'<<xi_11x11[17]<<'\t'<<xi_11x11[28]<<'\t'<<xi_11x11[39]<<'\t'<<xi_11x11[50]<<'\t'<<xi_11x11[61]<<std::endl ;
		outFile1<<xi_11x11[7]<<'\t'<<xi_11x11[18]<<'\t'<<xi_11x11[29]<<'\t'<<xi_11x11[40]<<'\t'<<xi_11x11[51]<<'\t'<<xi_11x11[62]<<std::endl ;
		outFile1<<xi_11x11[8]<<'\t'<<xi_11x11[19]<<'\t'<<xi_11x11[30]<<'\t'<<xi_11x11[41]<<'\t'<<xi_11x11[52]<<'\t'<<xi_11x11[63]<<std::endl ;
		outFile1<<xi_11x11[9]<<'\t'<<xi_11x11[20]<<'\t'<<xi_11x11[31]<<'\t'<<xi_11x11[42]<<'\t'<<xi_11x11[53]<<'\t'<<xi_11x11[64]<<std::endl ;
		outFile1<<xi_11x11[10]<<'\t'<<xi_11x11[21]<<'\t'<<xi_11x11[32]<<'\t'<<xi_11x11[43]<<'\t'<<xi_11x11[54]<<'\t'<<xi_11x11[65]<<std::endl ;
*/

		outFile1<<std::endl ;
		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[33]<<'\t'<<xi_11x11[44]<<'\t'<<xi_11x11[55]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[34]<<'\t'<<xi_11x11[45]<<'\t'<<xi_11x11[56]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[35]<<'\t'<<xi_11x11[46]<<'\t'<<xi_11x11[57]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[36]<<'\t'<<xi_11x11[47]<<'\t'<<xi_11x11[58]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[37]<<'\t'<<xi_11x11[48]<<'\t'<<xi_11x11[59]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[38]<<'\t'<<xi_11x11[49]<<'\t'<<xi_11x11[60]<<std::endl ;
  
	for (int a=0; a<3; a++)
		{
		for (int b=0; b<3; b++)
			{
			for (int g=0; g<3; g++)
				{
					h_clst_ijk[a][b][g] = 0.0;
					
				for (int p=0; p<5; p++)
					{
							h_clst_ijk[a][b][g]	+=		e[p][a][b]*Friction_Tnsr_rd.comp[g][p];		

					}													
				}
			}
		}

outFile1<<"h_ijk"<<endl;

for (int s=0; s<3; s++)
	{							
		outFile1 << setw(10) << h_clst_ijk[s][0][0] << "  " << setw(10) << h_clst_ijk[s][0][1] << "  " << setw(10) << h_clst_ijk[s][0][2] << endl;
		outFile1 << setw(10) << h_clst_ijk[s][1][0] << "  " << setw(10) << h_clst_ijk[s][1][1] << "  " << setw(10) << h_clst_ijk[s][1][2] << endl;
		outFile1 << setw(10) << h_clst_ijk[s][2][0] << "  " << setw(10) << h_clst_ijk[s][2][1] << "  " << setw(10) << h_clst_ijk[s][2][2] << endl;
    }
					
		
	for (int i=0;i<11;i++)
       { 
               for(int j=0;j<11;j++)
                       {
               outFile1<< std::setprecision(5) <<xi_11x11[i+11*j]<<'\t' ;
                       }
               outFile1<<std::endl; 
       }	
               outFile1<<std::endl; 
		
	inverse ( xi_11x11 , 11 )	 ; 			
	for (int i=0; i<121; i++)
		{
			xi_11x11[i]*=4.1419e-14;	// multiply by kbT in erg K-1
		} 	
	for (int i=0; i<6; i++)
		{
		for (int j=0; j<6; j++)
			{
				xi_6x6[i+j*6]=xi_11x11[i+j*11];	
			} 
		}

// center of diffusion calculation based on "Hydrodynamic properties of rigid particles: comparison of different modeling and computational procedures." Biophysical journal 76.6 (1999): 3044-3057.			
// Page 3046 equation 13.
vctr3D ctr_diff ; 
double temp_mat[3*3];
temp_mat[0] =  xi_6x6[28] + xi_6x6[35]	;
temp_mat[1] = -xi_6x6[27] 				; 
temp_mat[2] = -xi_6x6[33]				; 
temp_mat[3] = -xi_6x6[27]				;
temp_mat[4] =  xi_6x6[21] + xi_6x6[35]	;
temp_mat[5] = -xi_6x6[34]				;
temp_mat[6] = -xi_6x6[33]				;
temp_mat[7] = -xi_6x6[34]				;
temp_mat[8] =  xi_6x6[28] + xi_6x6[21]	;

inverse ( temp_mat , 3 )	 ; 	

ctr_diff.comp[0] = temp_mat[0]*(xi_6x6[16] - xi_6x6[11]) + temp_mat[3]*(xi_6x6[5] - xi_6x6[15]) + temp_mat[6]*(xi_6x6[9] - xi_6x6[4]) ;
ctr_diff.comp[1] = temp_mat[1]*(xi_6x6[16] - xi_6x6[11]) + temp_mat[4]*(xi_6x6[5] - xi_6x6[15]) + temp_mat[7]*(xi_6x6[9] - xi_6x6[4]) ;
ctr_diff.comp[2] = temp_mat[2]*(xi_6x6[16] - xi_6x6[11]) + temp_mat[5]*(xi_6x6[5] - xi_6x6[15]) + temp_mat[8]*(xi_6x6[9] - xi_6x6[4]) ;

mtrx3D D_tt, D_tr, D_rt , D_rr, U_OD;

D_tt.comp[0][0] = xi_6x6[0];
D_tt.comp[1][0] = xi_6x6[1];
D_tt.comp[2][0] = xi_6x6[2];
D_tt.comp[0][1] = xi_6x6[6];
D_tt.comp[1][1] = xi_6x6[7];
D_tt.comp[2][1] = xi_6x6[8];
D_tt.comp[0][2] = xi_6x6[12];
D_tt.comp[1][2] = xi_6x6[13];
D_tt.comp[2][2] = xi_6x6[14];

D_tr.comp[0][0] = xi_6x6[3];
D_tr.comp[1][0] = xi_6x6[4];
D_tr.comp[2][0] = xi_6x6[5];
D_tr.comp[0][1] = xi_6x6[9];
D_tr.comp[1][1] = xi_6x6[10];
D_tr.comp[2][1] = xi_6x6[11];
D_tr.comp[0][2] = xi_6x6[15];
D_tr.comp[1][2] = xi_6x6[16];
D_tr.comp[2][2] = xi_6x6[17];

D_rt.comp[0][0] = xi_6x6[18];
D_rt.comp[1][0] = xi_6x6[19];
D_rt.comp[2][0] = xi_6x6[20];
D_rt.comp[0][1] = xi_6x6[24];
D_rt.comp[1][1] = xi_6x6[25];
D_rt.comp[2][1] = xi_6x6[26];
D_rt.comp[0][2] = xi_6x6[30];
D_rt.comp[1][2] = xi_6x6[31];
D_rt.comp[2][2] = xi_6x6[32];

D_rr.comp[0][0] = xi_6x6[21];
D_rr.comp[1][0] = xi_6x6[22];
D_rr.comp[2][0] = xi_6x6[23];
D_rr.comp[0][1] = xi_6x6[27];
D_rr.comp[1][1] = xi_6x6[28];
D_rr.comp[2][1] = xi_6x6[29];
D_rr.comp[0][2] = xi_6x6[33];
D_rr.comp[1][2] = xi_6x6[34];
D_rr.comp[2][2] = xi_6x6[35];

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

xi_6x6[0]  = D_tt_CoD.comp[0][0];
xi_6x6[1]  = D_tt_CoD.comp[1][0];
xi_6x6[2]  = D_tt_CoD.comp[2][0];
xi_6x6[6]  = D_tt_CoD.comp[0][1];
xi_6x6[7]  = D_tt_CoD.comp[1][1];
xi_6x6[8]  = D_tt_CoD.comp[2][1];
xi_6x6[12] = D_tt_CoD.comp[0][2];
xi_6x6[13] = D_tt_CoD.comp[1][2];
xi_6x6[14] = D_tt_CoD.comp[2][2];

xi_6x6[3]  = D_tr_CoD.comp[0][0];
xi_6x6[4]  = D_tr_CoD.comp[1][0];
xi_6x6[5]  = D_tr_CoD.comp[2][0];
xi_6x6[9]  = D_tr_CoD.comp[0][1];
xi_6x6[10] = D_tr_CoD.comp[1][1];
xi_6x6[11] = D_tr_CoD.comp[2][1];
xi_6x6[15] = D_tr_CoD.comp[0][2];
xi_6x6[16] = D_tr_CoD.comp[1][2];
xi_6x6[17] = D_tr_CoD.comp[2][2];

xi_6x6[18]  = D_rt_CoD.comp[0][0];
xi_6x6[19]  = D_rt_CoD.comp[1][0];
xi_6x6[20]  = D_rt_CoD.comp[2][0];
xi_6x6[24]  = D_rt_CoD.comp[0][1];
xi_6x6[25]  = D_rt_CoD.comp[1][1];
xi_6x6[26]  = D_rt_CoD.comp[2][1];
xi_6x6[30]  = D_rt_CoD.comp[0][2];
xi_6x6[31]  = D_rt_CoD.comp[1][2];
xi_6x6[32]  = D_rt_CoD.comp[2][2];

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

		outFile1<<std::endl ;
		outFile1<<xi_11x11[0]<<'\t'<<xi_11x11[11]<<'\t'<<xi_11x11[22]<<'\t'<<xi_11x11[33]<<'\t'<<xi_11x11[44]<<'\t'<<xi_11x11[55]<<std::endl ;
		outFile1<<xi_11x11[1]<<'\t'<<xi_11x11[12]<<'\t'<<xi_11x11[23]<<'\t'<<xi_11x11[34]<<'\t'<<xi_11x11[45]<<'\t'<<xi_11x11[56]<<std::endl ;
		outFile1<<xi_11x11[2]<<'\t'<<xi_11x11[13]<<'\t'<<xi_11x11[24]<<'\t'<<xi_11x11[35]<<'\t'<<xi_11x11[46]<<'\t'<<xi_11x11[57]<<std::endl ;
		outFile1<<std::endl ;
		outFile1<<xi_11x11[3]<<'\t'<<xi_11x11[14]<<'\t'<<xi_11x11[25]<<'\t'<<xi_11x11[36]<<'\t'<<xi_11x11[47]<<'\t'<<xi_11x11[58]<<std::endl ;
		outFile1<<xi_11x11[4]<<'\t'<<xi_11x11[15]<<'\t'<<xi_11x11[26]<<'\t'<<xi_11x11[37]<<'\t'<<xi_11x11[48]<<'\t'<<xi_11x11[59]<<std::endl ;
		outFile1<<xi_11x11[5]<<'\t'<<xi_11x11[16]<<'\t'<<xi_11x11[27]<<'\t'<<xi_11x11[38]<<'\t'<<xi_11x11[49]<<'\t'<<xi_11x11[60]<<std::endl ;
		outFile1<<xi_11x11[72]<<'\t'<<xi_11x11[84]<<'\t'<<xi_11x11[96]<<'\t'<<xi_11x11[109]<<'\t'<<xi_11x11[120]<<std::endl ;
	//	outFile1<<ctr_diff.comp[0]<<'\t'<<ctr_diff.comp[1]<<'\t'<<ctr_diff.comp[2]<<std::endl ;

			for (int l=0; l<3; l++)
				{
					for (int k=0; k<3; k++)
						{				
							// column major format
							Friction_Tnsr_tt.comp[k][l]	=	xi_11x11[k	+	11*l					];
							Friction_Tnsr_tr.comp[k][l]	=	xi_11x11[k	+	11*l	+	33			];
							Friction_Tnsr_rt.comp[k][l]	=	xi_11x11[k	+	11*l	+	3			];
							Friction_Tnsr_rr.comp[k][l]	=	xi_11x11[k	+	11*l	+	33	+	3	];							
						}
				}
	/*			
	for (int i=0; i<11; i++)
		{
			outFile1<<xi_11x11[0+i]<<'\t'<<xi_11x11[11+i]<<'\t'<<xi_11x11[22+i]<<'\t'<<xi_11x11[33+i]<<'\t'<<xi_11x11[44+i]<<'\t'<<xi_11x11[55+i]<<'\t'<<xi_11x11[66+i]<<'\t'<<xi_11x11[77+i]<<'\t'<<xi_11x11[88+i]<<'\t'<<xi_11x11[99+i]<<'\t'<<xi_11x11[110+i]<<std::endl ;
		} 
	*/					
/*	Friction_Tnsr_tt.echo();
	Friction_Tnsr_tr.echo();
	Friction_Tnsr_rr.echo(); */
}
		
