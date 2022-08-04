
#include "vec.h"
#include "nondestructivetrimesh.h"
#include <vector>
#include <Eigen/Dense> // most of the vector functions I will need inside of an element
#include <Eigen/Sparse> // functions for solution of linear systems
#include <Eigen/OrderingMethods>
#include <math.h> 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>


//  Linlin Li 2019
//  update note: adding sizzled 

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
using namespace Eigen;



void advance_diffusion(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, std::vector<std::vector<double> > &species,double dt, double t, const std::vector<double> &diffusionParam, int &iterations, int &max_it ,int &force_postive_check)
{
    

//%% get parameters
int ne = mesh.num_triangles();// % get number of element
int nn = x.size(); //            % get total node number
//int nne = 3;//                   % get number of node per element 3 for triangle mesh

//int total_parameters = 23; // setup total number of the screening parameters



//std::map<std::string, int> parameters;
//set up the parmeters dicts------name should match with the sequential position with the input file -----------------------------
//string parameters_names[total_parameters] = {"DN","DBC","DBN","j1","j2","j3","k1","k2","decN","decBC","decBN","lambda_tld_C","lambda_tld_BC","DS","VS","n","decS","k4"};
//string species_name[total_species]={"B","C","BC","N","BN","S"}
//string parameters_names[total_parameters] = {"DN","DBC","DBN","j1","j2","j3","k1","k2","decN","decBC","decBN","lambda_tld_C","lambda_tld_BC"};
//string species_name[total_species]={"B","C","BC","N","BN"}


// --------------constant parameters for geomatry-----------------
//double r   = 350; // redius of embryo ( may need to change if the mesh changes)
// set up mesh number for expression map
int meshN = 101;
double pi = 3.14159265;
int nodeN = 3; // set up node number per mesh

// --------------initialize constant parameters-----------------
// ---------------diffusion rate---------------------
int total_species;
double DB  = 4.4; //            % diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
double DC  = 7; //            % diffusion rate of ligand (Chd)    (microns^2*s^-1)*60s*m^-1
double DN;//  = 10; //            % diffusion rate of ligand (Nog)    (microns^2*s^-1)*60s*m^-1
double DBC;// = 4.7475; //            % diffusion rate of ligand (BC)    (microns^2*s^-1)*60s*m^-1
double DBN;// = 7.4747; //            % diffusion rate of ligand (BC)    (microns^2*s^-1)*60s*m^-1
double DS = 10 ;
double DTld;
// --------------- expression constants ---------------------
double j1;// = 6.33; //          % production rate of BMP          nM*m^-1
double j2;// = 55.04; //          % production rate of Chd          nM*m^-1
double j3;// = 38.51; //          % production rate of Nog          nM*m^-1
double j4;// = 38.51; //          % production rate of Tld          nM*m^-1
// --------------- reaction  constants ---------------------
double k1;//  = 0.7539; //      % binding rates for BMP ligand and Chordin          nM^-1*m^-1
double k_1;//  = k2 ; //        % unbinding rates for BMP ligand and Chordin        m^-1
double k2;//   = 0.6044; //      % binding rates for BMP ligand and Noggin          nM^-1*m^-1
double k_2;//  = 0.1*k2; //     % unbinding rates for BMP ligand and Noggin        m^-1
// --------------- Nature decay rates ---------------------
double decB  = 8.7e-5; //  % decay rate of Ligand (BMP)    nM*m^-1
double decC  = 9.6e-5; //  % decay rate of Chd             nM*m^-1
double decN;//  = 9.3658e-4; //  % decay rate of Nog            nM*m^-1
double decBC;// = 4.7279e-4; //  % decay rate of BC             nM*m^-1
double decBN;// = 6.2045e-4; //  % decay rate of BN             nM*m^-1
double decS; // dycay rate for sizzled
double decTld; // dycay rate for Tld
// --------------- Tld processing constants ---------------------
double lambda_tld_C;//  = 6571;  //Tld processing rate of Chd nM^-1*m^-1
double lambda_tld_BC;// = 5353;  //Tld processing rate of BC nM^-1*m^-1


//double tld_parameter = 1;      // Tld concentration
//-------------------paramters related sizzled feedback------------
double kit; // tld feedback with sizzled 
double kmt; // tld feedback with C+BC
// double kia; // bmp1a feedback with sizzled
// double kma; // bmp1a feedback with C+BC
double Vs;
double sn;
double ks;

// --------------- sizzled parameter proparing constants ---------------------
double mutant_bmpmax;
double mutant_szlmax;
double szl_prameter_sampling;

//double bmpmax;
//double szlmax;


/// replace with the parameters from the script file
// diffusionParam[0] == species_n;
total_species = diffusionParam[0]; // setup total nuber of species (equations)
DN = diffusionParam[1];
DBC= diffusionParam[2];
DBN= diffusionParam[3];
j1 = diffusionParam[4];
j2 = diffusionParam[5];
j3 = diffusionParam[6];
k1 = diffusionParam[7];
k_1 = diffusionParam[8];
k2 = diffusionParam[9];
k_2 = diffusionParam[10];
decN = diffusionParam[11];
decBC = diffusionParam[12];
decBN = diffusionParam[13];
lambda_tld_C = diffusionParam[14];
lambda_tld_BC = diffusionParam[15];
decS = diffusionParam[16];
kit = diffusionParam[17]; // tld feedback with sizzled 
kmt = diffusionParam[18]; // tld feedback with C+BC
Vs  = diffusionParam[19];
sn = diffusionParam[20];
ks = diffusionParam[21];


DTld = diffusionParam[22];
decTld = diffusionParam[23];
j4 = diffusionParam[24];

mutant_bmpmax = diffusionParam[25];
mutant_szlmax = diffusionParam[26];
szl_prameter_sampling = diffusionParam[27];


//std::cout<<"total_species "<<total_species<<"\n";



double tol=1e-3;          //  % tolerance for Newton method
// double max_it=15;       // max iteration number

// -------------initial conditions---------------------------------------------
// initializing the residual vector, stiffness matrix, sparse linear solver
VectorXd RR(nn*total_species); 
VectorXd dU(nn*total_species);
SparseMatrix<double, ColMajor> KK(nn*total_species,nn*total_species);
std::vector<T> KK_triplets;
SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> >   solver;
 //gass value for the current concetration


// need to setup for the initial concentration for all the species--------------

//std::cout<<"nn="<<nn<<"\n";  
std::vector<double> c_B0(nn);
std::vector<double> c_C0(nn); 
std::vector<double> c_BC0(nn); 
std::vector<double> c_N0(nn); 
std::vector<double> c_BN0(nn); 
std::vector<double> c_S0(nn);
std::vector<double> c_Tld0(nn);

std::vector<double> c_B(nn);
std::vector<double> c_C(nn); 
std::vector<double> c_BC(nn); 
std::vector<double> c_N(nn); 
std::vector<double> c_BN(nn); 
std::vector<double> c_S(nn);
std::vector<double> c_Tld(nn);
//std::cout<<"c_B size = "<<c_B.size()<<"\n";  
//std::cout<<"species size = "<<species.size()<<"\n";
//std::cout<<"species 1 1 = "<<species[0][0]<<"\n";

// force concentration to positive if iteration over max number

for (int i=0; i<nn; ++i){

    if (force_postive_check == 1)
    {
    //std::cout << "force_postive_check turned on, inside diffusion\n"; 
    double zero = 0.0;
    c_B0[i]= std::max(species[i][0],zero);
    c_C0[i]= std::max(species[i][1],zero);
    c_BC0[i]= std::max(species[i][2],zero);
    c_N0[i]= std::max(species[i][3],zero);
    c_BN0[i]= std::max(species[i][4],zero);
    c_S0[i]= std::max(species[i][5],zero);
    c_Tld0[i]= std::max(species[i][6],zero);

    c_B[i]= std::max(species[i][0],zero);
    c_C[i]= std::max(species[i][1],zero);
    c_BC[i]= std::max(species[i][2],zero);
    c_N[i]= std::max(species[i][3],zero);
    c_BN[i]= std::max(species[i][4],zero);
    c_S[i]= std::max(species[i][5],zero);
    c_Tld[i]= std::max(species[i][6],zero);  
    }
    else
    {
    c_B0[i]= species[i][0];
    c_C0[i]= species[i][1];
    c_BC0[i]= species[i][2];
    c_N0[i]= species[i][3];
    c_BN0[i]= species[i][4];
    c_S0[i]= species[i][5];
    c_Tld0[i]= species[i][6];

    c_B[i]= species[i][0];
    c_C[i]= species[i][1];
    c_BC[i]= species[i][2];
    c_N[i]= species[i][3];
    c_BN[i]= species[i][4];
    c_S[i]= species[i][5];
    c_Tld[i]= species[i][6];   
    }
    
}


//std::cout<<"c_B="<<c_B[1]<<"\n";


//std::cout<<"c_B0="<<c_B0<<"\n";  
//------------------------------------------------------------------------------

// set up time scale for reae the expression map
double t35 = 0;
double t47 = 4320;
double t57 = 7920;
double t63 = 10080;
double t8 = 16020;
double t9 = 20000;

// std::cout<<"dt ="<<dt<<"\n";
// std::cout<<"t="<<t<<"\n";

int t1;
int t2;
std::string filename_bmp;
std::string filename_chd;
std::string filename_nog;
// std::string filename_szl;
// file for pervious time point
std::string filename_bmp_prep;
std::string filename_chd_prep;
std::string filename_nog_prep;

if (t<t47)
{ 
t1 = t35;  // data range between 35-47 
t2 = t47;
filename_bmp = "data/BMPmap_47.txt";
filename_chd = "data/Chdmap_47.txt";
filename_nog = "data/Nogmap_47.txt";

filename_bmp_prep = "data/BMPmap_35.txt";
filename_chd_prep = "data/Chdmap_35.txt";
filename_nog_prep = "data/Nogmap_35.txt";

//filename_szl = "data/Szlmap_47.txt";
}  
else if (t<t57)
{
t1 = t47;  // data range between 47-57 
t2 = t57;
filename_bmp="data/BMPmap_57.txt";
filename_chd = "data/Chdmap_57.txt";
filename_nog = "data/Nogmap_57.txt";
//filename_szl = "data/Szlmap_57.txt";

filename_bmp_prep = "data/BMPmap_47.txt";
filename_chd_prep = "data/Chdmap_47.txt";
filename_nog_prep = "data/Nogmap_47.txt";
}        
else if (t<t63)
{ 

t1 = t57;  // data range between 57-63 
t2 = t63; 
filename_bmp="data/BMPmap_63.txt"; 
filename_chd = "data/Chdmap_63.txt";
filename_nog = "data/Nogmap_63.txt";
//filename_szl = "data/Szlmap_63.txt";   

filename_bmp_prep = "data/BMPmap_57.txt";
filename_chd_prep = "data/Chdmap_57.txt";
filename_nog_prep = "data/Nogmap_57.txt";
}
else if (t<t8)
{ 
t1 = t63;  // data range between 57-63 
t2 = t8; 
filename_bmp = "data/BMPmap_8.txt"; 
filename_chd = "data/Chdmap_8.txt";
filename_nog = "data/Nogmap_8.txt";
//filename_szl = "data/Szlmap_8.txt"; 
filename_bmp_prep = "data/BMPmap_63.txt";
filename_chd_prep = "data/Chdmap_63.txt";
filename_nog_prep = "data/Nogmap_63.txt";  
}

else
{ 
t1 = t8;  // data range between 57-63 
t2 = t9; 
filename_bmp = "data/BMPmap_9.txt"; 
filename_chd = "data/Chdmap_9.txt";
filename_nog = "data/Nogmap_9.txt";
//filename_szl = "data/Szlmap_8.txt"; 
filename_bmp_prep = "data/BMPmap_8.txt";
filename_chd_prep = "data/Chdmap_8.txt";
filename_nog_prep = "data/Nogmap_8.txt";  
}    

// std::cout<<"frame "<<currframe<<"\n";
// std::cout<<"time "<<t<<"\n";
// std::cout<<"filename "<<filename<<"\n";

double bmpdata[10201];
double chddata[10201];
double nogdata[10201];

double bmpdata_pre[10201];
double chddata_pre[10201];
double nogdata_pre[10201];
//double szldata[10201];

// std::string line;
std::ifstream bmpfile(filename_bmp.c_str());
std::ifstream chdfile(filename_chd.c_str());
std::ifstream nogfile(filename_nog.c_str());

std::ifstream bmpfile_pre(filename_bmp_prep.c_str());
std::ifstream chdfile_pre(filename_chd_prep.c_str());
std::ifstream nogfile_pre(filename_nog_prep.c_str());
//std::ifstream szlfile(filename_szl.c_str());

for (int i = 0; i < 10201; i++) 
    {
    bmpfile >> bmpdata[i];
    chdfile >> chddata[i];
    nogfile >> nogdata[i];

    bmpfile_pre >> bmpdata_pre[i];
    chdfile_pre >> chddata_pre[i];
    nogfile_pre >> nogdata_pre[i];
//    szlfile >> szldata[i];
    }


//% iteration for  Newton-Raphson method
//l------------------loop for iteration-----------------------
 for (int it=0;it<max_it+3;it++)
{    
  // clear residual and stiffness for iteration
  KK_triplets.clear();
  RR.setZero();

// --------------loop for element start
// this code is for solve the diffusion reaction equatio
  for(int i=0;i<ne;i++)
  {
  // for triangle i   
// ---------------------loop for node start ---------
   // initialize node vectors
    Vector3d nodeindex;
    Vector3d nodex;
    Vector3d nodey;
    Vector3d nodez;

    Vector3d Uc_B0;
    Vector3d Uc_C0;
    Vector3d Uc_BC0;
    Vector3d Uc_N0;
    Vector3d Uc_BN0;
    Vector3d Uc_S0;
    Vector3d Uc_Tld0;

    Vector3d Uc_Bt;
    Vector3d Uc_Ct;
    Vector3d Uc_BCt;
    Vector3d Uc_Nt;
    Vector3d Uc_BNt;
    Vector3d Uc_St;
    Vector3d Uc_Tldt;

    Vector3d az;
    Vector3d el;
    //Vector3d index;

    // initialize source terms
    Vector3d bmpmap;
    Vector3d chdmap;
    Vector3d nogmap;
    Vector3d szlmap;
    Vector3d s_B;
    Vector3d s_C;
    Vector3d s_N;
    Vector3d s_Tld;

 
    for (int j=0; j<nodeN; j++)
    {
    
    nodeindex(j) = mesh.get_triangle(i)[j];   
    //std::cout <<"mesh nodes index = "<<nodeindex(j)<<"\n";

    // get coordinates
    // coords of node x    
    nodex(j) = x[nodeindex(j)][0];
    // coords of node x  
    nodey(j) = x[nodeindex(j)][1];
    // coords of node x
    nodez(j) = x[nodeindex(j)][2];

    // get concentrations at time 0
    Uc_B0(j) = c_B0[nodeindex[j]];    
    Uc_C0(j) = c_C0[nodeindex[j]];
    Uc_BC0(j) = c_BC0[nodeindex[j]];  
    Uc_N0(j) = c_N0[nodeindex[j]];   
    Uc_BN0(j) = c_BN0[nodeindex[j]];   
    Uc_S0(j) = c_S0[nodeindex[j]];
    Uc_Tld0(j) = c_Tld0[nodeindex[j]];

    // set concentration at time t
    Uc_Bt(j) = c_B[nodeindex[j]];    
    Uc_Ct(j) = c_C[nodeindex[j]];
    Uc_BCt(j) = c_BC[nodeindex[j]];  
    Uc_Nt(j) = c_N[nodeindex[j]];   
    Uc_BNt(j) = c_BN[nodeindex[j]];   
    Uc_St(j) = c_S[nodeindex[j]];
    Uc_Tldt(j) = c_Tld[nodeindex[j]];
    // get currtet position and index in az/el map
    az(j) = atan2(nodey[j],nodex[j]);
    el(j) = atan2(nodez[j],sqrt(nodex[j]*nodex[j] + nodey[j]*nodey[j]));
    int index = (round(((az[j]/pi+1)/2*(meshN-1)))*(meshN)+round((el[j]/(pi/2)+1)/2*(meshN-1)));
    //std::cout<<"index="<<index<<"\n";  
    // simple expression source
    bmpmap(j) = (bmpdata[index]*(t-t1)+bmpdata_pre[index]*(t2-t))/(t2-t1);
    chdmap(j) = (chddata[index]*(t-t1)+chddata_pre[index]*(t2-t))/(t2-t1);
    nogmap(j) = (nogdata[index]*(t-t1)+nogdata_pre[index]*(t2-t))/(t2-t1);
    // bmpmap(j) = bmpdata[index];
    // chdmap(j) = chddata[index];
    // nogmap(j) = nogdata[index];

    //szlmap(j) = szldata[index];
 
    s_B(j) = j1*bmpmap[j];//j1/(1+pow((nodex[j]+350)/350,50));//j1*bmpmap[j];
    s_C(j) = j2*chdmap[j];//pow((nodex[j]+r)/(2*r),50)/(pow(0.8,50)+pow((nodex[j]+r)/(2*r),50));
    s_N(j) = j3*nogmap[j];//pow((nodex[j]+r)/(2*r),50)/(pow(0.8,50)+pow((nodex[j]+r)/(2*r),50));
    s_Tld(j) = j4*bmpmap[j];
    }// end of node loop for initialization

    //std::cout<<s_B<<"\n";
    //double tld_conc=(t/(2000+t))*tld_parameter;
    //double tld_conc=tld_parameter;
      
    Vector3d dNxi;
    dNxi(0)=-1.0;
    dNxi(1)=1.0;
    dNxi(2)=0.0;
    Vector3d dNeta;
    dNeta(0)=-1.0;
    dNeta(1)=0.0;
    dNeta(2)=1.0;
    Vector3d N;
    N(0)=1.0/3.0;
    N(1)=1.0/3.0;
    N(2)=1.0/3.0;

    //double GradientN[2][3];
    MatrixXd GradientN(2,3);
    GradientN(0,0)=dNxi(0);
    GradientN(0,1)=dNxi(1);
    GradientN(0,2)=dNxi(2);
    GradientN(1,0)=dNeta(0);
    GradientN(1,1)=dNeta(1);
    GradientN(1,2)=dNeta(2);
    Vector3d g1;
    g1(0)=dNxi(0)*nodex(0)+dNxi(1)*nodex(1)+dNxi(2)*nodex(2);
    g1(1)=dNxi(0)*nodey(0)+dNxi(1)*nodey(1)+dNxi(2)*nodey(2);
    g1(2)=dNxi(0)*nodez(0)+dNxi(1)*nodez(1)+dNxi(2)*nodez(2);
    //std::cout<<"g1\n"<<g1<<"\n";
    Vector3d g2;
    g2(0)=dNeta(0)*nodex(0)+dNeta(1)*nodex(1)+dNeta(2)*nodex(2);
    g2(1)=dNeta(0)*nodey(0)+dNeta(1)*nodey(1)+dNeta(2)*nodey(2);
    g2(2)=dNeta(0)*nodez(0)+dNeta(1)*nodez(1)+dNeta(2)*nodez(2);
    //std::cout<<"g2\n"<<g2<<"\n";

    //inverse matrix 
    Vector3d e1=g1/g1.norm();
    Vector3d g1xg2 = g1.cross(g2);
    Vector3d n = g1xg2/g1xg2.norm();
    Vector3d e2=n.cross(e1);
    //std::cout<<"e2\n"<<e2<<"\n";
        
    Matrix3d A;
    A.col(0) = e1;
    A.col(1) = e2;
    A.col(2) = n;
    
    Matrix3d px;
    px(0,0)=nodex(0);
    px(1,0)=nodey(0);
    px(2,0)=nodez(0);
    px(0,1)=nodex(1);
    px(1,1)=nodey(1);
    px(2,1)=nodez(1);
    px(0,2)=nodex(2);
    px(1,2)=nodey(2);
    px(2,2)=nodez(2);

    // get new coordinate
    Matrix3d poiprime = A.inverse()*px;

    // took the first two row as the local coordinates for  surface diffusion
    MatrixXd poiprime2D(2,3);
    poiprime2D.row(0)=poiprime.row(0);
    poiprime2D.row(1)=poiprime.row(1);
    //std::cout<<"poiprime2D\n"<<poiprime2D<<"\n";
    // get the element concentraion for current time point
    double c_B0 = N.dot(Uc_B0);
    double c_C0 = N.dot(Uc_C0);
    double c_BC0 = N.dot(Uc_BC0);
    double c_N0 = N.dot(Uc_N0);
    double c_BN0 = N.dot(Uc_BN0);
    double c_S0 = N.dot(Uc_S0);
    double c_Tld0 = N.dot(Uc_Tld0);
    // get the element concentraion for test concentraion  
    double c_Bt = N.dot(Uc_Bt);
    double c_Ct = N.dot(Uc_Ct);
    double c_BCt = N.dot(Uc_BCt);
    double c_Nt = N.dot(Uc_Nt);
    double c_BNt = N.dot(Uc_BNt);
    double c_St = N.dot(Uc_St);  
    double c_Tldt = N.dot(Uc_Tldt);

    Matrix2d J = GradientN*poiprime2D.transpose(); //% Jacobian Matrix
    double DetJ = J.determinant();//% Jacobian determinate
    //std::cout<<"DetJ "<<DetJ<<"\n";
    MatrixXd B = J.inverse()*GradientN;
    //std::cout<<"B "<<B<<"\n";
       
    Matrix2d eye = MatrixXd::Identity(2,2);
    Matrix2d D_B = DB*eye;//                 % apply diffusion rate
    Matrix2d D_C = DC*eye;//                 % apply diffusion rate
    Matrix2d D_BC = DBC*eye;//                 % apply diffusion rate
    Matrix2d D_N = DN*eye;//                 % apply diffusion rate
    Matrix2d D_BN = DBN*eye;//                 % apply diffusion rate
    Matrix2d D_S = DS*eye;//                 % apply diffusion rate
    Matrix2d D_Tld = DTld*eye;//                 % apply diffusion rate

    Matrix3d ke_B = B.transpose()*D_B*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_C = B.transpose()*D_C*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_BC = B.transpose()*D_BC*B*DetJ;             //% element stiffness matrix
    Matrix3d ke_N = B.transpose()*D_N*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_BN = B.transpose()*D_BN*B*DetJ;             //% element stiffness matrix
    Matrix3d ke_S = B.transpose()*D_S*B*DetJ;  
    Matrix3d ke_Tld = B.transpose()*D_Tld*B*DetJ; 

    //std::cout<<"ke_BMB= "<<ke_B<<"\n";
    //std::cout<<"ke_CHD "<<ke_CHD<<"\n";
    //std::cout<<"ke_BC "<<ke_BC<<"\n";
    double se_B=N.dot(s_B);//                % element source term of BMP
    double se_C=N.dot(s_C);//                % element source term of CHD
    double se_N=N.dot(s_N);//                % element source term of BMP
    double se_Tld=N.dot(s_Tld);//                % element source term of BMP

    //std::cout<<"se_BMP="<<se_B<<"\n";  
    //std::cout<<"se_Tld="<<se_Tld<<"\n"; 

    // setup total source term for each equations----------------------------------------------------------
    // set up the extra term for tld and bmp1a first since we need them later------------------------------
    //double Tld = 1/(1+c_St/kit+(c_Ct+c_BCt)/kmt);
    //double bmp1a = /(1+c_St/kia+(c_Ct+c_BCt)/kma)


    // double bmp1a_dS = -bmp1a*bmp1a/kia;
    // double bmp1a_dC = -bmp1a*bmp1a/kma;
    // double bmp1a_dBC = -bmp1a*bmp1a/kma;

    // ---------------------eqaution for souce terms original---------------------------------------------
    // setup the feedback terms
    double feedback =1/(1+c_St/kit+(c_Ct+c_BCt)/kmt);
    double feedback_dC  = -1/(kmt*pow(1+c_St/kit+(c_Ct+c_BCt)/kmt,2));
    double feedback_dBC = -1/(kmt*pow(1+c_St/kit+(c_Ct+c_BCt)/kmt,2));
    double feedback_dS  = -1/(kit*pow(1+c_St/kit+(c_Ct+c_BCt)/kmt,2));

    double Esource_B   = se_B - k1*c_Bt*c_Ct + k_1*c_BCt - decB*c_Bt - k2*c_Bt*c_Nt + k_2*c_BNt + lambda_tld_BC*feedback*c_Tldt*c_BCt;
    double Esource_C   = se_C - k1*c_Bt*c_Ct + k_1*c_BCt - decC*c_Ct - lambda_tld_C*feedback*c_Tldt*c_Ct;
    double Esource_BC  = k1*c_Bt*c_Ct - k_1*c_BCt - decBC*c_BCt - lambda_tld_BC*feedback*c_Tldt*c_BCt;
    double Esource_N   = se_N - k2*c_Bt*c_Nt + k_2*c_BNt - decN*c_Nt;
    double Esource_BN  = k2*c_Bt*c_Nt - k_2*c_BNt - decBN*c_BNt;
    double Esource_S   = Vs*pow(c_Bt,sn)/(pow(ks,sn)+pow(c_Bt,sn))-decS*c_St;
    double Esource_Tld = se_Tld-decTld*c_Tldt;


 
    //--------------derivitives for each speices vs each other-------------------------------------------------
    // derivitives for BMP------------------------------------
    double dEB_dB = -k1*c_Ct-k2*c_Nt-decB;
    double dEB_dC = -k1*c_Bt+lambda_tld_BC*feedback_dC*c_Tldt*c_BCt;
    double dEB_dBC = k_1+lambda_tld_BC*feedback*c_Tldt+lambda_tld_BC*feedback_dBC*c_Tldt*c_BCt;
    double dEB_dN = -k2*c_Bt;
    double dEB_dBN = k_2;
    double dEB_dS = lambda_tld_BC*feedback_dS*c_Tldt*c_BCt;
    double dEB_dTld = lambda_tld_BC*c_BCt*feedback;

    // derivitives for CHD------------------------------------
    double dEC_dB = -k1*c_Ct;
    double dEC_dC = -k1*c_Bt-lambda_tld_C*feedback*c_Tldt-decC-lambda_tld_C*feedback_dC*c_Tldt*c_Ct;
    double dEC_dBC = k_1-lambda_tld_C*feedback_dBC*c_Tldt*c_Ct;
    double dEC_dN = 0;
    double dEC_dBN = 0;
    double dEC_dS = -lambda_tld_C*feedback_dS*c_Tldt*c_Ct;
    double dEC_dTld = -lambda_tld_C*c_Ct*feedback;

    // derivitive for BC--------------------------------------
    double dEBC_dB = k1*c_Ct;
    double dEBC_dC = k1*c_Bt-lambda_tld_BC*feedback_dC*c_Tldt*c_BCt;
    double dEBC_dBC = -k_1-lambda_tld_BC*feedback*c_Tldt-decBC-lambda_tld_BC*feedback_dBC*c_Tldt*c_BCt;
    double dEBC_dN = 0;
    double dEBC_dBN = 0;
    double dEBC_dS = -lambda_tld_BC*feedback_dS*c_Tldt*c_BCt;
    double dEBC_dTld = -lambda_tld_BC*c_BCt*feedback;

    // derivitive for N--------------------------------------
    double dEN_dB = -k2*c_Nt;
    double dEN_dC = 0;
    double dEN_dBC = 0;
    double dEN_dN = -k2*c_Bt-decN;
    double dEN_dBN = k_2;
    double dEN_dS = 0;
    double dEN_dTld = 0;

    // derivitive for BN--------------------------------------
    double dEBN_dB = k2*c_Nt;
    double dEBN_dC = 0;
    double dEBN_dBC = 0;
    double dEBN_dN = k2*c_Bt;
    double dEBN_dBN = -k_2-decBN;
    double dEBN_dS = 0;
    double dEBN_dTld = 0;


    // derivitive for S--------------------------------------
    double dES_dB = pow(ks,sn)*sn*Vs*pow(c_Bt,sn-1)/pow(pow(c_Bt,sn)+pow(ks,sn),2);
    double dES_dC = 0;
    double dES_dBC = 0;
    double dES_dN = 0;
    double dES_dBN = 0;
    double dES_dS = - decS;
    double dES_dTld = 0;


    // derivitive for Tld--------------------------------------
    double dETld_dB = 0;
    double dETld_dC =0;// - lambda_tld_BC*feedback_dC*c_Tldt*c_BCt- lambda_tld_C*feedback_dC*c_Tldt*c_Ct- lambda_tld_C*feedback*c_Tldt;
    double dETld_dBC =0;// -lambda_tld_BC*feedback_dBC*c_Tldt*c_BCt-lambda_tld_BC*feedback*c_Tldt-lambda_tld_C*feedback_dBC*c_Tldt*c_Ct;
    double dETld_dN = 0;
    double dETld_dBN = 0;
    double dETld_dS = 0;//-lambda_tld_BC*feedback_dS*c_Tldt*c_BCt-lambda_tld_C*feedback_dS*c_Tldt*c_Ct;
    double dETld_dTld = - decTld;//-lambda_tld_BC*feedback*c_BCt-lambda_tld_C*feedback*c_Ct;

    //std::cout<<"dETld_dTld="<<dETld_dTld<<"\n";

    //---------------------Set Residual and K matrix------------------------------------------------------------------------

    // initialize residuals
    Vector3d Res_B;
    Vector3d Res_C;
    Vector3d Res_BC;
    Vector3d Res_N;
    Vector3d Res_BN;
    Vector3d Res_S;
    Vector3d Res_Tld;

    // initialize de knn matrixs
    Matrix3d K11; K11.setZero();
    Matrix3d K22; K22.setZero();
    Matrix3d K33; K33.setZero(); 
    Matrix3d K44; K44.setZero();
    Matrix3d K55; K55.setZero();
    Matrix3d K66; K66.setZero();
    Matrix3d K77; K77.setZero();

    Matrix3d K12; K12.setZero(); 
    Matrix3d K13; K13.setZero();
    Matrix3d K14; K14.setZero();  
    Matrix3d K15; K15.setZero();
    Matrix3d K16; K16.setZero(); 
    Matrix3d K17; K17.setZero(); 
    
    Matrix3d K21; K21.setZero(); 
    Matrix3d K23; K23.setZero();
    Matrix3d K24; K24.setZero();  
    Matrix3d K25; K25.setZero();
    Matrix3d K26; K26.setZero(); 
    Matrix3d K27; K27.setZero();

    Matrix3d K31; K31.setZero(); 
    Matrix3d K32; K32.setZero();
    Matrix3d K34; K34.setZero();  
    Matrix3d K35; K35.setZero();
    Matrix3d K36; K36.setZero();
    Matrix3d K37; K37.setZero();

    Matrix3d K41; K41.setZero(); 
    Matrix3d K42; K42.setZero();
    Matrix3d K43; K43.setZero();  
    Matrix3d K45; K45.setZero();
    Matrix3d K46; K46.setZero();
    Matrix3d K47; K47.setZero();

    Matrix3d K51; K51.setZero(); 
    Matrix3d K52; K52.setZero();
    Matrix3d K53; K53.setZero();  
    Matrix3d K54; K54.setZero();
    Matrix3d K56; K56.setZero();
    Matrix3d K57; K57.setZero();

    Matrix3d K61; K61.setZero(); 
    Matrix3d K62; K62.setZero();
    Matrix3d K63; K63.setZero();  
    Matrix3d K64; K64.setZero();
    Matrix3d K65; K65.setZero();
    Matrix3d K67; K67.setZero();

    Matrix3d K71; K71.setZero(); 
    Matrix3d K72; K72.setZero();
    Matrix3d K73; K73.setZero();  
    Matrix3d K74; K74.setZero();
    Matrix3d K75; K75.setZero();
    Matrix3d K76; K76.setZero();


    // loop for residual for node 
    for (int j=0; j<3; j++)
    {

    // Set up residuals
    Res_B(j)   = -N(j)*c_B0/dt*DetJ   +N(j)*c_Bt/dt*DetJ   -N(j)*(Esource_B)*DetJ;   
    Res_C(j)   = -N(j)*c_C0/dt*DetJ   +N(j)*c_Ct/dt*DetJ   -N(j)*(Esource_C)*DetJ;      
    Res_BC(j)  = -N(j)*c_BC0/dt*DetJ  +N(j)*c_BCt/dt*DetJ  -N(j)*(Esource_BC)*DetJ;     
    Res_N(j)   = -N(j)*c_N0/dt*DetJ   +N(j)*c_Nt/dt*DetJ   -N(j)*(Esource_N)*DetJ;  
    Res_BN(j)  = -N(j)*c_BN0/dt*DetJ  +N(j)*c_BNt/dt*DetJ  -N(j)*(Esource_BN)*DetJ;
    Res_S(j)   = -N(j)*c_S0/dt*DetJ   +N(j)*c_St/dt*DetJ   -N(j)*(Esource_S)*DetJ;
    Res_Tld(j) = -N(j)*c_Tld0/dt*DetJ +N(j)*c_Tldt/dt*DetJ -N(j)*(Esource_Tld)*DetJ;

    // setup for Kmm matrix
        for (int k=0; k<3; k++)
        {
            // diaganal elements--------------------------------------------
            K11(k,j) = N(k)*N(j)*DetJ/dt-N(k)*(dEB_dB)*N(j)*DetJ+ke_B(k,j);
            K22(k,j) = N(k)*N(j)*DetJ/dt-N(k)*(dEC_dC)*N(j)*DetJ+ke_C(k,j);
            K33(k,j) = N(k)*N(j)*DetJ/dt-N(k)*(dEBC_dBC)*N(j)*DetJ+ke_BC(k,j);
            K44(k,j) = N(k)*N(j)*DetJ/dt-N(k)*(dEN_dN)*N(j)*DetJ+ke_N(k,j);
            K55(k,j) = N(k)*N(j)*DetJ/dt-N(k)*(dEBN_dBN)*N(j)*DetJ+ke_BN(k,j);
            K66(k,j) = N(k)*N(j)*DetJ/dt-N(k)*(dES_dS)*N(j)*DetJ+ke_S(k,j);
            K77(k,j) = N(k)*N(j)*DetJ/dt-N(k)*(dETld_dTld)*N(j)*DetJ+ke_Tld(k,j);
            
            // Coupled elements--------------------------------------------
            //-----------------BMP related---------------------------------
            K12(k,j) = - N(k)*(dEB_dC)*N(j)*DetJ;
            K13(k,j) = - N(k)*(dEB_dBC)*N(j)*DetJ;
            K14(k,j) = - N(k)*(dEB_dN)*N(j)*DetJ;
            K15(k,j) = - N(k)*(dEB_dBN)*N(j)*DetJ;
            K16(k,j) = - N(k)*(dEB_dS)*N(j)*DetJ;
            K17(k,j) = - N(k)*(dEB_dTld)*N(j)*DetJ;

            //-----------------Chd related---------------------------------
            K21(k,j) = - N(k)*(dEC_dB)*N(j)*DetJ;
            K23(k,j) = - N(k)*(dEC_dBC)*N(j)*DetJ;
            K24(k,j) = - N(k)*(dEC_dN)*N(j)*DetJ;
            K25(k,j) = - N(k)*(dEC_dBN)*N(j)*DetJ;
            K26(k,j) = - N(k)*(dEC_dS)*N(j)*DetJ;
            K27(k,j) = - N(k)*(dEC_dTld)*N(j)*DetJ;

            //-----------------BC related---------------------------------
            K31(k,j) = - N(k)*(dEBC_dB)*N(j)*DetJ;
            K32(k,j) = - N(k)*(dEBC_dC)*N(j)*DetJ;
            K34(k,j) = - N(k)*(dEBC_dN)*N(j)*DetJ;
            K35(k,j) = - N(k)*(dEBC_dBN)*N(j)*DetJ;
            K36(k,j) = - N(k)*(dEBC_dS)*N(j)*DetJ;
            K37(k,j) = - N(k)*(dEBC_dTld)*N(j)*DetJ;

            //-----------------Nog related---------------------------------
            K41(k,j) = - N(k)*(dEN_dB)*N(j)*DetJ;
            K42(k,j) = - N(k)*(dEN_dC)*N(j)*DetJ;
            K43(k,j) = - N(k)*(dEN_dBC)*N(j)*DetJ;
            K45(k,j) = - N(k)*(dEN_dBN)*N(j)*DetJ;
            K46(k,j) = - N(k)*(dEN_dS)*N(j)*DetJ;
            K47(k,j) = - N(k)*(dEN_dTld)*N(j)*DetJ;

            //-----------------BN related---------------------------------
            K51(k,j) = - N(k)*(dEBN_dB)*N(j)*DetJ;
            K52(k,j) = - N(k)*(dEBN_dC)*N(j)*DetJ;
            K53(k,j) = - N(k)*(dEBN_dBC)*N(j)*DetJ;
            K54(k,j) = - N(k)*(dEBN_dN)*N(j)*DetJ;
            K56(k,j) = - N(k)*(dEBN_dS)*N(j)*DetJ;
            K57(k,j) = - N(k)*(dEBN_dTld)*N(j)*DetJ;

            //-----------------Sizzled related---------------------------------
            K61(k,j) = - N(k)*(dES_dB)*N(j)*DetJ;
            K62(k,j) = - N(k)*(dES_dC)*N(j)*DetJ;
            K63(k,j) = - N(k)*(dES_dBC)*N(j)*DetJ;
            K64(k,j) = - N(k)*(dES_dN)*N(j)*DetJ;
            K65(k,j) = - N(k)*(dES_dBN)*N(j)*DetJ;
            K67(k,j) = - N(k)*(dES_dTld)*N(j)*DetJ;

            //-----------------Sizzled related---------------------------------
            K71(k,j) = - N(k)*(dETld_dB)*N(j)*DetJ;
            K72(k,j) = - N(k)*(dETld_dC)*N(j)*DetJ;
            K73(k,j) = - N(k)*(dETld_dBC)*N(j)*DetJ;
            K74(k,j) = - N(k)*(dETld_dN)*N(j)*DetJ;
            K75(k,j) = - N(k)*(dETld_dBN)*N(j)*DetJ;
            K76(k,j) = - N(k)*(dETld_dS)*N(j)*DetJ;
        }     
    }//end of loop for residuals

    Res_B = Res_B + ke_B*Uc_Bt;
    Res_C = Res_C+ke_C*Uc_Ct;
    Res_BC = Res_BC+ke_BC*Uc_BCt;
    Res_N = Res_N+ke_N*Uc_Nt;
    Res_BN =Res_BN+ke_BN*Uc_BNt;
    Res_S =Res_S+ke_S*Uc_St;
    Res_Tld =Res_Tld+ke_Tld*Uc_Tldt;

    //std::cout<<"Res_Tld="<<Res_Tld<<"\n";
    //std::cout<<"K76="<<K76<<"\n";
    //std::cout<<"total_species*nodeN="<<total_species*nodeN<<"\n";
    
    MatrixXd K_e(total_species*nodeN,total_species*nodeN); 
    K_e.block<3,3>(0,0) = K11;
    K_e.block<3,3>(0,3) = K12;
    K_e.block<3,3>(0,6) = K13;
    K_e.block<3,3>(0,9) = K14;
    K_e.block<3,3>(0,12) = K15;
    K_e.block<3,3>(0,15) = K16;
    K_e.block<3,3>(0,18) = K17;

    K_e.block<3,3>(3,0) = K21;
    K_e.block<3,3>(3,3) = K22;
    K_e.block<3,3>(3,6) = K23;
    K_e.block<3,3>(3,9) = K24;
    K_e.block<3,3>(3,12) = K25;
    K_e.block<3,3>(3,15) = K26;
    K_e.block<3,3>(3,18) = K27;

    K_e.block<3,3>(6,0) = K31;
    K_e.block<3,3>(6,3) = K32;
    K_e.block<3,3>(6,6) = K33; 
    K_e.block<3,3>(6,9) = K34;
    K_e.block<3,3>(6,12) = K35; 
    K_e.block<3,3>(6,15) = K36; 
    K_e.block<3,3>(6,18) = K37;

    K_e.block<3,3>(9,0) = K41;
    K_e.block<3,3>(9,3) = K42;
    K_e.block<3,3>(9,6) = K43; 
    K_e.block<3,3>(9,9) = K44;
    K_e.block<3,3>(9,12) = K45;
    K_e.block<3,3>(9,15) = K46;
    K_e.block<3,3>(9,18) = K47;

    K_e.block<3,3>(12,0) = K51;
    K_e.block<3,3>(12,3) = K52;
    K_e.block<3,3>(12,6) = K53; 
    K_e.block<3,3>(12,9) = K54;
    K_e.block<3,3>(12,12) = K55;
    K_e.block<3,3>(12,15) = K56;
    K_e.block<3,3>(12,18) = K57;

    K_e.block<3,3>(15,0) = K61;
    K_e.block<3,3>(15,3) = K62;
    K_e.block<3,3>(15,6) = K63; 
    K_e.block<3,3>(15,9) = K64;
    K_e.block<3,3>(15,12) = K65;
    K_e.block<3,3>(15,15) = K66;
    K_e.block<3,3>(15,18) = K67;

    K_e.block<3,3>(18,0) = K71;
    K_e.block<3,3>(18,3) = K72;
    K_e.block<3,3>(18,6) = K73; 
    K_e.block<3,3>(18,9) = K74;
    K_e.block<3,3>(18,12) = K75;
    K_e.block<3,3>(18,15) = K76;
    K_e.block<3,3>(18,18) = K77;




    //std::cout<<"ke11=\n"<<K11<<"\n";
    //std::cout<<"ke22=\n"<<K22<<"\n";
    //std::cout<<"ke33=\n"<<K33<<"\n";
    //std::cout<<"ke=\n"<<K_e<<"\n";
    
    // set up the total resu
    VectorXd Res_e(total_species*nodeN);
    Res_e.head<3>() = Res_B;
    Res_e.segment<3>(3) = Res_C;
    Res_e.segment<3>(6) = Res_BC;
    Res_e.segment<3>(9) = Res_N;
    Res_e.segment<3>(12) = Res_BN;
    Res_e.segment<3>(15) = Res_S;
    Res_e.tail<3>() = Res_Tld;
    //std::cout<<"Res_e=\n"<<Res_e<<"\n";  
    // ASSEMBLY
    // loop over species number
    for (int s = 0; s < total_species; s++)
    {
        for (int k = 0; k<3;k++)
        {
        RR(nodeindex(k)+s*nn) += Res_e((s*3)+k); // residual for cB for node 1 of element i
        }
    }


    //std::cout<<"RR=\n"<<RR<<"\n";  
    //int varArray[] = {&node1index, &node2index, &node3index};
    std::vector<int> varArray (3.0);
    varArray[0]=nodeindex(0);
    varArray[1]=nodeindex(1);
    varArray[2]=nodeindex(2);

    // assmbly form KK matrix 
    for(int a=0;a<3;a++)
    {
      for(int b=0;b<3;b++)
      {
         for (int s1 = 0; s1 < total_species; s1++)
        {
            for (int s2 = 0; s2 < total_species; s2++)
            {
            KK_triplets.push_back(T(varArray[a]+nn*s1,varArray[b]+nn*s2,K_e(a+3*s1,b+3*s2)));
            }
        }

       }

    } 
    //std::cout<<"KK_triplets=\n"<<KK_triplets<<"\n";     
}// closing the element loop

    double error;
    error=RR.norm();
    if (isdigit(error) == 0)
    {
        //std::cout<<"isdigit(error)= "<<isdigit(error)<<"\n";
        if(error>tol)
        {
          std::cout<<"error "<<error<<"\n";
          KK.setFromTriplets(KK_triplets.begin(), KK_triplets.end());
          KK.makeCompressed();
          solver.compute(KK);
          dU = solver.solve(-1.*RR);
          for(int i=0;i<nn;i++)
          {
            c_B[i] += dU(i);
            c_C[i] += dU(i+nn);
            c_BC[i] += dU(i+2*nn);
            c_N[i] += dU(i+3*nn);
            c_BN[i] += dU(i+4*nn);
            c_S[i] += dU(i+5*nn);
            c_Tld[i] += dU(i+6*nn);
            //std::cout<<"updated c_B"<<c_B[i]<<"\n";

          }
          //std::cout<<"updated c_B"<<c_B[1]<<"\n";
         }
        else
        {
           std::cout<<"yay!, error "<<error<<"\n";
           
           break;
        }
    }
    else
    {
        it = max_it+3;
    }

    // relocated concentrtion to species matrix

    
// count howmany iterations done

iterations = it;     
} 
if (iterations < max_it)
{
    for (int j=0; j<nn; ++j){
    species[j][0]=c_B[j];
    species[j][1]=c_C[j];
    species[j][2]=c_BC[j];
    species[j][3]=c_N[j];
    species[j][4]=c_BN[j];
    species[j][5]=c_S[j];
    species[j][6]=c_Tld[j];
 }
}else{
    std::cout<<"no convergence of diffusion\n";
}
// closes the Newton Raphson
} // closes the function 



