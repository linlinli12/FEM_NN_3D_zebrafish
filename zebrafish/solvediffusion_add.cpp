
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

//  Linlin Li 2019
//  update note: adding sizzled 

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
using namespace Eigen;

void advance_diffusion(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, std::vector<double> &c_B, std::vector<double> &c_C, std::vector<double> &c_BC, std::vector<double> &c_N,std::vector<double> &c_BN,std::vector<double> &c_Siz,double dt, double t, const std::vector<double> &diffusionParam)
{

//%% get parameters
int ne = mesh.num_triangles();// % get number of element
int nn = x.size(); //            % get total node number
//int nne = 3;//                   % get number of node per element 3 for triangle mesh

int total_parameters = 19;
int total_equations = 6;

std::map<std::string, int> parameters;

string parameters[4] = {"DN","DBC","DBN","j1","j2","j3","k2","k_2","k3","k_3","decN","decBC","decBN","lambda_tld_C","lambda_tld_BC","DS",""};

// --------------constant parameters-----------------
double r   = 350; // redius of embryo ( may need to change if the mesh changes)
double DB  = 4.4; //            % diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
double DC  = 7; //            % diffusion rate of ligand (Chd)    (microns^2*s^-1)*60s*m^-1
double DN;//  = 10; //            % diffusion rate of ligand (Nog)    (microns^2*s^-1)*60s*m^-1
double DBC;// = 4.7475; //            % diffusion rate of ligand (BC)    (microns^2*s^-1)*60s*m^-1
double DBN;// = 7.4747; //            % diffusion rate of ligand (BC)    (microns^2*s^-1)*60s*m^-1


double j1;// = 6.33; //          % production rate of BMP          nM*m^-1
double j2;// = 55.04; //          % production rate of Chd          nM*m^-1
double j3;// = 38.51; //          % production rate of Nog          nM*m^-1

double k2;//  = 0.7539; //      % binding rates for BMP ligand and Chordin          nM^-1*m^-1
double k_2;//  = k2 ; //        % unbinding rates for BMP ligand and Chordin        m^-1
double k3;//   = 0.6044; //      % binding rates for BMP ligand and Noggin          nM^-1*m^-1
double k_3;//  = 0.1*k3; //     % unbinding rates for BMP ligand and Noggin        m^-1


double decB  = 8.7e-5; //  % decay rate of Ligand (BMP)    nM*m^-1
double decC  = 9.6e-5; //  % decay rate of Chd             nM*m^-1
double decN;//  = 9.3658e-4; //  % decay rate of Nog             nM*m^-1
double decBC;// = 4.7279e-4; //  % decay rate of Ligand (BMP)    nM*m^-1
double decBN;// = 6.2045e-4; //  % decay rate of Chd             nM*m^-1

double lambda_tld_C;//  = 6571;  //Tld processing rate of Chd nM^-1*m^-1
double lambda_tld_BC;// = 5353;  //Tld processing rate of BC nM^-1*m^-1
double tld_parameter = 0.2;      // Tld concentration

// Sizzle related parameters
double DS;
double VS;
double n;
double decSiz;

/// replace with the parameters from the script file
DN = diffusionParam[0];
DBC= diffusionParam[1];
DBN= diffusionParam[2];
j1 = diffusionParam[3];
j2 = diffusionParam[4];
j3 = diffusionParam[5];
k2 = diffusionParam[6];
k_2 = diffusionParam[7];
k3 = diffusionParam[8];
k_3 = diffusionParam[9];
decN = diffusionParam[10];
decBC = diffusionParam[11];
decBN = diffusionParam[12];
lambda_tld_C = diffusionParam[13];
lambda_tld_BC = diffusionParam[14];
DS = diffusionParam[15]; 
VS = diffusionParam[16]; 
n = diffusionParam[17];
decSiz = diffusionParam[18];
k4 = diffusionParam[19];

double tol=1e-3;          //  % tolerance for Newton method
double max_it=20;       // max iteration number


// -------------initial conditions-------------------
// initializing the residual vector, stiffness matrix, sparse linear solver
VectorXd RR(nn*6); 
VectorXd dU(nn*6);
SparseMatrix<double, ColMajor> KK(nn*6,nn*6);
std::vector<T> KK_triplets;
SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> >   solver;
 //gass value for the current concetration
  std::vector<double> c_B0 = c_B;
  std::vector<double> c_C0 = c_C;
  std::vector<double> c_BC0 = c_BC;
  std::vector<double> c_N0 = c_N;
  std::vector<double> c_BN0 = c_BN;
  std::vector<double> c_Siz0 = c_Siz;

double t47 = 4320;
double t57 = 7920;
//double t63 = 10080;


int currframe;
std::string filename;

if (t<t47){ currframe = 47;   
filename="data/BMPmap_47.txt";
}  
else if (t<t57){currframe = 57;   
filename="data/BMPmap_57.txt";
}        
else{ currframe = 63;   
filename="data/BMPmap_63.txt";    
}

// std::cout<<"frame "<<currframe<<"\n";
// std::cout<<"time "<<t<<"\n";
// std::cout<<"filename "<<filename<<"\n";

double data[10201];

std::string line;
std::ifstream bmpfile(filename.c_str());

for (int i = 0; i < 10201; i++) {
    bmpfile >> data[i];

    }


//% iteration for  Newton-Raphson method
//l------------------oop for iteration-----------------------
 for (int it=0;it<max_it;it++)
{    
  // clear residual and stiffness for iteration
  KK_triplets.clear();
  RR.setZero();

 

// --------------loop for element start
// this code is for solve the diffusion reaction equatio
  for(int i=0;i<ne;i++){
  // for triangle i   
// ---------------------mesh information---------
    int node1index = mesh.get_triangle(i)[0];
    int node2index = mesh.get_triangle(i)[1];
    int node3index = mesh.get_triangle(i)[2];
   // std::cout <<"mesh nodes1 = "<<node1index<<"\n";
   // std::cout <<"mesh nodes2 = "<<node2index<<"\n";
    //std::cout <<"mesh nodes3 = "<<node3index<<"\n";
// get coordinates
// coords of node 1
    double node1x = x[node1index][0];
    double node1y = x[node1index][1];
    double node1z = x[node1index][2];

    // std::cout <<"node1x = "<<node1x<<"\n";
    // std::cout <<"node1y = "<<node1y<<"\n";
    // std::cout <<"node1z = "<<node1z<<"\n";
    // coords of node 2
    double node2x = x[node2index][0];
    double node2y = x[node2index][1];
    double node2z = x[node2index][2];
    // std::cout <<"node2x = "<<node2x<<"\n";
    // std::cout <<"node2y = "<<node2y<<"\n";
    // std::cout <<"node2z = "<<node2z<<"\n";

     // coords of node 3
    double node3x = x[node3index][0];
    double node3y = x[node3index][1];
    double node3z = x[node3index][2];

    // std::cout <<"node3x = "<<node3x<<"\n";
    // std::cout <<"node3y = "<<node3y<<"\n";
    // std::cout <<"node3z = "<<node3z<<"\n";
    // get concentrations
    // node 1 conc
    Vector3d Uc_BMP0;
    Uc_BMP0(0) = c_B0[node1index];
    Uc_BMP0(1) = c_B0[node2index];
    Uc_BMP0(2) = c_B0[node3index];
    
    Vector3d Uc_CHD0;
    Uc_CHD0(0) = c_C0[node1index];
    Uc_CHD0(1) = c_C0[node2index];
    Uc_CHD0(2) = c_C0[node3index];

    Vector3d Uc_BC0;
    Uc_BC0(0) = c_BC0[node1index];
    Uc_BC0(1) = c_BC0[node2index];
    Uc_BC0(2) = c_BC0[node3index];

    Vector3d Uc_NOG0;
    Uc_NOG0(0) = c_N0[node1index];
    Uc_NOG0(1) = c_N0[node2index];
    Uc_NOG0(2) = c_N0[node3index];

    Vector3d Uc_BN0;
    Uc_BN0(0) = c_BN0[node1index];
    Uc_BN0(1) = c_BN0[node2index];
    Uc_BN0(2) = c_BN0[node3index];

    Vector3d Uc_Siz0;
    Uc_Siz0(0) = c_Siz0[node1index];
    Uc_Siz0(1) = c_Siz0[node2index];
    Uc_Siz0(2) = c_Siz0[node3index];


    // node 2 conc

    Vector3d Uc_BMPt;
    Uc_BMPt(0) = c_B[node1index];
    Uc_BMPt(1) = c_B[node2index];
    Uc_BMPt(2) = c_B[node3index];
    
    Vector3d Uc_CHDt;
    Uc_CHDt(0) = c_C[node1index];
    Uc_CHDt(1) = c_C[node2index];
    Uc_CHDt(2) = c_C[node3index];

    Vector3d Uc_BCt;
    Uc_BCt(0) = c_BC[node1index];
    Uc_BCt(1) = c_BC[node2index];
    Uc_BCt(2) = c_BC[node3index];

    Vector3d Uc_NOGt;
    Uc_NOGt(0) = c_N[node1index];
    Uc_NOGt(1) = c_N[node2index];
    Uc_NOGt(2) = c_N[node3index];

    Vector3d Uc_BNt;
    Uc_BNt(0) = c_BN[node1index];
    Uc_BNt(1) = c_BN[node2index];
    Uc_BNt(2) = c_BN[node3index];
    
    Vector3d Uc_BNt;
    Uc_Sizt(0) = c_Siz[node1index];
    Uc_Sizt(1) = c_Siz[node2index];
    Uc_Sizt(2) = c_Siz[node3index];

    //std::cout<<"element "<<i<<"\n";
    //std::cout<<"Uc_BMP0=\n"<<Uc_BMP0<<"\n";
    //std::cout<<"Uc_BMPt=\n"<<Uc_BMPt<<"\n";
    //std::cout<<"Uc_CHD0=\n"<<Uc_CHD0<<"\n";
    //std::cout<<"Uc_CHDt=\n"<<Uc_CHDt<<"\n";
    //std::cout<<"Uc_BC0=\n"<<Uc_BC0<<"\n";
    //std::cout<<"Uc_BCt=\n"<<Uc_BCt<<"\n";
    //std::cout<<"current time="<<t<<"\n";

// // get angles for sphere coodinate (formate matching with matlab)



int meshN = 101;
double pi = 3.14159265;

double az1 = atan2(node1y,node1x);
double el1 = atan2(node1z,sqrt(node1x*node1x + node1y*node1y));
//std::cout<<"az1 "<<az1<<"\n";
//std::cout<<"el1 "<<el1<<"\n";

double az2 = atan2(node2y,node2x);
double el2 = atan2(node2z,sqrt(node2x*node2x + node2y*node2y));

double az3 = atan2(node3y,node3x);
double el3 = atan2(node3z,sqrt(node3x*node3x + node3y*node3y));

int index1 = round(((az1/pi+1)/2*(meshN-1)))*(meshN)+round((el1/(pi/2)+1)/2*(meshN-1));
int index2 = round(((az2/pi+1)/2*(meshN-1)))*(meshN)+round((el2/(pi/2)+1)/2*(meshN-1));
int index3 = round(((az3/pi+1)/2*(meshN-1)))*(meshN)+round((el3/(pi/2)+1)/2*(meshN-1));

// std::cout<<"index1 "<<index1<<"\n";
// std::cout<<"index2 "<<index2<<"\n";
// std::cout<<"index3 "<<index3<<"\n";


double bmpmap1=data[index1];
double bmpmap2=data[index2];
double bmpmap3=data[index3];

// std::cout<<  "bmp1 "<<bmpmap1<<std::endl;
// std::cout<<  "bmp2 "<<bmpmap2<<std::endl;
// std::cout<<  "bmp3 "<<bmpmap3<<std::endl;

// if(node1x<=0){bmpmap1=1;}
//     else {bmpmap1=0;}
// if(node2x<=0){bmpmap2=1;}
//     else {bmpmap2=0;}
// if(node3x<=0){bmpmap3=1;}
//     else {bmpmap3=0;}


    Vector3d s_BMP;
    // simple expression source
     s_BMP(0)=j1*bmpmap1;
     s_BMP(1)=j1*bmpmap2;
     s_BMP(2)=j1*bmpmap3;



    //std::cout<<s_BMP<<"\n";
//set source term for CHD
    Vector3d s_CHD;
    //s_CHD(0)=j2*pow((node1x+r)/(2*r),50)/(pow((0.6818*t+2477)/(t+2477),50)+pow((node1x+r)/(2*r),50));//%initialize source terms for chd !!!!!! modify function by chang of x
    //s_CHD(1)=j2*pow((node2x+r)/(2*r),50)/(pow((0.6818*t+2477)/(t+2477),50)+pow((node2x+r)/(2*r),50));//%initialize source terms for chd !!!!!! modify function by chang of x
    //s_CHD(2)=j2*pow((node3x+r)/(2*r),50)/(pow((0.6818*t+2477)/(t+2477),50)+pow((node3x+r)/(2*r),50));//%initialize source terms for chd !!!!!! modify function by chang of x
    //------------------modified for not changing by time
    s_CHD(0)=j2*pow((node1x+r)/(2*r),50)/(pow(0.8,50)+pow((node1x+r)/(2*r),50));
    s_CHD(1)=j2*pow((node2x+r)/(2*r),50)/(pow(0.8,50)+pow((node2x+r)/(2*r),50));
    s_CHD(2)=j2*pow((node3x+r)/(2*r),50)/(pow(0.8,50)+pow((node3x+r)/(2*r),50));
    
    //set source term for NOG
    Vector3d s_NOG;
    //s_NOG(0)=j3*pow((node1x+r)/(2*r),50)/(pow((0.72*t+3543)/(t+3541),50)+pow((node1x+r)/(2*r),50));//%initialize source terms for chd !!!!!! modify function by chang of x
    //s_NOG(1)=j3*pow((node2x+r)/(2*r),50)/(pow((0.72*t+3543)/(t+3541),50)+pow((node2x+r)/(2*r),50));//%initialize source terms for chd !!!!!! modify function by chang of x
    //s_NOG(2)=j3*pow((node3x+r)/(2*r),50)/(pow((0.72*t+3543)/(t+3541),50)+pow((node3x+r)/(2*r),50));//%initialize source terms for chd !!!!!! modify function by chang of x
    s_NOG(0)=j3*pow((node1x+r)/(2*r),50)/(pow(0.8,50)+pow((node1x+r)/(2*r),50));
    s_NOG(1)=j3*pow((node2x+r)/(2*r),50)/(pow(0.8,50)+pow((node2x+r)/(2*r),50));
    s_NOG(2)=j3*pow((node3x+r)/(2*r),50)/(pow(0.8,50)+pow((node3x+r)/(2*r),50));  



    //double tld_conc=(t/(2000+t))*tld_parameter;
    double tld_conc=tld_parameter;
      
    Vector3d  dNxi;
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
    g1(0)=dNxi(0)*node1x+dNxi(1)*node2x+dNxi(2)*node3x;
    g1(1)=dNxi(0)*node1y+dNxi(1)*node2y+dNxi(2)*node3y;
    g1(2)=dNxi(0)*node1z+dNxi(1)*node2z+dNxi(2)*node3z;
    //std::cout<<"g1\n"<<g1<<"\n";

    Vector3d g2;
    g2(0)=dNeta(0)*node1x+dNeta(1)*node2x+dNeta(2)*node3x;
    g2(1)=dNeta(0)*node1y+dNeta(1)*node2y+dNeta(2)*node3y;
    g2(2)=dNeta(0)*node1z+dNeta(1)*node2z+dNeta(2)*node3z;
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
    px(0,0)=node1x;
    px(1,0)=node1y;
    px(2,0)=node1z;
    px(0,1)=node2x;
    px(1,1)=node2y;
    px(2,1)=node2z;
    px(0,2)=node3x;
    px(1,2)=node3y;
    px(2,2)=node3z;

    // get new coordinate
    Matrix3d poiprime = A.inverse()*px;

    // took the first two row as the local coordinates for  surface diffusion
    MatrixXd poiprime2D(2,3);
    poiprime2D.row(0)=poiprime.row(0);
    poiprime2D.row(1)=poiprime.row(1);
    //std::cout<<"poiprime2D\n"<<poiprime2D<<"\n";
    // get the element concentraion for current time point
    double c_BMP0 = N.dot(Uc_BMP0);
    double c_CHD0 = N.dot(Uc_CHD0);
    double c_BC0 = N.dot(Uc_BC0);
    double c_NOG0 = N.dot(Uc_NOG0);
    double c_BN0 = N.dot(Uc_BN0);
    double c_Siz0 = N.dot(Uc_Siz0);
    // get the element concentraion for test concentraion  
    double c_BMPt = N.dot(Uc_BMPt);
    double c_CHDt = N.dot(Uc_CHDt);
    double c_BCt = N.dot(Uc_BCt);
    double c_NOGt = N.dot(Uc_NOGt);
    double c_BNt = N.dot(Uc_BNt);
    double c_Sizt = N.dot(Uc_Sizt);  

    Matrix2d J = GradientN*poiprime2D.transpose(); //% Jacobian Matrix
    double DetJ = J.determinant();//% Jacobian determinate
    //std::cout<<"DetJ "<<DetJ<<"\n";
    MatrixXd B = J.inverse()*GradientN;
    //std::cout<<"B "<<B<<"\n";
       
    Matrix2d eye = MatrixXd::Identity(2,2);
    Matrix2d D_BMP=DB*eye;//                 % apply diffusion rate
    Matrix2d D_CHD=DC*eye;//                 % apply diffusion rate
    Matrix2d D_BC=DBC*eye;//                 % apply diffusion rate
    Matrix2d D_NOG=DN*eye;//                 % apply diffusion rate
    Matrix2d D_BN=DBN*eye;//                 % apply diffusion rate
    Matrix2d D_Siz=DS*eye;//                 % apply diffusion rate

    Matrix3d ke_BMP=B.transpose()*D_BMP*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_CHD=B.transpose()*D_CHD*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_BC=B.transpose()*D_BC*B*DetJ;             //% element stiffness matrix
    Matrix3d ke_NOG=B.transpose()*D_NOG*B*DetJ;          //% element stiffness matrix
    Matrix3d ke_BN=B.transpose()*D_BN*B*DetJ;             //% element stiffness matrix
    Matrix3d ke_Siz=B.transpose()*D_Siz*B*DetJ;  

    //std::cout<<"ke_BMB "<<ke_BMP<<"\n";
    //std::cout<<"ke_CHD "<<ke_CHD<<"\n";
    //std::cout<<"ke_BC "<<ke_BC<<"\n";
    double se_BMP=N.dot(s_BMP);//                % element source term of BMP
    double se_CHD=N.dot(s_CHD);//                % element source term of CHD
    double se_NOG=N.dot(s_NOG);//                % element source term of BMP


    //std::cout<<"se_BMP="<<se_BMP<<"\n";  
    //std::cout<<"se_CHD="<<se_CHD<<"\n";     
    //---------------------Set Residual-------------------------------------------------------------------------

    Vector3d Res_BME;
    Res_BME(0) = -N(0)*c_BMP0/dt*DetJ +N(0)*c_BMPt/dt*DetJ -N(0)*(se_BMP-k2*c_BMPt*c_CHDt+k_2*c_BCt-decB*c_BMPt-k3*c_BMPt*c_NOGt+k_3*c_BNt+lambda_tld_BC*tld_conc*c_BCt)*DetJ;
    Res_BME(1) = -N(1)*c_BMP0/dt*DetJ +N(1)*c_BMPt/dt*DetJ -N(1)*(se_BMP-k2*c_BMPt*c_CHDt+k_2*c_BCt-decB*c_BMPt-k3*c_BMPt*c_NOGt+k_3*c_BNt+lambda_tld_BC*tld_conc*c_BCt)*DetJ;
    Res_BME(2) = -N(2)*c_BMP0/dt*DetJ +N(2)*c_BMPt/dt*DetJ -N(2)*(se_BMP-k2*c_BMPt*c_CHDt+k_2*c_BCt-decB*c_BMPt-k3*c_BMPt*c_NOGt+k_3*c_BNt+lambda_tld_BC*tld_conc*c_BCt)*DetJ;
    Res_BME = Res_BME + ke_BMP*Uc_BMPt;

    Vector3d Res_CHD;
    Res_CHD(0) = -N(0)*c_CHD0/dt*DetJ +N(0)*c_CHDt/dt*DetJ -N(0)*(se_CHD-k2*c_BMPt*c_CHDt+k_2*c_BCt-decC*c_CHDt-lambda_tld_C*tld_conc*c_CHDt)*DetJ;
    Res_CHD(1) = -N(1)*c_CHD0/dt*DetJ +N(1)*c_CHDt/dt*DetJ -N(1)*(se_CHD-k2*c_BMPt*c_CHDt+k_2*c_BCt-decC*c_CHDt-lambda_tld_C*tld_conc*c_CHDt)*DetJ;
    Res_CHD(2) = -N(2)*c_CHD0/dt*DetJ +N(2)*c_CHDt/dt*DetJ -N(2)*(se_CHD-k2*c_BMPt*c_CHDt+k_2*c_BCt-decC*c_CHDt-lambda_tld_C*tld_conc*c_CHDt)*DetJ;
    Res_CHD=Res_CHD+ke_CHD*Uc_CHDt;

    Vector3d Res_BC;
    Res_BC(0) = -N(0)*c_BC0/dt*DetJ +N(0)*c_BCt/dt*DetJ -N(0)*(k2*c_BMPt*c_CHDt-k_2*c_BCt-lambda_tld_BC*tld_conc*c_BCt-decBC*c_BCt)*DetJ;
    Res_BC(1) = -N(1)*c_BC0/dt*DetJ +N(1)*c_BCt/dt*DetJ -N(1)*(k2*c_BMPt*c_CHDt-k_2*c_BCt-lambda_tld_BC*tld_conc*c_BCt-decBC*c_BCt)*DetJ;
    Res_BC(2) = -N(2)*c_BC0/dt*DetJ +N(2)*c_BCt/dt*DetJ -N(2)*(k2*c_BMPt*c_CHDt-k_2*c_BCt-lambda_tld_BC*tld_conc*c_BCt-decBC*c_BCt)*DetJ;
    Res_BC =Res_BC+ke_BC*Uc_BCt;

    Vector3d Res_NOG;
    Res_NOG(0) = -N(0)*c_NOG0/dt*DetJ +N(0)*c_NOGt/dt*DetJ -N(0)*(se_NOG-k3*c_BMPt*c_NOGt+k_3*c_BNt-decN*c_NOGt)*DetJ;
    Res_NOG(1) = -N(1)*c_NOG0/dt*DetJ +N(1)*c_NOGt/dt*DetJ -N(1)*(se_NOG-k3*c_BMPt*c_NOGt+k_3*c_BNt-decN*c_NOGt)*DetJ;
    Res_NOG(2) = -N(2)*c_NOG0/dt*DetJ +N(2)*c_NOGt/dt*DetJ -N(2)*(se_NOG-k3*c_BMPt*c_NOGt+k_3*c_BNt-decN*c_NOGt)*DetJ;
    Res_NOG=Res_NOG+ke_NOG*Uc_NOGt;

    Vector3d Res_BN;
    Res_BN(0) = -N(0)*c_BN0/dt*DetJ +N(0)*c_BNt/dt*DetJ -N(0)*(k3*c_BMPt*c_NOGt-k_3*c_BNt-decBN*c_BNt)*DetJ;
    Res_BN(1) = -N(1)*c_BN0/dt*DetJ +N(1)*c_BNt/dt*DetJ -N(1)*(k3*c_BMPt*c_NOGt-k_3*c_BNt-decBN*c_BNt)*DetJ;
    Res_BN(2) = -N(2)*c_BN0/dt*DetJ +N(2)*c_BNt/dt*DetJ -N(2)*(k3*c_BMPt*c_NOGt-k_3*c_BNt-decBN*c_BNt)*DetJ;
    Res_BN =Res_BN+ke_BN*Uc_BNt;

    Vector3d Res_Siz;
    Res_Siz(0) = -N(0)*c_Siz0/dt*DetJ +N(0)*c_Sizt/dt*DetJ -N(0)*(VS*pow(c_BMPt,n)/(pow(k4,n)+pow(c_BMPt,n))-decSiz*c_Sizt)*DetJ;
    Res_Siz(1) = -N(1)*c_Siz0/dt*DetJ +N(1)*c_Sizt/dt*DetJ -N(1)*(VS*pow(c_BMPt,n)/(pow(k4,n)+pow(c_BMPt,n))-decSiz*c_Sizt)*DetJ;
    Res_Siz(2) = -N(2)*c_Siz0/dt*DetJ +N(2)*c_Sizt/dt*DetJ -N(2)*(VS*pow(c_BMPt,n)/(pow(k4,n)+pow(c_BMPt,n))-decSiz*c_Sizt)*DetJ;
    Res_Siz =Res_Siz+ke_Siz*Uc_Sizt;



    Matrix3d K11; K11.setZero();
    K11 << N(0)*N(0)*DetJ/dt-N(0)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(0)*DetJ+ke_BMP(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(1)*DetJ+ke_BMP(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(2)*DetJ+ke_BMP(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(0)*DetJ+ke_BMP(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(1)*DetJ+ke_BMP(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(2)*DetJ+ke_BMP(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(0)*DetJ+ke_BMP(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(1)*DetJ+ke_BMP(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k2*c_CHDt-k3*c_NOGt-decB)*N(2)*DetJ+ke_BMP(2,2);

    Matrix3d K22; K22.setZero();  
    K22 << N(0)*N(0)*DetJ/dt-N(0)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(0)*DetJ+ke_CHD(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(1)*DetJ+ke_CHD(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(2)*DetJ+ke_CHD(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(0)*DetJ+ke_CHD(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(1)*DetJ+ke_CHD(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(2)*DetJ+ke_CHD(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(0)*DetJ+ke_CHD(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(1)*DetJ+ke_CHD(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k2*c_BMPt-lambda_tld_C*tld_conc-decC)*N(2)*DetJ+ke_CHD(2,2);
    
    Matrix3d K33; K33.setZero();   
    K33 << N(0)*N(0)*DetJ/dt-N(0)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(0)*DetJ+ke_BC(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(1)*DetJ+ke_BC(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(2)*DetJ+ke_BC(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(0)*DetJ+ke_BC(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(1)*DetJ+ke_BC(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(2)*DetJ+ke_BC(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(0)*DetJ+ke_BC(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(2)*DetJ+ke_BC(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k_2-lambda_tld_BC*tld_conc-decBC)*N(2)*DetJ+ke_BC(2,2);

    Matrix3d K44; K44.setZero();
    K44 << N(0)*N(0)*DetJ/dt-N(0)*(-k3*c_BMPt-decN)*N(0)*DetJ+ke_NOG(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k3*c_BMPt-decN)*N(1)*DetJ+ke_NOG(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k3*c_BMPt-decN)*N(2)*DetJ+ke_NOG(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k3*c_BMPt-decN)*N(0)*DetJ+ke_NOG(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k3*c_BMPt-decN)*N(1)*DetJ+ke_NOG(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k3*c_BMPt-decN)*N(2)*DetJ+ke_NOG(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k3*c_BMPt-decN)*N(0)*DetJ+ke_NOG(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k3*c_BMPt-decN)*N(1)*DetJ+ke_NOG(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k3*c_BMPt-decN)*N(2)*DetJ+ke_NOG(2,2);

    Matrix3d K55; K55.setZero();   
    K55 << N(0)*N(0)*DetJ/dt-N(0)*(-k_3-decBN)*N(0)*DetJ+ke_BN(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-k_3-decBN)*N(1)*DetJ+ke_BN(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-k_3-decBN)*N(2)*DetJ+ke_BN(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-k_3-decBN)*N(0)*DetJ+ke_BN(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-k_3-decBN)*N(1)*DetJ+ke_BN(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-k_3-decBN)*N(2)*DetJ+ke_BN(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-k_3-decBN)*N(0)*DetJ+ke_BN(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-k_3-decBN)*N(2)*DetJ+ke_BN(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-k_3-decBN)*N(2)*DetJ+ke_BN(2,2);

    Matrix3d K66; K66.setZero();   
    K66 << N(0)*N(0)*DetJ/dt-N(0)*(-decSiz)*N(0)*DetJ+ke_Siz(0,0), N(0)*N(1)*DetJ/dt-N(0)*(-decSiz)*N(1)*DetJ+ke_Siz(0,1), N(0)*N(2)*DetJ/dt-N(0)*(-decSiz)*N(2)*DetJ+ke_Siz(0,2),
           N(1)*N(0)*DetJ/dt-N(1)*(-decSiz)*N(0)*DetJ+ke_Siz(1,0), N(1)*N(1)*DetJ/dt-N(1)*(-decSiz)*N(1)*DetJ+ke_Siz(1,1), N(1)*N(2)*DetJ/dt-N(1)*(-decSiz)*N(2)*DetJ+ke_Siz(1,2),
           N(2)*N(0)*DetJ/dt-N(2)*(-decSiz)*N(0)*DetJ+ke_Siz(2,0), N(2)*N(1)*DetJ/dt-N(2)*(-decSiz)*N(2)*DetJ+ke_Siz(2,1), N(2)*N(2)*DetJ/dt-N(2)*(-decSiz)*N(2)*DetJ+ke_Siz(2,2);




    //------------------------BMP equation
    Matrix3d K12; K12.setZero();  
    K12 << N(0)*(-k2*c_BMPt)*N(0)*DetJ,N(0)*(-k2*c_BMPt)*N(1)*DetJ,N(0)*(-k2*c_BMPt)*N(2)*DetJ,
           N(1)*(-k2*c_BMPt)*N(0)*DetJ,N(1)*(-k2*c_BMPt)*N(1)*DetJ,N(1)*(-k2*c_BMPt)*N(2)*DetJ,
           N(2)*(-k2*c_BMPt)*N(0)*DetJ,N(2)*(-k2*c_BMPt)*N(1)*DetJ,N(2)*(-k2*c_BMPt)*N(2)*DetJ;
    K12=-K12;

    Matrix3d K13; K13.setZero();
    K13 << N(0)*(k_2+lambda_tld_BC*tld_conc)*N(0)*DetJ,N(0)*(k_2+lambda_tld_BC*tld_conc)*N(1)*DetJ,N(0)*(k_2+lambda_tld_BC*tld_conc)*N(2)*DetJ,
           N(1)*(k_2+lambda_tld_BC*tld_conc)*N(0)*DetJ,N(1)*(k_2+lambda_tld_BC*tld_conc)*N(1)*DetJ,N(1)*(k_2+lambda_tld_BC*tld_conc)*N(2)*DetJ,
           N(2)*(k_2+lambda_tld_BC*tld_conc)*N(0)*DetJ,N(2)*(k_2+lambda_tld_BC*tld_conc)*N(1)*DetJ,N(2)*(k_2+lambda_tld_BC*tld_conc)*N(2)*DetJ;
    K13=-K13;

    Matrix3d K14; K14.setZero();  
    K14 << N(0)*(-k3*c_BMPt)*N(0)*DetJ,N(0)*(-k3*c_BMPt)*N(1)*DetJ,N(0)*(-k3*c_BMPt)*N(2)*DetJ,
           N(1)*(-k3*c_BMPt)*N(0)*DetJ,N(1)*(-k3*c_BMPt)*N(1)*DetJ,N(1)*(-k3*c_BMPt)*N(2)*DetJ,
           N(2)*(-k3*c_BMPt)*N(0)*DetJ,N(2)*(-k3*c_BMPt)*N(1)*DetJ,N(2)*(-k3*c_BMPt)*N(2)*DetJ;
    K14=-K14;

    Matrix3d K15; K15.setZero();
    K15 << N(0)*(k_3)*N(0)*DetJ,N(0)*(k_3)*N(1)*DetJ,N(0)*(k_3)*N(2)*DetJ,
           N(1)*(k_3)*N(0)*DetJ,N(1)*(k_3)*N(1)*DetJ,N(1)*(k_3)*N(2)*DetJ,
           N(2)*(k_3)*N(0)*DetJ,N(2)*(k_3)*N(1)*DetJ,N(2)*(k_3)*N(2)*DetJ;
    K15=-K15;
    

    //------------------------CHD equation
    Matrix3d K21; K21.setZero();   
    K21 << N(0)*(-k2*c_CHDt)*N(0)*DetJ,N(0)*(-k2*c_CHDt)*N(1)*DetJ,N(0)*(-k2*c_CHDt)*N(2)*DetJ,
           N(1)*(-k2*c_CHDt)*N(0)*DetJ,N(1)*(-k2*c_CHDt)*N(1)*DetJ,N(1)*(-k2*c_CHDt)*N(2)*DetJ,
           N(2)*(-k2*c_CHDt)*N(0)*DetJ,N(2)*(-k2*c_CHDt)*N(1)*DetJ,N(2)*(-k2*c_CHDt)*N(2)*DetJ;
    K21=-K21;

    Matrix3d K23; K23.setZero();
    K23 << N(0)*(k_2)*N(0)*DetJ,N(0)*(k_2)*N(1)*DetJ,N(0)*(k_2)*N(2)*DetJ,
           N(1)*(k_2)*N(0)*DetJ,N(1)*(k_2)*N(1)*DetJ,N(1)*(k_2)*N(2)*DetJ,
           N(2)*(k_2)*N(0)*DetJ,N(2)*(k_2)*N(1)*DetJ,N(2)*(k_2)*N(2)*DetJ;
    K23=-K23;   
    
    Matrix3d K24; K24.setZero();

    Matrix3d K25; K25.setZero();


    //------------------------BC equation
    Matrix3d K31;     
    K31 << N(0)*(k2*c_CHDt)*N(0)*DetJ,N(0)*(k2*c_CHDt)*N(1)*DetJ,N(0)*(k2*c_CHDt)*N(2)*DetJ,
           N(1)*(k2*c_CHDt)*N(0)*DetJ,N(1)*(k2*c_CHDt)*N(1)*DetJ,N(1)*(k2*c_CHDt)*N(2)*DetJ,
           N(2)*(k2*c_CHDt)*N(0)*DetJ,N(2)*(k2*c_CHDt)*N(1)*DetJ,N(2)*(k2*c_CHDt)*N(2)*DetJ; 
    K31=-K31; 
    
    Matrix3d K32;
    K32 << N(0)*(k2*c_BMPt)*N(0)*DetJ,N(0)*(k2*c_BMPt)*N(1)*DetJ,N(0)*(k2*c_BMPt)*N(2)*DetJ,
           N(1)*(k2*c_BMPt)*N(0)*DetJ,N(1)*(k2*c_BMPt)*N(1)*DetJ,N(1)*(k2*c_BMPt)*N(2)*DetJ,
           N(2)*(k2*c_BMPt)*N(0)*DetJ,N(2)*(k2*c_BMPt)*N(1)*DetJ,N(2)*(k2*c_BMPt)*N(2)*DetJ; 
    K32=-K32; 

    Matrix3d K34; K34.setZero();

    Matrix3d K35; K35.setZero();

   //------------------------Nog equation
    Matrix3d K41; K41.setZero();   
    K41 << N(0)*(-k3*c_NOGt)*N(0)*DetJ,N(0)*(-k3*c_NOGt)*N(1)*DetJ,N(0)*(-k3*c_NOGt)*N(2)*DetJ,
           N(1)*(-k3*c_NOGt)*N(0)*DetJ,N(1)*(-k3*c_NOGt)*N(1)*DetJ,N(1)*(-k3*c_NOGt)*N(2)*DetJ,
           N(2)*(-k3*c_NOGt)*N(0)*DetJ,N(2)*(-k3*c_NOGt)*N(1)*DetJ,N(2)*(-k3*c_NOGt)*N(2)*DetJ;
    K41=-K41;

    Matrix3d K42; K42.setZero();

    Matrix3d K43; K43.setZero();

    Matrix3d K45; K45.setZero();
    K45 << N(0)*(k_3)*N(0)*DetJ,N(0)*(k_3)*N(1)*DetJ,N(0)*(k_3)*N(2)*DetJ,
           N(1)*(k_3)*N(0)*DetJ,N(1)*(k_3)*N(1)*DetJ,N(1)*(k_3)*N(2)*DetJ,
           N(2)*(k_3)*N(0)*DetJ,N(2)*(k_3)*N(1)*DetJ,N(2)*(k_3)*N(2)*DetJ;
    K45=-K45; 

   //------------------------BN equation
    Matrix3d K51; K51.setZero();   
    K51 << N(0)*(k3*c_NOGt)*N(0)*DetJ,N(0)*(k3*c_NOGt)*N(1)*DetJ,N(0)*(k3*c_NOGt)*N(2)*DetJ,
           N(1)*(k3*c_NOGt)*N(0)*DetJ,N(1)*(k3*c_NOGt)*N(1)*DetJ,N(1)*(k3*c_NOGt)*N(2)*DetJ,
           N(2)*(k3*c_NOGt)*N(0)*DetJ,N(2)*(k3*c_NOGt)*N(1)*DetJ,N(2)*(k3*c_NOGt)*N(2)*DetJ;
    K51=-K51;

    Matrix3d K52; K52.setZero();

    Matrix3d K53; K53.setZero();

    Matrix3d K54; K54.setZero();
    K54 << N(0)*(k3*c_BMPt)*N(0)*DetJ,N(0)*(k3*c_BMPt)*N(1)*DetJ,N(0)*(k3*c_BMPt)*N(2)*DetJ,
           N(1)*(k3*c_BMPt)*N(0)*DetJ,N(1)*(k3*c_BMPt)*N(1)*DetJ,N(1)*(k3*c_BMPt)*N(2)*DetJ,
           N(2)*(k3*c_BMPt)*N(0)*DetJ,N(2)*(k3*c_BMPt)*N(1)*DetJ,N(2)*(k3*c_BMPt)*N(2)*DetJ;
    K54=-K54;

    //-------------------------Sizzled equation

    Matrix3d K61; K61.setZero();   
    K61 << 
    K61=-K61;

    Matrix3d K62; K62.setZero();

    Matrix3d K63; K63.setZero();

    Matrix3d K64; K64.setZero();
    K64 << 
    K64=-K64;

    Matrix3d K65; K65.setZero();

    // add Sizzled equaiton here 





      
    MatrixXd K_e(18,18); 
    K_e.block<3,3>(0,0) = K11;
    K_e.block<3,3>(0,3) = K12;
    K_e.block<3,3>(0,6) = K13;
    K_e.block<3,3>(0,9) = K14;
    K_e.block<3,3>(0,12) = K15;
    K_e.block<3,3>(0,15) = K16;

    K_e.block<3,3>(3,0) = K21;
    K_e.block<3,3>(3,3) = K22;
    K_e.block<3,3>(3,6) = K23;
    K_e.block<3,3>(3,9) = K24;
    K_e.block<3,3>(3,12) = K25;
    K_e.block<3,3>(3,15) = K26;

    K_e.block<3,3>(6,0) = K31;
    K_e.block<3,3>(6,3) = K32;
    K_e.block<3,3>(6,6) = K33; 
    K_e.block<3,3>(6,9) = K34;
    K_e.block<3,3>(6,12) = K35; 
    K_e.block<3,3>(6,15) = K36; 

    K_e.block<3,3>(9,0) = K41;
    K_e.block<3,3>(9,3) = K42;
    K_e.block<3,3>(9,6) = K43; 
    K_e.block<3,3>(9,9) = K44;
    K_e.block<3,3>(9,12) = K45;
    K_e.block<3,3>(9,15) = K46;

    K_e.block<3,3>(12,0) = K51;
    K_e.block<3,3>(12,3) = K52;
    K_e.block<3,3>(12,6) = K53; 
    K_e.block<3,3>(12,9) = K54;
    K_e.block<3,3>(12,12) = K55;
    K_e.block<3,3>(12,15) = K56;

    K_e.block<3,3>(15,0) = K61;
    K_e.block<3,3>(15,3) = K62;
    K_e.block<3,3>(15,6) = K63; 
    K_e.block<3,3>(15,9) = K64;
    K_e.block<3,3>(15,12) = K65;
    K_e.block<3,3>(15,15) = K66;


    //std::cout<<"ke11=\n"<<K11<<"\n";
    //std::cout<<"ke22=\n"<<K22<<"\n";
    //std::cout<<"ke33=\n"<<K33<<"\n";
    //std::cout<<"ke=\n"<<K_e<<"\n";
    

    VectorXd Res_e(18);
    Res_e.head<3>() = Res_BME;
    Res_e.segment<3>(3) = Res_CHD;
    Res_e.segment<3>(6) = Res_BC;
    Res_e.segment<3>(9) = Res_NOG;
    Res_e.segment<3>(12) = Res_BN;
    Res_e.tail<3>() = Res_Siz;
    //std::cout<<"Res_e=\n"<<Res_e<<"\n";  
// ASSEMBLY
    RR(node1index) += Res_e(0); // residual for cB for node 1 of element i
    RR(node2index) += Res_e(1); // residual for cB for node 2 of element i
    RR(node3index) += Res_e(2); // residual for cB for node 3 of element i

    RR(node1index+nn) += Res_e(3); // residual for cC for node 1 of element i
    RR(node2index+nn) += Res_e(4); // residual for cC for node 2 of element i
    RR(node3index+nn) += Res_e(5); // residual for cC for node 3 of element i

    RR(node1index+2*nn) += Res_e(6); // residual for cBC for node 1 of element i
    RR(node2index+2*nn) += Res_e(7); // residual for cBC for node 2 of element i
    RR(node3index+2*nn) += Res_e(8); // residual for cBC for node 3 of element i

    RR(node1index+3*nn) += Res_e(9); // residual for cN for node 1 of element i
    RR(node2index+3*nn) += Res_e(10); // residual for cN for node 2 of element i
    RR(node3index+3*nn) += Res_e(11); // residual for cN for node 3 of element i

    RR(node1index+4*nn) += Res_e(12); // residual for cBN for node 1 of element i
    RR(node2index+4*nn) += Res_e(13); // residual for cBN for node 2 of element i
    RR(node3index+4*nn) += Res_e(14); // residual for cBN for node 3 of element i

    RR(node1index+5*nn) += Res_e(15); // residual for cSiz for node 1 of element i
    RR(node2index+5*nn) += Res_e(16); // residual for cSiz for node 2 of element i
    RR(node3index+5*nn) += Res_e(17); // residual for cSiz for node 3 of element i
    // KK_triplets.push_back(row,col,value)  //?????????????????????? necessary?
    //int varArray[] = {&node1index, &node2index, &node3index};
    std::vector<int> varArray (3.0);
    varArray[0]=node1index;
    varArray[1]=node2index;
    varArray[2]=node3index;

    for(int a=0;a<3;a++)
    {
      for(int b=0;b<3;b++)
      {
            KK_triplets.push_back(T(varArray[a],varArray[b],K_e(a,b)));
            KK_triplets.push_back(T(varArray[a],varArray[b]+nn,K_e(a,b+3)));
            KK_triplets.push_back(T(varArray[a],varArray[b]+2*nn,K_e(a,b+6)));
            KK_triplets.push_back(T(varArray[a],varArray[b]+3*nn,K_e(a,b+9)));
            KK_triplets.push_back(T(varArray[a],varArray[b]+4*nn,K_e(a,b+12)));

            KK_triplets.push_back(T(varArray[a]+nn,varArray[b],K_e(a+3,b)));
            KK_triplets.push_back(T(varArray[a]+nn,varArray[b]+nn,K_e(a+3,b+3)));
            KK_triplets.push_back(T(varArray[a]+nn,varArray[b]+2*nn,K_e(a+3,b+6)));
            KK_triplets.push_back(T(varArray[a]+nn,varArray[b]+3*nn,K_e(a+3,b+9)));
            KK_triplets.push_back(T(varArray[a]+nn,varArray[b]+4*nn,K_e(a+3,b+12)));

            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b],K_e(a+6,b)));
            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b]+nn,K_e(a+6,b+3)));
            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b]+2*nn,K_e(a+6,b+6)));
            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b]+3*nn,K_e(a+6,b+9)));
            KK_triplets.push_back(T(varArray[a]+2*nn,varArray[b]+4*nn,K_e(a+6,b+12)));

            KK_triplets.push_back(T(varArray[a]+3*nn,varArray[b],K_e(a+9,b)));
            KK_triplets.push_back(T(varArray[a]+3*nn,varArray[b]+nn,K_e(a+9,b+3)));
            KK_triplets.push_back(T(varArray[a]+3*nn,varArray[b]+2*nn,K_e(a+9,b+6)));
            KK_triplets.push_back(T(varArray[a]+3*nn,varArray[b]+3*nn,K_e(a+9,b+9)));
            KK_triplets.push_back(T(varArray[a]+3*nn,varArray[b]+4*nn,K_e(a+9,b+12)));

            KK_triplets.push_back(T(varArray[a]+4*nn,varArray[b],K_e(a+12,b)));
            KK_triplets.push_back(T(varArray[a]+4*nn,varArray[b]+nn,K_e(a+12,b+3)));
            KK_triplets.push_back(T(varArray[a]+4*nn,varArray[b]+2*nn,K_e(a+12,b+6)));
            KK_triplets.push_back(T(varArray[a]+4*nn,varArray[b]+3*nn,K_e(a+12,b+9)));
            KK_triplets.push_back(T(varArray[a]+4*nn,varArray[b]+4*nn,K_e(a+12,b+12)));

            KK_triplets.push_back(T(varArray[a]+5*nn,varArray[b],K_e(a+15,b)));
            KK_triplets.push_back(T(varArray[a]+5*nn,varArray[b]+nn,K_e(a+15,b+3)));
            KK_triplets.push_back(T(varArray[a]+5*nn,varArray[b]+2*nn,K_e(a+15,b+6)));
            KK_triplets.push_back(T(varArray[a]+5*nn,varArray[b]+3*nn,K_e(a+15,b+9)));
            KK_triplets.push_back(T(varArray[a]+5*nn,varArray[b]+4*nn,K_e(a+15,b+12)));
      }
    }
  }// closing the element loop
    double error;
    error=RR.norm();

    if(error>tol)
    {
      std::cout<<"error "<<error<<"\n";
      KK.setFromTriplets(KK_triplets.begin(), KK_triplets.end());
      //std::cout<<"K=\n"<<KK<<"\n";
      KK.makeCompressed();
      
      //std::cout<<"KK2\n"<<KK2<<"\n";
      //solver2.analyzePattern(KK2);
      // Compute the numerical factorization
      //solver2.factorize(KK2);
      solver.compute(KK);
      //std::cout<<solver.info()<<\"n";
      //std::cout<<"here\n";
      //Use the factors to solve the linear system
      dU = solver.solve(-1.*RR);
      //std::cout<<"dU\n"<<dU<<"\n";
      //VectorXd dU(nn*3);
      //double dU=KK_triplets.inverse()*RR;
       //dU=-K\RES;
       //c_B0=c_B0+dU.head<nn>;
       //c_C0=c_B0+dU.segment<nn>(nn);
       //c_BC0=c_BC0+dU.tail<nn>;
      for(int i=0;i<nn;i++){
        c_B[i] += dU(i);
        c_C[i] += dU(i+nn);
        c_BC[i] += dU(i+2*nn);
        c_N[i] += dU(i+3*nn);
        c_BN[i] += dU(i+4*nn);
        c_Siz[i] += dU(i+5*nn);
        //std::cout<<"c_B"<<c_B[i]<<"\n";
      }
     }
    else
    {
       std::cout<<"yay!, error "<<error<<"\n";
       break;
    }
} // closes the Newton Raphson
} // closes the function 

