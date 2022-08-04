#include "vec.h"
#include "nondestructivetrimesh.h"
#include <vector>
#include <math.h>   
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
//#include "densitytable.h"

//using namespace Eigen;

// Updated lagrangian formulation for growing membrane 

// function declarations
//void driverVelocities_Explicit(std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt);
//void updateForces_Explicit(std::vector<Vec3d> &X_t, std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, double dt, double dt_real, const NonDestructiveTriMesh &mesh);
void calculateRes(const std::vector<Vec3d> &X_t,const std::vector<Vec3d> &x_t, std::vector<Eigen::VectorXd> &ea_t,const NonDestructiveTriMesh &mesh, double dt_real, Eigen::VectorXd &Res, Eigen::VectorXd &Res_ea, Eigen::VectorXd &MM);
void evalElementRe(Eigen::Vector3d &node1_X, Eigen::Vector3d &node2_X, Eigen::Vector3d &node3_X, Eigen::Vector3d &node1_x, Eigen::Vector3d &node2_x, Eigen::Vector3d &node3_x, Eigen::VectorXd &node1_ea, Eigen::VectorXd &node2_ea, Eigen::VectorXd &node3_ea, double dt_real, Eigen::VectorXd &Re, Eigen::VectorXd &Re_ea, Eigen::Vector3d &MMe);
void calculateConstraintForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t);
void calculatePressureForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t, const Eigen::VectorXd &MM, const NonDestructiveTriMesh &mesh, double t_tot);
void calculateBendingRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, Eigen::VectorXd &Res_bending);
double calculateBendingEnergy(const std::vector<Eigen::Vector3d> &v_vec_X , const std::vector<Eigen::Vector3d> &v_vec_x,bool print_info);

// Driver for Zebrafish based on reading from file
void driverVelocities_Zebrafish(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){
	// Open the corresponding file with the velocities and read into vectors
	int currframe = 4*int(t_tot/6)+130; // each frame is 1.5 min 
	std::cout<<"frame "<<currframe<<"\n";
	std::string filename="velmap_reflexed/velocity_at";
	filename.append(std::to_string(currframe));
	filename.append(".txt");
	std::ifstream velocityfile(filename.c_str());
	std::cout<<"file: "<<filename<<"\n";
    std::string line;
    std::vector<Vec3d> cell_positions; cell_positions.clear(); // positions in microns 
    std::vector<Vec3d> cell_u; cell_u.clear(); // displacement in microns 
    if (velocityfile.is_open())
    {
        while ( getline (velocityfile,line) )
        {
        	std::stringstream ss0(line);
        	std::vector<std::string> result;
			while( ss0.good() )
			{
    			std::string substr;
    			getline( ss0, substr, ',' );
    			result.push_back( substr );
			}
            Vec3d new_cell, new_u; 
            new_cell[0]= std::stod(result[0]);
            new_cell[1]= std::stod(result[1]);
            new_cell[2]= std::stod(result[2]);
            new_u[0]= std::stod(result[3]);
            new_u[1]= std::stod(result[4]);
            new_u[2]= std::stod(result[5]);
            cell_positions.push_back(new_cell);
            cell_u.push_back(new_u);
            //std::cout<<new_cell<<"\n";
            //std::cout<<new_u<<"\n\n";
        }
    }
    // Loop over the nodes of the mesh and find the closest cells and compute their average velocity
    int ncell = cell_positions.size();
    double radius = 15; // microns
    std::vector<int> neighbors;neighbors.clear();
    int nn = X_t.size();
	for(int ni=0;ni<nn;ni++){
		for(int nj=0;nj<ncell;nj++){
			Vec3d d_vec = X_t[ni]-cell_positions[nj];
			double d_mag = sqrt(d_vec[0]*d_vec[0]+d_vec[1]*d_vec[1]+d_vec[2]*d_vec[2]);
			if(d_mag<radius){
				neighbors.push_back(nj);
			}
		}
		// given the neighbor list get their average velocity 
		Vec3d vel_aux; vel_aux[0]=0.;vel_aux[1]=0.;vel_aux[2]=0.;
		for(int nj=0;nj<neighbors.size();nj++){
			vel_aux = vel_aux + cell_u[neighbors[nj]]/6.;
			// note: displacements are in microns every 6min 
			// velocity is in um/min
		}
		if(neighbors.size()>1){vel_aux = vel_aux/neighbors.size();}

		// Update the velocities of the node 
		V_t_hlf[ni] = vel_aux; // um/min
		//std::cout<<vel_aux<<"\n";
	}
}

void driverVelocities_Explicit(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar examples
		//V_t_hlf[ni][2]=0.0;
	}
	
	std::cout<<"done.\n";
}
void driverVelocities_Uniaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	int n_cont_v = 0;
	int n_cont_s = 0;
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		V_t_hlf[ni][2]=0.0;
		// For the left fix and for the right impose a deformation 
		if(X_t[ni][0]<-0.0499){
			//std::cout<<"constrain left.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][0]>0.9499+0.1*t_tot && t_tot<=2.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.1;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_v+=1;
		}
		if(X_t[ni][0]>0.9499+2.00*0.1 && t_tot>=2.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_s+=1;
		}
	}
	
	std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}
///
void driverVelocities_BulgeTest(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";

	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		//V_t_hlf[ni][2]=0.0;
		// For the left fix and for the right impose a deformation 
		if(X_t[ni][0]<0.001 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain left.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][0]>0.999 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]<0.001 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain bottom.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]>0.999 && X_t[ni][2]<0.001 && X_t[ni][2]>-0.001){
			//std::cout<<"constrain top.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
	}
	
	//std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}
/// 
void driverVelocities_OffBiaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	int n_cont_v = 0;
	int n_cont_s = 0;
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		V_t_hlf[ni][2]=0.0;
		// For the left fix and for the right impose a deformation 
		if(X_t[ni][0]<-0.04){
			//std::cout<<"constrain left.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][0]>0.94+0.1*t_tot && t_tot<1.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.1;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_v+=1;
		}
		if(X_t[ni][0]>0.94+0.99*0.1 && t_tot>=1.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_s+=1;
		}
		if(X_t[ni][1]>0.99){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]<0.01){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
	}
	
	//std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}
void driverVelocities_Biaxial(const std::vector<Vec3d> &X_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t,double dt, double t_tot){

	// Central difference following Belytshko Page 333

	// X_t contains the positions at time t
	// V_t contains the velocities at time t
	// A_t contains the velocities at time t
	// ea_t contains the euler almansi strain per node at time t (not needed in this function)
	// dt is the desired time step, note that Topo might take smaller time step

	// Velocities that we are going to give to Topo
	// These are velocities at the half time step
	int nn = X_t.size();
	std::cout<<"loop over nodes to get half velocities...\n";
	int n_cont_v = 0;
	int n_cont_s = 0;
	for(int ni=0;ni<nn;ni++){
		V_t_hlf[ni] = V_t[ni] + (dt/2.0)*A_t[ni] ;
		// turn off the normal component for the planar example
		V_t_hlf[ni][2]=0.0;
		// For the left fix and for the right impose a deformation 
		if(X_t[ni][0]<-0.04){
			//std::cout<<"constrain left.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][0]>0.94+0.1*t_tot && t_tot<5.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.1;
			V_t_hlf[ni][2]=0.0;
			n_cont_v+=1;
		}
		if(X_t[ni][0]>0.94+0.99*0.1 && t_tot>=5.0){
			//std::cout<<"constrain right.\n";
			V_t_hlf[ni][0]=0.0;
			V_t_hlf[ni][2]=0.0;
			n_cont_s+=1;
		}
		if(X_t[ni][1]>0.99+0.1*t_tot && t_tot<5.0){
			//std::cout<<"constrain top.\n";
			V_t_hlf[ni][1]=0.1;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]>0.99+0.99*0.1 && t_tot>=5.0){
			//std::cout<<"constrain top.\n";
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
		if(X_t[ni][1]<0.01){
			//std::cout<<"constrain bottom.\n";
			V_t_hlf[ni][1]=0.0;
			V_t_hlf[ni][2]=0.0;
		}
	}
	
	std::cout<<"done, imposed "<<n_cont_v<<", "<<n_cont_s<<" right constraints.\n";
}

void updateForces_Explicit(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Vec3d> &V_t, std::vector<Vec3d> &V_t_hlf, std::vector<Vec3d> &A_t, std::vector<Eigen::VectorXd> &ea_t, double dt, double dt_real, const NonDestructiveTriMesh &mesh, double t_tot){

	// Central difference following Belytshko Page 333 

	// at this point we passed the velocities of the mid time point to Topo
	// and Topo integrated the equation and gives us the actual current positions
	// and the actual time step, same mesh as before, the mesh changes are done
	// separately

	//** Update the force f_t  **//
	// calculate the internal residual
	int n_node = X_t.size();
	Eigen::VectorXd Res(n_node*3);
	// I also need the residual for the update of the strain
	Eigen::VectorXd Res_ea(n_node*6);
	// And here I can get the lumped mass 'matrix' easily as well
	Eigen::VectorXd MM(n_node);

	std::cout<<"going into calculateRes() \n";
	std::cout<<"number of nodes= "<<n_node<<"\n";
	//******************************************************//
	calculateRes(X_t,x_t,ea_t,mesh,dt_real, Res,Res_ea,MM); 
	//******************************************************//

	// Calculate the force constraint, forces needed to enforce that the
	// growing tissue wants to stick to the sphere
	Eigen::VectorXd F_constraint(n_node*3);F_constraint.setZero();
	//calculateConstraintForce(F_constraint,x_t);

	// Pressure force for the bulge test
	//calculatePressureForce(F_constraint,x_t,MM, mesh, t_tot);

	// Calculate bending residual, to prevent buckling, treat as regularization
	Eigen::VectorXd Res_bending(n_node*3);Res_bending.setZero();
	//calculateBendingRes(X_t,x_t,mesh,Res_bending);
	//std::cout<<"bending residual\n"<<Res_bending<<"\n";

	// net force
	Eigen::VectorXd f_t = F_constraint - Res - Res_bending;
	std::cout<<"net force norm "<<Res.norm()<<"\n";
	// Update the Acceleration for the next time step, and also the strain
	double damp = 2.0; // [(kg/s)] not sure this is good or bad? 1 was ok for uniaxial, 0.1 was ok for ZF
	//double rho = 1.0; // [ug/um^3] also not sure if this is good or not
	double rho = 1000.0; // [kg/m^3] 
	// Update acceleration and strain
	for(int ni=0;ni<n_node;ni++){
		// updating the acceleration
		A_t[ni][0] = (f_t(ni*3+0)-damp*MM(ni)*rho*V_t[ni][0])/(MM(ni)*rho);
		A_t[ni][1] = (f_t(ni*3+1)-damp*MM(ni)*rho*V_t[ni][1])/(MM(ni)*rho);
		A_t[ni][2] = (f_t(ni*3+2)-damp*MM(ni)*rho*V_t[ni][2])/(MM(ni)*rho);	
		// updating the strains at the nodes
		ea_t[ni](0) = Res_ea(ni*6+0)/MM(ni);
		ea_t[ni](1) = Res_ea(ni*6+1)/MM(ni);
		ea_t[ni](2) = Res_ea(ni*6+2)/MM(ni);
		ea_t[ni](3) = Res_ea(ni*6+3)/MM(ni);
		ea_t[ni](4) = Res_ea(ni*6+4)/MM(ni);
		ea_t[ni](5) = Res_ea(ni*6+5)/MM(ni);
		//std::cout<<"node "<<ni<<", f_t = "<<f_t(ni*3+0)<<", "<<f_t(ni*3+1)<<", "<<f_t(ni*3+2)<<"\n";
		//std::cout<<"node "<<ni<<", Res = "<<Res(ni*3+0)<<", "<<Res(ni*3+1)<<", "<<Res(ni*3+2)<<"\n";
	}
	
	
	// Update the velocities for next time step
	for(int ni=0;ni<n_node;ni++){
		if(dt_real>dt/2.0){
			V_t[ni] = V_t_hlf[ni] + (dt_real - dt/2.0)*A_t[ni];
		}else{
			V_t[ni] = V_t[ni] + dt_real*A_t[ni];
		}
	}
}

void calculateRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, std::vector<Eigen::VectorXd> &ea_t,const NonDestructiveTriMesh &mesh, double dt_real, Eigen::VectorXd &Res, Eigen::VectorXd &Res_ea, Eigen::VectorXd &MM){
	Res.setZero();
	Res_ea.setZero();
	MM.setZero();
	// START LOOP OVER ELEMENTS
	Eigen::VectorXd Re(9);
	Eigen::VectorXd Re_ea(18);
	Eigen::Vector3d MMe;
    int n_elem = mesh.num_triangles();
    for(int ei=0;ei<n_elem;ei++){
       	// connectivity of the linear elements
        int node1index = mesh.get_triangle(ei)[0];
    	int node2index = mesh.get_triangle(ei)[1];
    	int node3index = mesh.get_triangle(ei)[2];
        
        // previous time step nodal positions for this element
		Eigen::Vector3d node1_X; node1_X<<X_t[node1index][0],X_t[node1index][1],X_t[node1index][2];
		Eigen::Vector3d node2_X; node2_X<<X_t[node2index][0],X_t[node2index][1],X_t[node2index][2];
		Eigen::Vector3d node3_X; node3_X<<X_t[node3index][0],X_t[node3index][1],X_t[node3index][2];

		// previous time step euler almansi strains at the nodes
		Eigen::VectorXd node1_ea(6); 
		node1_ea<<ea_t[node1index](0),ea_t[node1index](1),ea_t[node1index](2),ea_t[node1index](3),ea_t[node1index](4),ea_t[node1index](5);
		Eigen::VectorXd node2_ea(6); 
		node2_ea<<ea_t[node2index](0),ea_t[node2index](1),ea_t[node2index](2),ea_t[node2index](3),ea_t[node2index](4),ea_t[node2index](5);
		Eigen::VectorXd node3_ea(6); 
		node3_ea<<ea_t[node3index](0),ea_t[node3index](1),ea_t[node3index](2),ea_t[node3index](3),ea_t[node3index](4),ea_t[node3index](5);

		// current nodal positions for this element
		Eigen::Vector3d node1_x; node1_x<<x_t[node1index][0],x_t[node1index][1],x_t[node1index][2];
		Eigen::Vector3d node2_x; node2_x<<x_t[node2index][0],x_t[node2index][1],x_t[node2index][2];
		Eigen::Vector3d node3_x; node3_x<<x_t[node3index][0],x_t[node3index][1],x_t[node3index][2];
						
        // and calculate the element Re
        Re.setZero();
        Re_ea.setZero();
        MMe.setZero();
				
        // subroutine to evaluate the element
        //std::cout<<"evaluating element "<<ei<<"\n";
        evalElementRe(node1_X,node2_X,node3_X,node1_x,node2_x,node3_x,node1_ea,node2_ea,node3_ea, dt_real, Re, Re_ea, MMe);
        //std::cout<<"element residual norm "<<Re.norm()<<"\n";
        // LOOP OVER NODES
		for(int nodei=0;nodei<3;nodei++){
			// ASSEMBLE DISPLACEMENT RESIDUAL
			for(int coordi=0;coordi<3;coordi++){
				// residual
				Res(mesh.get_triangle(ei)[nodei]*3+coordi) += Re(nodei*3+coordi);
			}
			for(int coordi=0;coordi<6;coordi++){
				Res_ea(mesh.get_triangle(ei)[nodei]*6+coordi) += Re_ea(nodei*6+coordi);
			}
			MM(mesh.get_triangle(ei)[nodei]) += MMe(nodei);
		}
	}
}

void evalElementRe(Eigen::Vector3d &node1_X, Eigen::Vector3d &node2_X, Eigen::Vector3d &node3_X, Eigen::Vector3d &node1_x, Eigen::Vector3d &node2_x, Eigen::Vector3d &node3_x, Eigen::VectorXd &node1_ea, Eigen::VectorXd &node2_ea, Eigen::VectorXd &node3_ea, double dt_real, Eigen::VectorXd &Re, Eigen::VectorXd &Re_ea, Eigen::Vector3d &MMe){

	// constants and parameters
	//double mu = 1.0e1; // [mili-Pa] material parameter
	double thetag_dot = 0.0; //0.0009; //0.005; // growth rate area_change/time_units, 0.01 is growth of 1% per second

	// example 1, biaxial test with no growth
	double mu = 75000.0; // Pa, from a value we are using for skin 
	//double thetag_dot = 0.0; 
	double thick = 0.005; // [meters] skin is ~0.5cm thick 

	// reference metric plus normal (reference meaning X_t, or configuration at time t)
	Eigen::Vector3d G1 = -1.0*node1_X + 1.0*node2_X;  // 3x1
	Eigen::Vector3d G2 = -1.0*node1_X + 1.0*node3_X; // 3x1
	Eigen::Vector3d G1xG2 = G1.cross(G2);
	Eigen::Vector3d N = G1xG2/G1xG2.norm(); // 3x1

	// interpolation of the euler almansi at the integration point 
	// this is ok for strains. 
	Eigen::VectorXd ea_t_vec = 1.0/3.0*node1_ea + 1.0/3.0*node2_ea + 1.0/3.0*node3_ea;
	// the euler almansi strain is stored as vector. Need to reconstruct in Matrix
	Eigen::Matrix3d ea_t; 
	ea_t<<ea_t_vec(0),ea_t_vec(3),ea_t_vec(4),
		  ea_t_vec(3),ea_t_vec(1),ea_t_vec(5),
		  ea_t_vec(4),ea_t_vec(5),ea_t_vec(2);

	// will need to go from ea_t to be_t to get the area change and evaluate stress
	Eigen::Matrix3d Id; Id<<1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0;
	Eigen::Matrix3d aux = Id-2.0*ea_t;
	Eigen::Matrix3d be_t = aux.inverse();
	// note that the euler almansi strains that are stored correspond to the surface only
	// Moreover, they are stored with respect to the surface at the previous step 't'
	// which now has normal N 
	// the identity matrix used is the full 3D identity so be_t will have the extra 
	// normal component with no stretch. To get the surface only component we have. 
	Eigen::Matrix3d NdyadN = N*N.transpose();
	Eigen::Matrix3d be_t_s = be_t - NdyadN ;
	// do i need to store the thickness strain? I don't think so because the new
	// thickness strain is calculated based only on the current elastic deformation
	// we should just remember that the real be_t has some thickness strain to guarantee
	// incompressibility, just like the bew be has to be adjusted for incompressibility 
	
	// deformed metric and normal (meaning the predicted state at time t+dt)
	Eigen::Vector3d g1 = -1.0*node1_x + 1.0*node2_x; // 3x1
	Eigen::Vector3d g2 = -1.0*node1_x + 1.0*node3_x; // 3x1
	Eigen::Vector3d g1crossg2 = g1.cross(g2);
	Eigen::Vector3d n = g1crossg2/g1crossg2.norm(); // 3x1
	Eigen::Matrix3d ndyadn = n*n.transpose();

	// put Gi in matrix G and get inverse
	Eigen::Matrix3d G; 
	G << G1(0),G2(0),N(0),
		 G1(1),G2(1),N(1),
		 G1(2),G2(2),N(2);
	Eigen::Matrix3d g; 
	g << g1(0),g2(0),n(0),
		 g1(1),g2(1),n(1),
		 g1(2),g2(2),n(2);
	//std::cout<<"G = "<<G<<"\ng="<<g<<"\n";
	Eigen::Matrix3d Ginv = G.inverse();
	Eigen::Vector3d G1inv; G1inv<<Ginv(0,0),Ginv(0,1),Ginv(0,2);
	Eigen::Vector3d G2inv; G2inv<<Ginv(1,0),Ginv(1,1),Ginv(1,2);
	// full deformation gradient at IP before static condensation 
	// this deformation gradient is the one corresponding to the time step, F_delta
	Eigen::Matrix3d F = g1*G1inv.transpose() + g2*G2inv.transpose()+n*N.transpose();
	//std::cout<<"F="<<F<<"\n";
	// not sure I need the following, shouldn't I be able to recast in terms of b?  
	Eigen::Matrix3d C = F.transpose()*F;
	Eigen::Matrix3d Cinv = C.inverse();
	Eigen::Matrix3d Cinv_s = Cinv - NdyadN;

	// Area changes to get the current normal strain 
	double theta_e_t = sqrt(be_t.determinant()); // use the total be_t with 1*NdyadN
	double theta_g = 1.0+thetag_dot*dt_real; // the growth rate is a constant  

	double theta = F.determinant(); // total area change from the deformation of the current step
	// the area changes satisfy theta*theta_e_t = theta_e*theta_g
	double theta_e = theta*theta_e_t/theta_g;
	// what about the volume change? Je_t = Je = 1, so Jdelta = Jg = theta_g

	// ABT 08/02/2018. Growth in response to elastic stretch
	//double K_thetae = 1.5;
	//double p_thetag = 0.0; // [1/sec]
	//if(theta_e>1.1){
	//	theta_g = 1.0+p_thetag*(theta_e/(K_thetae+theta_e));
	//}

	double lamda_N = 1.0/theta_e; // normal strain from the current Fe needed for incompressible
	double p = -mu*(lamda_N*lamda_N); // pressure for sigma_n = 0
	
	// STRESS refered to the 'current' configuration t, in my notes this is 'Sstar'
	Eigen::Matrix3d S = ((mu/theta_g)*be_t_s + theta_g*p*Cinv_s);
	//std::cout<<"S_norm "<<S.norm()<<"\n";
	// UPDATE the strain 
	// well, the strain doesn't get updated right now, here we get the IP value of the updated
	// strain to interpolate it back to the nodes 
	Eigen::Matrix3d be = (1.0/(theta_g*theta_g))*F*be_t*F.transpose();
	Eigen::Matrix3d ea_s = 0.5*(Id-be.inverse());
	Eigen::VectorXd ea_s_vec(6);
	ea_s_vec<<ea_s(0,0),ea_s(1,1),ea_s(2,2),ea_s(0,1),ea_s(0,2),ea_s(1,2);

	// Given the stress fill out the element residual, lets leave the tangent alone
	Re.setZero();
	Re_ea.setZero();
	MMe.setZero();
	for(int i=0;i<3;i++){
	    for(int p=0;p<3;p++){
	        Eigen::Matrix3d temp; temp.setZero();
	        temp(p,i) = 1.0;
	        Eigen::Vector3d varg1 = -1.0*temp.col(0) + temp.col(1);
	        Eigen::Vector3d varg2 = -1.0*temp.col(0) + temp.col(2);
	        // varn =  cross(varg1,varg2)/norm(cross(varg1,varg2)) ??
	        Eigen::Matrix3d varF = varg1*G1inv.transpose() + varg2*G2inv.transpose(); //+outer(varn,Ginv[2,:])?
	        Eigen::Matrix3d varC = varF.transpose()*F + F.transpose()*varF;
	        Eigen::Matrix3d varE = 0.5*varC;
	        double trvarE = varE(0,0)+varE(1,1)+varE(2,2);
	        Eigen::VectorXd varE_voigt(6);
	        varE_voigt<<varE(0,0),varE(1,1),varE(2,2),2*varE(0,1),2*varE(0,2),2*varE(1,2);
	        Eigen::VectorXd S_voigt(6);
	        S_voigt<< S(0,0),S(1,1),S(2,2),S(0,1),S(0,2),S(1,2);
	        Re(i*3+p) = S_voigt.dot(varE_voigt)*G1xG2.norm()/2.0*thick;
	    }
	    for(int q=0;q<6;q++){
	    	Re_ea(i*6+q) = (1./3.)*ea_s_vec(q)*G1xG2.norm()/2.0*thick;
	    }
	    MMe(i) = (1./3.)*G1xG2.norm()/2.0*thick;
	}
}

void calculateConstraintForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t){
	// Force per node to guarantee it is close to the sphere 
	double dist;
	double radius = 350.0; // [m]
	double penalty_dist = 4000.0; // 1000.0; 1000 works
	F_constraint.setZero();
	int n_node = x_t.size();
	for(int ni=0;ni<n_node;ni++){
		dist = sqrt(x_t[ni][0]*x_t[ni][0] + x_t[ni][1]*x_t[ni][1] + x_t[ni][2]*x_t[ni][2]);
		// force should be proportional to the distance to the sphere of radius r 
		// and directed along the vector connecting the center of the sphere to the current x
		Vec3d u_x = x_t[ni]/dist;
		F_constraint(ni*3+0) = penalty_dist*(radius-dist)*u_x[0];
		F_constraint(ni*3+1) = penalty_dist*(radius-dist)*u_x[1];
		F_constraint(ni*3+2) = penalty_dist*(radius-dist)*u_x[2];
	}
}

void calculatePressureForce(Eigen::VectorXd &F_constraint, const std::vector<Vec3d> &x_t, const Eigen::VectorXd &MM,const NonDestructiveTriMesh &mesh, double t_tot ){
	F_constraint.setZero();
	double pres = 2000.0; // [Pa] 
	if(t_tot<1){
		pres=2000.0*t_tot;
		//pres = 0.0;
	}
	if(t_tot>2){
		pres=2000.0*(3-t_tot);
		//pres = 0.0;
	}
	if(t_tot>3){
		pres=0.0;
	}
	
	// loop over nodes and get the normal vector
	double thick = 0.1; //[m] this should be passed actually, I just know it is 0.1 for everywhere 
	for(int ni=0;ni<x_t.size();ni++){
		Eigen::Vector3d normal; normal.setZero();
		// Get the triangles inciden in this vertex
		double Asum = 0.0;
		std::vector<size_t> onering = mesh.m_vertex_to_triangle_map[ni];
		for(int ti=0;ti<onering.size();ti++){
			// get the normal and area of each tringle
        	const Vec3st& tri = mesh.get_triangle(onering[ti]);
        	Eigen::Vector3d tri_n0; tri_n0<<x_t[tri[0]][0],x_t[tri[0]][1],x_t[tri[0]][2] ;
        	Eigen::Vector3d tri_n1; tri_n1<<x_t[tri[1]][0],x_t[tri[1]][1],x_t[tri[1]][2] ;
        	Eigen::Vector3d tri_n2; tri_n2<<x_t[tri[2]][0],x_t[tri[2]][1],x_t[tri[2]][2] ;
        	Eigen::Vector3d E1 = tri_n1 - tri_n0; 
        	Eigen::Vector3d E2 = tri_n2 - tri_n0; 
			Eigen::Vector3d normal_ti = E1.cross(E2);
			double area_ti = 0.5*normal_ti.norm();
			normal_ti = normal_ti/normal_ti.norm();
			
			// add to the average calculation
			normal += area_ti*normal_ti;
			Asum += area_ti; 
		}
		normal = normal/Asum;
		normal = normal/normal.norm();
		
		//std::cout<<"vertex "<<ni<<", normal = "<<normal(0)<<", "<<normal(1)<<", "<<normal(2)<<", area = "<<Asum<<"\n";
		// pressure is constant traction, needs to be scaled by area of each node
		F_constraint(ni*3+0) = -pres*(MM(ni)/thick)*normal(0);
		F_constraint(ni*3+1) = -pres*(MM(ni)/thick)*normal(1);
		F_constraint(ni*3+2) = -pres*(MM(ni)/thick)*normal(2);
		//F_constraint(ni*3+0) = 0.0;
		//F_constraint(ni*3+1) = 0.0;
		//F_constraint(ni*3+2) = -1.0;
		//std::cout<<"F_constraint = "<<F_constraint(ni*3+0)<<", "<<F_constraint(ni*3+1)<<", "<<F_constraint(ni*3+2)<<"\n";  
	}
}

void calculateBendingRes(const std::vector<Vec3d> &X_t, const std::vector<Vec3d> &x_t, const NonDestructiveTriMesh &mesh, Eigen::VectorXd &Res_bending){
	// Given the coordinates of the nodes at the end of time step t, and the new positions at t+dt, get the
	// forces for the next time step, and that includes some bending opposing force. 
	// this is mostly to regularize the deformation and avoid the bad crumpling effect

	Res_bending.setZero();
	// loop over edges
	int n_edges = mesh.m_edges.size();
	//std::cout<<"loop over "<<n_edges<<"edges\n";
	for(int e=0;e<n_edges;e++){
		// if it is an interior triangle then get the triangles incident on it
		if(!mesh.m_is_boundary_edge[e]){
			// get triangles incident to the edge
        	int t0 = mesh.m_edge_to_triangle_map[e][0];
        	int t1 = mesh.m_edge_to_triangle_map[e][1];
        	const Vec3st& tri0 = mesh.get_triangle(t0);
        	const Vec3st& tri1 = mesh.get_triangle(t1);
			// get vertex opposite the edge for each triangle
        	int v2_t0 = mesh.get_third_vertex( e, tri0 );
        	int v2_t1 = mesh.get_third_vertex( e, tri1 );
			// get also the two vertices that form this edge
			int v0 = mesh.m_edges[e][0];      
			int v1 = mesh.m_edges[e][1];

			// get the reference coordinates of the vertices
			Eigen::Vector3d v0_X, v1_X, v2t0_X, v2t1_X;
			v0_X<<X_t[v0][0],X_t[v0][1],X_t[v0][2];
			v1_X<<X_t[v1][0],X_t[v1][1],X_t[v1][2];
			v2t0_X<<X_t[v2_t0][0],X_t[v2_t0][1],X_t[v2_t0][2];
			v2t1_X<<X_t[v2_t1][0],X_t[v2_t1][1],X_t[v2_t1][2];
			// get the deformed coordinates of the vertices
			Eigen::Vector3d v0_x, v1_x, v2t0_x, v2t1_x;
			v0_x<<x_t[v0][0],x_t[v0][1],x_t[v0][2];
			v1_x<<x_t[v1][0],x_t[v1][1],x_t[v1][2];
			v2t0_x<<x_t[v2_t0][0],x_t[v2_t0][1],x_t[v2_t0][2];
			v2t1_x<<x_t[v2_t1][0],x_t[v2_t1][1],x_t[v2_t1][2];

			// put the indices and coordinates of the vertices in vectors
			std::vector<int> v_vec_ind(4,0);
			v_vec_ind[0]= v0; v_vec_ind[1]= v1;
			v_vec_ind[2]= v2_t0; v_vec_ind[3]= v2_t1;
			std::vector<Eigen::Vector3d> v_vec_X(4,Eigen::Vector3d());
			v_vec_X[0] = v0_X; v_vec_X[1] = v1_X; 
			v_vec_X[2] = v2t0_X; v_vec_X[3] = v2t1_X; 
			std::vector<Eigen::Vector3d> v_vec_x(4,Eigen::Vector3d());
			v_vec_x[0] = v0_x; v_vec_x[1] = v1_x; 
			v_vec_x[2] = v2t0_x; v_vec_x[3] = v2t1_x; 

			// ENERGY from bending
			double Energy_b = calculateBendingEnergy(v_vec_X,v_vec_x,false); 
			//std::cout<<"Energy_b = "<<Energy_b<<"\n";
			// variation of the bending in terms of the displacement to update force
			// loop over the four nodes that affect this edge
			double delta = 1e-6;
			double Eb_plus;
			for(int ni=0;ni<4;ni++){
				for(int ci=0;ci<3;ci++){
					// variation of vj, first numerically for simplicity
					v_vec_x[ni](ci) += delta; 
					Eb_plus = calculateBendingEnergy(v_vec_X,v_vec_x,false); 
					Res_bending(v_vec_ind[ni]*3+ci)+=(Eb_plus-Energy_b)/delta;
					v_vec_x[ni](ci) -= delta;
				}
			}
		}
	}

}

double calculateBendingEnergy(const std::vector<Eigen::Vector3d> &v_vec_X , const std::vector<Eigen::Vector3d> &v_vec_x, bool print_info){
	
	// beding energy is controlled by a material parameter
	//double mu_b = 10000.0;
	double mu_b = 10.0; //[Pa]
	double h = 0.1; // [m], thickness

	// unpack
	Eigen::Vector3d v0_X, v1_X, v2t0_X, v2t1_X;
	Eigen::Vector3d v0_x, v1_x, v2t0_x, v2t1_x;
	v0_X = v_vec_X[0];
	v1_X = v_vec_X[1];
	v2t0_X = v_vec_X[2];
	v2t1_X = v_vec_X[3];
	v0_x = v_vec_x[0];
	v1_x = v_vec_x[1];
	v2t0_x = v_vec_x[2];
	v2t1_x = v_vec_x[3];

	// get the angle between the two triangles before, based on X
	// normal of triangles
	Eigen::Vector3d edge_X = v1_X-v0_X;
	Eigen::Vector3d edget0_X = v2t0_X-v0_X;
	Eigen::Vector3d Nt0_X = edge_X.cross(edget0_X);
	// area
	double At0 = 0.5*Nt0_X.norm();
	// normal
	Nt0_X = Nt0_X/Nt0_X.norm();
	//
	Eigen::Vector3d edget1_X = v2t1_X-v0_X;
	Eigen::Vector3d Nt1_X = edget1_X.cross(edge_X);
	// area
	double At1 = 0.5*Nt1_X.norm();
	Nt1_X = Nt1_X/Nt1_X.norm();
	// get the angle directly from normals
	double N0dotN1 = Nt0_X.dot(Nt1_X)>=1 ? 0.9999999999 : Nt0_X.dot(Nt1_X);
	double thetae_X = acos(N0dotN1);
	if(print_info){
	std::cout<<"initial areas, At0="<<At0<<", At1="<<At1<<"\n";
	std::cout<<"reference normals, N0=("<<Nt0_X(0)<<","<<Nt0_X(1)<<","<<Nt0_X(2)<<"), N1=("<<Nt1_X(0)<<","<<Nt0_X(1)<<","<<Nt0_X(2)<<")\n";
	std::cout<<"dot product = "<<Nt0_X.dot(Nt1_X)<<"\n";
	std::cout<<"reference angle, theta_X = "<<thetae_X<<"\n";}
	// do the same for the deformed triangle
	Eigen::Vector3d edge_x = v1_x-v0_x;
	Eigen::Vector3d edget0_x = v2t0_x-v0_x;
	Eigen::Vector3d Nt0_x = edge_x.cross(edget0_x);
	Nt0_x = Nt0_x/Nt0_x.norm();
	//
	Eigen::Vector3d edget1_x = v2t1_x-v0_x;
	Eigen::Vector3d Nt1_x = edget1_x.cross(edge_x);
	Nt1_x = Nt1_x/Nt1_x.norm();
	// deformed angle
	double n0dotn1 = Nt0_x.dot(Nt1_x)>=1 ? 0.9999999999 : Nt0_x.dot(Nt1_x);
	double thetae_x = acos(n0dotn1);
	if(print_info){
	std::cout<<"deformed normals, n0=("<<Nt0_x(0)<<","<<Nt0_x(1)<<","<<Nt0_x(2)<<"), n1=("<<Nt1_x(0)<<","<<Nt0_x(1)<<","<<Nt0_x(2)<<")\n";
	std::cout<<"dot product = "<<Nt0_x.dot(Nt1_x)<<"\n";
	std::cout<<"deformed angle, theta_x = "<<thetae_x<<"\n";}
	// length of the edge
	double lbar_edge = edge_X.norm();
	// barycentric length
	double lstar_edge = (At0+At1)/(3*lbar_edge);

	// ENERGY from bending [N*m]
	double Energy_b = mu_b*h*h*h*(thetae_x*thetae_x)*lbar_edge/lstar_edge; 


	return Energy_b;
}
