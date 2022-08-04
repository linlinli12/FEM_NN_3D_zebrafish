// ---------------------------------------------------------
//
//  ScriptInit.cpp
//  Tyson Brochu 2011
//
//  Parse a script text file to initialize the simulation and mesh obects.
//
// ---------------------------------------------------------


#include <scriptinit.h>

#include <invagination.h>
#include <fstream>
#include <geometryinit.h>
#include <iomesh.h>
#include <subdivisionscheme.h>

// add by linlin Jun24,2020
#include <string>
#include <iostream>
#include <sstream>
#include <dirent.h>
#include <stdlib.h>                                                             
#include <stdio.h> 

#include <random>

#ifdef _MSC_VER
//this is kind of hacky, but seems to do the trick for now (on Windows, the name is prefaced by an underscore)
//so we'll just rename is here
#define snprintf _snprintf 
#endif
// ---------------------------------------------------------

void ScriptInit::parse_surftrack_parameters( const ParseTree& surftrack_branch )
{
    
    int use_fraction, perform_improvement, topology_changes, collision_safety, non_manifold;
    std::string subdivision_scheme;
    
    surftrack_branch.get_int( "use_fraction", use_fraction );
    surftrack_branch.get_number( "min_edge_length", surf_track_params.m_min_edge_length );
    surftrack_branch.get_number( "max_edge_length", surf_track_params.m_max_edge_length  );
    surftrack_branch.get_number( "max_volume_change", surf_track_params.m_max_volume_change );
    surftrack_branch.get_number( "min_triangle_angle", surf_track_params.m_min_triangle_angle );   
    surftrack_branch.get_number( "max_triangle_angle", surf_track_params.m_max_triangle_angle );      
    surftrack_branch.get_number( "min_triangle_area", surf_track_params.m_min_triangle_area );
    
    int use_curvature_when_splitting;
    if ( surftrack_branch.get_int( "use_curvature_when_splitting", use_curvature_when_splitting ) )
    {
        surf_track_params.m_use_curvature_when_splitting = ( use_curvature_when_splitting != 0 );
    }
    
    int use_curvature_when_collapsing;
    if ( surftrack_branch.get_int( "use_curvature_when_collapsing", use_curvature_when_collapsing ) )
    {
        surf_track_params.m_use_curvature_when_collapsing = ( use_curvature_when_collapsing != 0 );
    }
    
    surftrack_branch.get_number( "min_curvature_multiplier", surf_track_params.m_min_curvature_multiplier );
    surftrack_branch.get_number( "max_curvature_multiplier", surf_track_params.m_max_curvature_multiplier );
    surftrack_branch.get_number( "merge_proximity", surf_track_params.m_merge_proximity_epsilon );
    surftrack_branch.get_number( "repulsion_proximity", surf_track_params.m_proximity_epsilon );
    surftrack_branch.get_number( "friction_coefficient", surf_track_params.m_friction_coefficient );
    surftrack_branch.get_int( "perform_improvement", perform_improvement );
    surftrack_branch.get_int( "allow_topology_changes", topology_changes );   
    surftrack_branch.get_int( "collision_safety", collision_safety );   
    surftrack_branch.get_string( "subdivision_scheme", subdivision_scheme );
    
    surf_track_params.m_use_fraction = ( use_fraction != 0 );   
    surf_track_params.m_perform_improvement = (perform_improvement != 0);
    surf_track_params.m_allow_topology_changes = (topology_changes != 0);
    surf_track_params.m_collision_safety = (collision_safety != 0);   
    
    if ( surftrack_branch.get_int( "allow_non_manifold", non_manifold ) )
    {
      surf_track_params.m_allow_non_manifold = (non_manifold != 0);
    }
    
    if ( strcmp( subdivision_scheme.c_str(), "butterfly" ) == 0 )
    {
        surf_track_params.m_subdivision_scheme = new ButterflyScheme();
    }
    else
    {
        surf_track_params.m_subdivision_scheme = new MidpointScheme();
    }
    
    int allow_vertex_movement;
    if ( surftrack_branch.get_int( "allow_vertex_movement", allow_vertex_movement ) )
    {
        surf_track_params.m_allow_vertex_movement = ( allow_vertex_movement != 0 );
    }
    
}

// ---------------------------------------------------------


void ScriptInit::parse_invagination( const ParseTree& invagination_sim_branch )
{
    double speed;
    invagination_sim_branch.get_number( "speed", speed );
    
    
    driver = new InvaginationDriver( speed );
}

// ---------------------------------------------------------


void ScriptInit::parse_camera( const ParseTree& camera_branch )
{
    camera_branch.get_vec3d( "target", camera_target );
    camera_branch.get_number( "distance", camera_distance );
    camera_branch.get_number( "heading", camera_heading );   
    camera_branch.get_number( "pitch", camera_pitch );
}

// ---------------------------------------------------------

void ScriptInit::parse_sheet( const ParseTree& sheet_branch )
{
    Vec3d lower_corner;
    sheet_branch.get_vec3d( "corner", lower_corner );
    
    double sheet_dx;
    sheet_branch.get_number( "dx", sheet_dx );
    
    int sheet_ni, sheet_nj;
    sheet_branch.get_int( "ni", sheet_ni );
    sheet_branch.get_int( "nj", sheet_nj );
    
    std::vector<Vec3d> sheet_verts;
    std::vector<Vec3st> sheet_tris;
    
    create_sheet( sheet_verts, sheet_tris, lower_corner, sheet_dx, sheet_ni, sheet_nj );      
    
    
    Vec3d rotate_axis;
    double rotate_radians;
    if ( sheet_branch.get_vec3d( "rotate_axis", rotate_axis ) )
    {
        if ( sheet_branch.get_number( "rotate_radians", rotate_radians ) )
        {
            Vec3d centre(0,0,0);
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                centre += sheet_verts[i];
            }
            centre /= static_cast<double>( sheet_verts.size() );
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] -= centre;
            }
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] = rotate( sheet_verts[i], rotate_radians, rotate_axis );
                sheet_verts[i] += centre;
            }
            
        }
    }
    
    
    int is_solid = 0;
    sheet_branch.get_int( "is_solid", is_solid );  
    
    std::vector<double> sheet_masses;
    
    if ( is_solid )
    {
        sheet_masses.resize( sheet_verts.size(), std::numeric_limits<double>::infinity() );      
    }
    else
    {
        sheet_masses.resize( sheet_verts.size(), 1.0 );
    }
    
    append_mesh( triangles, vertices, masses, sheet_tris, sheet_verts, sheet_masses );
    
}


// ---------------------------------------------------------

void ScriptInit::parse_curved_sheet( const ParseTree& sheet_branch )
{
    Vec3d lower_corner;
    sheet_branch.get_vec3d( "corner", lower_corner );
    
    double sheet_dx;
    sheet_branch.get_number( "dx", sheet_dx );
    
    int sheet_ni, sheet_nj;
    sheet_branch.get_int( "ni", sheet_ni );
    sheet_branch.get_int( "nj", sheet_nj );
    
    std::vector<Vec3d> sheet_verts;
    std::vector<Vec3st> sheet_tris;
    
    create_curved_sheet( sheet_verts, sheet_tris, lower_corner, sheet_dx, sheet_ni, sheet_nj );      
    
    Vec3d rotate_axis;
    double rotate_radians;
    if ( sheet_branch.get_vec3d( "rotate_axis", rotate_axis ) )
    {
        if ( sheet_branch.get_number( "rotate_radians", rotate_radians ) )
        {
            Vec3d centre(0,0,0);
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                centre += sheet_verts[i];
            }
            centre /= static_cast<double>( sheet_verts.size() );
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] -= centre;
            }
            
            for ( unsigned int i = 0; i < sheet_verts.size(); ++i )
            {
                sheet_verts[i] = rotate( sheet_verts[i], rotate_radians, rotate_axis );
                sheet_verts[i] += centre;
            }
            
        }
    }
    
    
    std::vector<double> sheet_masses( sheet_verts.size(), 1.0 );
    
    append_mesh( triangles, vertices, masses, sheet_tris, sheet_verts, sheet_masses );
    
}

// ---------------------------------------------------------

void ScriptInit::parse_sphere( const ParseTree& sphere_branch )
{
    Vec3d sphere_center;
    sphere_branch.get_vec3d( "sphere_center", sphere_center );
    double sphere_radius;
    sphere_branch.get_number( "sphere_radius", sphere_radius );
    
    double dx;
    sphere_branch.get_number( "sphere_dx", dx );
    
    int is_solid = 0;
    sphere_branch.get_int( "is_solid", is_solid );
    
    std::vector<Vec3d> sphere_vertices;
    std::vector<Vec3st> sphere_triangles;
    
    create_sphere( sphere_center, sphere_radius, dx, sphere_vertices, sphere_triangles );
    
    std::vector<Vec3d> sphere_velocities( sphere_vertices.size(), Vec3d(0) );
    
    std::vector<double> sphere_masses;
    if ( is_solid == 0 )
    {
        sphere_masses.resize( sphere_vertices.size(), 1.0 );
    }
    else
    {
        sphere_masses.resize( sphere_vertices.size(), std::numeric_limits<double>::infinity() );         
    }
    
    append_mesh( triangles, vertices, masses, sphere_triangles, sphere_vertices, sphere_masses );
}

// ---------------------------------------------------------


void ScriptInit::parse_dumbbell( const ParseTree& dumbbell_branch )
{
    double domain_dx;
    dumbbell_branch.get_number( "domain_dx", domain_dx );
    
    Vec3d centre_a;
    dumbbell_branch.get_vec3d( "centre_a", centre_a );
    
    Vec3d centre_b;
    dumbbell_branch.get_vec3d( "centre_b", centre_b );
    
    double sphere_radius;
    dumbbell_branch.get_number( "sphere_radius", sphere_radius );
    
    double handle_width;
    dumbbell_branch.get_number( "handle_width", handle_width );   
    
    Vec3d domain_low = min_union( centre_a, centre_b ) - 2.0 * Vec3d( sphere_radius );
    Vec3d domain_high = max_union( centre_a, centre_b ) + 2.0 * Vec3d( sphere_radius );
    
    Array3d phi;   
    create_dumbbell_signed_distance( centre_a, centre_b, sphere_radius, handle_width, domain_dx, domain_low, domain_high, phi );      
    
    std::vector<Vec3st> new_tris;
    std::vector<Vec3d> new_verts;
    contour_phi( domain_low, domain_dx, phi, new_tris, new_verts );
    project_to_exact_dumbbell( new_verts, centre_a, centre_b, sphere_radius, handle_width );
    
    std::vector<double> new_masses( new_verts.size(), 1.0 );
    std::vector<Vec3d> new_velocities( new_verts.size(), Vec3d(0,0,0) );
    append_mesh( triangles, vertices, masses, new_tris, new_verts, new_masses );
}


// ---------------------------------------------------------


// funcion for get the last line of vtk file by linlin 062420
// std::string getLastLineInFile(std::ifstream& fin)
// {

//         fin.seekg(-2,std::ios_base::end);                // go to one spot before the EOF
//         // jump the last line of the data since it is empty 

//         bool keepLooping = true;
//         while(keepLooping) {
//             char ch;
//             fin.get(ch);                            // Get current byte's data

//             if((int)fin.tellg() <= 1) {             // If the data was at or before the 0th byte
//                 fin.seekg(0);                       // The first line is the last line
//                 keepLooping = false;                // So stop there
//             }
//             else if(ch == '\n') {                   // If the data was a newline
//                 keepLooping = false;                // Stop at the current position.
//             }
//             else {                                  // If the data was neither a newline nor at the 0 byte
//                 fin.seekg(-2,std::ios_base::cur);        // Move to the front of that data, then to the front of the data before it
//             }
//         }

//         std::string lastLine;            
//         getline(fin,lastLine);                      // Read the current line        
//         fin.close();

//         std::cout << "get last line from file " << lastLine << std::endl;


//         return lastLine;
// }

void ScriptInit::parse_script( const char* filename )
{
    
    std::ifstream filestream( filename );
    if ( !filestream.good() )
    {
        std::cerr << "Could not open script file" << std::endl;
        exit(1);
    }
    
    std::cout << "script file: " << filename << std::endl;
    
    ParseTree tree;
    parse_stream( filestream, tree );
    
    
    //
    // Frame stepper
    //
    
    bool ok = tree.get_number( "frame_dt", frame_dt );
    assert( ok );
    
    int num_substeps;
    bool substeps_specified = tree.get_int( "num_substeps", num_substeps );
    if ( substeps_specified )
    {
        sim_dt = frame_dt / (double) num_substeps;
    }
    
    double read_sim_dt;
    if ( tree.get_number( "sim_dt", read_sim_dt ) )
    {
        if ( substeps_specified )
        {
            std::cerr << "WARNING: Both sim_dt and num_substeps specified in config script.  Going with sim_dt." << std::endl;
        }
        
        sim_dt = read_sim_dt;
    }
    
    tree.get_number( "end_sim_t", end_sim_t );
    
    curr_t_specified = tree.get_number( "curr_t", curr_t );
    
    
    //
    // File output
    //
    
    if ( tree.get_string( "output_path", output_path ) )
    {
        output_path_is_relative = false;
    }
    else if ( tree.get_string( "relative_output_path", output_path ) )
    {
        output_path_is_relative = true;
    }         
    else
    {
        // no path specified
        output_path_is_relative = false;
        output_path = std::string( "./" );
    }
    
    //
    // OpenGL camera
    //
    
    const ParseTree* camera_branch = tree.get_branch( "camera" );
    parse_camera( *camera_branch );
    
    
    //
    // Surface geometry
    //
    
    {
        unsigned int curved_sheet_n = 0;
        char curved_sheet_name[256];
        snprintf( curved_sheet_name, 256, "curved_sheet%d", curved_sheet_n );
        const ParseTree* curved_sheet_branch = tree.get_branch( curved_sheet_name );
        
        while ( curved_sheet_branch != NULL )
        {
            parse_curved_sheet( *curved_sheet_branch );
            curved_sheet_branch = NULL;
            ++curved_sheet_n;
            snprintf( curved_sheet_name, 256, "curved_sheet%d", curved_sheet_n );
            curved_sheet_branch = tree.get_branch( curved_sheet_name );
        }
    }
    
    {
        unsigned int sheet_n = 0;   
        char sheet_name[256];
        snprintf( sheet_name, 256, "sheet%d", sheet_n );
        const ParseTree* sheet_branch = tree.get_branch( sheet_name );
        
        while ( sheet_branch != NULL )
        {
            parse_sheet( *sheet_branch );      
            sheet_branch = NULL;      
            ++sheet_n;
            snprintf( sheet_name, 256, "sheet%d", sheet_n );
            sheet_branch = tree.get_branch( sheet_name );
        }
    }
    
    {
        unsigned int sphere_n = 0;   
        char sphere_name[256];
        snprintf( sphere_name, 256, "sphere%d", sphere_n );
        const ParseTree* sphere_branch = tree.get_branch( sphere_name );
        
        while ( sphere_branch != NULL )
        {
            parse_sphere( *sphere_branch );
            sphere_branch = NULL;
            ++sphere_n;
            snprintf( sphere_name, 256, "sphere%d", sphere_n );
            sphere_branch = tree.get_branch( sphere_name );
        }
    }
    
    
    const ParseTree* sphere_branch = tree.get_branch( "sphere" );
    if ( sphere_branch != NULL )
    {
        parse_sphere( *sphere_branch );
    }      
    
    const ParseTree* dumbbell_branch = tree.get_branch( "dumbbell" );
    if ( dumbbell_branch != NULL )
    {
        parse_dumbbell( *dumbbell_branch );
    }      
    
    const ParseTree* trimesh_branch = tree.get_branch( "trimesh" );
    if ( trimesh_branch != NULL )
    {
        printf("Found trimesh branch\n");
        
        std::string meshpath;
        trimesh_branch->get_string("filepath", meshpath);
        printf("Got path: %s\n", meshpath.c_str());
        
        NonDestructiveTriMesh trimesh;
        
        printf("Reading file\n");
        
        std::vector<Vec3d> input_vertices;
        std::vector<double> in_masses;
        read_binary_file( trimesh, input_vertices, in_masses, curr_t, meshpath.c_str() );
        curr_t_specified = true;
        
        int is_solid = 0;
        if ( trimesh_branch->get_int( "is_solid", is_solid ) )
        {
            in_masses.clear();
            if ( is_solid )
            {
                in_masses.resize( input_vertices.size(), std::numeric_limits<double>::infinity() );
            }
            else
            {
                in_masses.resize( input_vertices.size(), 1.0 );
            }
        }
        
        Vec3d translate;
        if ( trimesh_branch->get_vec3d("translate", translate) )
        {
            for ( size_t i = 0; i < input_vertices.size(); ++i )
            {
                input_vertices[i] += translate;
            }
        }
                
        append_mesh( triangles, vertices, masses, trimesh.get_triangles(), input_vertices, in_masses );        
        
        printf("loaded file %s", meshpath.c_str());
    }
    
    
    const ParseTree* obj_branch = tree.get_branch( "objfile" );
    if ( obj_branch != NULL )
    {
        printf("Found obj branch\n");
        
        std::string meshpath;
        obj_branch->get_string("filepath", meshpath);
        printf("Got path: %s\n", meshpath.c_str());
        
        NonDestructiveTriMesh trimesh;
        
        printf("Reading file\n");
        
        std::vector<Vec3d> obj_vertices;
        read_objfile( trimesh, obj_vertices, meshpath.c_str() );
        
        std::vector<Vec3st> obj_triangles = trimesh.get_triangles();
        
        Vec3d translate;
        if ( obj_branch->get_vec3d("translate", translate) )
        {
            for ( size_t i = 0; i < obj_vertices.size(); ++i )
            {
                obj_vertices[i] += translate;
            }
        }
        
        std::vector<double> obj_masses(0);
        int is_solid = 0;
        obj_branch->get_int( "is_solid", is_solid );
        
        if ( is_solid )
        {
            obj_masses.resize( obj_vertices.size(), std::numeric_limits<double>::infinity() );
        }
        else
        {
            obj_masses.resize( obj_vertices.size(), 1.0 );
        }
                
        append_mesh( triangles, vertices, masses, obj_triangles, obj_vertices, obj_masses );
    }

    const ParseTree* ecomsol_branch = tree.get_branch( "ecomsolfile" );
    if ( ecomsol_branch != NULL )
    {
        printf("Found ecomsol branch\n");
        
        std::string meshpath;
        ecomsol_branch->get_string("filepath", meshpath);
        printf("Got path: %s\n", meshpath.c_str());
        
        NonDestructiveTriMesh trimesh;
        
        printf("Reading file\n");
        
        std::vector<Vec3d> obj_vertices;
        //read_ecomsol( trimesh, obj_vertices, meshpath.c_str() );

        // NOTE changed to include concentrations
        // std::vector<double> in_c_B;
        // std::vector<double> in_c_C;
        // std::vector<double> in_c_BC;
        // std::vector<double> in_c_N;
        // std::vector<double> in_c_BN;
        // std::vector<double> in_c_S;
        // std::vector<double> in_c_Tld;
        //read_ecomsol( trimesh, obj_vertices,in_c_B,in_c_C,in_c_BC,in_c_N,in_c_BN,in_c_S,in_c_Tld, meshpath.c_str() );        
        std::vector<std::vector<double> > in_c_species;
        
        read_ecomsol(trimesh, obj_vertices,in_c_species, meshpath.c_str() );

        std::vector<Vec3st> obj_triangles = trimesh.get_triangles();
        
        
        
        std::vector<double> obj_masses(0);
        obj_masses.resize( obj_vertices.size(), 1.0 );
        
        // NOTE changed to account for concentrations
        //append_mesh( triangles, vertices, masses, c_B, c_C, c_BC, c_N, c_BN,c_S,c_Tld, obj_triangles, obj_vertices, obj_masses, in_c_B, in_c_C, in_c_BC ,in_c_N, in_c_BN,in_c_S,in_c_Tld );
        append_mesh( triangles, vertices, masses,c_species,obj_triangles, obj_vertices, obj_masses, in_c_species);
    }
    
    //
    // SurfTrack parameters
    //
    
    const ParseTree* surftrack_branch = tree.get_branch( "surftrack_parameters" );
    parse_surftrack_parameters( *surftrack_branch );
    
    
    //
    // Mesh drivers
    //
    
    
    const ParseTree* invagination_sim_branch = tree.get_branch( "invagination_simulation" );
    if ( invagination_sim_branch != NULL )
    {
        parse_invagination( *invagination_sim_branch );
    }
    
    const ParseTree* zebrafish_diffusion_branch = tree.get_branch("zebrafish_parameters");
    if(zebrafish_diffusion_branch != NULL)
    {
        double DN;//= 10;
        double DBC;//= 4.7475;
        double DBN;//= 7.4747;
        double j1;// = 6.33;
        double j2;// = 55.04;
        double j3 ;//= 38.51;
        double k2 ;//= 0.7539;
        double k_2 ;//= 0.7539;
        double k3  ;// = 0.6044;
        double k_3 ;//= 0.06044;
        double decN;//= 9.3658e-4;
        double decBC;// = 4.7279e-4;
        double decBN;// = 6.2045e-4;
        double lambda_tld_C;// = 6571;
        double lambda_tld_BC;// = 5353;

        double decS; // siz decay
        double kit; // tld feedback with sizzled 
        double kmt;// tld feedback with C+BC
        double Vs;  // siz feedback
        double n;   // siz feedback
        double ks; // siz feedback

        double DTld; //Tld diffusion
        double decTld; //Tld decay
        double j4;    //Tld expression
        double species_n;
        double mutant_bmpmax;
        double mutant_szlmax;
        double szl_prameter_sampling;

        zebrafish_diffusion_branch->get_number( "species_n", species_n );
        zebrafish_diffusion_branch->get_number( "DN", DN );
        zebrafish_diffusion_branch->get_number( "DBC", DBC );
        zebrafish_diffusion_branch->get_number( "DBN", DBN );
        zebrafish_diffusion_branch->get_number( "j1", j1 );
        zebrafish_diffusion_branch->get_number( "j2", j2 );
        zebrafish_diffusion_branch->get_number( "j3", j3 );
        zebrafish_diffusion_branch->get_number( "k2", k2 );
        zebrafish_diffusion_branch->get_number( "k_2", k_2 );
        zebrafish_diffusion_branch->get_number( "k3", k3 );
        zebrafish_diffusion_branch->get_number( "k_3", k_3 );
        zebrafish_diffusion_branch->get_number( "decN", decN );
        zebrafish_diffusion_branch->get_number( "decBC", decBC );
        zebrafish_diffusion_branch->get_number( "decBN", decBN );
        zebrafish_diffusion_branch->get_number( "lambda_tld_C", lambda_tld_C );
        zebrafish_diffusion_branch->get_number( "lambda_tld_BC", lambda_tld_BC );

        zebrafish_diffusion_branch->get_number( "decS", decS );
        zebrafish_diffusion_branch->get_number( "kit", kit );
        zebrafish_diffusion_branch->get_number( "kmt", kmt );
        zebrafish_diffusion_branch->get_number( "Vs", Vs );
        zebrafish_diffusion_branch->get_number( "n", n );
        zebrafish_diffusion_branch->get_number( "ks", ks );
        zebrafish_diffusion_branch->get_number( "DTld", DTld );
        zebrafish_diffusion_branch->get_number( "decTld", decTld );
        zebrafish_diffusion_branch->get_number( "j4", j4 );
        zebrafish_diffusion_branch->get_number( "mutant_bmpmax", mutant_bmpmax );
        zebrafish_diffusion_branch->get_number( "mutant_szlmax", mutant_szlmax );
        zebrafish_diffusion_branch->get_number( "szl_prameter_sampling", szl_prameter_sampling );

        //std::cout<<"mutant_szlmax = "<<mutant_szlmax<<"\n";
        int max_dir = 0; // parameter for get the current max folder
        // find the last previous result 
        DIR *dp;
        struct dirent *ep;   
        dp = opendir ( output_path.c_str() );
                    
        
        if (dp != NULL)
        {      
            while ( (ep = readdir (dp)) != 0 )
            {
                
                std::string name_string( ep->d_name );
                
                //if ( ep->d_type == DT_DIR )    // doesn't exist on all platforms
                {
                    int i;
                    int num_variables_read = sscanf (name_string.c_str(), "run%04d", &i );
                    if ( num_variables_read > 0 )
                    {
                        max_dir = max( max_dir, i );
                    }
                }
            }
            
            (void) closedir (dp);
        }
        else
        {
            std::cout << "Couldn't open the specified output directory" << std::endl;
        }

        //------------------- chack mutant parameter------------------------------------
        if (mutant_bmpmax==1)
        { 
            j2 = 0;
            Vs = 0;
            j4 = 0;

        }  

        if (mutant_szlmax==1)
        { 
            j2 = 0;
            j4 = 0;

            // set up the number for last frame
            char filename_max[1024] ;
            sprintf( filename_max, "%s/run%04d/max_bmpszl.txt", output_path.c_str(), max_dir);

            std::cout << "file path :  " <<filename_max<< std::endl;
            std::ifstream fs;
            fs.open(filename_max, std::fstream::in);

            // check if file opened
            if(fs.is_open())
            {
                std::cout << "Could open file " <<filename_max<< std::endl;
            }
            else{
                std::cout << "Could not open file" << std::endl;
            }            
            std::cout<<"----start reading max bmp/szl level from previous simulaion---"<<"\n"; // Prints our STRING.
            std::string STRING1;
            std::string STRING2;
            char *stopstring1, *stopstring2;

            getline(fs,STRING1); // Saves the line in STRING.
            std::cout<<"STRING1"<<STRING1<<"\n"; // Prints our STRING.
            std::stringstream ss1(STRING1);
            std::string s1; ss1 >> s1;
            double max_bmp_in; ss1 >> max_bmp_in;
               
            getline(fs,STRING2); // Saves the line in STRING.
            std::cout<<"STRING2"<<STRING2<<"\n"; // Prints our STRING.
            std::stringstream ss2(STRING2);
            std::string s2; ss2 >> s2;
            double max_szl_in; ss2 >> max_szl_in;

            std::cout << "s1=" << s1<<"\n";
            std::cout << "max_bmp_in=" << max_bmp_in<<"\n";
            std::cout << "s2=" << s2<<"\n";
            std::cout << "max_szl_in=" << max_szl_in<<"\n";

            std::cout<<"ks read in = "<<ks<<"\n";

            ks = ks * max_bmp_in;
            std::cout<<"ks * max_bmp_in= "<<ks<<"\n";
            std::cout<<"mutant_szlmax = "<<mutant_szlmax<<"\n";
        }
         

        if (szl_prameter_sampling!=0)
        { 

            // set up the number for first chd mutant for bmp max  
            char filename_max[1024] ;
            sprintf( filename_max, "%s/run%04d/max_bmpszl.txt", output_path.c_str(), max_dir-(int)szl_prameter_sampling);

            std::cout << "szl_prameter_sampling=" << szl_prameter_sampling<<"\n";
            std::cout << "file path :  " <<filename_max<< std::endl;
            std::ifstream fs;
            fs.open(filename_max, std::fstream::in);

            // check if file opened
            if(fs.is_open())
            {
                std::cout << "Could open file " <<filename_max<< std::endl;
            }
            else{
                std::cout << "Could not open file" << std::endl;
            }
            
            std::cout<<"----start reading max bmp level from previous simulaion---"<<"\n"; // Prints our STRING.
            std::string STRING1;
            std::string STRING2;
            char *stopstring1, *stopstring2;

            getline(fs,STRING1); // Saves the line in STRING.
            std::cout<<"STRING1"<<STRING1<<"\n"; // Prints our STRING.
            std::stringstream ss1(STRING1);
            std::string s1; ss1 >> s1;
            double max_bmp_in; ss1 >> max_bmp_in;
            //max_bmp_in = strtod(STRING1.c_str(), &stopstring1); 
            fs.close();

            // getline(fs,STRING2); // Saves the line in STRING.
            // std::cout<<"STRING2"<<STRING2<<"\n"; // Prints our STRING.
            // std::stringstream ss2(STRING2);
            // std::string s2; ss2 >> s2;
            // double max_szl_in; ss2 >> max_szl_in;
            // //max_szl_in = strtod(STRING2.c_str(), &stopstring2); 

            std::cout << "s1=" << s1<<"\n";
            std::cout << "max_bmp_in=" << max_bmp_in<<"\n";
            // std::cout << "s2=" << s2<<"\n";
            // std::cout << "max_szl_in=" << max_szl_in<<"\n";

            std::cout<<"ks read in = "<<ks<<"\n";
            ks = ks * max_bmp_in;
            std::cout<<"ks * bmpmax = "<<ks<<"\n";

            

            // set up the number for first chd mutant for szl max  
            // char filename[1024] ;
            sprintf( filename_max, "%s/run%04d/max_bmpszl.txt", output_path.c_str(), max_dir-(int)szl_prameter_sampling+1);

            std::cout << "file path :  " <<filename_max<< std::endl;
            std::ifstream fs1;
            fs1.open(filename_max, std::fstream::in);

            // check if file opened
            if(fs1.is_open())
            {
                std::cout << "Could open file " <<filename_max<< std::endl;
            }
            else{
                std::cout << "Could not open file" << std::endl;
            }
            
            std::cout<<"----start reading max szl level from previous simulaion---"<<"\n"; // Prints our STRING.
            // std::string STRING1;
            // std::string STRING2;
            // char *stopstring1, *stopstring2;

            getline(fs1,STRING1); // Saves the line in STRING.
            std::cout<<"STRING1"<<STRING1<<"\n"; // Prints our STRING.
            std::stringstream ss11(STRING1);
            std::string s11;
            ss11 >> s1;
            double max_bmp_notuse; ss11 >> max_bmp_notuse;
 
             

            getline(fs1,STRING2); // Saves the line in STRING.
            std::cout<<"STRING2"<<STRING2<<"\n"; // Prints our STRING.
            std::stringstream ss2(STRING2);
            std::string s2; 
            ss2 >> s2;
            double max_szl_in; ss2 >> max_szl_in;
            fs1.close();

            std::cout << "s1=" << s1<<"\n";
            std::cout << "max_bmp_notuse=" << max_bmp_notuse<<"\n";
            std::cout << "s2=" << s2<<"\n";
            std::cout << "max_szl_in=" << max_szl_in<<"\n";


            std::random_device rd;
            std::default_random_engine generator(rd()); // rd() provides a random seed
            std::uniform_real_distribution<double> distribution(0.1,10);
            //double number = distribution(generator);
            //std::cout << "rand_numbere =  " << number <<"\n";
            kit = max_szl_in * kit; 
            std::cout << "kit = max_szl_in * rand_numbere = " << kit <<"\n";
        }

        diffusionParam.clear();
        diffusionParam.push_back(species_n);
        diffusionParam.push_back(DN);
        diffusionParam.push_back(DBC);
        diffusionParam.push_back(DBN);
        diffusionParam.push_back(j1);
        diffusionParam.push_back(j2);
        diffusionParam.push_back(j3);
        diffusionParam.push_back(k2);
        diffusionParam.push_back(k_2);
        diffusionParam.push_back(k3);
        diffusionParam.push_back(k_3);
        diffusionParam.push_back(decN);
        diffusionParam.push_back(decBC);
        diffusionParam.push_back(decBN);
        diffusionParam.push_back(lambda_tld_C);
        diffusionParam.push_back(lambda_tld_BC);

        diffusionParam.push_back(decS);
        diffusionParam.push_back(kit);
        diffusionParam.push_back(kmt);
        diffusionParam.push_back(Vs);
        diffusionParam.push_back(n);
        diffusionParam.push_back(ks);
        diffusionParam.push_back(DTld);
        diffusionParam.push_back(decTld);
        diffusionParam.push_back(j4);
        diffusionParam.push_back(mutant_bmpmax);
        diffusionParam.push_back(mutant_szlmax);
        diffusionParam.push_back(szl_prameter_sampling);    

        // save the parameter used in the simulation (Linlin Li 070620)
        char data_filename[256];
        sprintf(data_filename, "%s/parameter_list_run%04d.txt", output_path.c_str(), max_dir+1); 
        std::ofstream outfile(data_filename,std::ios::trunc);
        outfile<<"submit_filename "<< filename <<"\n";
        outfile<<"species_n "<< species_n <<"\n";
        outfile<<"DN "<< DN <<"\n";
        outfile<<"DBC "<< DBC <<"\n";
        outfile<<"DBN "<< DBN <<"\n";
        outfile<<"j1 "<< j1 <<"\n";

        outfile<<"j2 "<< j2 <<"\n";
        outfile<<"j3 "<< j3 <<"\n";
        outfile<<"k2 "<< k2 <<"\n";
        outfile<<"k_2 "<< k_2 <<"\n";
        outfile<<"k3 "<< k3 <<"\n";

        outfile<<"k_3 "<< k_3 <<"\n";
        outfile<<"decN "<< decN <<"\n";
        outfile<<"decBC "<< decBC <<"\n";
        outfile<<"decBN "<< decBN <<"\n";
        outfile<<"lambda_tld_C "<< lambda_tld_C <<"\n";

        outfile<<"lambda_tld_BC "<< lambda_tld_BC <<"\n";
        outfile<<"decS "<< decS <<"\n";
        outfile<<"kit "<< kit <<"\n";
        outfile<<"kmt "<< kmt <<"\n";
        outfile<<"Vs "<< Vs <<"\n";

        outfile<<"n "<< n <<"\n";
        outfile<<"ks "<< ks <<"\n";
        outfile<<"DTld "<< DTld <<"\n";
        outfile<<"decTld "<< decTld <<"\n";
        outfile<<"j4 "<< j4 <<"\n";

        outfile<<"mutant_bmpmax "<< mutant_bmpmax <<"\n";
        outfile<<"mutant_szlmax "<< mutant_szlmax <<"\n";
        outfile<<"szl_prameter_sampling "<< szl_prameter_sampling <<"\n";



        outfile.close ();    
        

    }
        
}

