// ---------------------------------------------------------
//
//  main.cpp
//  Tyson Brochu 2008
//  Adrian Buganza 2017
//  Linlin Li 2019
//  update note: adding sizzled 
//
//  Functions for setting up and running an explicit surface simulation, 
//  as well as visualizing it with gluvi.
//  And then driving the deformation of a sphere with prescribed velocities.
//
// ---------------------------------------------------------


// ---------------------------------------------------------
// Defines
// ---------------------------------------------------------

// Whether to use an OpenGL GUI, or run command-line style
// (Define NO_GUI in your build to suppress the GUI.)

//#define NO_GUI

#ifndef NO_GUI
#define USE_GUI
#endif

// Whether to run rendering and simulation on separate threads
#define RUN_ASYNC

// Whether to use the El Topo C API.  If not defined, uses the El Topo classes directly.
//#define USE_C_API


// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

// std
#include <cstdio>
#include <fenv_include.h>
#include <fstream>
#include <vector>
#include <queue>

// common
#include <array2.h>
#include <ccd_wrapper.h>
#include <collisionqueries.h>
#include <expansion.h>
#include <marching_tiles_hires.h>
#include <util.h>
#include <vec.h>
#include <wallclocktime.h>

// el topo
#include <collisionpipeline.h>
#include <eltopo.h>
#include <iomesh.h>
#include <meshrenderer.h>
#include <runstats.h>
#include <surftrack.h>
#include <trianglequality.h>

// zebrafish
#include <framestepper.h>
#include <meshdriver.h>
#include <invagination.h>
#include <scriptinit.h>
#include <simulation.h>

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>    //for windows, you'll need http://www.softagalleria.net/dirent.php or something similar
#include <sys/stat.h>
// for count the execution time
#include <time.h>
// for stop the programs 

// for find the max number (add by Linlin Li, Jun 27,2020)
#include <algorithm>

#ifdef _MSC_VER

#include <pthread.h> //requires use pthreads-win from http://sources.redhat.com/pthreads-win32/, or something like it

//this is kind of hacky, but seems to do the trick for now (on Windows, the name is prefaced by an underscore)
//so we'll just rename is here
#define snprintf _snprintf 

#endif

// conditional includes
#ifdef USE_GUI 
#include <gluvi.h> 
#endif

void advance_diffusion(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, std::vector<std::vector<double> > &species, double dt,double t, const std::vector<double> &diffusionParam, int &iteration, int &max_it,int &force_postive_check);

// ---------------------------------------------------------
// Global interface declarations
// ---------------------------------------------------------

extern "C" {
    void exactinit();    // in predicates.c
}

int main(int argc, char **argv); // defined in ~1550

// ---------------------------------------------------------
// Global variable
// ---------------------------------------------------------

RunStats g_stats;

char g_output_path[256];  // Where to write output data
char g_base_output_path[256];  // Where to write output data
bool g_render_to_texture;

#ifdef USE_GUI

float camera_target[] = {0.0f, 0.0f, 0.0f};    // OpenGL camera points here
Gluvi::Target3D gluvi_cam( camera_target, 10.0f, 0.0f, 0.0f, 45.0f, 0.1f, 20.0f );   // OpenGL camera
bool g_take_screenshots = false;
#else

struct hack_camera
{
    float target[3];
    float dist;
    float pitch;
    float heading;
} gluvi_cam;

#endif



SurfTrack* g_surf = NULL;        // The dynamic surface

// global variable for diffusion parameters
std::vector<double> diffusionParam;

// global variable for check if calculate diffusion edit linlin 091620
int diffusion_off = 0;

std::vector<Vec3d> renderable_vertices;
std::vector<Vec3d> renderable_vertex_normals;
std::vector<Vec3st> renderable_solid_triangles;
std::vector<Vec3st> renderable_non_solid_triangles;
std::vector<Vec2st> renderable_edges;


// Local to main.cpp
namespace {    

    MeshDriver* driver = NULL;     // The thing that makes the dynamic surface go

    // ---------------------------------------------------------
    // Interface declarations for functions local to this file
    // ---------------------------------------------------------
    
    void advance_sim( double &dt ); // defined in ~434
    void advance_frame(); // defined in ~616
    void run_simulation( ); // defined in ~650
    void set_up_output_path( ); // defined in ~1390
    
#ifdef USE_GUI
    enum { RAY_HIT_VERTEX, RAY_HIT_EDGE, RAY_HIT_TRIANGLE, RAY_HIT_NOTHING };
    
    void update_renderable_objects();
    void* advance_frame_async( void* nothing );
    void start_advance_frame();
    void advance_frame_done();
    void timer_advance_frame( int junk );
    unsigned int ray_cast( const DynamicSurface& surface, const Vec3f& ray_origin, const Vec3f& ray_direction, unsigned int primitive_type, size_t& hit_index );          
    void keyboard(unsigned char key, int x, int y);
    void mouse(int button, int state, int x, int y);
    void mouseMotion(int x, int y);
    void display(void);
    void init_gui( const ScriptInit& script_init );
#endif
    
    void init_simulation( char *script_filename ); // defined in ~1450
    
    // ---------------------------------------------------------
    // Variables local to this file
    // ---------------------------------------------------------
    
    
    Simulation* sim = NULL;        // Timekeeping
    FrameStepper* frame_stepper = NULL;
    
    
#ifdef USE_GUI
    
    bool display_driver = false;
    MeshRenderer mesh_renderer;
    const MeshDriver::MeshDriverRenderer* driver_renderer;
    Gluvi::DynamicText* status_text_widget;
    
    // Mouse click / ray casting data
    ptrdiff_t selected_vertex = -1;
    ptrdiff_t selected_edge = -1;
    ptrdiff_t selected_triangle = -1;
    
    //
    // async stuff
    //
    
    bool advance_single_frame = false;
    bool display_status_text = true;
        
    pthread_mutex_t thread_is_running_mutex = PTHREAD_MUTEX_INITIALIZER;
    bool thread_is_running = false;
    pthread_t advance_frame_thread;
    bool waiting_for_thread_to_finish = false;
    pthread_mutex_t surf_mutex = PTHREAD_MUTEX_INITIALIZER;
    
#endif  // def USE_GUI
    
    
    // ---------------------------------------------------------
    // Function definitions
    // ---------------------------------------------------------
    
    
    
#ifdef USE_C_API
    
    // ---------------------------------------------------------
    ///
    /// Wrapper for El Topo API to mesh improvement/topology change functions.
    ///
    /// This function is terribly inefficient: it copies the current state of g_surf to buffers, then calls the C API, which internally
    /// creates a new SurfTrack object.  This function is only included for demonstrating and testing the C API.
    ///
    // ---------------------------------------------------------
    
    void el_topo_wrapper_static_operations()
    {
        
        // Build input & options structures
        
        ElTopoMesh inputs;
        inputs.num_vertices = g_surf->get_num_vertices();
        inputs.vertex_locations = new double[ 3*inputs.num_vertices ];
        memcpy( inputs.vertex_locations, &(g_surf->get_position(0)), 3*inputs.num_vertices*sizeof(double) );   
        
        inputs.num_triangles = g_surf->m_mesh.num_triangles();
        inputs.triangles = new int[ 3*inputs.num_triangles ];
        for ( int i = 0; i < inputs.num_triangles; ++i )
        {
            inputs.triangles[3*i+0] = g_surf->m_mesh.get_triangle(i)[0];
            inputs.triangles[3*i+1] = g_surf->m_mesh.get_triangle(i)[1];
            inputs.triangles[3*i+2] = g_surf->m_mesh.get_triangle(i)[2];
        }
        
        inputs.vertex_masses = new double[ inputs.num_vertices ];
        memcpy( inputs.vertex_masses, &(g_surf->m_masses[0]), inputs.num_vertices*sizeof(double) );
        
        ElTopoGeneralOptions general_options;
        general_options.m_verbose = g_surf->m_verbose;
        general_options.m_collision_safety = g_surf->m_collision_safety;
        general_options.m_proximity_epsilon = g_surf->m_proximity_epsilon;
        
        ElTopoStaticOperationsOptions options;
        options.m_perform_improvement = g_surf->m_perform_improvement;
        options.m_allow_topology_changes = g_surf->m_allow_topology_changes;
        options.m_max_volume_change = g_surf->m_max_volume_change;
        options.m_min_edge_length = g_surf->m_min_edge_length;
        options.m_max_edge_length = g_surf->m_max_edge_length;
        options.m_min_triangle_area = g_surf->m_min_triangle_area;
        options.m_min_triangle_angle = g_surf->m_min_triangle_angle;
        options.m_max_triangle_angle = g_surf->m_max_triangle_angle;
        options.m_use_curvature_when_splitting = g_surf->m_splitter.m_use_curvature;
        options.m_use_curvature_when_collapsing = g_surf->m_collapser.m_use_curvature;
        options.m_min_curvature_multiplier = g_surf->m_collapser.m_min_curvature_multiplier;
        options.m_max_curvature_multiplier = g_surf->m_splitter.m_max_curvature_multiplier;
        options.m_allow_vertex_movement = g_surf->m_allow_vertex_movement;
        options.m_edge_flip_min_length_change = g_surf->m_edge_flip_min_length_change;
        options.m_merge_proximity_epsilon = g_surf->m_merge_proximity_epsilon;
        options.m_subdivision_scheme = g_surf->m_subdivision_scheme;
        
        ElTopoMesh outputs;
        ElTopoDefragInformation defrag_info;
        
        el_topo_static_operations( &inputs, &general_options, &options, &defrag_info, &outputs );
        
        {
            g_surf->improve_mesh();
            g_surf->topology_changes();
            g_surf->defrag_mesh();
            
            assert( g_surf->m_mesh.num_triangles() == outputs.num_triangles );
            assert( g_surf->get_num_vertices() == outputs.num_vertices );      
        }
        
        
        //
        // Update surface info
        //
        
        g_surf->set_all_positions( outputs.num_vertices, outputs.vertex_locations );
        
        g_surf->m_masses.resize( outputs.num_vertices );
        for( unsigned int i = 0; i < outputs.num_vertices; ++i )
        {
            g_surf->m_masses[i] = outputs.vertex_masses[i];
        }
        
        std::vector<Vec3ui> new_triangles( outputs.num_triangles );
        for ( size_t i = 0; i < new_triangles.size(); ++i )
        {
            new_triangles[i][0] = outputs.triangles[3*i+0];
            new_triangles[i][1] = outputs.triangles[3*i+1];
            new_triangles[i][2] = outputs.triangles[3*i+2];      
        }
        
        g_surf->m_mesh.replace_all_triangles( new_triangles );
        
        //
        // Send change history to mesh driver
        //
        
        g_surf->m_vertex_change_history.resize( defrag_info.num_vertex_changes );
        
        for ( unsigned int i = 0; i < defrag_info.num_vertex_changes; ++i )
        {
            g_surf->m_vertex_change_history[i].is_remove = defrag_info.vertex_is_remove[i];
            g_surf->m_vertex_change_history[i].vertex_index = defrag_info.vertex_index[i];
            g_surf->m_vertex_change_history[i].split_edge[0] = defrag_info.split_edge[2*i+0];
            g_surf->m_vertex_change_history[i].split_edge[1] = defrag_info.split_edge[2*i+1];
        }
        
        g_surf->m_triangle_change_history.resize( defrag_info.num_triangle_changes );
        
        for ( unsigned int i = 0; i < defrag_info.num_triangle_changes; ++i )
        {
            g_surf->m_triangle_change_history[i].is_remove = defrag_info.triangle_is_remove[i];
            g_surf->m_triangle_change_history[i].triangle_index = defrag_info.triangle_index[i];
            g_surf->m_triangle_change_history[i].tri[0] = defrag_info.new_tri[3*i+0];
            g_surf->m_triangle_change_history[i].tri[1] = defrag_info.new_tri[3*i+1];
            g_surf->m_triangle_change_history[i].tri[2] = defrag_info.new_tri[3*i+2];
        }
        
        g_surf->m_defragged_triangle_map.resize( defrag_info.defragged_triangle_map_size );
        
        for ( unsigned int i = 0; i < defrag_info.defragged_triangle_map_size; ++i )
        {
            g_surf->m_defragged_triangle_map[i][0] = defrag_info.defragged_triangle_map[2*i+0];
            g_surf->m_defragged_triangle_map[i][1] = defrag_info.defragged_triangle_map[2*i+1];
        }
        
        g_surf->m_defragged_vertex_map.resize( defrag_info.defragged_vertex_map_size );
        
        for ( unsigned int i = 0; i < defrag_info.defragged_vertex_map_size; ++i )
        {
            g_surf->m_defragged_vertex_map[i][0] = defrag_info.defragged_vertex_map[2*i+0];
            g_surf->m_defragged_vertex_map[i][1] = defrag_info.defragged_vertex_map[2*i+1];
        }
        
        driver->update_simulation_elements( *g_surf );
        
        
        // Free allocated memory
        
        delete[] inputs.triangles;
        delete[] inputs.vertex_locations;
        delete[] inputs.vertex_masses;
        
        el_topo_free_static_operations_results( &outputs, &defrag_info );
        
    }
    
    // ---------------------------------------------------------
    ///
    /// Wrapper for El Topo API to mesh integration/collision detection.
    ///
    /// This function is terribly inefficient: it copies the current state of g_surf to buffers, then calls the C API, which internally
    /// creates a new SurfTrack object.  This function is only included for demonstrating and testing the C API.
    ///
    // ---------------------------------------------------------
    
    void el_topo_wrapper_integrate( double dt, double& out_dt )
    {
        
        ElTopoMesh inputs;
        inputs.num_vertices = g_surf->get_num_vertices();
        inputs.vertex_locations = new double[ 3*inputs.num_vertices ];
        memcpy( inputs.vertex_locations, &(g_surf->get_position(0)), 3*inputs.num_vertices*sizeof(double) );   
        
        inputs.num_triangles = g_surf->m_mesh.num_triangles();
        inputs.triangles = new int[ 3*inputs.num_triangles ];
        memcpy( inputs.triangles, &(g_surf->m_mesh.get_triangles()[0]), 3*inputs.num_triangles*sizeof(int) );
        
        inputs.vertex_masses = new double[ inputs.num_vertices ];
        memcpy( inputs.vertex_masses, &(g_surf->m_masses[0]), inputs.num_vertices*sizeof(double) );
        
        assert( g_surf->m_masses.size() == inputs.num_vertices );
        for ( unsigned int i = 0; i < inputs.num_vertices; ++i )
        {
            assert( g_surf->m_masses[i] == g_surf->m_masses[i] );
            assert( g_surf->m_masses[i] != 0.0 );
            assert( inputs.vertex_masses[i] == inputs.vertex_masses[i] );
            assert( inputs.vertex_masses[i] != 0.0  );
        }
        
        ElTopoGeneralOptions general_options;
        general_options.m_verbose = g_surf->m_verbose;
        general_options.m_collision_safety = g_surf->m_collision_safety;
        general_options.m_proximity_epsilon = g_surf->m_proximity_epsilon;
        
        ElTopoIntegrationOptions options;
        options.m_dt = dt;
        options.m_friction_coefficient = g_surf->m_collision_pipeline.m_friction_coefficient;
        
        double* new_vertex_locations = new double[ 3 * g_surf->get_num_vertices() ];
        
        memcpy( new_vertex_locations, &(g_surf->get_newposition(0)), 3 * g_surf->get_num_vertices() * sizeof(double) );
        
        double* out_vertex_locations;
        
        el_topo_integrate( &inputs, new_vertex_locations, &general_options, &options, &out_vertex_locations, &out_dt );
        
        g_surf->set_all_positions( g_surf->get_num_vertices(), out_vertex_locations );
        g_surf->set_all_newpositions( g_surf->get_num_vertices(), out_vertex_locations );
        
        delete[] inputs.triangles;
        delete[] inputs.vertex_locations;
        delete[] inputs.vertex_masses;
        
        delete[] new_vertex_locations;
        
        el_topo_free_integrate_results( out_vertex_locations );
        
        
    }
    
#endif   // USE_C_API
    
    
    // ---------------------------------------------------------
    ///
    /// Advance simulation, GUI-ignorant
    ///
    // ---------------------------------------------------------
    
    void advance_sim( double &sim_dt )
    {
        
        double start_time = get_time_in_seconds();
        
        double accum_dt = 0;


                // CHECK which time step 'dt' was actually advanced in time 
        // this function now is just the same as the matlab, takes in concentrations at 
        // previous time step and gives back new concentrations
        // advance_diffusion(mesh, dt);
        // modified by linlin 
        double time_c = sim->m_curr_t;
        int iteration;
        int diffusion_check;

        if (diffusion_off == 1)
        {
            diffusion_check = 0;
            sim_dt = 250;
        }
        else 
        {
          diffusion_check = 1;  
        }

        double dt_change = 0;
        int force_postive_check = 0;

        int count_iter = 0;

        while (diffusion_check == 1)
        {
            if (sim_dt > 0.00001)
            {
                if (count_iter > 4)
                    {
                        diffusion_check == 0;
                        diffusion_off = 1; // tuned of diffusion calulation for next time step

                          // force to pass the diffusion with error to keep the program running  
                        char data_filename[256];
                        sprintf( data_filename, "%s/diffusion_not_converge.txt", g_output_path);  
                        std::ofstream outfile(data_filename,std::ios::trunc);
                        outfile<<"diffusion not converged ignor the result\n";
                        outfile.close ();

                        break;

                    }
                else if (count_iter > 2)
                    {
                      force_postive_check = 1; 
                      std::cout << "force_postive_check turned on \n"; 

                    }


                std::cout << "count_iter" <<count_iter << "\n";
                int max_iteration = 19;
                //advance_diffusion(g_surf->m_mesh, g_surf->get_positions(), g_surf->get_cB(),g_surf->get_cC(),g_surf->get_cBC(),g_surf->get_cN(),g_surf->get_cBN(),g_surf->get_cS(),g_surf->get_cTld(), sim_dt,time_c,diffusionParam);
                advance_diffusion(g_surf->m_mesh, g_surf->get_positions(), g_surf->get_cspecies(),sim_dt,time_c,diffusionParam,iteration,max_iteration,force_postive_check);

                std::cout << "iteration" <<iteration << "\n";
                std::cout << "sim_dt in diffusion= " <<sim_dt<< "\n";
                if (iteration<4)
                {
                    dt_change = 0.5*sim_dt;
                    diffusion_check = 0;
                }
                else if(iteration >= max_iteration)
                {
                    sim_dt = 0.01*sim_dt;
                    diffusion_check = 1 ;
                    count_iter = count_iter+1;

                }
                else if (iteration > 10 )
                {
                    dt_change = -0.9*sim_dt;
                    diffusion_check = 0 ;
                }

                else{
                    dt_change = 0.0;
                    diffusion_check = 0; 
                }
            }
            else 
            {
                diffusion_check == 0;
                diffusion_off = 1; // tuned of diffusion calulation for next time step

                  // force to pass the diffusion with error to keep the program running  
                char data_filename[256];
                sprintf( data_filename, "%s/diffusion_not_converge.txt", g_output_path);  
                std::ofstream outfile(data_filename,std::ios::trunc);
                outfile<<"sim_dt is too small --------diffusion not converged ignor the result\n";
                outfile.close ();

                break;  
            }    

            
            // if keep looping force  negtive concetration to positive


        }
        // update timestep if iteration <4 or >10 change made by linlin  








        
        while ( (accum_dt < 0.99 * sim_dt) && (sim->m_curr_t + accum_dt < sim->m_max_t) )
        {       
            
            std::cout << "\n\n ------- sim step ------- \n\n";
            
            std::cout << "curr t: " << sim->m_curr_t + accum_dt << " / " << sim->m_max_t << std::endl;
            
            // ----------
            // collect stats
            // ----------
            
            // comment out all g_stats 

            //g_stats.add_to_int( "total_sim_steps", 1 );
            //g_stats.update_min_int( "min_num_triangles", (int) g_surf->m_mesh.num_triangles() ); 
            //g_stats.update_max_int( "max_num_triangles", (int) g_surf->m_mesh.num_triangles() );
            //g_stats.add_to_int( "total_num_vertices", g_surf->get_num_vertices() );
            //g_stats.update_min_int( "min_num_vertices", (int) g_surf->get_num_vertices() );
            //g_stats.update_max_int( "max_num_vertices", (int) g_surf->get_num_vertices() );
            //g_stats.add_to_int( "total_num_triangles", g_surf->m_mesh.num_triangles() );
            
            double start_time_sub_step = get_time_in_seconds();


            // ---------- 
            // mesh maintenance & topology changes
            // ---------- 
            
            
#ifdef USE_C_API
            el_topo_wrapper_static_operations();
#else
            {         
                // Improve
                
                double pre_improve_time = get_time_in_seconds();
                g_surf->improve_mesh();
                double post_improve_time = get_time_in_seconds();
                //g_stats.add_to_double( "total_improve_time", post_improve_time - pre_improve_time );
                
                //g_stats.add_per_frame_double( "frame_improve_time", frame_stepper->get_frame(), post_improve_time - pre_improve_time );

                // NEED to have an array, same way the nodes or mass array is defined, just add one more
                // but instead of node coordinates or masses keep track of concentrations
                // -> added to dynamic surface class, both concentration field and in the constructor

                // NEED to go to the source code of "improve_mesh()" and interpolate B, C and BC values
                // I think is done but might need to double check 

                // Topology changes
                // NOTE: no topology changes for zebrafish

                double pre_topology_time = get_time_in_seconds();
                g_surf->topology_changes();
                double post_topology_time = get_time_in_seconds();
                //g_stats.add_to_double( "total_topology_time", post_topology_time - pre_topology_time );
                
                //g_stats.add_per_frame_double( "frame_topo_time", frame_stepper->get_frame(), post_topology_time - pre_topology_time );
                
                double pre_defrag_time = get_time_in_seconds();
                g_surf->defrag_mesh();
                double post_defrag_time = get_time_in_seconds();
                
                //g_stats.add_per_frame_double( "frame_defrag_time", frame_stepper->get_frame(), post_defrag_time - pre_defrag_time );
                
                // Update driver
                
                driver->update_simulation_elements( *g_surf );
                
            }
#endif
            
            
            double min_angle = min_triangle_angle( *g_surf );   
            //g_stats.add_per_frame_double( "min_angle", frame_stepper->get_frame(), rad2deg(min_angle) );
            double max_angle = max_triangle_angle( *g_surf );
            //g_stats.add_per_frame_double( "max_angle", frame_stepper->get_frame(), rad2deg(max_angle) );
            
            char imp_stats_filename[256];
            sprintf( imp_stats_filename, "%s/aaa-imp-stats.txt", g_output_path );      
            //g_stats.write_to_file( imp_stats_filename );
            
            // ---------- 
            // advance underlying simulation
            // ----------
            
            double curr_dt = sim_dt - accum_dt;
            curr_dt = min( curr_dt, sim->m_max_t - sim->m_curr_t - accum_dt );
            
            std::cout << "curr_dt: " << curr_dt << std::endl;
            
            double time_before_driver = get_time_in_seconds();
            
            std::vector<Vec3d> new_positions( g_surf->get_num_vertices() );
            driver->set_predicted_vertex_positions( *g_surf, new_positions, sim->m_curr_t + accum_dt, curr_dt );
            
            g_surf->set_all_newpositions( new_positions );
            
            std::cout << "curr_dt in topo: " << curr_dt << std::endl;
            
            double time_after_driver = get_time_in_seconds();
            //g_stats.add_to_double( "total_driver_time", time_after_driver - time_before_driver );
            
            // ----------
            // move & handle collision detection
            // ----------
            
            double time_before_integration = get_time_in_seconds();             
            double actual_dt;
            
            std::vector<Vec3d> initial_positions = g_surf->get_positions();
            
#ifdef USE_C_API
            el_topo_wrapper_integrate( curr_dt, actual_dt );
#else
            g_surf->integrate( curr_dt, actual_dt );
#endif
            
            curr_dt = actual_dt;    // the time step actually taken by el topo
            
            std::cout << "actual_dt: " << actual_dt << std::endl;
            
            std::vector<Vec3d> final_positions = g_surf->get_positions();
            
            double curr_integration_time = get_time_in_seconds() - time_before_integration;
            //g_stats.add_to_double( "total_integration_time", curr_integration_time );
            
            accum_dt += curr_dt;
            
            //
            // Need to inform the driver what the actual vertex motion was (e.g. for velocity update)
            //
            
            driver->notify_done_integration( initial_positions, final_positions, actual_dt );
            
            driver->compute_error( *g_surf, sim->m_curr_t + accum_dt );
            

            double end_time_sub_step = get_time_in_seconds();
            //g_stats.add_per_frame_double( "sim_time", frame_stepper->get_frame(), end_time_sub_step - start_time_sub_step );
            
            //g_stats.add_per_frame_int( "num_vertices", frame_stepper->get_frame(), g_surf->get_num_vertices() );
            //g_stats.add_per_frame_int( "num_triangles", frame_stepper->get_frame(), g_surf->m_mesh.num_triangles() );

            //
            // file output
            //
            std::vector<std::vector<double> > c_output = g_surf-> get_cspecies();

            int node_n = c_output.size();
            //std::cout<<"coutput"<< c_output.size()<<"\n";


            std::vector<double> c_B(node_n);
            std::vector<double> c_C(node_n);
            std::vector<double> c_BC(node_n);
            std::vector<double> c_N(node_n);
            std::vector<double> c_BN(node_n);
            std::vector<double> c_S(node_n);
            std::vector<double> c_Tld(node_n);

            for (int i=0; i<node_n;++i){
                c_B[i] = c_output[i][0];
                c_C[i] = c_output[i][1];
                c_BC[i] = c_output[i][2];
                c_N[i] = c_output[i][3];
                c_BN[i] = c_output[i][4];
                c_S[i] = c_output[i][5];
                c_Tld[i] = c_output[i][6];
                }
            //std::cout<<"updated c_B in main"<<c_B[1]<<"\n";
            char binary_filename[256];
            sprintf( binary_filename, "%s/frame%04d.bin", g_output_path, frame_stepper->get_frame() );      
            //write_binary_file( g_surf->m_mesh, g_surf->get_positions(), g_surf->m_masses, sim->m_curr_t, binary_filename );   

            char vtk_filename[256];
            sprintf( vtk_filename, "%s/Aframe%04d.vtk", g_output_path, frame_stepper->get_frame() );      
            // if there are concentrations... 
            write_vtk( g_surf->m_mesh, g_surf->get_positions(), c_B,c_C,c_BC,c_N,sim->m_curr_t, vtk_filename ); 

            char vtk_filename1[256];
            sprintf( vtk_filename1, "%s/Bframe%04d.vtk", g_output_path, frame_stepper->get_frame() );      
            // if there are concentrations... 
            write_vtk( g_surf->m_mesh, g_surf->get_positions(),c_B,c_BN,c_S,c_Tld,sim->m_curr_t, vtk_filename1 ); 
          
            char stats_filename[256];
            sprintf( stats_filename, "%s/aaa-stats.txt", g_output_path );      
            //g_stats.write_to_file( stats_filename );

           
        }



        
        std::cout << " -------------- end sim step -------------- \n" << std::endl;
        
        double sim_step_time = get_time_in_seconds() - start_time;


        //g_stats.add_to_double( "total_sim_time", sim_step_time );
        
        sim->m_curr_t += accum_dt;
        //        
        
        // ----------
        // check if max time is reached
        // ----------
        
        if ( sim->done_simulation() )
        {         
            sim->m_running = false;
            std::cout << "total time steps: " << g_stats.get_int( "total_sim_steps" ) << std::endl;
        }
        
        // change dt based on interation of diffusion _ edit by linlin 
        sim_dt = sim_dt + dt_change;
        



    }
    
    
    // ---------------------------------------------------------
    ///
    /// Advance simulation by one frame
    ///
    // ---------------------------------------------------------
    
    void advance_frame()
    {
        
        if( !sim->m_currently_advancing_simulation )
        {
            sim->m_currently_advancing_simulation = true;
            
            //
            // Advance frame
            //
            
            std::cout << "\n --------------- frame " << frame_stepper->get_frame() << " --------------- \n" << std::endl;
            
            while ( !(frame_stepper->done_frame()) )
            {
                double dt = frame_stepper->get_step_length( sim->m_dt );
                frame_stepper->advance_step( dt ); 
                advance_sim( dt );
                sim->m_dt = dt;
                        
            }
            
            std::cout << " --------------- end frame " << frame_stepper->get_frame() << " --------------- \n" << std::endl;
            
            // update frame stepper      
            frame_stepper->next_frame();
            
            sim->m_currently_advancing_simulation = false;
            
        }
        
    }
    
    
    
    // ---------------------------------------------------------
    ///
    /// Run an entire simulation without GUI.  No threading.
    ///
    // ---------------------------------------------------------
    void run_simulation( )
    {

        sim->m_running = true;
        while ( sim->m_running )
        {
            advance_frame();
        }
        sim->m_running = false;
    }
    
    // =========================================================
    // GUI FUNCTIONS
    // =========================================================
    
#ifdef USE_GUI
    
    // ---------------------------------------------------------
    // ---------------------------------------------------------
    
    void update_renderable_objects()
    {
        // aquire mutex on fluid sim object
#ifdef RUN_ASYNC
        pthread_mutex_lock( &surf_mutex );   
#endif
        
        renderable_vertices = g_surf->get_positions();
        
        const std::vector<Vec3st>& mesh_triangles = g_surf->m_mesh.get_triangles();
        
        renderable_solid_triangles.clear();
        renderable_non_solid_triangles.clear();
        for ( size_t i = 0; i < mesh_triangles.size(); ++i )
        {
            if ( g_surf->triangle_is_solid(i) )
            {
                renderable_solid_triangles.push_back( mesh_triangles[i] );
            }
            else
            {
                renderable_non_solid_triangles.push_back( mesh_triangles[i] );
            }
        }
        
        renderable_edges.clear();
        for ( size_t i = 0; i < g_surf->m_mesh.m_edges.size(); ++i )
        {
            renderable_edges.push_back( g_surf->m_mesh.m_edges[i] );
        }
        
        // compute vertex normals
        g_surf->get_all_vertex_normals( renderable_vertex_normals );
        
        // get cached driver visualizer
        driver_renderer = driver->get_renderer();
        
#ifdef RUN_ASYNC   
        // release mutex
        int ok = pthread_mutex_unlock( &surf_mutex );
        assert( ok == 0 );
#endif
        
    }
    
    // ---------------------------------------------------------
    
    // runs on a separate thread:
    void* advance_frame_async( void* )
    {      
        // advance one frame
        pthread_mutex_lock( &surf_mutex );
        advance_frame();   
        int ok = pthread_mutex_unlock( &surf_mutex );   
        assert( ok == 0 );
        
        // and signal we're done
        pthread_mutex_lock( &thread_is_running_mutex );
        thread_is_running = false;
        ok = pthread_mutex_unlock( &thread_is_running_mutex );
        assert( ok == 0 );
        
        return NULL;
    }
    
    
    // ---------------------------------------------------------
    
    void start_advance_frame()
    {
        
        // make sure no frame is currently running
        
        pthread_mutex_lock( &thread_is_running_mutex );
        
        if ( !thread_is_running && !waiting_for_thread_to_finish )
        {
            thread_is_running = true;
            
            int ok = pthread_mutex_unlock( &thread_is_running_mutex );
            assert( ok == 0 );
            
            if ( g_take_screenshots )
            {
                // If first frame, output initial screen cap
                if ( frame_stepper->get_frame() == 0 )
                {
                    // output initial conditions
                    char sgi_filename[256];
                    sprintf( sgi_filename, "%s/frame%04d.sgi", g_output_path, frame_stepper->get_frame() );      
                    Gluvi::sgi_screenshot( sgi_filename );         
                }
            }
            
            waiting_for_thread_to_finish = true;
            
#ifdef RUN_ASYNC
            // kick off advance_frame_async  
            pthread_create( &advance_frame_thread, NULL, advance_frame_async, NULL );
#else
            advance_frame_async( 0 );
#endif
            
            pthread_mutex_lock( &thread_is_running_mutex );
            
        }
        
        int ok = pthread_mutex_unlock( &thread_is_running_mutex );
        assert( ok == 0 );
        
    }
    
    // ---------------------------------------------------------
    
    void advance_frame_done()
    {
        // copy sim data into render buffers
        update_renderable_objects();
        
        pthread_mutex_lock( &surf_mutex );
        
        if ( g_take_screenshots )
        {
            char sgi_filename[256];
            sprintf( sgi_filename, "%s/frame%04d.sgi", g_output_path, frame_stepper->get_frame() );      
            Gluvi::sgi_screenshot( sgi_filename );
        }
        
        // allow the driver to write to disk (e.g. for caching simulation data)
        driver->write_to_disk( g_output_path, frame_stepper->get_frame() );
        
        int ok = pthread_mutex_unlock( &surf_mutex );
        assert( ok == 0 );
        
    }
    
    
    // ---------------------------------------------------------
    
    void timer_advance_frame( int )
    {
        
        if ( sim->m_running || advance_single_frame )
        {
            start_advance_frame();
            advance_single_frame = false;
        }
        
        glutPostRedisplay();
        
    }
    
    // ---------------------------------------------------------
    ///
    /// Handle keyboard input
    ///
    // ---------------------------------------------------------
    
    void keyboard(unsigned char key, int, int )
    {
        
        if(key == 'q')
        {
            pthread_mutex_lock( &surf_mutex );   
            delete g_surf;
            pthread_mutex_unlock( &surf_mutex );   
            delete driver;
            exit(0);
        }
        
        if ( key == 's' )
        {
            display_driver = !display_driver;
        }
        
        //   if ( key == 'S' )
        //   {
        //      g_fluid_render_type = ( g_fluid_render_type + 1 ) % 7;
        //   }
        
        // run one frame
        // 
        if(key == 'n')
        {
            advance_single_frame = true;
        }
        
        if ( key == 'e' )
        {
            mesh_renderer.render_edges = !mesh_renderer.render_edges;
        }
        
        if ( key == 't' )
        {
            mesh_renderer.render_fill_triangles = !mesh_renderer.render_fill_triangles;
        }
        
        if ( key == 'v' )
        {
            mesh_renderer.render_vertex_rank = !mesh_renderer.render_vertex_rank;
        }
        
        
        // return to default camera
        //
        if ( key == 'r' )
        {
            gluvi_cam.return_to_default();
        }
        
        // define default camera
        //
        if ( key == 'd' )
        {
            gluvi_cam.default_target[0] = gluvi_cam.target[0];
            gluvi_cam.default_target[1] = gluvi_cam.target[1];
            gluvi_cam.default_target[2] = gluvi_cam.target[2];      
            gluvi_cam.default_dist = gluvi_cam.dist;
            gluvi_cam.default_heading = gluvi_cam.heading;
            gluvi_cam.default_pitch = gluvi_cam.pitch;   
            
            std::cout << "cam_target " << gluvi_cam.default_target[0] << " " <<  gluvi_cam.default_target[1] << " " << gluvi_cam.default_target[2] << std::endl;
            std::cout << "cam_distance " << gluvi_cam.default_dist << std::endl;
            std::cout << "cam_heading " << gluvi_cam.default_heading << std::endl;
            std::cout << "cam_pitch " << gluvi_cam.default_pitch << std::endl;
        }
        
        // output binary file
        //
        if(key == 'b')
        {
            pthread_mutex_lock( &surf_mutex );   
            char binary_filename[256];
            sprintf( binary_filename, "%s/frame%04d.bin", g_output_path, frame_stepper->get_frame() );      
            write_binary_file( g_surf->m_mesh, g_surf->get_positions(), g_surf->m_masses, sim->m_curr_t, binary_filename );   
            std::cout << "binary file written" << std::endl;      
            pthread_mutex_unlock( &surf_mutex );   
        }
        
        // output OBJ file
        //
        if(key == 'o')
        {
            pthread_mutex_lock( &surf_mutex );   
            char obj_filename[256];
            sprintf( obj_filename, "%s/frame%04d.obj", g_output_path, frame_stepper->get_frame() );
            write_objfile( g_surf->m_mesh, g_surf->get_positions(), obj_filename );
            std::cout << "obj file written" << std::endl;
            pthread_mutex_unlock( &surf_mutex );   
        }
        
        // SGI screenshot
        //
        if(key == 'p')
        {
            pthread_mutex_lock( &surf_mutex );   
            char sgi_filename[256];
            sprintf( sgi_filename, "%s/frame%04d.sgi", g_output_path, frame_stepper->get_frame() );      
            Gluvi::sgi_screenshot( sgi_filename );
            std::cout << "sgi screen shot taken" << std::endl;
            pthread_mutex_unlock( &surf_mutex );   
        }
        
        // toggle simulation running
        //
        if(key == ' ')
        {
            advance_single_frame = false;
            sim->m_running = !sim->m_running;
            std::cout << "running: " << (sim->m_running ? "yes" : "no") << std::endl;
        }
        
        //   if (key == 'h')
        //   {
        //      pthread_mutex_lock( &surf_mutex );   
        //      char v_binary_filename[256];
        //      sprintf( v_binary_filename, "%s/pre-improve.bin", g_output_path );      
        //      read_binary_file( g_surf->m_mesh, , g_surf->m_masses, sim->m_curr_t, v_binary_filename );   
        //      g_surf->m_newpositions.resize( g_surf->get_num_vertices() );      
        //      g_surf->improve_mesh();
        //      pthread_mutex_unlock( &surf_mutex );   
        //   }   
        
        
        if ( key == 'Z' )
        {
            run_all_sisc_examples( "./sisc-scripts/" );
        }
        
        if ( key == 'c' )
        {
            std::cout << "loading meshes" << std::endl;
            
            //g_surf->m_positions.clear();
            
            double dt;
            std::vector<Vec3d> xs;
            read_binary_file( g_surf->m_mesh, xs, g_surf->m_masses, dt, "/Users/tyson/scratch/tbrochu/collisiondebug/current.bin" );
            g_surf->set_all_positions(xs);
            
            NonDestructiveTriMesh temp_mesh;
            std::vector<Vec3d> new_xs;
            read_binary_file( temp_mesh, new_xs, g_surf->m_masses, dt, "/Users/tyson/scratch/tbrochu/collisiondebug/predicted.bin" );
            g_surf->set_all_newpositions(new_xs);
            
            // TEMP: unique-ify
            std::vector< Vec3st > tris;
            for ( unsigned int i = 0; i < temp_mesh.get_triangles().size(); ++i )
            {
                bool found = false;
                for ( unsigned int j = 0; j < tris.size(); ++j )
                {
                    if ( tris[j] == temp_mesh.get_triangle(i) )
                    {
                        found = true;
                        break;
                    }
                }
                
                if ( !found )
                {
                    tris.push_back( temp_mesh.get_triangle(i) );
                }
            }
            
            if ( dt == 0.0 ) { dt = 1.0; }
            
            g_surf->m_mesh.clear();
            g_surf->m_mesh.replace_all_triangles( tris );
            g_surf->defrag_mesh();      
            
            std::cout << "integrating" << std::endl;
            
            double actual_dt;
            g_surf->integrate( dt, actual_dt );      
        }
                
        timer_advance_frame(0);
        
        glutPostRedisplay();
        
    }
    
    // ---------------------------------------------------------
    // Rendering/raycasting
    
    /// Ray cast into the scene, return index of the closest primitive of the type specified
    ///
    unsigned int ray_cast(const DynamicSurface& surface,
                          const Vec3f& ray_origin, 
                          const Vec3f& ray_direction, 
                          unsigned int primitive_type, 
                          size_t& hit_index )
    {
        
        const double RAY_CAST_HIT_DISTANCE = 2e-2;
        
        Vec3d rorigin( ray_origin[0], ray_origin[1], ray_origin[2] );
        Vec3d rhead = rorigin + 100.0 * Vec3d( ray_direction[0], ray_direction[1], ray_direction[2] );
        
        
        Vec3d aabb_low, aabb_high;
        minmax( rorigin, rhead, aabb_low, aabb_high );
        
        if ( primitive_type == 0)
        {
            // -----------------------------------------------------------------
            // Vertices
            // -----------------------------------------------------------------
            
            if ( surface.m_verbose ) { std::cout << "looking for vertices..." << std::endl; }
            
            double vertex_min_parameter = BIG_DOUBLE;
            size_t vertex_min_index = (size_t)~0;
            
            std::vector<size_t> overlapping_vertices;
            surface.m_broad_phase->get_potential_vertex_collisions( aabb_low, aabb_high, true, true, overlapping_vertices );
            
            if ( surface.m_verbose ) { std::cout << "potential vertices: " << overlapping_vertices.size() << std::endl; }
            
            for ( size_t i = 0; i < overlapping_vertices.size(); ++i )
            {
                const Vec3d& v = surface.get_position( overlapping_vertices[i] );
                
                // initialized to silence compiler warnings
                double distance = UNINITIALIZED_DOUBLE;       
                double ray_parameter = UNINITIALIZED_DOUBLE;
                Vec3d normal(UNINITIALIZED_DOUBLE, UNINITIALIZED_DOUBLE, UNINITIALIZED_DOUBLE);
                double normal_multiplier = UNINITIALIZED_DOUBLE;
                
                check_point_edge_proximity( false, v, rhead, rorigin, distance, ray_parameter, normal, normal_multiplier );
                
                if ( distance < RAY_CAST_HIT_DISTANCE )
                {
                    if ( ray_parameter < vertex_min_parameter )
                    {
                        vertex_min_parameter = ray_parameter - RAY_CAST_HIT_DISTANCE;
                        vertex_min_index = overlapping_vertices[i];
                    }
                }
            }
            
            if ( vertex_min_parameter < BIG_DOUBLE )
            {
                hit_index = vertex_min_index;
                return RAY_HIT_VERTEX;
            }
            else
            {
                return RAY_HIT_NOTHING;
            }
            
        }
        else if ( primitive_type == 1)
        {
            
            // -----------------------------------------------------------------
            // Edges
            // -----------------------------------------------------------------
            
            double edge_min_distance = BIG_DOUBLE;
            size_t edge_min_index = (size_t)~0;
            
            std::vector<size_t> overlapping_edges;
            surface.m_broad_phase->get_potential_edge_collisions( aabb_low, aabb_high, true, true, overlapping_edges );
            
            for ( size_t i = 0; i < overlapping_edges.size(); ++i )
            {
                const Vec2st& e = surface.m_mesh.m_edges[ overlapping_edges[i] ];
                const Vec3d& v0 = surface.get_position( e[0] );
                const Vec3d& v1 = surface.get_position( e[1] );
                
                double distance;
                double ray_parameter, edge_parameter;
                Vec3d normal;
                
                check_edge_edge_proximity( rhead,  
                                          rorigin,
                                          v0,
                                          v1,
                                          distance, ray_parameter, edge_parameter, normal );
                
                if ( distance < RAY_CAST_HIT_DISTANCE && edge_parameter > 0.0 && edge_parameter < 1.0 )
                {
                    if ( distance < edge_min_distance )
                    {
                        edge_min_distance = distance; // - RAY_CAST_HIT_DISTANCE;
                        edge_min_index = overlapping_edges[i];
                    }
                }         
                
            }
            
            if ( edge_min_distance < BIG_DOUBLE )
            {
                hit_index = edge_min_index;
                return RAY_HIT_EDGE;
            }
            else
            {
                return RAY_HIT_NOTHING;
            }
            
            
        }
        else
        {
            // -----------------------------------------------------------------
            // Triangles
            // -----------------------------------------------------------------
            
            double triangle_min_parameter = BIG_DOUBLE;
            size_t triangle_min_index = (size_t)~0;
            
            std::vector<size_t> overlapping_triangles;
            surface.m_broad_phase->get_potential_triangle_collisions( aabb_low, aabb_high, true, true, overlapping_triangles );
            
            for ( size_t i = 0; i < overlapping_triangles.size(); ++i )
            {
                const Vec3st& tri = surface.m_mesh.get_triangle( overlapping_triangles[i] );
                
                Vec3st t = sort_triangle( tri );
                assert( t[0] < t[1] && t[0] < t[2] && t[1] < t[2] );
                
                const Vec3d& v0 = surface.get_position( t[0] );
                const Vec3d& v1 = surface.get_position( t[1] );
                const Vec3d& v2 = surface.get_position( t[2] );      
                
                double origin_parameter, ray_parameter, s0, s1, s2;
                Vec3d normal;
                size_t dummy_index = surface.get_num_vertices();
                
                
                bool ray_hit = segment_triangle_intersection( rorigin, dummy_index,
                                                             rhead, dummy_index + 1,                                                      
                                                             v0, t[0],
                                                             v1, t[1],
                                                             v2, t[2],   
                                                             origin_parameter, ray_parameter, s0, s1, s2,
                                                             false );
                
                if ( ray_hit )
                {
                    if ( ray_parameter < triangle_min_parameter )
                    {
                        triangle_min_parameter = ray_parameter;
                        triangle_min_index = overlapping_triangles[i];
                    }
                }         
            }
            
            if ( triangle_min_parameter < BIG_DOUBLE )
            {
                hit_index = triangle_min_index;
                return RAY_HIT_TRIANGLE;
            }
            else
            {
                return RAY_HIT_NOTHING;
            }
            
        }
        
    }


    // ---------------------------------------------------------
    ///
    /// Handle mouse click
    ///
    // ---------------------------------------------------------
    
    void mouse(int button, int state, int x, int y)
    {
        if ( glutGetModifiers() == GLUT_ACTIVE_ALT )
        {
            selected_vertex = selected_edge = selected_triangle = -1;
            
            float ray_origin[3], ray_direction[3];
            
            gluvi_cam.transform_mouse( x, y, ray_origin, ray_direction );
            
            size_t hit_index;
            
            unsigned int ray_cast_type;
            if ( button==GLUT_LEFT_BUTTON ) { ray_cast_type = 0; }
            else if ( button==GLUT_MIDDLE_BUTTON ) { ray_cast_type = 1; }
            else { ray_cast_type = 2; }
            
            pthread_mutex_lock( &surf_mutex );   
            unsigned int ray_hit_type = ray_cast( *g_surf, Vec3f( ray_origin ), Vec3f( ray_direction ), ray_cast_type, hit_index );
            pthread_mutex_unlock( &surf_mutex );   
            
            switch (ray_hit_type)
            {
                case RAY_HIT_VERTEX:
                    selected_vertex = hit_index;
                    break;
                case RAY_HIT_EDGE:
                    selected_edge = hit_index;
                    break;
                case RAY_HIT_TRIANGLE:
                    selected_triangle = hit_index;      
                    break;
                default:
                    assert( ray_hit_type == RAY_HIT_NOTHING );
                    break;
            }
            
            glutPostRedisplay();
            return;
        }
        
        gluvi_cam.click(button, state, x, y);
    }
    
    // ---------------------------------------------------------
    ///
    /// Handle mouse motion
    ///
    // ---------------------------------------------------------
    
    void mouseMotion(int x, int y)
    {
        gluvi_cam.drag(x, y);
    }
    
    
    // ---------------------------------------------------------
    ///
    /// Render to the OpenGL GUI.  Rendering the triangle surface is handled by SurfTrack::render
    ///
    // ---------------------------------------------------------
    
    void display(void)
    {
        
        if ( ! Gluvi::taking_screenshot )
        {
            
            if ( display_status_text )
            {
                pthread_mutex_lock( &thread_is_running_mutex );
                if ( thread_is_running )
                {
                    if ( sim->m_running )
                    {
                        // running async
                        status_text_widget->set_color( 0.0f, 0.0f, 0.0f );
                        status_text_widget->text = "Running --- path:" + std::string( g_output_path );
                    }
                    else
                    {
                        // stopping after the current frame
                        status_text_widget->set_color( 0.8f, 0.2f, 0.2f );         
                        status_text_widget->text = "Stopping...";
                    }
                }
                else
                {
                    if ( waiting_for_thread_to_finish )
                    {
                        
                        // thread has signalled that it's done
                        // wait for the thread to finish and be destroyed
                        pthread_join( advance_frame_thread, NULL );
                        
                        waiting_for_thread_to_finish = false;
                        
                        advance_frame_done();
                        
                        status_text_widget->set_color( 0.0f, 0.0f, 0.0f );
                        status_text_widget->text = "Ready";
                        
                    }
                }
                int ok = pthread_mutex_unlock( &thread_is_running_mutex );
                assert( ok == 0 );
            }
            else
            {
                status_text_widget->set_color( 1.0f, 1.0f, 1.0f );
                status_text_widget->text = "";
            }
            
            // now should immediately check if we should launch another thread
            
            timer_advance_frame(0);
            
        }
        
        
        //
        // Render the mesh
        //
                
        mesh_renderer.render(renderable_vertices, 
                             renderable_vertex_normals,
                             renderable_non_solid_triangles,
                             renderable_edges );
        
        //
        // Render the simulation
        //
        
        if ( display_driver && driver_renderer != NULL )
        {
            driver_renderer->render();
        }
        
    }
    
    
    // ---------------------------------------------------------
    ///
    /// Initialize the OpenGL GUI
    ///
    // ---------------------------------------------------------
    
    void init_gui( const ScriptInit& script_init )
    {
        
        glClearColor(1,1,1,0);
        Gluvi::camera = &gluvi_cam;
        
        gluvi_cam.target[0] = static_cast<float>(script_init.camera_target[0]);
        gluvi_cam.target[1] = static_cast<float>(script_init.camera_target[1]);
        gluvi_cam.target[2] = static_cast<float>(script_init.camera_target[2]);   
        gluvi_cam.dist = static_cast<float>(script_init.camera_distance);
        gluvi_cam.heading = static_cast<float>(script_init.camera_heading);
        gluvi_cam.pitch = static_cast<float>(script_init.camera_pitch);
        
        // set to default
        gluvi_cam.default_target[0] = gluvi_cam.target[0];
        gluvi_cam.default_target[1] = gluvi_cam.target[1];
        gluvi_cam.default_target[2] = gluvi_cam.target[2];      
        gluvi_cam.default_dist = gluvi_cam.dist;
        gluvi_cam.default_heading = gluvi_cam.heading;
        gluvi_cam.default_pitch = gluvi_cam.pitch;   
        
        
        Gluvi::userDisplayFunc=display;
        
        glutKeyboardFunc(keyboard);
        glutMouseFunc(mouse);
        glutMotionFunc(mouseMotion);
        
        status_text_widget = new Gluvi::DynamicText( "Ready" );
        Gluvi::root.list.push_back( status_text_widget );
        
        gluvi_cam.gl_transform();
        
        update_renderable_objects();
    }
    
#endif  // #ifdef USE_GUI
    
    
    // ---------------------------------------------------------
    
    void set_up_output_path( )
    {
        
        DIR *dp;
        struct dirent *ep;     
        dp = opendir ( g_output_path );
        
        int max_dir = 0;
        
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
        
        sprintf( g_output_path, "%s/run%04d", g_output_path, max_dir + 1 );
        std::cout << "Output path: " << g_output_path << std::endl;
        
        char mkdir_command[1024];
        sprintf( mkdir_command, "mkdir -p %s", g_output_path );
        system(mkdir_command);
                
    }
    
    
    // ---------------------------------------------------------
    ///
    /// Initialize the simulation.  Set the simulation parameters, initial geometry, etc.
    ///
    // ---------------------------------------------------------
    
    void init_simulation( char *script_filename )
    {
        
        SurfTrackInitializationParameters init_params;
        
        std::vector<Vec3ui> tris;
        std::vector<Vec3d> verts;
        std::vector<Vec3d> velocities;
        std::vector<double> masses;

        // NOTE: changing also for concentrations
        std::vector<std::vector<double> > c_species;  
  
        // std::vector<double> c_C;
        // std::vector<double> c_BC;
        // std::vector<double> c_N;
        // std::vector<double> c_BN;
        // std::vector<double> c_S;  
        // std::vector<double> c_Tld;               
        
        ScriptInit script_init;
        script_init.parse_script( script_filename );

        // set the global variable of diffusion param to whatever was read in the 'parsing script function'
        diffusionParam = script_init.diffusionParam;
        std::cout<<"read the following diffusion parameters for this simulation\n";
        for(int i=0;i<diffusionParam.size()-1;i++){
            std::cout<<diffusionParam[i]<<",    ";
        }
        std::cout<<diffusionParam[diffusionParam.size()-1]<<"\n";

        if ( script_init.output_path_is_relative )
        {
            snprintf( g_output_path, 256, "%s/%s", g_base_output_path, script_init.output_path.c_str() );
        }
        else
        {
            snprintf( g_output_path, 256, "%s", script_init.output_path.c_str() );
        }
        
        // Init frame stepper
        
        frame_stepper = new FrameStepper( script_init.frame_dt );
        
        // init simulation
        
        if ( script_init.end_sim_t != UNINITIALIZED_DOUBLE )
        {
            sim = new Simulation( script_init.sim_dt, script_init.end_sim_t );
        }
        else
        {
            sim = new Simulation( script_init.sim_dt );
        }
        
        if ( script_init.curr_t_specified )
        {
            sim->m_curr_t = script_init.curr_t;
            
            unsigned int curr_frame = static_cast<unsigned int>(script_init.curr_t / script_init.frame_dt); 
            std::cout << "curr_t: " << sim->m_curr_t << std::endl;
            std::cout << "curr_frame: " << curr_frame << std::endl; 
            
            frame_stepper->frame_count = curr_frame;
        }
        
        
        // init SurfTrack
        
        // g_surf = new SurfTrack( script_init.vertices, script_init.triangles, script_init.masses, script_init.surf_track_params );   

        // NOTE: changing also for concentrations
        //g_surf = new SurfTrack( script_init.vertices, script_init.triangles, script_init.masses,script_init.c_B,script_init.c_C,script_init.c_BC,script_init.c_N,script_init.c_BN,script_init.c_S,script_init.c_Tld, script_init.surf_track_params );   
        g_surf = new SurfTrack( script_init.vertices, script_init.triangles, script_init.masses,script_init.c_species, script_init.surf_track_params );   
        
        g_surf->defrag_mesh();
        
        // init driver
        
        driver = script_init.driver;   
        driver->initialize( *g_surf );
        
        if ( g_surf->m_collision_safety )
        {
            g_surf->m_collision_pipeline.assert_mesh_is_intersection_free( false );      
        }
        
        //g_stats.set_double( "max_edge_length", g_surf->m_splitter.m_max_edge_length );
        //g_stats.set_double( "min_edge_length", g_surf->m_collapser.m_min_edge_length );
        //g_stats.set_double( "max_volume_change", g_surf->m_max_volume_change );
        
        set_time_base();
        
#ifdef USE_GUI
        
        init_gui( script_init );
        
#endif
        
        
    }
    
}  // unnamed namespace

// ---------------------------------------------------------
///
/// MAIN
///
// ---------------------------------------------------------

int main(int argc, char **argv)
{   
    clock_t tStart = clock();
    std::cout << "\nzebrafish: use of El Topo to simulate invagination" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl << std::endl;
    
    if( argc < 2 )
    {
        std::cout << "Usage: <executable> <scriptfile> <outputbasedirectory>" << std::endl;
        std::cout << "e.g.: " << std::endl;
        std::cout << "$ ./zebrafish_release simulation_script.txt ./ \n" << std::endl;
        return 0;
    }
    
    // set path for outputting obj, bin files, etc.
    
    if ( argc > 2 )
    {
        strncpy( g_base_output_path, argv[2], sizeof(g_base_output_path) );
    }
    else
    {
        // argc == 2
        strncpy( g_base_output_path, "./", sizeof(g_base_output_path) );
    }
    
#ifdef USE_GUI
    pthread_mutexattr_t mta;
    int rc = pthread_mutexattr_init(&mta);
    rc = pthread_mutexattr_settype(&mta, PTHREAD_MUTEX_ERRORCHECK );
    rc = pthread_mutex_init(&thread_is_running_mutex, &mta);
    
    Gluvi::winwidth = 800;
    Gluvi::winheight = 600;
    
    Gluvi::init( "Talpa", &argc, argv );
#endif
    
    //
    // Initialize the simulation using the script file
    //
    
    init_simulation( argv[1] );
    
    
    //
    // Make a new directory for output
    //
    
    set_up_output_path();
    
    //
    // Make a copy of the input script in the output directory
    //
    
    char* script_filename = argv[1];
    std::ifstream original_file_stream( script_filename );
    assert( original_file_stream.good() );
    
    char script_copy_filename[256];
    snprintf( script_copy_filename, 256, "%s/aaa-script.txt", g_output_path );
    
    char command[1024];
    sprintf( command, "cp %s %s", script_filename, script_copy_filename );
    system(command);
    
    srand( 1 );

    int node_n = g_surf->get_num_vertices(); 
    std::vector<std::vector<double> > c_output;
    c_output = g_surf->get_cspecies();
   

    std::vector<double> c_B(node_n);
    std::vector<double> c_C(node_n);
    std::vector<double> c_BC(node_n);
    std::vector<double> c_N(node_n);
    std::vector<double> c_BN(node_n);
    std::vector<double> c_S(node_n);
    std::vector<double> c_Tld(node_n);

    for (int i=0; i<node_n;++i){
        c_B[i] = c_output[i][0];
        c_C[i] = c_output[i][1];
        c_BC[i] = c_output[i][2];
        c_N[i] = c_output[i][3];
        c_BN[i] = c_output[i][4];
        c_S[i] = c_output[i][5];
        c_Tld[i] = c_output[i][6];
    }

    
    // write frame 0
    char binary_filename[256];
    sprintf( binary_filename, "%s/frame%04d.bin", g_output_path, frame_stepper->get_frame() );      
    //write_binary_file( g_surf->m_mesh, g_surf->get_positions(), g_surf->m_masses, sim->m_curr_t, binary_filename );   
    
    char vtk_filename[256];
    sprintf( vtk_filename, "%s/Aframe%04d.vtk", g_output_path, frame_stepper->get_frame() );      
    write_vtk( g_surf->m_mesh, g_surf->get_positions(),c_B,c_C,c_BC,c_N,sim->m_curr_t, vtk_filename );

    char vtk_filename1[256];
    sprintf( vtk_filename1, "%s/Bframe%04d.vtk", g_output_path, frame_stepper->get_frame() );      
    write_vtk( g_surf->m_mesh, g_surf->get_positions(),c_B,c_BN,c_S,c_Tld,sim->m_curr_t, vtk_filename1 );


    driver->write_to_disk( g_output_path, frame_stepper->get_frame() );
    
    //
    // Now start
    //
    
#ifdef USE_GUI
    // start the GUI
    Gluvi::run();
#else
    
    // start the simulation (hands off)
    run_simulation();
    
    // uncomment to run all examples in our SISC paper:
    //run_all_sisc_examples( "./sisc-scripts/" );

    // get the last frame of the simulaiton result _ add by linlin Li063020

    std::vector<std::vector<double> > c_output_last = g_surf-> get_cspecies();

    node_n = c_output_last.size();
    //std::cout<<"coutput"<< c_output.size()<<"\n";


    std::vector<double> c_B_last(node_n);
    std::vector<double> c_C_last(node_n);
    std::vector<double> c_BC_last(node_n);
    std::vector<double> c_N_last(node_n);
    std::vector<double> c_BN_last(node_n);
    std::vector<double> c_S_last(node_n);
    std::vector<double> c_Tld_last(node_n);

    for (int i=0; i<node_n;++i){
        c_B_last[i] = c_output_last[i][0];
        c_C_last[i] = c_output_last[i][1];
        c_BC_last[i] = c_output_last[i][2];
        c_N_last[i] = c_output_last[i][3];
        c_BN_last[i] = c_output_last[i][4];
        c_S_last[i] = c_output_last[i][5];
        c_Tld_last[i] = c_output_last[i][6];
    }
    auto maxPosition_b = max_element(c_B_last.begin(), c_B_last.end());
    double max_bmp = *maxPosition_b;
    auto maxPosition_s = max_element(c_S_last.begin(), c_S_last.end());
    double max_szl = *maxPosition_s;

        // find max value of BMP and Szl for the very last frame ( add by Linlin Li, Jun 27,2020)
    std::cout << " max_bmp \n" << max_bmp <<std::endl;
    std::cout << " max_szl \n" << max_szl <<std::endl;


    char data_filename[256];
    sprintf(data_filename, "%s/max_bmpszl.txt", g_output_path ); 
    std::ofstream outfile(data_filename,std::ios::trunc);

    outfile<<"max_bmp "<<max_bmp<<"\n";
    outfile<<"max_szl "<<max_szl<<"\n";

    outfile.close ();

    
#endif

// print out the execution time    
printf("execution time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
}



