//********************************************************
//********************************************************
//****************VERTX2D*********************************
//********************************************************
//********************************************************
//2D vertex model code, originally written by Matej Krajnc et al. (please thank him in Ackowledgemnts)


//****************GLOBAL PARAMETERS************************
//this part sets global variables - best not to change for now
//AREA-ENERGY-RELATED PARAMETERS
double c_kAb=1.; //area modulus
double A0=1., A01=1, A02=0.5; //A0, A01 are the target area

//PERIMETER-ENERGY-RELATED PARAMETERS
double kPer = 0.01;
double P0 = 3;

//nematic activities
double zeta_inter = 0; //an old activity, doesn't matter by default (leave = 0)
double zeta_1 = 0; //nematic activity
double zeta_2 = 0;

//chiral activities
double chiral_alpha_1 = 0;
double chiral_alpha_2 = 0;

//FLUCTUATIONS 
double gamma0=0.; //mean tension, leave =0
double tauM = 0; //timescale
double Sig = 0; //amplitude
int step_count = 0; //irrelevant

//T1
double lth=0.01; //threshold for a T1 transition

//MISC
double tMAX=1000; //simulation time
int transition_type = 0; //transition type for disk initial condition (0: sharp, 1: smooth)
double box_ratio = 1.; //leave 1
int Nx = 12; //size of simulation box in cells along x

//None of these matter, but the code wants them here:
int N_grid = 22;//60;//4;
int N_grid_max = 1000;//60;//4;
int T1_angle_count[362];
int N_defects = 0;

int T1_count_pop_11 = 0;
int T1_count_pop_12 = 0;
int T1_count_pop_22 = 0;


#include<random>
#include <Eigen/Sparse>
std::mt19937 gen_sigma; //RNG generator, doesn't get used
std::normal_distribution<double> dis_sigma(0,1);

//**************************************************
//*******************MATRIX*************************
//**************************************************
//for wet dynamics
using Tripl   = Eigen::Triplet<double>;
using SpMat   = Eigen::SparseMatrix<double>;
using SolType = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;

SpMat _M;   // Connectivity matrix
SolType _solver; // Sparse solver

double _gamma_wet = 1.; //external friction (should not be 0)
double _zeta_wet = 0.; //internal friction
    
//**************************************************
//*******************DECLARATIONS*******************
//**************************************************
#include "fstream"
#include "vertX2D//_functions.h"
#include "vertX2D//_output.h"
#include <filesystem>
#include <string>
#include <iostream>

//********************************************************
//********************************************************
//********************MAIN********************************
//********************************************************
//********************************************************
//start of the actual code - everything above here is best left unchaged
int main(int argc, char *argv[]){

    //********************************************************
    //********************SIMULATION PARAMETERS***************
    //********************************************************
    //this sets the simulation parameters - reads from the command line when you run the code

    //********************READ PARAMETERS********************
    //this parameters are read for each simulation when you run it from the command line, e.g.: 
    // ./main.opp 32 0.5 0.01 3.95 ...
    //would mean Nx = 32; c_kAb = 0.5; kPer = 0.01; P0 = 3.95; and so on
    if(argc != 24) { std::cerr << "Error! Wrong number of input elements!" << std::endl; return -1; }
    
    //general vertex model parameters
    Nx = std::atoi(argv[1]); //size of initial system along X (good value: 12 for simple tests, 32 for serious simulations)
    c_kAb = std::stod(argv[2]);;//target area modulus (good value: 0.5)
    kPer = std::stod(argv[3]);;//target perimeter modulus (good value: 0.01)
    P0 = std::stod(argv[4]);//target periemter (good values between 3.6 and 3.95)
    double disorder = std::stod(argv[5]); //this will determine how disordered the initial tissue is - 0.03 for very disordered, 0.02 for moderately disordered, 0.001 for only hexagons
    tMAX = std::stod(argv[6]); //how long the simulation runs (good value: 1000)
    
    //write to file parameters
    double delta_write_time = std::stod(argv[7]); //time interval between writing to file (good value: 10)
    double delta_image_time = std::stod(argv[8]); //time interval between outputting tissue (good value: 100 or 1000, depending on tMAX)
    
    //chiral activity parameters
    chiral_alpha_1 = std::stod(argv[9]); //chiral activity of the first population (inside of droplet if there is a droplet)
    chiral_alpha_2 = std::stod(argv[10]); //chiral activity of the second population
    
    //initial condition activity
    int initial_condition = std::atoi(argv[11]); //initial condition type (0: drop, 1: random mixture, 2: split along x, 3: split along y)
    double initial_parameter = std::stod(argv[12]); //parameter for generating the initial condition (for the drop, radius, for the random mixture, fraction of second population); the radius is given as a fraction of 1/2 the y dimension of the box, so r = 1 means the diameter of the patch matches the size of the box
    int initial_seed = std::atoi(argv[13]); //seed for chosing the random parts of the initial condition 
    int tiling_seed = std::atoi(argv[14]); //seed for chosing the random parts of the initial condition 

    int dynamics_type = std::atoi(argv[15]); //this chooses if the model is edges-wet (dynamics_type == 2), vertex-vertex wet (dynamics_type == 1) or dry (dynamics_type == 0)
    _gamma_wet = std::stod(argv[16]); //external friction (should not be 0)
    _zeta_wet = std::stod(argv[17]); //internal friction (if you are doing wet dynamics, I would set this to one and change _gamma_wet; in general, things seem to be interesting if _gamma_wet is at least an order of magnitude less than _zeta_wet)
    transition_type = std::atoi(argv[18]); //transition type for disk initial condition (0: sharp, 1: smooth)

    double streak_length   = std::stod(argv[19]); //length of streak
    double streak_width    = std::stod(argv[20]); //width of streak
    double streak_center_y = std::stod(argv[21]); //y coordinate of streak center
    double node_radius     = std::stod(argv[22]); //radius of node
    double node_center_y   = std::stod(argv[23]); //y coordinate of node center

    

    //********************FIXED PARAMETERS********************
    //These are some other parameteirs - they don't really matter for now, so they are just fixed
    zeta_1 = 0;
    zeta_2 = 0;
    tauM =  1;
    Sig = 0;
    gamma0 =  0;  
    double pert_scale = 0.1;
    int pert_seed = 1;
    array_max = 300000; //size of arrays, this is a good value
    seed  = 1; //random seed, doesn't get used
    srand(seed);

    //********************************************************
    //********************FILES*******************************
    //********************************************************
    
    //This creates the folder into which the output files will go
    char filename_description[300]; snprintf(filename_description, sizeof(char) * 300, "./output/out_Nx_%d_kA_%g_kP_%g_P0_%g_disorder_%g_tMAX_%g_dw_%g_di_%g_ca1_%g_ca2_%g_iC_%d_iP_%g_iS_%d_iT_%d_dynamics_%d_gammaW_%g_zetaW_%g_tT_%d", Nx, c_kAb, kPer, P0, disorder, tMAX, delta_write_time, delta_image_time, chiral_alpha_1, chiral_alpha_2, initial_condition, initial_parameter, initial_seed, tiling_seed, dynamics_type, _gamma_wet, _zeta_wet, transition_type);
    filesystem::create_directory(filename_description);
    cout << filename_description << endl; //writes the output folder to command line

    //Creates the various output files
    char filename_energies[500]; snprintf(filename_energies, sizeof(char) * 500, "%s/energies.dat", filename_description);
    char filename_sorting[500]; snprintf(filename_sorting,   sizeof(char) * 500, "%s/sorting.dat", filename_description);
    char filename_msd[500];      snprintf(filename_msd,      sizeof(char) * 500, "%s/msd.dat",      filename_description);
    char filename_moves[500];    snprintf(filename_moves,    sizeof(char) * 500, "%s/max_move.dat", filename_description);
    char filename_chiral_area[500]; snprintf(filename_chiral_area,  sizeof(char) *500, "%s/chiral_area.dat", filename_description);
    char filename_T1_transitions[500]; snprintf(filename_T1_transitions, sizeof(char) * 500, "%s/T1_transitions.dat", filename_description);
    
    
    FILE *file_energies; file_energies = fopen(filename_energies, "wt");
    FILE *file_sorting;  file_sorting  = fopen(filename_sorting, "wt");
    FILE *file_msd;      file_msd      = fopen(filename_msd, "wt");
    FILE *file_moves;    file_moves    = fopen(filename_moves, "wt");
    FILE *file_chiral_area; file_chiral_area    =   fopen(filename_chiral_area, "wt");
    FILE *file_T1_transitions; file_T1_transitions = fopen(filename_T1_transitions, "wt");

   

    //********************************************************
    //********************INITIALIZE**************************
    //********************************************************
    //sets up the model tissue & does some background stuff

    //*******************SETUP FUNCTIONS***********************
    //functions that need to run at the start of the simulation
    allocate(); //background stuff
    set_initial_regHex(Nx); //makes Nx by Nx hexagonal lattice
    
    //*******************PREPARE*******************************
    //Sets some additional variables that are necessary
    Time=0; //sets time to 0
    h=0.01; //time step dt
    double write_time = 0.9* tMAX; //time of next output (=delta_write_time, originally)
    double image_time = 0.9* tMAX; //time of next image output (=delta_image_time, originally)

    //*******************DISORDER******************************
    //creates a disordered tissue
    make_disordered(disorder, tiling_seed); //makes the initial tissue disordered

    if (initial_condition == 0) disk_initial(initial_parameter, initial_seed, Nx, transition_type);
    if (initial_condition == 1) mixed_initial(initial_parameter, initial_seed, Nx);
    if (initial_condition == 2) split_initial_x(initial_parameter, initial_seed, Nx);
    if (initial_condition == 3) split_initial_y(initial_parameter, initial_seed, Nx);
    if (initial_condition == 4) chick_initial(streak_length, streak_width, streak_center_y, node_radius, node_center_y);

    outWholeTissue(filename_description, 0); //outputs the tissue initial condition to file
    if (dynamics_type == 1) update_M(); //this creates the viscosity matrix
    //********************************************************
    //********************RUN*********************************
    //********************************************************
    //This is the main part of the code - it will solve the equation of motion until t = tMAX

    while(Time < tMAX - h /2.){  


        //MOTION
        eqOfMotion(filename_description, image_time, dynamics_type);//calculates forces and moves vertices, then increments time by h

        //WRITE //this is the part that writes to files, if simulation time is greater than write time, it outputs to files, then increments write time by delta_write_time (same for image time)
        if (Time > write_time - h/2.) write_time = print_output(write_time, delta_write_time, file_energies, file_sorting, file_msd, file_chiral_area, file_T1_transitions, false);
        if (Time > image_time - h/2.) image_time = image_output(image_time, delta_image_time, filename_description);
        
        //TOPOLOGY
        update_topology(0.011, dynamics_type); //T1 transitions
    }


    //********************************************************
    //********************COMPLETE****************************
    //********************************************************
    //completes the simulation, closing the files & deallocating all arrays
    fprintf(file_moves, "%g\t%g\n", Time,        max_move);
    fclose(file_energies); fclose(file_sorting); fclose(file_msd); fclose(file_moves); fclose(file_chiral_area); fclose(file_T1_transitions);
    deallocate(); 
    
    return 0;
}
//********************************************************
