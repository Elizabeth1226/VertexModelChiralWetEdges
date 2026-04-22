//DECLARATIONS
#include <stdio.h>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include "_freeid.h"
#include <time.h>
#include <algorithm>
#define PI 3.141592654

using namespace std;

//GEOMETRIC ELEMENTS
size_t Nv, Nv_pass, Ne, Nf, Nc;
double **v, **v_pass;
int **e, **f;
int **basal_edges;
int **basal_vertices, **basal_facets;

//ATTRIBUTES
//vertices
double **v_F    ;
double **v_vel    ;
int **v_cells   ;
int **v_cell_neighbors   ;
int **v_edges   ;
int *v_T1dir    ;
int *v_T1dir2    ;
int *v_vertT1   ;
int *v_edgeT1   ;
double *v_clock  ;
int *v_type  ;
double *v_charge  ;
double *defect_charge  ;
double *defect_x  ;
double *defect_y  ;

//passive vertices
int *v_pass_cell;

//edges
int *e_cell1    ;
int *e_cell2    ;
double *e_length;
double *e_dldt;
double *e_g;
double *e_target;
double **e_l;
int *e_bin_categories;


//facets
int *f_cell;
    
//cells
int *c_cent;
double *c_A0;
double *c_dA;
int *c_type;
int *c_fixed;
double *c_activity;
double **c_initPos;
double **c_prevPos;
double **c_DxDy;

int *c_hexagons;
int *c_bin_categories;


double *c_MSD_initial_x;
double *c_MSD_initial_y;
double *c_MSD_previous_x;
double *c_MSD_previous_y;

double *c_hexatic_initial_re;
double *c_hexatic_initial_im;

double *d_pair_sum_re;
double *d_pair_sum_im;
double *d_pair_sum_count;

//grid
double **g_Q_xx;
double **g_Q_xy;
double **g_beta;
double **g_count;
double **g_winding;
double **g_v_x;
double **g_v_y;
double **g_v_count;
double **g_direction;

//MISC
double *perioXYZ;
double *chiral_cell_center;
double *v_azimuthal_ave;
double *vr_azimuthal_ave;
double *vphi_azimuthal_ave;
double *Q_type1_azimuthal_ave;
double *Q_type2_azimuthal_ave;
double *Qphi_type1_azimuthal_ave;
double *Qphi_type2_azimuthal_ave;
double *Kp_type1_azimuthal_ave;
double *Kp_type2_azimuthal_ave;
double *bin_positions;
double *fine_bin_positions;
int *azimuthal_bin_T1_count_pop_11;
int *azimuthal_bin_T1_count_pop_12;
int *azimuthal_bin_T1_count_pop_22;
double *azimuthal_edge_length;
int *azimuthal_edge_count;
double h, Time, max_move;
double wA, wP, wl, A_tot, A0tot;
//int nrOfEdges;
//int v_free_id, v_pass_free_id, e_free_id, f_free_id, c_free_id;
FreeId v_freeId, v_pass_freeId, e_freeId, f_freeId, c_freeId;
size_t array_max, seed;
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
#include "_allocate.h"
#include "_rndom.h"
#include "_torus.h"
#include "_distances.h"
#include "_vert_edg_fac.h"
#include "_basal_network.h"
#include "_basal_side.h"
#include "_force_area.h"
#include "_force_length.h"
#include "_force_nematic.h"
#include "_force_chiral.h"
#include "_force_perimeter.h"
//#include "_output.h"
#include "_equation_of_motion.h"
#include "_initial_structure.h"
#include "_dissolve.h"
#include "_list_manipulation.h"
#include "_T1_transformation.h"
#include "_cell_extrusion.h"
#include "_MSD.h"
#include "contractionResolution.h"
#include "_disordered.h"
#include "_azimuthal_average.h"

//****************************************************************************
//****************************************************************************
