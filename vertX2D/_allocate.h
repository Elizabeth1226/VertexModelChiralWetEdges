//****************************************************************************
//****************************************************************************
//******************************ALLOCATE**************************************
//****************************************************************************
//****************************************************************************
void reset_arrays(){
    
    for(int i = 1; i<=array_max; i++){
        
        //***********************
        //GEOMETRIC ELEMENTS*****
        //***********************
        for(int j = 0; j < 3; j++){
            v[i][j]=0;
            v_pass[i][j]=0;
            f[i][j]=0;
            c_initPos[i][j]=0;
            c_prevPos[i][j]=0;
            c_DxDy[i][j]=0;
        }
        for(int j = 0; j < 3; j++){
            e[i][j]=0;
            e_l[i][j]=0;
        }
        //***********************
        //CELL ARRAYS************
        //***********************
        //basal_edges
        basal_edges[i][1]=0;
        basal_edges[i][2]=0;
        //basal_vertices
        basal_vertices[i][1]=0;
        basal_vertices[i][2]=0;
        //basal_facets
        basal_facets[i][1]=0;
        basal_facets[i][2]=0;
        //v_cells
        v_cells[i][1]=0;
        v_cells[i][2]=0;
        v_edges[i][1]=0;
        v_edges[i][2]=0;

        v_cell_neighbors[i][0]=0;
        v_cell_neighbors[i][1]=0;
        v_cell_neighbors[i][2]=0;
        v_cell_neighbors[i][3]=0;
        v_cell_neighbors[i][4]=0;
        v_cell_neighbors[i][5]=0;
        v_cell_neighbors[i][6]=0;
        v_cell_neighbors[i][7]=0;
        v_cell_neighbors[i][8]=0;
        //***********************
        
        
        //***********************
        //ATTRIBUTES*************
        //***********************
        for(int j = 1; j < 3; j++){
            v_F[i][j]=0;
            v_vel[i][j]=0;
        }
        v_T1dir[i]=0;
        v_T1dir2[i]=0;
        v_vertT1[i]=0;
        v_edgeT1[i]=0;
        v_clock[i]=0;
        v_type[i]=0;
        v_charge[i]=0;
        defect_charge[i]=0;
        defect_x[i]=0;
        defect_y[i]=0;
        //***********************
        //passive vertices
        v_pass_cell[i]=0;
        //***********************
        //edges
        e_cell1[i]=0;
        e_cell2[i]=0;
        e_length[i]=0;
        e_dldt[i]=0;
        e_g[i]=0;
        e_target[i]=0;
        e_bin_categories[i]=-1;
        //***********************
        //facets
        f_cell[i]=0;
        //***********************
        //cells
        c_cent[i]=0;
        c_A0[i]=0;
        c_dA[i]=0;
        c_type[i]=0;
        c_fixed[i]=0;
        c_activity[i]=0;
        

        c_hexagons[i]=0;
        c_bin_categories[i]=-1;


        c_MSD_initial_x[i] = 0.;
        c_MSD_initial_y[i] = 0.;
        c_MSD_previous_x[i] = 0.;
        c_MSD_previous_y[i] = 0.;

        c_hexatic_initial_re[i] = 0.;
        c_hexatic_initial_im[i] = 0.;


        //***********************
    }
    for (int i = 0; i < 80; i ++ )
    {
        d_pair_sum_re[i] = 0.;
        d_pair_sum_im[i] = 0.;
        d_pair_sum_count[i] = 0.;
    }
    for (int i = 0; i < N_grid_max; i ++ )
    {
        for (int j = 0; j < N_grid_max; j ++ )
        {
            g_Q_xx[i][j] = 0;
            g_Q_xy[i][j] = 0;
            g_beta[i][j] = 0;
            g_count[i][j] = 0;
            g_winding[i][j] = 0;
            g_v_x[i][j] = 0;
            g_v_y[i][j] = 0;
            g_v_count[i][j] = 0;
            g_direction[i][j] = 0;
        }
        azimuthal_bin_T1_count_pop_11 [i] = 0;
        azimuthal_bin_T1_count_pop_12 [i] = 0;
        azimuthal_bin_T1_count_pop_22 [i] = 0;
        azimuthal_edge_length [i] = 0;
        azimuthal_edge_count [i] = 0;
    }
}
//****************************************************************************
void allocate(){// allocate space
    
    int arrayMax=array_max+1;
    //std::clog << "Max Array Size: " << arrayMax << std::endl;
    
    //GEOMETRIC ELEMENTS
    Nv      = 0;
    Nv_pass = 0;
    Ne      = 0; 
    Nf      = 0; 
    Nc      = 0;
    v        = new double*[arrayMax];
    v_pass   = new double*[arrayMax];
    e = new int*[arrayMax];
    f = new int*[arrayMax];
    c_initPos      = new double*[arrayMax];
    c_prevPos      = new double*[arrayMax];
    c_DxDy         = new double*[arrayMax];
    basal_edges    = new int*[arrayMax];
    basal_vertices = new int*[arrayMax];
    basal_facets   = new int*[arrayMax];
    v_cells   = new int*[arrayMax];
    v_cell_neighbors  = new int*[arrayMax];
    v_edges   = new int*[arrayMax];
    e_l = new double*[arrayMax];
    for(int i=0; i < arrayMax; i++){
        v[i]        = new double[3];
        v_pass[i]   = new double[3];
        e[i]   = new int[3];
        f[i]   = new int[4];
        c_initPos[i]        = new double[3];
        c_prevPos[i]        = new double[3];
        c_DxDy[i]        = new double[3];
        basal_edges[i]    = new int[25];
        basal_vertices[i] = new int[25];
        basal_facets[i]   = new int[25];
        v_cells[i]   = new int[25];
        v_cell_neighbors[i]   = new int[25];
        v_edges[i]   = new int[25];
        e_l[i]   = new double[3];
    }
    
    //ATTRIBUTES
    //vertices
    v_F      = new double*[arrayMax];
    v_vel    = new double*[arrayMax];
    for(int i=0; i < arrayMax; i++){
        v_F[i]      = new double[3];
        v_vel[i]      = new double[3];
    }
    v_T1dir    = new int[arrayMax];
    v_T1dir2    = new int[arrayMax];
    v_vertT1   = new int[arrayMax];
    v_edgeT1   = new int[arrayMax];
    v_clock   = new double[arrayMax];
    v_type    = new int[arrayMax];
    v_charge   = new double[arrayMax];
    defect_charge   = new double[arrayMax];
    defect_x   = new double[arrayMax];
    defect_y   = new double[arrayMax];
    
    //passive vertices
    v_pass_cell = new int[arrayMax];

    //edges
    e_cell1    = new int[arrayMax];
    e_cell2    = new int[arrayMax];
    e_length = new double[arrayMax];
    e_dldt = new double[arrayMax];
    e_g = new double[arrayMax];
    e_target = new double[arrayMax];
    e_bin_categories = new int [arrayMax];

    //facets
    f_cell = new int[arrayMax];

    //cells
    c_cent = new int[arrayMax];
    c_A0 = new double[arrayMax];
    c_dA = new double[arrayMax];
    c_type = new int[arrayMax];
    c_fixed = new int[arrayMax];
    c_activity = new double[arrayMax];

    c_hexagons = new int[arrayMax];
    c_bin_categories = new int[arrayMax];
 

    c_MSD_initial_x  = new double[arrayMax];
    c_MSD_initial_y  = new double[arrayMax];
    c_MSD_previous_x = new double[arrayMax];
    c_MSD_previous_y = new double[arrayMax];

    c_hexatic_initial_re  = new double[arrayMax];
    c_hexatic_initial_im  = new double[arrayMax];

    //grid
    g_Q_xx      = new double*[N_grid_max];
    g_Q_xy      = new double*[N_grid_max];
    g_beta      = new double*[N_grid_max];
    g_count     = new double*[N_grid_max];
    g_winding   = new double*[N_grid_max];
    g_v_x       = new double*[N_grid_max];
    g_v_y       = new double*[N_grid_max];
    g_v_count   = new double*[N_grid_max];
    g_direction = new double*[N_grid_max];

    for(int i=0; i < N_grid_max; i++){
        g_Q_xx[i]      = new double[N_grid_max];
        g_Q_xy[i]      = new double[N_grid_max];
        g_beta[i]      = new double[N_grid_max];
        g_count[i]     = new double[N_grid_max];
        g_winding[i]   = new double[N_grid_max];
        g_v_x[i]       = new double[N_grid_max];
        g_v_y[i]       = new double[N_grid_max];
        g_v_count[i]   = new double[N_grid_max];
        g_direction[i] = new double[N_grid_max];
    }
    
    //MISC
    perioXYZ = new double[2];
    bin_positions           = new double[N_grid_max];
    fine_bin_positions      = new double[10*N_grid_max];
    chiral_cell_center      = new double[2];
    v_azimuthal_ave         = new double[N_grid_max];
    vr_azimuthal_ave        = new double[N_grid_max];
    vphi_azimuthal_ave      = new double[N_grid_max];
    Q_type1_azimuthal_ave    = new double[N_grid_max];
    Qphi_type1_azimuthal_ave = new double[N_grid_max];
    Q_type2_azimuthal_ave    = new double[N_grid_max];
    Qphi_type2_azimuthal_ave = new double[N_grid_max];
    Kp_type1_azimuthal_ave   = new double[N_grid_max];
    Kp_type2_azimuthal_ave   = new double[N_grid_max];
    azimuthal_bin_T1_count_pop_11 = new int[10*N_grid_max];
    azimuthal_bin_T1_count_pop_12 = new int[10*N_grid_max];
    azimuthal_bin_T1_count_pop_22 = new int[10*N_grid_max];
    azimuthal_edge_length = new double[10*N_grid_max];
    azimuthal_edge_count = new int[10*N_grid_max];
    h              = 0;
    Time           = 0;
    max_move       = 0;
    /*v_free_id      = 1;
    v_pass_free_id = 1;
    e_free_id      = 1;
    f_free_id      = 1;
    c_free_id      = 1;*/
    wA             = 0;
    wP             = 0;
    wl             = 0;
    A_tot          = 0;
    A0tot          = 0;

    d_pair_sum_re = new double[80];
    d_pair_sum_im = new double[80];
    d_pair_sum_count = new double[80];
    
    
    //RESET ARRAYS
    reset_arrays();

    return ;
}
//****************************************************************************
void deallocate(){// Deallocate space
    
    int arrayMax=array_max+1;

    //GEOMETRIC ELEMENTS
    for(int i=0; i < arrayMax; i++){
        delete [] v[i]      ;
        delete [] v_pass[i];
        delete [] e[i];
        delete [] f[i];
        delete [] c_initPos[i]      ;
        delete [] c_prevPos[i]      ;
        delete [] c_DxDy[i]      ;
        delete [] basal_edges[i]   ;
        delete [] basal_vertices[i];
        delete [] basal_facets[i]  ;
        delete [] v_cells[i]  ;
        delete [] v_edges[i]  ;
        delete [] e_l[i];
        delete [] v_cell_neighbors[i]  ;
    }
    delete [] v;
    delete [] v_pass;
    delete [] e;
    delete [] f;
    delete [] c_initPos;
    delete [] c_prevPos;
    delete [] c_DxDy;
    delete [] basal_edges   ;
    delete [] basal_vertices;
    delete [] basal_facets  ;
    delete [] v_cells  ;
    delete [] v_cell_neighbors  ;
    delete [] v_edges  ;
    delete [] e_l;

    //vertices
    for(int i=0; i < arrayMax; i++){
        delete [] v_F[i]    ;
        delete [] v_vel[i]    ;
    }
    delete [] v_vel    ;
    delete [] v_T1dir   ;
    delete [] v_T1dir2   ;
    delete [] v_vertT1  ;
    delete [] v_edgeT1  ;
    delete [] v_clock  ;
    delete [] v_type  ;
    delete [] v_charge  ;
    delete [] defect_charge  ;
    delete [] defect_x  ;
    delete [] defect_y  ;

    //passive vertices
    delete [] v_pass_cell;
    
    //edges
    delete [] e_cell1   ;
    delete [] e_cell2   ;
    delete [] e_length;
    delete [] e_dldt;
    delete [] e_g;
    delete [] e_target;
    delete [] e_bin_categories;
    
    //facets
    delete [] f_cell;
    
    //cells
    delete [] c_cent;
    delete [] c_A0;
    delete [] c_dA;
    delete [] c_type;
    delete [] c_fixed;
    delete [] c_activity;

    delete[] c_hexagons;
    delete[] c_bin_categories;


    delete [] c_MSD_initial_x ;
    delete [] c_MSD_initial_y ;
    delete [] c_MSD_previous_x;
    delete [] c_MSD_previous_y;

    delete [] c_hexatic_initial_re ;
    delete [] c_hexatic_initial_im ;

    //grid
     for(int i=0; i < N_grid_max; i++){
        delete [] g_Q_xx[i]      ;
        delete [] g_Q_xy[i]      ;
        delete [] g_beta[i]      ;
        delete [] g_count[i]     ;
        delete [] g_winding[i]   ;
        delete [] g_v_x[i]       ;
        delete [] g_v_y[i]       ;
        delete [] g_v_count[i]   ;
        delete [] g_direction[i] ;
    }
    delete [] g_Q_xx;
    delete [] g_Q_xy;
    delete [] g_beta;
    delete [] g_count;
    delete [] g_winding;
    delete [] g_v_x;
    delete [] g_v_y;
    delete [] g_v_count;
    delete [] g_direction;
    
    //MISC
    delete [] perioXYZ;
    delete [] bin_positions;
    delete [] fine_bin_positions;
    delete [] chiral_cell_center;
    delete [] v_azimuthal_ave;
    delete [] vr_azimuthal_ave;
    delete [] vphi_azimuthal_ave;
    delete [] Q_type1_azimuthal_ave;
    delete [] Qphi_type1_azimuthal_ave;
    delete [] Q_type2_azimuthal_ave;
    delete [] Qphi_type2_azimuthal_ave;
    delete [] Kp_type1_azimuthal_ave;
    delete [] Kp_type2_azimuthal_ave;
    delete [] azimuthal_bin_T1_count_pop_11;
    delete [] azimuthal_bin_T1_count_pop_12;
    delete [] azimuthal_bin_T1_count_pop_22;
    delete [] azimuthal_edge_length;
    delete [] azimuthal_edge_count;

    delete [] d_pair_sum_re;
    delete [] d_pair_sum_im;
    delete [] d_pair_sum_count;

    return ;
}
//****************************************************************************
