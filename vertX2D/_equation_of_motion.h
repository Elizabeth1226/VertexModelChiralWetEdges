//****************************************************************************
//****************************************************************************
//***********************EQUATION OF MOTION***********************************
//****************************************************************************
//****************************************************************************
void reset_line_tensions(){
    for(int i = 1; i <= Ne; i++) if(e[i][0]!=0){
        //TENSION DYNAMICS
        e_g[i] = gamma0;
        
    }
}
//****************************************************************************
void set_target_tensions_gamma0(){
    for(int i = 1; i <= Ne; i++) if(e[i][0]!=0){
        //TENSION DYNAMICS
        e_target[i] = gamma0;
        
    }
}
//****************************************************************************
void set_nematic_tensions(){
    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {
        c_nematic_edge_tension(i);
    }
}
//****************************************************************************
void line_tension_fluctuations(){
    for(int i = 1; i <= Ne; i++) if(e[i][0]!=0){
        
        //TENSION DYNAMICS
        e_g[i]+=-h*(1/tauM)*(e_g[i]-e_target[i])+sqrt(2.*Sig*Sig*h/tauM)*dis_sigma(gen_sigma);
        
    }
}
//****************************************************************************
void calc_forces(){
    

    //AREA TERM
    wA=0;
    wP=0;
    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {


        wA += c_AreaCompressibility_force_New(i);
        wP += c_perimeter_force(i);
        c_chiral_force(i);


    }
    


    // //PERIMETER TERM & TENSION DYNAMICS
    // wl=0;
    // for(int i = 1; i <= Ne; i++) if(e[i][0]!=0){
        
    //     //FORCES
    //     wl += e_BasalLength_force(i);
        
    // }
}
//****************************************************************************
void calc_forces_output(FILE *file_cell_force, FILE *file_cell_perimeter_force, FILE *file_cell_area_force, FILE *file_cell_elastic_force, FILE *file_cell_chiral_force, double image_time){
    
    double *cell_force_old_x = new double[10];
    double *cell_force_old_y = new double[10];
    double *cell_perimeter_force_old_x = new double[10];
    double *cell_perimeter_force_old_y = new double[10];
    double *cell_chiral_force_old_x = new double[10];
    double *cell_chiral_force_old_y = new double[10];
    double *cell_force_x = new double[10];
    double *cell_force_y = new double[10];
    double *cell_perimeter_force_x = new double[10];
    double *cell_perimeter_force_y = new double[10];
    double *cell_area_force_x = new double[10];
    double *cell_area_force_y = new double[10];
    double *cell_elastic_force_x = new double[10];
    double *cell_elastic_force_y = new double[10];
    double *cell_chiral_force_x = new double[10];
    double *cell_chiral_force_y = new double[10];

    //AREA TERM
    wA=0;
    wP=0;
    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {

        for (int j = 1; j <= basal_edges[i][2]; j ++){
            cell_force_old_x[j-1]= v_F[basal_vertices[i][j+2]][1];
            cell_force_old_y[j-1]= v_F[basal_vertices[i][j+2]][2];
        }
        


        wA += c_AreaCompressibility_force_New(i);


        fprintf(file_cell_area_force, "%d\t", i);
        for (int j = 1; j <= basal_edges[i][2]; j++){
            cell_area_force_x[j-1]= v_F[basal_vertices[i][j+2]][1]-cell_force_old_x[j-1];
            cell_area_force_y[j-1]= v_F[basal_vertices[i][j+2]][2]-cell_force_old_y[j-1];
            cell_perimeter_force_old_x[j-1]= v_F[basal_vertices[i][j+2]][1];
            cell_perimeter_force_old_y[j-1]= v_F[basal_vertices[i][j+2]][2];

            fprintf(file_cell_area_force, "%.10lf\t%.10lf\t", cell_area_force_x[j-1],cell_area_force_y[j-1]);
        }
        fprintf(file_cell_area_force, "\n");




        wP += c_perimeter_force(i);

        fprintf(file_cell_perimeter_force, "%d\t", i);
        for (int j = 1; j <= basal_edges[i][2]; j++){
            cell_perimeter_force_x[j-1]= v_F[basal_vertices[i][j+2]][1]-cell_perimeter_force_old_x[j-1];
            cell_perimeter_force_y[j-1]= v_F[basal_vertices[i][j+2]][2]-cell_perimeter_force_old_y[j-1];

            fprintf(file_cell_perimeter_force, "%.10lf\t%.10lf\t", cell_perimeter_force_x[j-1],cell_perimeter_force_y[j-1]);
        }
        fprintf(file_cell_perimeter_force, "\n");




        fprintf(file_cell_elastic_force, "%d\t", i);

        for (int j = 1; j <= basal_edges[i][2]; j++){
            cell_elastic_force_x[j-1]= v_F[basal_vertices[i][j+2]][1]-cell_force_old_x[j-1];
            cell_elastic_force_y[j-1]= v_F[basal_vertices[i][j+2]][2]-cell_force_old_y[j-1];
            cell_chiral_force_old_x[j-1]= v_F[basal_vertices[i][j+2]][1];
            cell_chiral_force_old_y[j-1]= v_F[basal_vertices[i][j+2]][2];

            fprintf(file_cell_elastic_force, "%.10lf\t%.10lf\t", cell_elastic_force_x[j-1],cell_elastic_force_y[j-1]);
        }
        fprintf(file_cell_elastic_force, "\n");




        c_chiral_force(i);

        


        fprintf(file_cell_force, "%d\t", i);
        fprintf(file_cell_chiral_force, "%d\t", i);

        for (int j = 1; j <= basal_edges[i][2]; j ++){
            cell_force_x[j-1]= v_F[basal_vertices[i][j+2]][1]-cell_force_old_x[j-1];
            cell_force_y[j-1]= v_F[basal_vertices[i][j+2]][2]-cell_force_old_y[j-1];
            cell_chiral_force_x[j-1] = v_F[basal_vertices[i][j+2]][1]-cell_chiral_force_old_x[j-1];
            cell_chiral_force_y[j-1] = v_F[basal_vertices[i][j+2]][2]-cell_chiral_force_old_y[j-1];

            fprintf(file_cell_force, "%.10lf\t%.10lf\t", cell_force_x[j-1],cell_force_y[j-1]);
            fprintf(file_cell_chiral_force, "%.10lf\t%.10lf\t", cell_chiral_force_x[j-1],cell_chiral_force_y[j-1]);
        }
        fprintf(file_cell_force, "\n");
        fprintf(file_cell_chiral_force, "\n");

    }
    


    // //PERIMETER TERM & TENSION DYNAMICS
    // wl=0;
    // for(int i = 1; i <= Ne; i++) if(e[i][0]!=0){
        
    //     //FORCES
    //     wl += e_BasalLength_force(i);
        
    // }
}
//****************************************************************************
void reset_forces(){
    for(int i = 1; i <= Nv; i++) if(v[i][0]>0.5) for(int j = 1; j<=2; j++) {v_F[i][j] = 0; v_vel[i][j] = 0;}
}
//****************************************************************************
void calc_central_vertices(){
    for(int i = 1; i<=Nc; i++) if(basal_edges[i][1]!=0) basal_center(i,0);
}
//****************************************************************************
double eqOfMotion_dry(char filename_description[],double image_time){
    

    if (Time > image_time -h*1.5){

        char filename_cell_force[500];
        char filename_cell_perimeter_force[500];
        char filename_cell_area_force[500];
        char filename_cell_elastic_force[500];
        char filename_cell_chiral_force[500];

        snprintf(filename_cell_force, sizeof(char)*500, "%s/cell_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_perimeter_force, sizeof(char)*500, "%s/cell_perimeter_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_area_force, sizeof(char)*500, "%s/cell_area_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_elastic_force, sizeof(char)*500, "%s/cell_elastic_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_chiral_force, sizeof(char)*500, "%s/cell_chiral_force_%g.dat",filename_description, image_time);

        FILE *file_cell_force; file_cell_force = fopen(filename_cell_force, "wt");
        FILE *file_cell_perimeter_force; file_cell_perimeter_force = fopen(filename_cell_perimeter_force, "wt");
        FILE *file_cell_area_force; file_cell_area_force = fopen(filename_cell_area_force, "wt");
        FILE *file_cell_elastic_force; file_cell_elastic_force = fopen(filename_cell_elastic_force, "wt");
        FILE *file_cell_chiral_force; file_cell_chiral_force = fopen(filename_cell_chiral_force, "wt");
        
        reset_forces();
        calc_forces_output(file_cell_force,file_cell_perimeter_force, file_cell_area_force, file_cell_elastic_force,file_cell_chiral_force,image_time);

        fclose(file_cell_force);
        fclose(file_cell_perimeter_force);
        fclose(file_cell_area_force);
        fclose(file_cell_elastic_force);
        fclose(file_cell_chiral_force);
        }


    //CALCULATES FORCES
    reset_forces();
    calc_forces();
    
    //EQUATION OF MOTION
    for(int i = 1; i <= Nv; i++) if(v[i][0]>0.5){
            
            //MOVES VERTICES
            double pm_val = 0;

            if (v_type[i] != 1) for(int j = 1; j<=2; j++) 
            {
                v[i][j] += h*v_F[i][j];
                pm_val += pow(v_F[i][j], 2.);
                v_vel[i][j] = v_F[i][j];
            }
            pm_val = h*sqrt(pm_val);
            if (pm_val > max_move) max_move = pm_val;
            
            //PERIODIC BOUNDARY CONDITIONS
            torus_vertex(i);
    }
    
    //FIXES PASSIVE VERTICES
    calc_central_vertices();
    
    //OUTPUT
    //printf("%g \t\t wA=%.20g \t\t wP=%.20g \t\t w=%.20g \t\t %g\n", Time, wA, wP, wA+wP, max_move);
    
    //t=t+dt
    Time+=h;
    
    return 0;
}
//****************************************************************************
double eqOfMotion_wet(char filename_description[],double image_time){
    double pm_val;
    //CALCULATES FORCES

    if (Time > image_time -h*1.5){

        char filename_cell_force[500];
        char filename_cell_perimeter_force[500];
        char filename_cell_area_force[500];
        char filename_cell_elastic_force[500];
        char filename_cell_chiral_force[500];

        snprintf(filename_cell_force, sizeof(char)*500, "%s/cell_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_perimeter_force, sizeof(char)*500, "%s/cell_perimeter_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_area_force, sizeof(char)*500, "%s/cell_area_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_elastic_force, sizeof(char)*500, "%s/cell_elastic_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_chiral_force, sizeof(char)*500, "%s/cell_chiral_force_%g.dat",filename_description, image_time);

        FILE *file_cell_force; file_cell_force = fopen(filename_cell_force, "wt");
        FILE *file_cell_perimeter_force; file_cell_perimeter_force = fopen(filename_cell_perimeter_force, "wt");
        FILE *file_cell_area_force; file_cell_area_force = fopen(filename_cell_area_force, "wt");
        FILE *file_cell_elastic_force; file_cell_elastic_force = fopen(filename_cell_elastic_force, "wt");
        FILE *file_cell_chiral_force; file_cell_chiral_force = fopen(filename_cell_chiral_force, "wt");
        
        reset_forces();
        calc_forces_output(file_cell_force,file_cell_perimeter_force, file_cell_area_force, file_cell_elastic_force,file_cell_chiral_force,image_time);

        fclose(file_cell_force);
        fclose(file_cell_perimeter_force);
        fclose(file_cell_area_force);
        fclose(file_cell_elastic_force);
        fclose(file_cell_chiral_force);
    }



    reset_forces();
    calc_forces();

    //INTEGRATOR
    Eigen::VectorXd rhs_x(_M.rows()), rhs_y(_M.rows());
    Eigen::VectorXd fx(_M.rows()), fy(_M.rows());
    Eigen::VectorXd x(_M.rows()), y(_M.rows());
    Eigen::VectorXd x_sol(_M.rows()), y_sol(_M.rows());

    // Populate the right hand side of the system
    for (int i = 1; i <= Nv; i++)
    {
      if (v[i][0]>0.5 && v_type[i] != 1)
      {
        x[i - 1] = v[i][1];
        y[i - 1] = v[i][2];

        fx[i - 1] = v_F[i][1];
        fy[i - 1] = v_F[i][2];
      }
      else
      {
        x[i - 1] = 0.;
        y[i - 1] = 0.;

        fx[i - 1] = 0.;
        fy[i - 1] = 0.;
      }
    }

    rhs_x = _M*x + h * fx;
    rhs_y = _M*y + h * fy;

    x_sol = _solver.solve(rhs_x);
    y_sol = _solver.solve(rhs_y);
    
    // Update positions
    for (int i = 1; i <= Nv; i++)
    {
      if (v[i][0]>0.5 && v_type[i] != 1)
      {
        
        
        pm_val = (x_sol[i - 1] - x[i - 1])*(x_sol[i - 1] - x[i - 1]) + (y_sol[i - 1] - y[i - 1])*(y_sol[i - 1] - y[i - 1]);
        
        pm_val = sqrt(pm_val);
        if (pm_val > max_move) max_move = pm_val;

        v[i][1] = x_sol[i - 1];
        v[i][2] = y_sol[i - 1];

        v_vel[i][1] =  (x_sol[i - 1] - x[i - 1])/h;
        v_vel[i][2] =  (y_sol[i - 1] - y[i - 1])/h;

        torus_vertex(i);
      }
    }
    
    
    //FIXES PASSIVE VERTICES
    calc_central_vertices();
    //t=t+dt
    Time+=h;
    
    return 0;
}
//****************************************************************************
double update_M(){

    double vertex_coordination;

    _M.resize(Nv, Nv);
    _M.reserve(5 * Nv);

    vector<Tripl> elems;

    
    for (int i = 1; i <= Nv; i++)
    {
        if (v[i][0] > 0.5 && v_type[i] != 1)
        {
            vertex_coordination = 3.;
            if (v_edges[i][2]>3) vertex_coordination = 4;
            elems.push_back(Tripl(i - 1,i - 1, _gamma_wet + _zeta_wet*vertex_coordination));
        }
        else
        {
            elems.push_back(Tripl(i - 1, i - 1, 1.0));
        }
    }

    int v1, v2;
    for (int i = 1; i <= Ne; i++)
    {
        if (e[i][0] != 0)
        {
            v1 = e[i][1];
            v2 = e[i][2];
            if (v_type[v1] != 1) elems.push_back(Tripl(v1 - 1, v2 - 1, -_zeta_wet));
            if (v_type[v2] != 1) elems.push_back(Tripl(v2 - 1, v1 - 1, -_zeta_wet));
        }
    }
    
    _M.setFromTriplets(elems.begin(), elems.end());
    _solver.analyzePattern(_M); 
    // Compute the numerical factorization 
    _solver.factorize(_M); 
    return 0;
}

//****************************************************************************
double update_edges_M(){
    //does not account for fixed edges
    double vertex_coordination;

    _M.resize(2 * Nv, 2* Nv);
    _M.reserve(10 * Nv);

    vector<Tripl> elems;

    
    for (int i = 1; i <= Nv; i++)
    {
        if (v[i][0] > 0.5 && v_type[i] != 1)
        {
            //add dry viscosity _gamma_wet to x and y row of vertex i
            elems.push_back(Tripl(2 * (i - 1),     2 * (i - 1)    , _gamma_wet));
            elems.push_back(Tripl(2 * (i - 1) + 1, 2 * (i - 1) + 1, _gamma_wet));
        }
        else
        {
            //add dry viscosity 1. to x and y row of vertex i
            elems.push_back(Tripl(2 * (i - 1),     2 * (i - 1)    , 1.));
            elems.push_back(Tripl(2 * (i - 1) + 1, 2 * (i - 1) + 1, 1.));
        }
    }

    int v1, v2;

    int coord_1;
    int coord_2;

    double t_x;
    double t_y;

    for (int i = 1; i <= Ne; i++)
    {
        if (e[i][0] != 0)
        {
            v1 = e[i][1];
            v2 = e[i][2];

            //coordinates of first rows of vertices v1 and v2
            coord_1 = 2 * (v1 - 1);
            coord_2 = 2 * (v2 - 1);

            //calculates edge direction vectors
            compute_edge_direction(v1, v2, t_x, t_y);

            //if v1 isn't fixed
            if (v_type[v1] != 1) 
            {
                //***********************************************
                //****************  X terms    ******************
                //***********************************************

                //self x-x term
                elems.push_back(Tripl(coord_1    , coord_1    ,       _zeta_wet * t_x * t_x));

                //self x-y term
                elems.push_back(Tripl(coord_1    , coord_1 + 1,       _zeta_wet * t_x * t_y));

                //ij x-x term
                elems.push_back(Tripl(coord_1    , coord_2    , -1. * _zeta_wet * t_x * t_x));
                
                //ij x-y term
                elems.push_back(Tripl(coord_1    , coord_2 + 1, -1. * _zeta_wet * t_x * t_y));

                //***********************************************
                //****************  Y terms    ******************
                //***********************************************

                //self y-y term
                elems.push_back(Tripl(coord_1 + 1, coord_1 + 1,       _zeta_wet * t_y * t_y));

                //self y-x term
                elems.push_back(Tripl(coord_1 + 1, coord_1    ,       _zeta_wet * t_x * t_y));

                //ij y-y term
                elems.push_back(Tripl(coord_1 + 1, coord_2 + 1, -1. * _zeta_wet * t_y * t_y));
                
                //ij y-x term
                elems.push_back(Tripl(coord_1 + 1, coord_2    , -1. * _zeta_wet * t_x * t_y));

            }
            //if v2 isn't fixed
            if (v_type[v2] != 1)
            {
                //***********************************************
                //****************  X terms    ******************
                //***********************************************

                //self x-x term
                elems.push_back(Tripl(coord_2    , coord_2    ,       _zeta_wet * t_x * t_x));

                //self x-y term
                elems.push_back(Tripl(coord_2    , coord_2 + 1,       _zeta_wet * t_x * t_y));

                //ij x-x term
                elems.push_back(Tripl(coord_2    , coord_1    , -1. * _zeta_wet * t_x * t_x));
                
                //ij x-y term
                elems.push_back(Tripl(coord_2    , coord_1 + 1, -1. * _zeta_wet * t_x * t_y));

                //***********************************************
                //****************  Y terms    ******************
                //***********************************************

                //self y-y term
                elems.push_back(Tripl(coord_2 + 1, coord_2 + 1,       _zeta_wet * t_y * t_y));

                //self y-x term
                elems.push_back(Tripl(coord_2 + 1, coord_2    ,       _zeta_wet * t_x * t_y));

                //ij y-y term
                elems.push_back(Tripl(coord_2 + 1, coord_1 + 1, -1. * _zeta_wet * t_y * t_y));
                
                //ij y-x term
                elems.push_back(Tripl(coord_2 + 1, coord_1    , -1. * _zeta_wet * t_x * t_y));
            }
        }
    }
    
    _M.setFromTriplets(elems.begin(), elems.end());
    _solver.analyzePattern(_M); 
    // Compute the numerical factorization 
    _solver.factorize(_M); 
    return 0;
}
//****************************************************************************
double eqOfMotion_edges(char filename_description[],double image_time){
    double pm_val;
    //CALCULATES FORCES

    if (Time > image_time -h*1.5){

        char filename_cell_force[500];
        char filename_cell_perimeter_force[500];
        char filename_cell_area_force[500];
        char filename_cell_elastic_force[500];
        char filename_cell_chiral_force[500];

        snprintf(filename_cell_force, sizeof(char)*500, "%s/cell_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_perimeter_force, sizeof(char)*500, "%s/cell_perimeter_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_area_force, sizeof(char)*500, "%s/cell_area_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_elastic_force, sizeof(char)*500, "%s/cell_elastic_force_%g.dat",filename_description, image_time);
        snprintf(filename_cell_chiral_force, sizeof(char)*500, "%s/cell_chiral_force_%g.dat",filename_description, image_time);

        FILE *file_cell_force; file_cell_force = fopen(filename_cell_force, "wt");
        FILE *file_cell_perimeter_force; file_cell_perimeter_force = fopen(filename_cell_perimeter_force, "wt");
        FILE *file_cell_area_force; file_cell_area_force = fopen(filename_cell_area_force, "wt");
        FILE *file_cell_elastic_force; file_cell_elastic_force = fopen(filename_cell_elastic_force, "wt");
        FILE *file_cell_chiral_force; file_cell_chiral_force = fopen(filename_cell_chiral_force, "wt");
        
        reset_forces();
        calc_forces_output(file_cell_force,file_cell_perimeter_force, file_cell_area_force, file_cell_elastic_force,file_cell_chiral_force,image_time);

        fclose(file_cell_force);
        fclose(file_cell_perimeter_force);
        fclose(file_cell_area_force);
        fclose(file_cell_elastic_force);
        fclose(file_cell_chiral_force);
    }


    reset_forces();
    calc_forces();

    //recalculates the dissipation matrix based on current geometry
    update_edges_M();

    //INTEGRATOR
    Eigen::VectorXd rhs_xy(_M.rows());
    Eigen::VectorXd fxy(_M.rows());
    Eigen::VectorXd xy(_M.rows());
    Eigen::VectorXd xy_sol(_M.rows());

    // Populate the right hand side of the system
    for (int i = 1; i <= Nv; i++)
    {
      if (v[i][0]>0.5 && v_type[i] != 1)
      {
        xy[2 * (i - 1)]     = v[i][1];
        xy[2 * (i - 1) + 1] = v[i][2];

        fxy[2 * (i - 1)]     = v_F[i][1];
        fxy[2 * (i - 1) + 1] = v_F[i][2];
      }
      else
      {
        xy[2 * (i - 1)]     = 0.;
        xy[2 * (i - 1) + 1] = 0.;

        fxy[2 * (i - 1)]     = 0.;
        fxy[2 * (i - 1) + 1] = 0.;
      }
    }

    //solves the coupled equations of motion
    rhs_xy = _M*xy + h * fxy;
    xy_sol = _solver.solve(rhs_xy);
    
    // Update positions
    for (int i = 1; i <= Nv; i++)
    {
      if (v[i][0]>0.5 && v_type[i] != 1)
      {
        
        
        //only relevant for storing max move
        pm_val = (xy_sol[2 * (i - 1)    ] - xy[2 * (i - 1)    ])*(xy_sol[2 * (i - 1)    ] - xy[2 * (i - 1)    ]) + (xy_sol[2 * (i - 1) + 1] - xy[2 * (i - 1) + 1])*(xy_sol[2 * (i - 1) + 1] - xy[2 * (i - 1) + 1]);
        pm_val = sqrt(pm_val);
        if (pm_val > max_move) max_move = pm_val;

        //sets the new x and y of vertex i, which are found at 2 * (i - 1) and 2 * (i - 1) + 1 of xy_sol
        v[i][1] = xy_sol[2 * (i - 1)    ];
        v[i][2] = xy_sol[2 * (i - 1) + 1];

        //calculates the new velocities of vertex i based on the change between xy and xy_sol
        v_vel[i][1] =  (xy_sol[2 * (i - 1)    ] - xy[2 * (i - 1)    ])/h;
        v_vel[i][2] =  (xy_sol[2 * (i - 1) + 1] - xy[2 * (i - 1) + 1])/h;

        //moves vertex back into periodic box if it escaped
        torus_vertex(i);
      }
    }
    
    
    //FIXES PASSIVE VERTICES
    calc_central_vertices();
    //t=t+dt
    Time+=h;
    
    return 0;
}
//****************************************************************************
void eqOfMotion(char filename_description[],double image_time,int dynamics_type){
    if (dynamics_type == 1) {eqOfMotion_wet(filename_description, image_time);}
    else if (dynamics_type == 2) {eqOfMotion_edges(filename_description, image_time);}
    else {eqOfMotion_dry(filename_description, image_time);}
}
//****************************************************************************
