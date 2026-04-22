//****************************************************************************
//****************************************************************************
//***********************INITIAL STRUCTURE************************************
//****************************************************************************
//****************************************************************************
void set_initial_regHex(int _Nx){
    
    //INTERNAL VARIABLES
    double *v_stitching_edge;//only useful for building regular hexagonal network
    v_stitching_edge = new double[array_max];
    int Ny=_Nx/2;
    double dd=1.2408064788027995;
    double ddx=0.08;
    
    //GLOBAL VARIABLES
    perioXYZ[0] = _Nx*sqrt(3)*dd/2.;
    perioXYZ[1] = Ny*(dd+dd/2);
    
    for(int i=1; i<=Ny; i++){
        for (int j =1; j<=_Nx; j++){
            double x0=(j-1)*sqrt(3)*dd/2+ddx;
            double y0=(3*dd/2)*(i-1)+ddx;
            make_vertex(
                        x0,
                        y0-dd/2+dd/2
                        );
            make_vertex(
                        x0+sqrt(3)*dd/4,
                        y0-dd/4+dd/2
                        );
            make_vertex(
                        x0+sqrt(3)*dd/4,
                        y0+dd/4+dd/2
                        );
            make_vertex(
                        x0,
                        y0+dd/2+dd/2
                        );
        }
        
        make_edge(i*_Nx*4-2, (i-1)*_Nx*4+1);
        make_edge((i-1)*_Nx*4+1, (i-1)*_Nx*4+2);
        make_edge((i-1)*_Nx*4+2, (i-1)*_Nx*4+3);
        make_edge((i-1)*_Nx*4+3, (i-1)*_Nx*4+4);
        make_edge((i-1)*_Nx*4+4, i*_Nx*4-1);
        
        for (int j =2; j<=_Nx; j++){
            make_edge((i-1)*_Nx*4+(j-2)*4+2, (i-1)*_Nx*4+(j-2)*4+5);
            make_edge((i-1)*_Nx*4+(j-2)*4+5, (i-1)*_Nx*4+(j-2)*4+6);
            make_edge((i-1)*_Nx*4+(j-2)*4+6, (i-1)*_Nx*4+(j-2)*4+7);
            make_edge((i-1)*_Nx*4+(j-2)*4+7, (i-1)*_Nx*4+(j-2)*4+8);
            make_edge((i-1)*_Nx*4+(j-2)*4+8, (i-1)*_Nx*4+(j-2)*4+3);
        }
    }
    
    for (int i =1; i<=Ny-1; i++){
        for (int j =1; j<=_Nx; j++){
            int edgid=make_edge(_Nx*4*(i-1)+4+(j-1)*4, _Nx*4*i+1+(j-1)*4);
            v_stitching_edge[_Nx*4*(i-1)+4+(j-1)*4]=edgid;
        }
    }
    
    int i = Ny;
    for (int j=1; j<=_Nx; j++){
        int edgid=make_edge(_Nx*4*(i-1)+4+(j-1)*4, 1+(j-1)*4);
        v_stitching_edge[_Nx*4*(i-1)+4+(j-1)*4]=edgid;
    }
    
    
    
    //*************************************
    //*************************************
    //*************************************
    //BASAL EDGES
    //*************************************x
    //LIHE VRSTICE
    for (int i=1; i<=Ny; i++){
        
        make_cell(
                  6,
                  (i-1)*_Nx*5+1,
                  (i-1)*_Nx*5+2,
                  (i-1)*_Nx*5+3,
                  (i-1)*_Nx*5+4,
                  (i-1)*_Nx*5+5,
                  -((i-1)*_Nx*5+(_Nx*5-2)),
                  0,0,0,0,0,0
                  );
        
        for (int j=6; j<=_Nx*5-4; j+=5){
            make_cell(
                      6,
                      (i-1)*_Nx*5+j,
                      (i-1)*_Nx*5+j+1,
                      (i-1)*_Nx*5+j+2,
                      (i-1)*_Nx*5+j+3,
                      (i-1)*_Nx*5+j+4,
                      -((i-1)*_Nx*5+j-3),
                      0,0,0,0,0,0
                      );
        }
    }
    
    //SODE VRSTICE
    for (int i=1; i<=Ny-1; i++){
        
        for (int j=1; j<=_Nx-1; j++){
            make_cell(
                      6,
                      -(4+(j-1)*5+(i-1)*_Nx*5),
                      -(10+(j-1)*5+(i-1)*_Nx*5),
                      v_stitching_edge[e[(10+(j-1)*5)+(i-1)*_Nx*5][1]],
                      -(5*_Nx+6+(j-1)*5+(i-1)*_Nx*5),
                      -(5*_Nx+6+(j-1)*5-4+(i-1)*_Nx*5),
                      -(v_stitching_edge[e[(10+(j-1)*5)+(i-1)*_Nx*5][1]]-1),
                      0,0,0,0,0,0
                      );
        }
        
        make_cell(
                  6,
                  -(4+(_Nx-1)*5+(i-1)*_Nx*5),
                  -(5+(i-1)*_Nx*5),
                  v_stitching_edge[e[10+(i-1)*_Nx*5][1]]-1,
                  -(5*_Nx+6-4-1+(i-1)*_Nx*5),
                  -(5*_Nx+6+(_Nx-1)*5-4+(i-1)*_Nx*5),
                  -(v_stitching_edge[e[(10+(_Nx-1-1)*5+(i-1)*_Nx*5)][1]]),
                  0,0,0,0,0,0
                  );
    }
    
    //ZADNJA VRSTICA
    for (int j=1; j<=_Nx-1; j++){
        
        make_cell(
                  6,
                  -(4+(j-1)*5+(Ny-1)*_Nx*5),
                  -(10+(j-1)*5+(Ny-1)*_Nx*5),
                  v_stitching_edge[e[(10+(j-1)*5)+(Ny-1)*_Nx*5][1]],
                  -(6+(j-1)*5),
                  -(2+(j-1)*5),
                  -(v_stitching_edge[e[(10+(j-1)*5)+(Ny-1)*_Nx*5][1]]-1),
                  0,0,0,0,0,0
                  );
    }
    
    make_cell(
              6,
              -(4+(_Nx-1)*5+(Ny-1)*_Nx*5),
              -(5+(Ny-1)*_Nx*5),
              v_stitching_edge[e[10+(Ny-1)*_Nx*5][1]]-1,
              -1,
              -(_Nx*5-3),
              -(v_stitching_edge[e[(10+(_Nx-1-1)*5+(Ny-1)*_Nx*5)][1]]),
              0,0,0,0,0,0
              );
    
    
    delete [] v_stitching_edge;
    
    //CREATE 3D TISSUE
    make_basal_side_ALL();
    
    Time=0;
    
    std::clog << "Initialization completed." << std::endl;
    printf("%g  %g\n", perioXYZ[0], perioXYZ[1]);
}
void set_initial_regHex_channel(int _Nx, int _Ny){
    
    //INTERNAL VARIABLES
    double *v_stitching_edge;//only useful for building regular hexagonal network
    v_stitching_edge = new double[array_max];
    int Ny=_Ny/2;
    double dd=1.2408064788027995;
    double ddx=0.08;
    
    //GLOBAL VARIABLES
    perioXYZ[0] = _Nx*sqrt(3)*dd/2.;
    perioXYZ[1] = Ny*(dd+dd/2);
    
    for(int i=1; i<=Ny; i++){
        for (int j =1; j<=_Nx; j++){
            double x0=(j-1)*sqrt(3)*dd/2+ddx;
            double y0=(3*dd/2)*(i-1)+ddx;
            make_vertex(
                        x0,
                        y0-dd/2+dd/2
                        );
            make_vertex(
                        x0+sqrt(3)*dd/4,
                        y0-dd/4+dd/2
                        );
            make_vertex(
                        x0+sqrt(3)*dd/4,
                        y0+dd/4+dd/2
                        );
            make_vertex(
                        x0,
                        y0+dd/2+dd/2
                        );
        }
        
        make_edge(i*_Nx*4-2, (i-1)*_Nx*4+1);
        make_edge((i-1)*_Nx*4+1, (i-1)*_Nx*4+2);
        make_edge((i-1)*_Nx*4+2, (i-1)*_Nx*4+3);
        make_edge((i-1)*_Nx*4+3, (i-1)*_Nx*4+4);
        make_edge((i-1)*_Nx*4+4, i*_Nx*4-1);
        
        for (int j =2; j<=_Nx; j++){
            make_edge((i-1)*_Nx*4+(j-2)*4+2, (i-1)*_Nx*4+(j-2)*4+5);
            make_edge((i-1)*_Nx*4+(j-2)*4+5, (i-1)*_Nx*4+(j-2)*4+6);
            make_edge((i-1)*_Nx*4+(j-2)*4+6, (i-1)*_Nx*4+(j-2)*4+7);
            make_edge((i-1)*_Nx*4+(j-2)*4+7, (i-1)*_Nx*4+(j-2)*4+8);
            make_edge((i-1)*_Nx*4+(j-2)*4+8, (i-1)*_Nx*4+(j-2)*4+3);
        }
    }
    
    for (int i =1; i<=Ny-1; i++){
        for (int j =1; j<=_Nx; j++){
            int edgid=make_edge(_Nx*4*(i-1)+4+(j-1)*4, _Nx*4*i+1+(j-1)*4);
            v_stitching_edge[_Nx*4*(i-1)+4+(j-1)*4]=edgid;
        }
    }
    
    int i = Ny;
    for (int j=1; j<=_Nx; j++){
        int edgid=make_edge(_Nx*4*(i-1)+4+(j-1)*4, 1+(j-1)*4);
        v_stitching_edge[_Nx*4*(i-1)+4+(j-1)*4]=edgid;
    }
    
    
    
    //*************************************
    //*************************************
    //*************************************
    //BASAL EDGES
    //*************************************x
    //LIHE VRSTICE
    for (int i=1; i<=Ny; i++){
        
        make_cell(
                  6,
                  (i-1)*_Nx*5+1,
                  (i-1)*_Nx*5+2,
                  (i-1)*_Nx*5+3,
                  (i-1)*_Nx*5+4,
                  (i-1)*_Nx*5+5,
                  -((i-1)*_Nx*5+(_Nx*5-2)),
                  0,0,0,0,0,0
                  );
        
        for (int j=6; j<=_Nx*5-4; j+=5){
            make_cell(
                      6,
                      (i-1)*_Nx*5+j,
                      (i-1)*_Nx*5+j+1,
                      (i-1)*_Nx*5+j+2,
                      (i-1)*_Nx*5+j+3,
                      (i-1)*_Nx*5+j+4,
                      -((i-1)*_Nx*5+j-3),
                      0,0,0,0,0,0
                      );
        }
    }
    
    //SODE VRSTICE
    for (int i=1; i<=Ny-1; i++){
        
        for (int j=1; j<=_Nx-1; j++){
            make_cell(
                      6,
                      -(4+(j-1)*5+(i-1)*_Nx*5),
                      -(10+(j-1)*5+(i-1)*_Nx*5),
                      v_stitching_edge[e[(10+(j-1)*5)+(i-1)*_Nx*5][1]],
                      -(5*_Nx+6+(j-1)*5+(i-1)*_Nx*5),
                      -(5*_Nx+6+(j-1)*5-4+(i-1)*_Nx*5),
                      -(v_stitching_edge[e[(10+(j-1)*5)+(i-1)*_Nx*5][1]]-1),
                      0,0,0,0,0,0
                      );
        }
        
        make_cell(
                  6,
                  -(4+(_Nx-1)*5+(i-1)*_Nx*5),
                  -(5+(i-1)*_Nx*5),
                  v_stitching_edge[e[10+(i-1)*_Nx*5][1]]-1,
                  -(5*_Nx+6-4-1+(i-1)*_Nx*5),
                  -(5*_Nx+6+(_Nx-1)*5-4+(i-1)*_Nx*5),
                  -(v_stitching_edge[e[(10+(_Nx-1-1)*5+(i-1)*_Nx*5)][1]]),
                  0,0,0,0,0,0
                  );
    }
    
    //ZADNJA VRSTICA
    for (int j=1; j<=_Nx-1; j++){
        
        make_cell(
                  6,
                  -(4+(j-1)*5+(Ny-1)*_Nx*5),
                  -(10+(j-1)*5+(Ny-1)*_Nx*5),
                  v_stitching_edge[e[(10+(j-1)*5)+(Ny-1)*_Nx*5][1]],
                  -(6+(j-1)*5),
                  -(2+(j-1)*5),
                  -(v_stitching_edge[e[(10+(j-1)*5)+(Ny-1)*_Nx*5][1]]-1),
                  0,0,0,0,0,0
                  );
    }
    
    make_cell(
              6,
              -(4+(_Nx-1)*5+(Ny-1)*_Nx*5),
              -(5+(Ny-1)*_Nx*5),
              v_stitching_edge[e[10+(Ny-1)*_Nx*5][1]]-1,
              -1,
              -(_Nx*5-3),
              -(v_stitching_edge[e[(10+(_Nx-1-1)*5+(Ny-1)*_Nx*5)][1]]),
              0,0,0,0,0,0
              );
    
    
    delete [] v_stitching_edge;
    
    //CREATE 3D TISSUE
    make_basal_side_ALL();
    
    Time=0;
    
    std::clog << "Initialization completed." << std::endl;
    printf("%g  %g\n", perioXYZ[0], perioXYZ[1]);
}
//****************************************************************************
void set_initial_fromFile2D(){
    
    char filename2[50];
    snprintf(filename2, sizeof(char) * 50, "./initial/filename.vt2d");
    FILE *file1; file1 = fopen(filename2, "rt");
    
    //V,E,C
    int nrV=0, nrE=0, nrC=0;
    int in_val = fscanf(file1, "%d  %d  %d\n", &nrV, &nrE, &nrC);
    
    //perioXYZ
    in_val = fscanf(file1, "%lf  %lf\n", &perioXYZ[0], &perioXYZ[1]);
    
    //VERTICES
    double xx, yy, vclock, t1vdir;
    int vID;
    for(int i=1; i<=nrV; i++){
        in_val = fscanf(file1, "%lf  %lf  %lf  %lf\n", &xx, &yy, &vclock, &t1vdir);
        vID=make_vertex(xx,yy);
        v_clock[vID]=vclock;
        v_T1dir[vID]=t1vdir;
    }
    
    //EDGES
    int v1, v2; double eG;
    for(int i=1; i<=nrE; i++){
        in_val = fscanf(file1, "%d  %d  %lf\n", &v1, &v2, &eG);
        int EIDd=make_edge(v1,v2);
        e_g[EIDd]=eG;
    }
    
    //CELLS
    int b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13;
    for(int i=1; i<=nrC; i++){
        in_val = fscanf(file1, "%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n", &b1, &b2, &b3, &b4, &b5, &b6, &b7, &b8, &b9, &b10, &b11, &b12, &b13);
        make_cell(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13);
    }
    fclose(file1);
    
    //CREATE 3D TISSUE
    make_basal_side_ALL();
    
    Time=0;
    
    //std::clog << "Initialization completed." << std::endl;
}
//****************************************************************************
void rescale_box(double desiredL){
    
    for(int i=1; i<=Nv; i++){
        v[i][1]*=desiredL/perioXYZ[0];
        v[i][2]*=desiredL/perioXYZ[1];
    }
    
    perioXYZ[0] = desiredL;
    perioXYZ[1] = desiredL;
    
    for(int i=1; i<=Nv; i++) torus_vertex(i);
    
    calc_central_vertices();
}
//****************************************************************************
void rescale_box_ratio(double box_ratio_val){
    
    for(int i=1; i<=Nv; i++){
        v[i][1]*=box_ratio_val;
        v[i][2]*=1/box_ratio_val;
    }
    
    perioXYZ[0] = perioXYZ[0] * box_ratio_val;
    perioXYZ[1] = perioXYZ[1] / box_ratio_val;
    
    for(int i=1; i<=Nv; i++) torus_vertex(i);
    
    calc_central_vertices();
}
//****************************************************************************
//****************************************************************************
void select_activities(double a2_probability, int pert_seed)
{   
    std::mt19937 gen(pert_seed);
    std::uniform_real_distribution<double> dis(0,1);

    for (int i = 1; i <= Nc; i++) 
    {
        c_activity[i] = 1;
        if (dis(gen) < a2_probability) c_activity[i] = 2;
    } 

    return;
}
//****************************************************************************
void one_activities()
{   
    for (int i = 1; i <= Nc; i++) 
    {
        c_activity[i] = 1;
    } 

    return;
}
//****************************************************************************
void line_activities(double a2_probability)
{   
    double xc_val;

    for (int i = 1; i <= Nc; i++) 
    {
        xc_val = center_coord(i, 1);
        c_activity[i] = 1;
        if (xc_val > perioXYZ[0] * (0.5 - a2_probability / 2.) && xc_val < perioXYZ[0] * (0.5 + a2_probability / 2.)) c_activity[i] = 2;
    } 

    return;
}
//****************************************************************************
void gap_activities(double a2_probabilityX, double a2_probabilityY, int ac1, int ac2)
{   
    double xc_val, yc_val;

    for (int i = 1; i <= Nc; i++) 
    {
        xc_val = center_coord(i, 1);
        yc_val = center_coord(i, 2);
        c_activity[i] = ac1;
        if (
            (
                xc_val < perioXYZ[0] * (0.5 - a2_probabilityX / 2.) || 
                xc_val > perioXYZ[0] * (0.5 + a2_probabilityX / 2.)
            ) &&
            yc_val > perioXYZ[1] * (a2_probabilityY) && 
            yc_val < perioXYZ[1] * (1. - a2_probabilityY) 
            
            ) c_activity[i] = ac2;
    } 

    return;
}
//****************************************************************************
void half_activities(double a2_probability1, double a2_probability2, int ac1, int ac2, int pert_seed)
{   
    std::mt19937 gen(pert_seed);
    std::uniform_int_distribution<> dis(1,(int) Nc);

    for (int i = 1; i <= Nc; i++) 
    {
        c_activity[i] = ac1;
    } 

    double selected_cells = 0;
    int rnd_cell_id;
    
    while (selected_cells < a2_probability1 * ((double) Nc) - 0.00001)
    {
        rnd_cell_id = dis(gen);
        if (c_activity[rnd_cell_id] == ac1)
        {
            c_activity[rnd_cell_id] = ac2;
            selected_cells += 1.;
        }
    }

    return;
}
//****************************************************************************
void dot_activities(double a2_probability1, double a2_probability2, int ac1, int ac2, int pert_seed)
{   
    double xc_val, yc_val;

    for (int i = 1; i <= Nc; i++) 
    {
        xc_val = center_coord(i, 1);
        yc_val = center_coord(i, 2);
        c_activity[i] = ac1;
        if (
            (
		xc_val > perioXYZ[0] * (0.5 - a2_probability1 / 2.) && 
		xc_val < perioXYZ[0] * (0.5 + a2_probability1 / 2.)
	    ) &&
            yc_val > perioXYZ[1] * (0.5 - a2_probability2 / 2.) && 
            yc_val < perioXYZ[1] * (0.5 + a2_probability2 / 2.) 
            
            ) c_activity[i] = ac2;
    } 

    return;
}
//****************************************************************************
void separate_activities(double a2_probability1, double a2_probability2, int ac1, int ac2, int pert_seed)
{   
    double xc_val, yc_val;

    for (int i = 1; i <= Nc; i++) 
    {
        xc_val = center_coord(i, 1);
        c_activity[i] = ac1;
        if (xc_val > perioXYZ[0] * (0.5 - a2_probability1 / 2.) && xc_val < perioXYZ[0] * (0.5 + a2_probability1 / 2.)) c_activity[i] = ac2;
    } 

    return;
}
//****************************************************************************
void pure_activities(double a2_probability1, double a2_probability2, int ac1, int ac2, int pert_seed)
{   
    std::mt19937 gen(pert_seed);
    std::uniform_int_distribution<> dis(1,(int) Nc);

    for (int i = 1; i <= Nc; i++) 
    {
        c_activity[i] = ac1;
    } 

    return;
}
void perturbe_vertices(double scale_pert, int pert_seed) //moves every vertex a bit
{   
    std::mt19937 gen(pert_seed);
    std::uniform_real_distribution<double> dis(0,1);
    
    double r_pert;
    double theta_pert;

    for (int i = 1; i <= Nv; i++) if (v_type[i] != 1)
    {
        r_pert = scale_pert * dis(gen);
        theta_pert = 2 * 3.141592654 * dis(gen);
        v[i][1] += r_pert * cos(theta_pert);
        v[i][2] += r_pert * sin(theta_pert);
        torus_vertex(i);
    } 
    calc_central_vertices();

    return;
}


