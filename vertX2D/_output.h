

//****************************************************************************
//****************************************************************************
//******************************OUTPUT****************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double out_Tissue(char filename_in[], double image_time){
    
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char) * 500, "%s/X_%g.dat", filename_in, image_time);
    snprintf(filename3, sizeof(char) * 500, "%s/Y_%g.dat", filename_in, image_time);
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    int vert_ref_id, vert_id;
    FILE *file2; file2 = fopen(filename2, "wt");
    FILE *file3; file3 = fopen(filename3, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            fprintf(file2, "%d ", basal_vertices[i][2]);
            fprintf(file3, "%d ", basal_vertices[i][2]);
            
            vert_ref_id=basal_vertices[i][3];
            fprintf(file2, "%.10lf ", v[vert_ref_id][1]);
            fprintf(file3, "%.10lf ", v[vert_ref_id][2]);
            
            for(int j = 4; j <= 2+basal_edges[i][2]; j++){
                vert_id=basal_vertices[i][j];
                torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
                fprintf(file2, "%.10lf ", v[vert_id][1] + dxdydz[0]);
                fprintf(file3, "%.10lf ", v[vert_id][2] + dxdydz[1]);
            }
            
            for(int j = basal_edges[i][2]; j <= 15; j++){
                fprintf(file2, "0 ");
                fprintf(file3, "0 ");
            }
            
            fprintf(file2, "\n");
            fprintf(file3, "\n");
        }
    }
    fclose(file2);
    fclose(file3);
    
    delete []dxdydz;
    
    return 0;
}
//****************************************************************************
double out_director(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/n_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            double nem_x, nem_y;
            compute_director_from_Q(i, nem_x, nem_y);

            double Q_xx, Q_xy, Q_yy;
            get_Q(i, Q_xx, Q_xy, Q_yy);

            fprintf(file2, "%d %.10lf %.10lf %.10lf %.10lf %.10lf", basal_vertices[i][2], nem_x, nem_y, Q_xx, Q_xy, Q_yy);            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
// double out_director_test(char filename_in[], double image_time){
    
    
//     char filename2[500]; 
//     snprintf(filename2, sizeof(char) * 500, "%s/n_test_%g.dat", filename_in, image_time);
    
//     FILE *file2; file2 = fopen(filename2, "wt");
    
//     for(int i=1; i<=Nc; i++){
//         if(basal_edges[i][1]!=0){
            
//             double nem_x, nem_y;
//             compute_director_from_Q_test(i, nem_x, nem_y);

//             double Q_xx, Q_xy, Q_yy;
//             get_Q(i, Q_xx, Q_xy, Q_yy);

//             fprintf(file2, "%d %.10lf %.10lf %.10lf %.10lf %.10lf", basal_vertices[i][2], nem_x, nem_y, Q_xx, Q_xy, Q_yy);            
//             fprintf(file2, "\n");
//         }
//     }
//     fclose(file2);
    
    
//     return 0;
// }



//****************************************************************************
double out_perio(char filename_in[], double image_time){
    
    
    char filename2[500];
    snprintf(filename2, sizeof(char) * 500, "%s/perio_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");

    
    fprintf(file2, "%lf %lf\n", perioXYZ[0], perioXYZ[1]);
    
    fclose(file2);
    
    
    return 0;
}

//****************************************************************************
double out_edgs(char filename_in[], double image_time){
    
    
    char filename2[500];
    snprintf(filename2, sizeof(char) * 500, "%s/EDG_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");

    int v1;
    int v2;
    
    for(int i=1; i<=Ne; i++){
        if(e[i][0]!=0){
            
            v1 = e[i][1];
            v2 = e[i][2];
            //cout << i << " " << v1 << " " << v2 << endl;
            //cout << v[v1][1] << " " << v[v1][2] << " " << v[v2][1] << " " << v[v2][2] << endl;
            fprintf(file2, "%d %d %d %d %lf %lf %lf %lf %lf\n", i, e[i][0], e[i][1], e[i][2], v[v1][1], v[v1][2], v[v2][1], v[v2][2], e_g[i]);
        }
        else 
        {
            fprintf(file2, "%d %d %d %d %lf %lf %lf %lf %lf\n", i, 0, 0, 0, 0., 0., 0., 0., 0.);
        }
    }
    fclose(file2);
    
    
    return 0;
}

//****************************************************************************
double out_vertex_positions(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/v_xy_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nv; i++){
        if(v[i][0]>0.5){
            fprintf(file2, "%d %g %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %d", i, v[i][0], v[i][1], v[i][2], v_vel[i][1], v_vel[i][2], v_F[i][1], v_F[i][2], v_type[i]);            
            fprintf(file2, "\n");
        } else {
            fprintf(file2, "%d %g %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %d", i, 0., 0., 0., 0., 0., 0., 0., 0);            
            fprintf(file2, "\n");

        }
    }
    fclose(file2);
    
    
    return 0;
}

//****************************************************************************
void compute_cell_velocity(int i, double& vel_x, double& vel_y){
    vel_x = 0;
    vel_y = 0;
    double vert_count = 0;
    for (int j = 1; j<= basal_vertices[i][2]; j++)
    {
        if (v_type[basal_vertices[i][j + 2]] != 1) {
            vel_x += v_vel[basal_vertices[i][j + 2]][1];
            vel_y += v_vel[basal_vertices[i][j + 2]][2];
            vert_count += 1.;
        }
    }
    vel_x/=vert_count;
    vel_y/=vert_count;

    return;
}
//****************************************************************************
void compute_cell_force(int i, double& vel_x, double& vel_y){
    vel_x = 0;
    vel_y = 0;
    double vert_count = 0;
    for (int j = 1; j<= basal_vertices[i][2]; j++)
    {
        if (v_type[basal_vertices[i][j + 2]] != 1) {
            vel_x += v_F[basal_vertices[i][j + 2]][1];
            vel_y += v_F[basal_vertices[i][j + 2]][2];
            vert_count += 1.;
        }
    }
    vel_x/=vert_count;
    vel_y/=vert_count;

    return;
}
//****************************************************************************
double out_cell_velocities(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/c_vel_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            double vel_x, vel_y, for_x, for_y;
            compute_cell_velocity(i, vel_x, vel_y);
            compute_cell_force(i, for_x, for_y);

            fprintf(file2, "%d %.10lf %.10lf %.10lf %.10lf %d", basal_vertices[i][2], vel_x, vel_y, for_x, for_y, c_fixed[i]);            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
double out_cell_vertices(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/c_vertices_%g.dat", filename_in, image_time);
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    int vert_ref_id, vert_id;
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            fprintf(file2, "%d %d", i, basal_vertices[i][2]);
            
            for(int j = 1; j <= basal_edges[i][2]; j++){
                fprintf(file2, " %d", basal_vertices[i][j + 2]);
            }
            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    delete []dxdydz;
    
    return 0;
}
//****************************************************************************
double out_cell_activities(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/c_activity_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            fprintf(file2, "%g\n", c_activity[i]);
        }
        else fprintf(file2, "0\n");   
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
double out_cell_neighbors(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/c_neighbors_%g.dat", filename_in, image_time);
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    int vert_ref_id, vert_id;
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            fprintf(file2, "%d %d", i, basal_vertices[i][2]);
            
            for(int j = 1; j <= basal_edges[i][2]; j++){
                fprintf(file2, " %d", other_cell(i, abs(basal_edges[i][j + 2])));
            }
            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    delete []dxdydz;
    
    return 0;
}
//****************************************************************************
double out_cell_edges(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/c_edges_%g.dat", filename_in, image_time);
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    int vert_ref_id, vert_id;
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            fprintf(file2, "%d %d", i, basal_vertices[i][2]);
            
            for(int j = 1; j <= basal_edges[i][2]; j++){
                fprintf(file2, " %d", basal_edges[i][j + 2]);
            }
            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    delete []dxdydz;
    
    return 0;
}


//****************************************************************************
double out_vertex_velocities(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/v_vel_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nv; i++){
        if(v[i][0]>0.5){

            fprintf(file2, "%.10lf %.10lf %.10lf %.10lf", v[i][1], v[i][2], v_F[i][1], v_F[i][2]);            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
void zero_vertex_neighbors(){
   for (int i = 1; i <= Nv; i++)
   {
      for (int j = 0; j<=8;j++) v_cell_neighbors[i][j] = 0;
   }
}
//****************************************************************************
void store_vertex_neighbors(){
   for (int i = 1; i <= Nc; i++) if (basal_edges[i][1] != 0)
   {
      for (int j = 1; j<=basal_edges[i][2];j++) 
      {
         int v1 = basal_vertices[i][j + 2];
         int current_found = v_cell_neighbors[v1][0];
         v_cell_neighbors[v1][current_found + 1] = i;
         v_cell_neighbors[v1][0] += 1;
      }
   }
}
//****************************************************************************
double out_vertex_neighbors(char filename_in[], double image_time){
    
    zero_vertex_neighbors();
    store_vertex_neighbors();
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/v_neighbors_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nv; i++){
        if(v[i][0]>0.5){

            fprintf(file2, "%d %d %d %d %d %d %d %d %d %d", i, v_cell_neighbors[i][0], v_cell_neighbors[i][1], v_cell_neighbors[i][2], v_cell_neighbors[i][3], v_cell_neighbors[i][4], v_cell_neighbors[i][5], v_cell_neighbors[i][6], v_cell_neighbors[i][7], v_cell_neighbors[i][8]);            
            fprintf(file2, "\n");
        }
        else {

            fprintf(file2, "%d 0 0 0 0 0 0 0 0 0", i);            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
double out_edge_neighbors(char filename_in[], double image_time){
    
    
    char filename2[500];
    snprintf(filename2, sizeof(char) * 500, "%s/e_neighbors_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");

    int v1;
    int v2;
    
    for(int i=1; i<=Ne; i++){
        if(e[i][0]!=0){
            
            v1 = e[i][1];
            v2 = e[i][2];
            //cout << i << " " << v1 << " " << v2 << endl;
            //cout << v[v1][1] << " " << v[v1][2] << " " << v[v2][1] << " " << v[v2][2] << endl;
            fprintf(file2, "%d %d %d\n", i, v1, v2);
        }
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
double out_defects(char filename_in[], double image_time){
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/defects_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for (int i = 0; i < N_defects; i ++ )
    {
        fprintf(file2, "%g %.10lf %.10lf\n", defect_charge[i], defect_x[i], defect_y[i]);            
    }
    fclose(file2);
    
    return 0;
}
//****************************************************************************
double out_T1(char filename_in[], double image_time){
    
    
    char filename2[500];
    snprintf(filename2, sizeof(char) * 500, "%s/T1_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");

    
    for(int i=0; i<=361; i++){
        fprintf(file2, "%d %d\n", i, T1_angle_count[i]);
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
double out_hexatic(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/hexatic_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){

            fprintf(file2, "%d %.10lf %.10lf %.10lf %.10lf", basal_vertices[i][2], re_psi_6(i), im_psi_6(i), abs_psi_6(i), cell_shape_index(i));            
            fprintf(file2, "\n");
        }
    }
    fclose(file2);
    
    
    return 0;
}
//****************************************************************************
double out_cell_angles(char filename_in[], double image_time){
    
    
    char filename2[500]; 
    snprintf(filename2, sizeof(char) * 500, "%s/cell_angles_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            for (int j = 1; j<=basal_edges[i][2]; j++)
            {
                int other_cell_val = other_cell(i, abs(basal_edges[i][j + 2]));
                int zero_count = 0;
                double angle_val = cell_cell_angle(i, other_cell_val, zero_count);
                fprintf(file2, "%d %d %.10lf %d %g %g %g %g\n", i, other_cell_val, angle_val, zero_count, c_activity[i], c_activity[other_cell_val], cell_edge_angle(i, abs(basal_edges[i][j + 2])), cell_edge_angle(other_cell_val, abs(basal_edges[i][j + 2])));            

            }
        }
    }
    fclose(file2);
    
    
    return 0;
}

//****************************************************************************
double out_avm_vertices(char filename_in[], double image_time){
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char) * 500, "%s/vertices_%g.dat", filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nv; i++){
        if (v[i][0] > 0.5){
            fprintf(file2, "%d %.22lf %.22lf\n", i - 1, v[i][1] - perioXYZ[0]/2., v[i][2] - perioXYZ[1]/2.);
        }
    }
    fclose(file2);
    
    return 0;
}
//****************************************************************************
double out_avm_faces(char filename_in[], double image_time){
    
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char) * 500, "%s/faces_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            fprintf(file2, "%d %d ", i - 1, basal_vertices[i][2]);
            
            for(int j = 3; j <= 2+basal_edges[i][2]; j++){
                fprintf(file2, "%d ", basal_vertices[i][j] - 1);
            }
            
            fprintf(file2, "0 0 0\n");
        }
    }
    fclose(file2);
    
    return 0;
}
//****************************************************************************
double out_avm_box(char filename_in[], double image_time){
    
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char) * 500, "%s/box_%g.dat", filename_in, image_time);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    
    fprintf(file2, "%.22lf %.22lf\n", perioXYZ[0], perioXYZ[1]);

    fclose(file2);
    
    return 0;
}
//****************************************************************************

double out_chiral_disk_area(){
    double current_area;

    double chiral_area=0;
    for (double i = 1; i <= Nc; i ++){
        current_area = CellArea_new(i);
        chiral_area += current_area;
    }

    
    return chiral_area;
}


//****************************************************************************** */
double out_c_spin(char filename_in[], double image_time){

    double Cspin_current;
    int v1;

    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char) * 500, "%s/cell_spin_%g.dat", filename_in, image_time);
    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 1; i <= Nc; i++){
        Cspin_current = 0;
        int vert_ref_id=basal_vertices[i][3];
        double vrefx = v[vert_ref_id][1];
        double vrefy = v[vert_ref_id][2];
        double cvx;
        double cvy;
        compute_cell_force(i,cvx,cvy);

        for (int j=3; j<=2+basal_vertices[i][2]; j++ ){

            v1=basal_vertices[i][j];
            //v1
            double v1x=v[v1][1];
            double v1y=v[v1][2];
            //x
            if(fabs(v1x-vrefx)>perioXYZ[0]/2.){
                if(v1x<vrefx) v1x+=perioXYZ[0];
                else if(v1x>vrefx) v1x-=perioXYZ[0];
            }
            //y
            if(fabs(v1y-vrefy)>perioXYZ[1]/2.){
                if(v1y<vrefy) v1y+=perioXYZ[1];
                else if(v1y>vrefy) v1y-=perioXYZ[1];
            }

            double u1x=v_F[v1][1]-cvx;
            double u1y=v_F[v1][2]-cvy;

            Cspin_current=Cspin_current+(v1x*u1y-v1y*u1x)/CellArea_new(i);
        }
        
        fprintf(file2, "%d %.22lf\n", i, Cspin_current);
    }

    fclose(file2);

    return 0;
}





//****************************************************************************
double out_azimuthal_average_v(char filename_in[], double image_time){

    double v_current;

    azimuthal_averaged_v();
    

    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_v_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        v_current = v_azimuthal_ave[i];
        if (v_current >=0.){
            fprintf(file2, "%.10lf %.10lf\n", bin_positions[i], v_current);
        }
    }

    fclose(file2);
    
    return 0;
}
//****************************************************************************
double out_azimuthal_average_vr(char filename_in[], double image_time){

    double vr_current;

    azimuthal_averaged_vr();

    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_vr_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        vr_current = vr_azimuthal_ave[i];
        if (vr_current != 99999.){
            fprintf(file2, "%.10lf %.10lf\n", bin_positions[i], vr_current);
        }
    }

    fclose(file2);

    return 0;
}
//****************************************************************************
double out_azimuthal_average_vphi(char filename_in[], double image_time){

    double vphi_current;

    azimuthal_averaged_vphi();

    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_vphi_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        vphi_current = vphi_azimuthal_ave[i];
        if (vphi_current != 99999.){
            fprintf(file2, "%.10lf %.10lf\n", bin_positions[i], vphi_current);
        }
    }

    fclose(file2);

    return 0;
}
//****************************************************************************
// double out_azimuthal_average_Q(char filename_in[], double image_time){

//     double grid_size, N_grid, Q_current;

//     grid_size = 0.5*perioXYZ[1]/sqrt(Nc);
//     N_grid = 0.5*perioXYZ[1]/grid_size;


//     char filename2[500]; char filename3[500];
//     snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_Q_%g.dat",filename_in, image_time);

//     FILE *file2; file2 = fopen(filename2, "wt");

//     for (double i = 0; i <= 0.5*perioXYZ[1]; i += grid_size){
//         Q_current = azimuthal_averaged_Q(i);
//         if (Q_current != 99999.){
//             fprintf(file2, "%.10lf %.10lf\n", i, Q_current);
//         }
//     }

//     fclose(file2);

//     return 0;
// }
//****************************************************************************
double out_azimuthal_average_Q_type(char filename_in[], double image_time){

    double Q_current;

    azimuthal_averaged_Q_type1();
    azimuthal_averaged_Q_type2();

    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_Q_type1_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        Q_current = Q_type1_azimuthal_ave[i];
        if (Q_current >=0.){
            fprintf(file2, "%.10lf %.10lf\n", bin_positions[i], Q_current);
        }
    }

    fclose(file2);


    char filename4[500]; char filename5[500];
    snprintf(filename4, sizeof(char)*500, "%s/azimuthal_ave_Q_type2_%g.dat",filename_in, image_time);

    FILE *file4; file4 = fopen(filename4, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        Q_current = Q_type2_azimuthal_ave[i];
        if (Q_current >=0.){
            fprintf(file4, "%.10lf %.10lf\n", bin_positions[i], Q_current);
        }
    }

    fclose(file4);

    return 0;
}
//****************************************************************************
// double out_azimuthal_average_Qphi(char filename_in[], double image_time){

//     double grid_size, N_grid, Qphi_current;

//     grid_size = 0.5*perioXYZ[1]/sqrt(Nc);
//     N_grid = 0.5*perioXYZ[1]/grid_size;


//     char filename2[500]; char filename3[500];
//     snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_Qphi_%g.dat",filename_in, image_time);

//     FILE *file2; file2 = fopen(filename2, "wt");

//     for (double i = 0; i <= 0.5*perioXYZ[1]; i += grid_size){
//         Qphi_current = azimuthal_averaged_Qphi(i);
//         fprintf(file2, "%.10lf %.10lf\n", i, Qphi_current);
//     }

//     fclose(file2);

//     return 0;
// }
//****************************************************************************
double out_azimuthal_average_Qphi_type(char filename_in[], double image_time){

    double Qphi_current;

    azimuthal_averaged_Qphi_type1();
    azimuthal_averaged_Qphi_type2();
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_Qphi_type1_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        Qphi_current = Qphi_type1_azimuthal_ave[i];
        if (Qphi_current != 99999.){
            fprintf(file2, "%.10lf %.10lf\n", bin_positions[i], Qphi_current);
        }
    }

    fclose(file2);



    char filename4[500]; char filename5[500];
    snprintf(filename4, sizeof(char)*500, "%s/azimuthal_ave_Qphi_type2_%g.dat",filename_in, image_time);

    FILE *file4; file4 = fopen(filename4, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        Qphi_current = Qphi_type2_azimuthal_ave[i];
        if (Qphi_current != 99999.){
            fprintf(file4, "%.10lf %.10lf\n", bin_positions[i], Qphi_current);
        }
    }

    fclose(file4);

    return 0;
}
//****************************************************************************
double out_azimuthal_average_Kp_type(char filename_in[], double image_time){

    double Kp_current;

    azimuthal_averaged_cell_shape_index_type1();
    azimuthal_averaged_cell_shape_index_type2();

    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_ave_Kp_type1_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        Kp_current = Kp_type1_azimuthal_ave[i];
        if (Kp_current >=0.){
            fprintf(file2, "%.10lf %.10lf\n", bin_positions[i], Kp_current);
        } 
    }

    fclose(file2);



    char filename4[500]; char filename5[500];
    snprintf(filename4, sizeof(char)*500, "%s/azimuthal_ave_Kp_type2_%g.dat",filename_in, image_time);

    FILE *file4; file4 = fopen(filename4, "wt");

    for (int i = 0; i <= sqrt(Nc); i ++){
        Kp_current = Kp_type2_azimuthal_ave[i];
        if (Kp_current >=0.){
            fprintf(file4, "%.10lf %.10lf\n", bin_positions[i], Kp_current);
        }
    }

    fclose(file4);

    return 0;
}
//****************************************************************************
double out_azimuthal_T1_count(char filename_in[], double image_time){



    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_T1_count_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= 5*sqrt(Nc); i ++){
        fprintf(file2, "%.10lf\t%d\t%d\t%d\n", fine_bin_positions[i], azimuthal_bin_T1_count_pop_11[i], azimuthal_bin_T1_count_pop_12[i], azimuthal_bin_T1_count_pop_22[i]);
    }

    fclose(file2);

    return 0;
}
//****************************************************************************
double out_azimuthal_edge_length(char filename_in[], double image_time){
    
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/azimuthal_edge_length_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= 5*sqrt(Nc); i ++){
        fprintf(file2, "%.10lf\t%.10lf\n", fine_bin_positions[i], azimuthal_edge_length[i]);
    }

    fclose(file2);

    return 0;

}

//****************************************************************************
double out_cell_bin_categories(char filename_in[], double image_time){
    
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/cell_bin_categories_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= Nc; i ++){
        fprintf(file2, "%d\n", c_bin_categories[i]);
    }

    fclose(file2);

    return 0;

}

//****************************************************************************
double out_edge_bin_categories(char filename_in[], double image_time){
    
    
    char filename2[500]; char filename3[500];
    snprintf(filename2, sizeof(char)*500, "%s/edge_bin_categories_%g.dat",filename_in, image_time);

    FILE *file2; file2 = fopen(filename2, "wt");

    for (int i = 0; i <= Ne; i ++){
        fprintf(file2, "%d\n", e_bin_categories[i]);
    }

    fclose(file2);

    return 0;

}


//****************************************************************************
void outWholeTissue(char filename_description[], double image_time){

    
    //tissue
    out_Tissue(filename_description, image_time);
    out_perio(filename_description, image_time);
    out_director(filename_description, image_time);
    // out_director_test(filename_description, image_time);
    out_edgs(filename_description, image_time);
    
    out_vertex_positions(filename_description, image_time);
    out_vertex_neighbors(filename_description, image_time);
    out_vertex_velocities(filename_description, image_time);

    out_cell_vertices(filename_description, image_time);
    out_cell_velocities(filename_description, image_time);
    out_cell_activities(filename_description, image_time); 
    out_cell_neighbors(filename_description, image_time); 
    out_cell_edges(filename_description, image_time); 
    out_c_spin(filename_description, image_time);

    azimuthal_ave_setup();
    out_azimuthal_average_Kp_type(filename_description, image_time);
    // out_azimuthal_average_Q(filename_description, image_time);
    out_azimuthal_average_Q_type(filename_description, image_time);
    // out_azimuthal_average_Qphi(filename_description, image_time);
    out_azimuthal_average_Qphi_type(filename_description, image_time);
    out_azimuthal_average_v(filename_description, image_time);
    out_azimuthal_average_vphi(filename_description, image_time);
    out_azimuthal_average_vr(filename_description, image_time);
    out_azimuthal_T1_count(filename_description, image_time);
    out_azimuthal_edge_length(filename_description, image_time);
    out_cell_bin_categories(filename_description, image_time);
    out_edge_bin_categories(filename_description, image_time);

    

    return;
}
//****************************************************************************
double sorting_index()
{
    double SI_sum = 0.;
    double SI_count = 0.;

    double n_c_same = 0;
    double n_c_other = 0;

    int other_cell_val;
    
    for(int i=1; i<=Nc; i++)
    {
        if(basal_edges[i][1]!=0)
        {
            n_c_same = 0;
            n_c_other = 0;
            for (int j = 1; j<=basal_edges[i][2]; j++)
            {
                other_cell_val = other_cell(i, abs(basal_edges[i][j + 2]));
                if (c_activity[i] == c_activity[other_cell_val]) {n_c_same += 1;}
                else {n_c_other += 1;}
            }   
            SI_sum += (n_c_same / (n_c_same + n_c_other));
            SI_count += 1;  
        }
    }
    return (SI_sum / SI_count);
}
//****************************************************************************
double print_output(double write_time, double delta_write_time, FILE *file_energies, FILE *file_sorting, FILE *file_msd, FILE *file_chiral_area, FILE *file_T1_transitions, bool print_MSD)
{
    
    fprintf(file_energies, "%g\t%g\t%g\t%g\t%g\t%g\n", Time, wA + wP + wl, wA, wP, wl, max_move);
    fprintf(file_sorting,  "%g\t%g\n", Time, sorting_index());
    fprintf(file_chiral_area, "%g\t%.15lf\n",Time, out_chiral_disk_area());
    fprintf(file_T1_transitions, "%g\t%d\t%d\t%d\n", Time, T1_count_pop_11, T1_count_pop_12, T1_count_pop_22);

    if (print_MSD) {
        fprintf(file_msd, "%g\t%g\t%g\t%g\n", Time, compute_MSD()/(1.*Nc), compute_normalised_MSD(1), compute_normalised_MSD(2));
    }

    cout << write_time << " " << endl;

    return write_time + delta_write_time;
}
//****************************************************************************
double image_output(double image_time, double delta_image_time, char filename_in[])
{
    //OUTPUTS TISSUE DATA

    outWholeTissue(filename_in, image_time);

    return image_time + delta_image_time;
}
//****************************************************************************
