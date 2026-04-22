//****************************************************************************
//****************************************************************************
//***********************RESOL CONTR******************************************
//****************************************************************************
//****************************************************************************
int rosetteResolution(double lfin){
    
    int total_resolutions = 0;
    int j;
    
    for(int i=1; i<=Nv; i++) if(v[i][0]>0.5){
            
        if(v_edges[i][2]>3){
            j=valence_reduction(i,lfin,v_T1dir[i],v_T1dir2[i]);
            e_length[j]=edge_length(j);
            e_g[j]=gamma0; //if(e_g[j]<0) e_g[j]=0;
            total_resolutions += 1;
        }
        
    }
    
    return total_resolutions;
}

//****************************************************************************
double T1_mean_angle(int i, bool& both_cells_nematic){

   int cell1_val = e_cell1[i];
   int cell2_val = e_cell2[i];

   //cell 1 elongation
   double nem_x_1 = 0, nem_y_1 = 0, nem_x_2 = 0, nem_y_2 = 0;
   
   compute_director(cell1_val, nem_x_1, nem_y_1);
   compute_director(cell2_val, nem_x_2, nem_y_2);

   double l_nem_1 = sqrt(nem_x_1*nem_x_1 + nem_y_1*nem_y_1);

   nem_x_1/=l_nem_1;
   nem_y_1/=l_nem_1;

   double l_nem_2 = sqrt(nem_x_2*nem_x_2 + nem_y_2*nem_y_2);

   nem_x_2/=l_nem_2;
   nem_y_2/=l_nem_2;

   both_cells_nematic = false;
   if (nonzero_vector(nem_x_1, nem_y_1) && nonzero_vector(nem_x_2, nem_y_2)) both_cells_nematic = true;

   int v1=e[i][1], v2=e[i][2];

   double *dxdydz = new double[2];
   dxdydz[0]=0; dxdydz[1]=0;

   double ddx,ddy,ddz;
   torus_dx_dy_dz(dxdydz,v1,v2);
   ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
   ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
   double lb=sqrt(ddx*ddx+ddy*ddy);

   delete []dxdydz;

   ddx/=lb;
   ddy/=lb;

   double nem_prod_1 = (nem_x_1*ddx + nem_y_1*ddy);
   double nem_theta_1;

   if (nem_prod_1 > 1.)  {nem_prod_1 = 1.;}
   if (nem_prod_1 < -1.) {nem_prod_1 = -1.;}
   nem_theta_1 = acos(nem_prod_1);
   if (nem_theta_1 > PI / 2.) nem_theta_1 = PI - nem_theta_1;

   double nem_prod_2 = (nem_x_2*ddx + nem_y_2*ddy);
   double nem_theta_2;

   if (nem_prod_2 > 1.)  {nem_prod_2 = 1.;}
   if (nem_prod_2 < -1.) {nem_prod_2 = -1.;}
   nem_theta_2 = acos(nem_prod_2);
   if (nem_theta_2 > PI / 2.) nem_theta_2 = PI - nem_theta_2;



   double mean_degree_angle = (fabs(nem_theta_1) + fabs(nem_theta_2))/2. * 360 /2. / PI;

   if (mean_degree_angle < 0 || mean_degree_angle > 360) cout << "angle out of bounds" << endl;

   return (fabs(nem_theta_1) + fabs(nem_theta_2))/2.;
}
//****************************************************************************
int edgeContraction(){

    int cell_edge_id_1;
    int cell_edge_id_2;
    
    int total_contractions = 0;
    for(int i=1; i<=Ne; i++) if(e[i][0]!=0){
        double elen=edge_length(i);
        e_dldt[i]=elen-e_length[i];
        e_length[i]=elen;
    }
    
    int vID;
    for(int i=1; i<=Ne; i++) if(e[i][0]!=0){
            
        if(  basal_edges[e_cell1[i]][2]==3 || (e_cell2[i]!=0 && basal_edges[e_cell2[i]][2]==3) || v_type[e[i][1]] == 1 || v_type[e[i][2]] == 1 ) continue;
                
        else if( e_length[i]<lth && v_edges[e[i][1]][2]<=3 && v_edges[e[i][2]][2]<=3 ){
   
            cell_edge_id_1 = e_cell1[i];
            cell_edge_id_2 = e_cell2[i];

            if (c_activity[cell_edge_id_1] != c_activity[cell_edge_id_2]) {
                T1_count_pop_12 += 1;
                azimuthal_bin_T1_count_pop_12[e_bin_categories[i]] +=1;
                }
            if (c_activity[cell_edge_id_1] == 1 && c_activity[cell_edge_id_2] == 1) {
                T1_count_pop_11 += 1;
                azimuthal_bin_T1_count_pop_11[e_bin_categories[i]] +=1;
                }
            if (c_activity[cell_edge_id_1] == 2 && c_activity[cell_edge_id_2] == 2) {
                T1_count_pop_22 += 1;
                azimuthal_bin_T1_count_pop_22[e_bin_categories[i]] +=1;
                }


            vID=merge_vertices(i);
            if(v_edges[vID][2]>3) v_clock[vID]=0;
            total_contractions += 1;
        }
            
    }
    
    return total_contractions;
}
//****************************************************************************
void update_topology(double lfin, int dynamics_type){
    
    int total_topology_changes = 0;
    total_topology_changes += edgeContraction();
    total_topology_changes += rosetteResolution(lfin);
    
    if (total_topology_changes != 0 && dynamics_type == 1) update_M();

    return;
}
