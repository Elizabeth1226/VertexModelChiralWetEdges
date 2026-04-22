//****************************************************************************
//****************************************************************************
//*****************************FLUIDIZATION***********************************
//****************************************************************************
//****************************************************************************
void initialize_cell_positions(){
 
    for(int i=1; i<=Nc; i++){
        c_initPos[i][1]=v_pass[c_cent[i]][1];
        c_initPos[i][2]=v_pass[c_cent[i]][2];
    }   
    
    return;
}
//****************************************************************************
void initialize_previous_positions(){
 
    int pass_id;
    for(int i=1; i<=Nc; i++){
        pass_id=c_cent[i];
        c_prevPos[i][1] = v_pass[pass_id][1];
        c_prevPos[i][2] = v_pass[pass_id][2];
    }
    
    return;
}
//****************************************************************************
void update_boundary_for_MSD(){
 
    int pass_id;
    for(int i = 1; i<=Nc; i++)
    {   
        pass_id = c_cent[i];
        //x 
        if(fabs(v_pass[pass_id][1]-c_prevPos[i][1])>0.5*perioXYZ[0]){
            if      (v_pass[pass_id][1]<c_prevPos[i][1]) c_DxDy[i][1]+=perioXYZ[0];
            else if (v_pass[pass_id][1]>c_prevPos[i][1]) c_DxDy[i][1]-=perioXYZ[0];
         }   

        //y 
        if(fabs(v_pass[pass_id][2]-c_prevPos[i][2])>0.5*perioXYZ[1]){
            if      (v_pass[pass_id][2]<c_prevPos[i][2]) c_DxDy[i][2]+=perioXYZ[1];
            else if (v_pass[pass_id][2]>c_prevPos[i][2]) c_DxDy[i][2]-=perioXYZ[1];
        }   
        c_prevPos[i][1] = v_pass[pass_id][1];
        c_prevPos[i][2] = v_pass[pass_id][2];
    }
    
    return;
}
//****************************************************************************
double square_displacementtCORR(int i){ 
    
    double ddx=c_initPos[i][1]-(v_pass[c_cent[i]][1]+c_DxDy[i][1]);
    double ddy=c_initPos[i][2]-(v_pass[c_cent[i]][2]+c_DxDy[i][2]);
    
    return ddx*ddx+ddy*ddy;
}
//****************************************************************************
double lin_displacementtCORR(int i){ 
    
    double ddx=c_initPos[i][1]-(v_pass[c_cent[i]][1]+c_DxDy[i][1]);
    
    return -ddx;
}
//****************************************************************************
double compute_MSD(){ 
    
    double SDispl = 0;
    double MSD = 0;

    for (int i=1; i<=Nc; i++)
    {   
        SDispl=square_displacementtCORR(i);
        MSD += SDispl;
    }  
    return MSD;
}
//****************************************************************************
double compute_MD_X(){ 
    
    double SDispl = 0;
    double MSD = 0;

    for (int i=1; i<=Nc; i++)
    {   
        SDispl=lin_displacementtCORR(i);
        MSD += SDispl;
    }  
    return MSD;
}
//****************************************************************************
double compute_normalised_MSD(int activity_type){ 
    
    double SDispl = 0;
    double MSD = 0;
    double N_type = 0;

    for (int i=1; i<=Nc; i++)
    {   
        if (c_activity[i] == activity_type)
        {
            SDispl = square_displacementtCORR(i);
            MSD += SDispl;
            N_type += 1;
        }
    }  

    return MSD / N_type;
}
//****************************************************************************

//****************************************************************************
//****************************************************************************
//*****************************HEXATIC PHASE**********************************
//****************************************************************************
//****************************************************************************
int second_edg_cell(int edg_id, int first_cell_id)
{
    int second_cell_id = -1;
    if (e_cell1[edg_id] != first_cell_id) second_cell_id = e_cell1[edg_id];
    else if (e_cell2[edg_id] != first_cell_id) second_cell_id = e_cell2[edg_id];
    else cout << "no cell found" << endl;

    return second_cell_id;
}
//****************************************************************************
double cell_distance_coord(int c1, int cRef, int coord) //cRef stays in same place, gives distance from cRef to c1
{
    //cout << c1 << " " << cRef << " " << c_cent[c1] << " " << c_cent[cRef] << " ";

    double v1_coord = v_pass[c_cent[c1]][coord];
    double vRef_coord = v_pass[c_cent[cRef]][coord];

    //cout << v1_coord << " " << vRef_coord << endl;

    if(fabs(v1_coord-vRef_coord)>0.5*perioXYZ[coord - 1]){
        if      (v1_coord < vRef_coord)   v1_coord += perioXYZ[coord - 1];
        else if (v1_coord > vRef_coord)   v1_coord -= perioXYZ[coord - 1];
    }

    return v1_coord - vRef_coord;
}
//****************************************************************************
double cell_distance_abs(int c1, int cRef) //cRef stays in same place, gives distance from cRef to c1
{
    double x_distance = cell_distance_coord(c1, cRef, 1);
    double y_distance = cell_distance_coord(c1, cRef, 2);

    return sqrt(x_distance * x_distance + y_distance * y_distance);
}
//****************************************************************************

//****************************************************************************
//*****************************PSI_6******************************************
//****************************************************************************
double re_psi_6(int i){
    double out_val = 0;
    double n_neighbors = (double) basal_facets[i][2];

    double theta_val = 0.;

    for (int j = 1; j <= basal_facets[i][2]; j++)
    {
        int neighbor_cell = second_edg_cell(abs(basal_edges[i][j + 2]), i);

        double x_distance = cell_distance_coord(neighbor_cell, i, 1);
        double y_distance = cell_distance_coord(neighbor_cell, i, 2);

        theta_val = atan2(y_distance, x_distance);

        out_val += cos(6 * theta_val);
    }

    return (out_val / ((double) n_neighbors));
}
//****************************************************************************
double im_psi_6(int i){
    double out_val = 0;
    double n_neighbors = (double) basal_facets[i][2];

    double theta_val = 0.;

    for (int j = 1; j <= basal_facets[i][2]; j++)
    {
        int neighbor_cell = second_edg_cell(abs(basal_edges[i][j + 2]), i);

        double x_distance = cell_distance_coord(neighbor_cell, i, 1);
        double y_distance = cell_distance_coord(neighbor_cell, i, 2);

        theta_val = atan2(y_distance, x_distance);

        //cout << x_distance << " " << y_distance << " " << theta_val << endl;

        out_val += sin(6 * theta_val);
    }

    return (out_val / ((double) n_neighbors));
}
//****************************************************************************
double abs_psi_6(int i){
    double out_val_re = 0;
    double out_val_im = 0;
    double n_neighbors = (double) basal_facets[i][2];

    double theta_val = 0.;

    for (int j = 1; j <= basal_facets[i][2]; j++)
    {
        int neighbor_cell = second_edg_cell(abs(basal_edges[i][j + 2]), i);

        double x_distance = cell_distance_coord(neighbor_cell, i, 1);
        double y_distance = cell_distance_coord(neighbor_cell, i, 2);

        theta_val = atan2(y_distance, x_distance);

        out_val_re += cos(6 * theta_val);
        out_val_im += sin(6 * theta_val);

        //cout << x_distance << " " << y_distance << " " << theta_val << endl;

        //out_val += sqrt(cos(6 * theta_val)*cos(6 * theta_val) + sin(6 * theta_val)*sin(6 * theta_val));
    }

    return (sqrt(out_val_re*out_val_re+out_val_im*out_val_im) / ((double) n_neighbors));
}
//****************************************************************************
double tissue_psi_6(){
    double total_re_val = 0.;
    double total_im_val = 0.;
    double total_count = 0;
    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {
        total_re_val += re_psi_6(i);
        total_im_val += im_psi_6(i);
        total_count += 1;
    }
    return sqrt(total_re_val * total_re_val + total_im_val * total_im_val) / total_count;
}
//****************************************************************************
double tissue_psi_square_6(){
    double total_val = 0.;
    double total_count = 0;
    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {
        total_val += pow(re_psi_6(i), 2.) + pow(im_psi_6(i), 2.);
        total_count += 1;
    }
    return total_val / total_count;
}
//****************************************************************************
//*****************************DEFECTS****************************************
//****************************************************************************
double defect_density(){
    double n6 =0.;
    double nDef = 0.;

    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {
        if(basal_facets[i][2] == 6) {n6 += 1.;} else {nDef += 1.;}
    }
    return nDef / (n6 + nDef);
}
//****************************************************************************
int defect_neighbors(int i){

    int defect_count = 0;
    
    for (int j = 1; j <= basal_facets[i][2]; j++)
    {
        int neighbor_cell = second_edg_cell(abs(basal_edges[i][j + 2]), i);
        if (basal_facets[neighbor_cell][2] != 6) defect_count += 1;
    }

    return defect_count;

}
//****************************************************************************
int check_disclination(int i){
    int is_disclination = 0;
    if (basal_facets[i][2] == 5 || basal_facets[i][2] == 7)
    {   
        if (defect_neighbors(i) == 0) is_disclination = 1;
    }
    return is_disclination;
}
//****************************************************************************
int dislocation_heptagon(int i){

    if (basal_facets[i][2] == 7)
    {
        if (defect_neighbors(i) == 1) {return 1;} else return 0;
    }
    else return 0;
}
//****************************************************************************
int check_dislocation(int i){
    if (basal_facets[i][2] == 5)
    {   
        int total_hexagon_neighbors = 0;
        int total_dislocation_heptagons = 0;
        for (int j = 1; j <= basal_facets[i][2]; j++)
        {
            int neighbor_cell = second_edg_cell(abs(basal_edges[i][j + 2]), i);

            if (basal_facets[neighbor_cell][2] == 6) {total_hexagon_neighbors += 1;}
            else {total_dislocation_heptagons += dislocation_heptagon(neighbor_cell);}

        }
        if (total_hexagon_neighbors == 4 && total_dislocation_heptagons == 1) {return 1;} else return 0;

    }
    else return 0;
}
//****************************************************************************
double dislocation_density(){
    double nAll =0.;
    double nDef = 0.;

    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {
        nAll += 1;
        nDef += check_dislocation(i);
    }
    return nDef / (nAll);
}
//****************************************************************************
double disclination_density(){
    double nAll =0.;
    double nDef = 0.;

    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {
        nAll += 1;
        nDef += check_disclination(i);
    }
    return nDef / (nAll);
}

//****************************************************************************
//*****************************SHAPE INDEX************************************
//****************************************************************************
double cell_shape_index(int cell_i)
{
    double current_A = CellArea_new(cell_i);
    double current_P = CellPerimeter(cell_i);
    return current_P / sqrt(current_A);
}
//****************************************************************************
double tissue_shape_index()
{
    double shape_index_sum = 0;
    double shape_index_count = 0;
    for (int i = 1; i <= Nc; i ++)
    {
        if (basal_edges[i][1] > 0 )
        {
            shape_index_sum += cell_shape_index(i);
            shape_index_count += 1;
        }
    }
    return shape_index_sum / shape_index_count;
}
//****************************************************************************
double tissue_shape_index_susceptibility()
{
    double shape_index_sum = 0;
    double shape_index_square = 0;
    double shape_index_count = 0;
    double current_shape_index;
    for (int i = 1; i <= Nc; i ++)
    {
        if (basal_edges[i][1] > 0 )
        {
            current_shape_index = cell_shape_index(i);
            shape_index_sum += current_shape_index;
            shape_index_square += current_shape_index * current_shape_index;
            shape_index_count += 1;
        }
    }
    return shape_index_square / shape_index_count - pow(shape_index_sum / shape_index_count, 2.);
}

