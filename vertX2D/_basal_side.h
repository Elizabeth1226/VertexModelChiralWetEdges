//****************************************************************************
//****************************************************************************
//***************************BASAL SIDE***************************************
//****************************************************************************
//****************************************************************************
int basal_center(int i, int make_or_not){
    
    int vert_id, vert_ref_id, pass_id;
    double x_cent=0, y_cent=0;
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    //vert_ref_id
    vert_ref_id=basal_vertices[i][3];
    x_cent+=v[vert_ref_id][1];
    y_cent+=v[vert_ref_id][2];
    
    //vert_id
    for(int j=4; j <= 2+basal_edges[i][2]; j++){
        vert_id=basal_vertices[i][j];
        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        x_cent += v[vert_id][1] + dxdydz[0];
        y_cent += v[vert_id][2] + dxdydz[1];
    }

    //divided by number of edges
    x_cent/=(1.*basal_edges[i][2]);
    y_cent/=(1.*basal_edges[i][2]);
    
    //MAKE PASSIVE CENTER VERTEX
    if(make_or_not==1){
        pass_id=make_vertex_pass(x_cent,y_cent);
        torus_vertex_pass(pass_id);
        v_pass_cell[pass_id]=i;
        c_cent[i]=pass_id;
    }
    //PLACE PASSIVE CENTER VERTEX
    else{
        pass_id=c_cent[i];
        v_pass[pass_id][1]=x_cent;
        v_pass[pass_id][2]=y_cent;
        torus_vertex_pass(pass_id);
    }
    
    delete []dxdydz;
    
    return pass_id;
    
}
//****************************************************************************
void set_e_cell(int edge_id, int cell_id){
    
    if      (e_cell1[edge_id]==0 && e_cell2[edge_id]!=cell_id) e_cell1[edge_id]=cell_id;
    else if (e_cell2[edge_id]==0 && e_cell1[edge_id]!=cell_id) e_cell2[edge_id]=cell_id;
    
}
//****************************************************************************
void set_v_cell(int vert_id, int i){
    
    int flg=0;
    for(int j=3; j<=2+v_cells[vert_id][2]; j++){
        if(v_cells[vert_id][j]==i){
            flg++;
            break;
        }
    }
    
    if(flg==0) {
        v_cells[vert_id][2]++;
        v_cells[vert_id][2+v_cells[vert_id][2]]=i;
    }
    
}
//****************************************************************************
void make_basal_vertices(int i){
    
    int vert_id;
    
    basal_vertices[i][1]=basal_edges[i][1];
    basal_vertices[i][2]=basal_edges[i][2];
    
    for(int j=3; j <= 2+basal_edges[i][2]; j++){
        
        if(basal_edges[i][j]>0) basal_vertices[i][j]=e[basal_edges[i][j]][1];
        if(basal_edges[i][j]<0) basal_vertices[i][j]=e[-basal_edges[i][j]][2];
        
    }
}
//****************************************************************************
void make_basal_side(int i){
    
    make_basal_vertices(i);
    
    //CALCULATE CELL CENTER
    int pass_id=basal_center(i,1);
    
    //MAKE TRIANGULAR FACETS
    basal_facets[i][1]=basal_edges[i][1];
    basal_facets[i][2]=basal_edges[i][2];
    
    for(int j=3; j < 2+basal_edges[i][2]; j++){
        basal_facets[i][j]=make_facet(
                                      basal_vertices[i][j],
                                      basal_vertices[i][j+1],
                                      pass_id
                                      );
        f_cell[basal_facets[i][j]]=i;//f_cell
        set_e_cell(abs(basal_edges[i][j]),i);//e_cell
    }
    int j = 2+basal_vertices[i][2];
    basal_facets[i][j]=make_facet(
                                  basal_vertices[i][j],
                                  basal_vertices[i][3],
                                  pass_id
                                  );
    f_cell[basal_facets[i][j]]=i;//f_cell
    set_e_cell(abs(basal_edges[i][j]),i);//e_cell
    
    //v_cell
    for(int j=3; j <= 2+basal_edges[i][2]; j++) set_v_cell(basal_vertices[i][j],i);
    
}
//****************************************************************************
void make_basal_side_ALL(){
    for(int i=1; i<=Nc; i++) make_basal_side(i);
}
//****************************************************************************










