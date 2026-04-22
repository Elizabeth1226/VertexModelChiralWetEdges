//****************************************************************************
//****************************************************************************
//***************************DISSOLVE*****************************************
//****************************************************************************
//****************************************************************************
void dissolve_vertex(int i){
    
    //printf("vertex %d dissolved\n", i);
    
    //VERTEX
    v[i][0]=0;
    v[i][1]=0;
    v[i][2]=0;
    
    //ATTRIBUTES
    v_F[i][1]=0; v_F[i][2]=0;
    v_cells[i][1]=0; v_cells[i][2]=0;
    v_edges[i][1]=0; v_edges[i][2]=0;
    v_T1dir[i]=0; v_T1dir2[i]=0; v_vertT1[i]=0; v_edgeT1[i]=0;
    
    v_freeId.add(i);
}
//****************************************************************************
void dissolve_vertex_pass(int i){
    
    //printf("passive vertex %d dissolved\n", i);
    
    //PASSIVE VERTEX
    v_pass[i][0]=0;
    v_pass[i][1]=0;
    v_pass[i][2]=0;
    
    //ATTRIBUTES OF OTHER ELEMENTS
    if(v_pass_cell[i]!=0) c_cent[v_pass_cell[i]]=0;
    
    //ATTRIBUTES
    v_pass_cell[i]=0;
    
    v_pass_freeId.add(i);
}
//****************************************************************************
void dissolve_edge(int i){
    
    //printf("edge %d dissolved\n", i);
    
    //ATTRIBUTES OF OTHER ELEMENTS
    for(int j=1; j<=2; j++){
        int vert_id=e[i][j];
        for(int k=3; k<=2+v_edges[vert_id][2]; k++){
            if(v_edges[vert_id][k]==i){
                v_edges[vert_id][k]=v_edges[vert_id][2+v_edges[vert_id][2]];
                v_edges[vert_id][2]--;
            }
        }
    }
    
    //EDGE
    e[i][0]=0;
    e[i][1]=0;
    e[i][2]=0;
    
    //ATTRIBUTES
    e_cell1[i]=0; e_cell2[i]=0;
    e_length[i]=0;
    e_dldt[i]=0;
    e_g[i]=0;
    
    e_freeId.add(i);
}
//****************************************************************************
void dissolve_facet(int i){
    
    //printf("facet %d dissolved\n", i);
    
    //FACET
    f[i][0]=0;
    f[i][1]=0;
    f[i][2]=0;
    f[i][3]=0;
    
    //ATTRIBUTES
    f_cell[i]=0;
    
    f_freeId.add(i);
}
//****************************************************************************
void dissolve_cell(int i){
    
    int pass_vert_basal=f[basal_facets[i][3]][3];
    
    //DISSOLVE FACETS
    for(size_t j=3; j<=2+basal_facets[i][2]; j++) dissolve_facet(basal_facets[i][j]);
    
    //DELETES FACETS FROM LISTS
    for(size_t j=0; j<=14; j++) basal_facets[i][j]=0;
    
    //PASSIVE VERTICES
    dissolve_vertex_pass(pass_vert_basal);
    
    int vert_id, edge_id;
    
    //DELETES LIST basal_vertices ( is recreated by make_basal_side() )
    for(size_t j=3; j<=2+basal_edges[i][2]; j++){
        
        //v_cell
        vert_id=basal_vertices[i][j];
        for(int k=3; k<=2+v_cells[vert_id][2]; k++){
            if(v_cells[vert_id][k]==i){
                v_cells[vert_id][k]=v_cells[vert_id][2+v_cells[vert_id][2]];
                v_cells[vert_id][2]--;
            }
        }
        
        //basal_vertices
        basal_vertices[i][j]=0;
        
        //e_cell
        edge_id=abs(basal_edges[i][j]);
        if(e_cell1[edge_id]==i) e_cell1[edge_id]=0;
        if(e_cell2[edge_id]==i) e_cell2[edge_id]=0;
    }
    basal_vertices[i][1]=0;
    basal_vertices[i][2]=0;
    
    //ATTRIBUTES
    c_cent[i]=0;
}
//****************************************************************************
void delete_cell_from_basal_list(int i){
    
    //printf("cell %d dissolved\n", i);
    
    for(size_t j=0; j<=14; j++) basal_edges[i][j]=0;
    
    c_freeId.add(i);
}
//****************************************************************************
