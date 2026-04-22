//****************************************************************************
//****************************************************************************
//******************************CELL EXTRUSION********************************
//****************************************************************************
//****************************************************************************
void ReadFromLists_T2(int i, int coordNum, int *edgesE, int *cellsC, int *edgesEE, int *verticesV){

    //IDENTIFIES EDGES e
    for(size_t j=1; j<=coordNum; j++) edgesE[j]=abs(basal_edges[i][2+j]);
    
    //IDENTIFIES VERTICES v
    for(size_t j=1; j<=coordNum; j++) verticesV[j]=basal_vertices[i][2+j];
    
    //IDENTIFIES CELLS c
    for(size_t j=1; j<=coordNum; j++) {
        if(e_cell1[edgesE[j]]!=i) cellsC[j]=e_cell1[edgesE[j]];
        else cellsC[j]=e_cell2[edgesE[j]];
    }
    
    //IDENTIFIES EDGES ee
    int edg_id, cell_id1, cell_id2, vert_id;
    for(size_t j=1; j<=coordNum; j++) {
        vert_id=verticesV[j];
        if(j==1) cell_id1=cellsC[coordNum];
        else cell_id1=cellsC[j-1];
        cell_id2=cellsC[j];
        for(int k=3; k<=2+basal_edges[cell_id2][2]; k++){
            edg_id=abs(basal_edges[cell_id2][k]);
            if( edg_id!=edgesE[j]){
                if (e[edg_id][1]==vert_id || e[edg_id][2]==vert_id){
                    if (e_cell1[edg_id]==cell_id1 && e_cell2[edg_id]==cell_id2) edgesEE[j]=edg_id;
                    else if (e_cell1[edg_id]==cell_id2 && e_cell2[edg_id]==cell_id1) edgesEE[j]=edg_id;
                }
            }
        }
    }
    
}
//****************************************************************************
int make_rosette(int i){
    
    int coordNum=basal_edges[i][2];
    
    printf("rosette created by extrusion of cell %d\n", i);
    
    //FINDS ELEMENTS ASSOCIATED WITH TRANSFORMATION
    int *edgesE; edgesE = new int[coordNum+1]; int *cellsC; cellsC = new int[coordNum+1]; int *edgesEE; edgesEE = new int[coordNum+1]; int *verticesV; verticesV = new int[coordNum+1];
    ReadFromLists_T2(i,coordNum,edgesE,cellsC,edgesEE,verticesV);
    
    //READS CENTER OF EXTRUDING CELL
    double x_cent=v_pass[c_cent[i]][1], y_cent=v_pass[c_cent[i]][2], z_cent=v_pass[c_cent[i]][3];
    
    //DISSOLVES CELLS i, c1, c2 & c3
    dissolve_cell(i);
    for(size_t j=1; j<=coordNum; j++) dissolve_cell(cellsC[j]);
    
    //DISSOLVES EDGES e1,e2 & e3
    for(size_t j=1; j<=coordNum; j++) dissolve_edge(edgesE[j]);
    
    //DISSOLVES VERTICES v2 & v3 AND THEIR APICAL PARTNERS
    for(size_t j=2; j<=coordNum; j++) dissolve_vertex(verticesV[j]);
    
    //MOVES VERTEX v1 TO THE CENTER
    v[verticesV[1]][1]=x_cent; v[verticesV[1]][2]=y_cent; v[verticesV[1]][3]=z_cent;
    
    //RESTITCHES EDGES ee
    int v1, v2;
    for(size_t j=1; j<=coordNum; j++) {
        if(e[edgesEE[j]][1]==verticesV[j]){
            v1=verticesV[1];
            v2=e[edgesEE[j]][2];
        }
        else{
            v1=e[edgesEE[j]][1];
            v2=verticesV[1];
        }
        dissolve_edge(edgesEE[j]);
        make_edge(v1,v2);
    }
    
    //REMOVES EDGES e1, e2 & e3 FROM basal_edges
    for(size_t j=1; j<=coordNum; j++) remove_edge_from_basal_edges(edgesE[j],cellsC[j]);
    
    //DELETES CELL i FROM LIST basal_edges
    delete_cell_from_basal_list(i);
    
    //REMAKES CELLS
    for(size_t j=1; j<=coordNum; j++) make_basal_side(cellsC[j]);
    
    //DELETE ELEMENTS
    int remVert=verticesV[1];
    delete [] edgesE; delete [] cellsC; delete [] edgesEE; delete [] verticesV;
    
    return remVert;
}
//****************************************************************************
