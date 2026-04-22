//****************************************************************************
//****************************************************************************
//*********************TOPOLOGICAL TRANSITIONS********************************
//****************************************************************************
//****************************************************************************
void shrink_edge(int v1, int v2){
    
    //CALCULATES CENTER POINT (NEW POSITION OF VERTEX v1)
    double x_cent, y_cent;
    //x
    if(fabs(dist_x(v1,v2))>perioXYZ[0]/2.) x_cent=0.5*(v[v1][1]+v[v2][1]+perioXYZ[0]);
    else x_cent=0.5*(v[v1][1]+v[v2][1]);
    //y
    if(fabs(dist_y(v1,v2))>perioXYZ[1]/2.) y_cent=0.5*(v[v1][2]+v[v2][2]+perioXYZ[1]);
    else y_cent=0.5*(v[v1][2]+v[v2][2]);
    
    //PUTS VERTEX v1 TO (x_cent,y_cent)
    v[v1][1]=x_cent;
    v[v1][2]=y_cent;
    torus_vertex(v1);
    
}
//****************************************************************************
void ReadFromLists_merge_vertices(int i, int v1, int v2, int c1, int c2, int *cntCells, int *edgesE1, int *edgesE2, int *cellsC1, int *cellsC2){
    
    //FIGURES OUT CELLS c1 & c2
    cellsC1[1]=c1;
    if(c2!=0){
        cellsC1[2]=c2;
        cntCells[1]=2;
    }
    else cntCells[1]=1;
    
    //FIGURES OUT CELLS
    //v1
    for(int j=3; j<=2+v_cells[v1][2]; j++){
        if( v_cells[v1][j]!=cellsC1[1] && v_cells[v1][j]!=cellsC1[2] ){
            cntCells[1]++;
            cellsC1[cntCells[1]]=v_cells[v1][j];
        }
    }
    
    //FIGURES OUT CELLS c1 & c2
    cellsC2[1]=c1;
    if(c2!=0){
        cellsC2[2]=c2;
        cntCells[2]=2;
    }
    else cntCells[2]=1;
    
    //FIGURES OUT CELLS
    //v2
    for(int j=3; j<=2+v_cells[v2][2]; j++){
        if( v_cells[v2][j]!=cellsC2[1] && v_cells[v2][j]!=cellsC2[2] ){
            cntCells[2]++;
            cellsC2[cntCells[2]]=v_cells[v2][j];
        }
    }
    
    //FIGURES OUT EDGES
    //v1
    int cntEdges=1;
    edgesE1[cntEdges]=i;
    for(int j=3; j<=2+v_edges[v1][2]; j++){
        if(v_edges[v1][j]!=i){
            cntEdges++;
            edgesE1[cntEdges]=v_edges[v1][j];
        }
    }
    //v2
    cntEdges=1;
    edgesE2[cntEdges]=i;
    for(int j=3; j<=2+v_edges[v2][2]; j++){
        if(v_edges[v2][j]!=i){
            cntEdges++;
            edgesE2[cntEdges]=v_edges[v2][j];
        }
    }
    
    /*for(int j=1; j<=v_edges[v1][2]; j++) printf("%d ", edgesE1[j]);
    printf("\n");
    for(int j=1; j<=v_edges[v2][2]; j++) printf("%d ", edgesE2[j]);
    printf("\n");
    
    for(int j=1; j<=cntCells[1]; j++) printf("%d ", cellsC1[j]);
    printf("\n");
    for(int j=1; j<=cntCells[2]; j++) printf("%d ", cellsC2[j]);
    printf("\n");*/
    
}
//****************************************************************************
int merge_vertices(int i){ //vertex v1 stays
    
    //FIGURES OUT VERTICES v1 & v2
    int v1=e[i][1], v2=e[i][2];
    
    //FIGURES OUT CELLS c1 & c2
    int c1=e_cell1[i], c2=e_cell2[i];
    if(c1==0 && c2!=0){
        c1=c2;
        c2=0;
    }
    
    int *cntEdges; cntEdges = new int[3];
    cntEdges[1]=v_edges[v1][2]; cntEdges[2]=v_edges[v2][2];
    
    //printf("coordNum1=%d  coordNum2=%d  v1=%d  v2=%d  c1=%d  c2=%d\n", cntEdges[1], cntEdges[2], v1, v2, c1, c2);
    
    //FINDS ELEMENTS ASSOCIATED WITH TRANSFORMATION
    int *edgesE1; edgesE1 = new int[15]; int *edgesE2; edgesE2 = new int[15];
    int *cellsC1; cellsC1 = new int[15]; int *cellsC2; cellsC2 = new int[15];
    int *cntCells; cntCells = new int[3];
    ReadFromLists_merge_vertices(i,v1,v2,c1,c2,cntCells,edgesE1,edgesE2,cellsC1,cellsC2);
    
    //printf("%d  %d\n", cntCells[1], cntCells[2]);
    
    //CHECKS IF VALENCE REDUCTION IS NEEDED AFTER merge_vertices
    if (cntEdges[1]==3 && cntEdges[2]==3) {v_T1dir[v1]=c1; v_T1dir2[v1]=c2;}
    
    //DISSOLVES CELLS
    dissolve_cell(c1);
    if(c2!=0){
        dissolve_cell(c2);
        for(int j=3; j<=cntCells[1]; j++) dissolve_cell(cellsC1[j]);
        for(int j=3; j<=cntCells[2]; j++) dissolve_cell(cellsC2[j]);
    }
    else {
        for(int j=2; j<=cntCells[1]; j++) dissolve_cell(cellsC1[j]);
        for(int j=2; j<=cntCells[2]; j++) dissolve_cell(cellsC2[j]);
    }
    
    //DISSOLVES EDGE i
    dissolve_edge(i);
    
    //SHRINKS EDGE
    shrink_edge(v1,v2);
    
    //DISSOLVES VERTEX v2
    dissolve_vertex(v2);
    
    //RESTITCHES EDGES
    int edg_id, vv1, vv2, eID;
    double eTENS;
    for(int j=2; j<=cntEdges[2]; j++){
        edg_id=edgesE2[j];
        if(e[edg_id][1]==v2){
            vv1=v1;
            vv2=e[edg_id][2];
        }
        else if(e[edg_id][2]==v2){
            vv1=e[edg_id][1];
            vv2=v1;
        }
        eTENS=e_g[edg_id];
        dissolve_edge(edg_id);
        eID=make_edge(vv1,vv2);
        e_g[eID]=eTENS;
    }
    
    //REMOVES EDGE i FROM basal_edges
    remove_edge_from_basal_edges(i,c1);
    if(c2!=0) remove_edge_from_basal_edges(i,c2);
    
    //REMAKES CELLS
    make_basal_side(c1);
    if(c2!=0){
        make_basal_side(c2);
        for(int j=3; j<=cntCells[1]; j++) make_basal_side(cellsC1[j]);
        for(int j=3; j<=cntCells[2]; j++) make_basal_side(cellsC2[j]);
    }
    else{
        for(int j=2; j<=cntCells[1]; j++) make_basal_side(cellsC1[j]);
        for(int j=2; j<=cntCells[2]; j++) make_basal_side(cellsC2[j]);
    }
    
    //DELETE ELEMENTS
    delete [] edgesE1; delete [] edgesE2; delete [] cellsC1; delete [] cellsC2; delete [] cntEdges; delete [] cntCells;
    
    //printf("%d\n", v_edges[v1][2]);
    //printf("%d\n", v_cells[v1][2]);
    
    return v1;
    
}
//****************************************************************************
void ReadFromLists_valence_reduction(int v1, int cres, int *_cntCells, int *cntAdjCells, int *edgesE, int *cellsC){
    
    //FIGURES OUT EDGES
    int cntEdges=0;
    for(int j=3; j<=2+v_edges[v1][2]; j++){
        if(e_cell1[v_edges[v1][j]]==cres || e_cell2[v_edges[v1][j]]==cres){
            cntEdges++;
            edgesE[cntEdges]=v_edges[v1][j];
        }
    }
    for(int j=3; j<=2+v_edges[v1][2]; j++){
        if(e_cell1[v_edges[v1][j]]!=cres && e_cell2[v_edges[v1][j]]!=cres){
            cntEdges++;
            edgesE[cntEdges]=v_edges[v1][j];
        }
    }
    
    //for(int j=1; j<=v_edges[v1][2]; j++) printf("%d ", edgesE[j]);
    //printf("\n");
    
    
    //FIGURES OUT CELLS
    int cntCells=0;
    cntCells++; cellsC[cntCells]=cres;
    for(int j=1; j<=cntEdges; j++){
        if(e_cell1[edgesE[j]]==cres && e_cell2[edgesE[j]]!=0){
            cntCells++; cellsC[cntCells]=e_cell2[edgesE[j]];
        }
        if(e_cell2[edgesE[j]]==cres && e_cell1[edgesE[j]]!=0){
            cntCells++; cellsC[cntCells]=e_cell1[edgesE[j]];
        }
    }
    cntAdjCells[1]=cntCells-1;
    
    //printf("%d\n", cntCells);
    
    if(cntCells==3){
        
        for(int j=3; j<=2+v_cells[v1][2]; j++){
            if( v_cells[v1][j]!=cellsC[1] && v_cells[v1][j]!=cellsC[2] && v_cells[v1][j]!=cellsC[3] ){
                cntCells++; cellsC[cntCells]=v_cells[v1][j];
            }
        }
        
    }
    
    else{
        for(int j=3; j<=2+v_cells[v1][2]; j++){
            if( v_cells[v1][j]!=cellsC[1] && v_cells[v1][j]!=cellsC[2]){
                cntCells++; cellsC[cntCells]=v_cells[v1][j];
            }
        }
    }
    
    _cntCells[1]=cntCells;
    
    //for(int j=1; j<=cntCells; j++) printf("%d ", cellsC[j]);
    //printf("\n");
    
}
//****************************************************************************
int valence_reduction(int v1, double fin_len, int cres1, int cres2){
    
    int nrosette=v_edges[v1][2];

    //FINDS ELEMENTS ASSOCIATED WITH TRANSFORMATION
    int *edgesE; edgesE = new int[nrosette+1];
    int *cellsC; cellsC = new int[nrosette+1];
    int *cntCells; cntCells = new int[2];
    int *cntAdjCells; cntAdjCells = new int[2];
    ReadFromLists_valence_reduction(v1,cres1,cntCells,cntAdjCells,edgesE,cellsC);
    
    //CALCULATES DIRECITONS OF CELL cres
    double *cdir; cdir=new double[4];
    for(int j=1; j<=2; j++){
        if(fabs(v_pass[c_cent[cres1]][j]-v_pass[c_cent[cres2]][j])>perioXYZ[j-1]/2.){
            if(v_pass[c_cent[cres1]][j]>v_pass[c_cent[cres2]][j]) cdir[j]=v_pass[c_cent[cres1]][j]-(v_pass[c_cent[cres2]][j]+perioXYZ[j-1]);
            if(v_pass[c_cent[cres1]][j]<v_pass[c_cent[cres2]][j]) cdir[j]=v_pass[c_cent[cres1]][j]-(v_pass[c_cent[cres2]][j]-perioXYZ[j-1]);
        }
        else cdir[j]=v_pass[c_cent[cres1]][j]-v_pass[c_cent[cres2]][j];
    }
    double norm=sqrt(cdir[1]*cdir[1]+cdir[2]*cdir[2]);
    
    //DISSOLVES CELLS
    for(int j=1; j<=cntCells[1]; j++) dissolve_cell(cellsC[j]);
    
    //MAKES VERTEX v2
    int v2=make_vertex(
                v[v1][1]+fin_len*0.5 * cdir[1]/norm,
                v[v1][2]+fin_len*0.5 * cdir[2]/norm
                );

    v[v1][1] -= fin_len*0.5 * cdir[1]/norm,
    v[v1][2] -= fin_len*0.5 * cdir[2]/norm,
    torus_vertex(v1);

    delete [] cdir;
    
    //RESTITCHES EDGES
    double eTENS;
    int eID;
    int vv1, vv2;
    for(int j=1; j<=2; j++){
        if(e[edgesE[j]][1]==v1){
            vv1=v2;
            vv2=e[edgesE[j]][2];
        }
        else{
            vv1=e[edgesE[j]][1];
            vv2=v2;
        }
        eTENS=e_g[edgesE[j]];
        dissolve_edge(edgesE[j]);
        eID=make_edge(vv1,vv2);
        e_g[eID]=eTENS;
    }
    
    //MAKES EDGE i
    int i=make_edge(v1,v2);
    
    //printf("%d\n", cntAdjCells[1]);
    
    //INSERTS EDGE i IN CELLS
    for(int m=2; m<=1+cntAdjCells[1]; m++){
        int e1=0, e2=0;
        for(int j=3; j<=2+basal_edges[cellsC[m]][2]; j++){
            for(int k=1; k<=nrosette; k++){
                if( edgesE[k]==abs(basal_edges[cellsC[m]][j]) ){
                    if(e1==0) e1=edgesE[k];
                    else e2=edgesE[k];
                }
            }
        }
        insert_edge_to_basal_edges(i,cellsC[m],e1,e2);
    }
    
    //REMAKES CELLS
    for(int j=1; j<=cntCells[1]; j++) make_basal_side(cellsC[j]);
    
    //DELETE ELEMENTS
    delete [] edgesE; delete [] cellsC; delete [] cntCells; delete [] cntAdjCells;
    
    return i;
}
//****************************************************************************
