//****************************************************************************
//****************************************************************************
//*****************************DISTANCES**************************************
//****************************************************************************
//****************************************************************************
double edge_length(int i){
    
    double ddx,ddy;
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    //basal
    int v1=e[i][1], v2=e[i][2];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
    double lb=sqrt(ddx*ddx+ddy*ddy);
    
    delete []dxdydz;
    
    return lb;
}
//****************************************************************************
void output_edge_lengths(char *fileName){
    
    FILE *file1;
    file1 = fopen(fileName, "wt");
    for(size_t i = 1; i<=Ne; ++i){
        if(e[i][0]!=0){
            fprintf(file1, "%g\n", edge_length(i));
        }
    }
    fclose(file1);
    
}
//****************************************************************************
double dist(int i, int j){
    double dx=v[i][1]-v[j][1];
    double dy=v[i][2]-v[j][2];
    return sqrt(dx*dx+dy*dy);
}
//****************************************************************************
double cell_dist(int i, int j){
    double dx=v_pass[c_cent[i]][1]-v_pass[c_cent[j]][1];
    double dy=v_pass[c_cent[i]][2]-v_pass[c_cent[j]][2];
    return sqrt(dx*dx+dy*dy);
}
//****************************************************************************
double dist_x(int i, int j){
    return v[j][1]-v[i][1];
}
//****************************************************************************
double dist_y(int i, int j){
    return v[j][2]-v[i][2];
}
//****************************************************************************
double cell_perimeter(int i){
    
    double _bperimeter=0.;
    
    for(size_t j = 3; j <= 2+basal_edges[i][2]; ++j){
        int edgID=abs(basal_edges[i][j]);
        _bperimeter += edge_length(edgID);
    }
    
    return _bperimeter;
}
//****************************************************************************

