//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double EdgeLen(int i){
    

    double ddx,ddy;
    double *dxdydz = new double[2]; dxdydz[0]=0; dxdydz[1]=0;
    
    int v1=e[i][1], v2=e[i][2];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
    double lb=sqrt(ddx*ddx+ddy*ddy);
    
    delete []dxdydz;
    
    
    
    return lb;
}
//****************************************************************************
double CellPerimeter(int i){
    

   double total_per_len = 0;
   for (int j = 3; j <= basal_edges[i][2] + 2; j ++)
   {
      total_per_len += EdgeLen(abs(basal_edges[i][j]));
   }
   return total_per_len;
            
}
//****************************************************************************
double e_BasalLength_force(int i){
    

    double ddx,ddy;
    double *dxdydz = new double[2]; dxdydz[0]=0; dxdydz[1]=0;
    
    int v1=e[i][1], v2=e[i][2];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
    double lb=sqrt(ddx*ddx+ddy*ddy);
    
    delete []dxdydz;
    
    
    double c0;
    if(lb<0.00000001) c0=0;
    else c0 = -(e_g[i])/lb;
    
    v_F[v1][1] += -c0*ddx;
    v_F[v1][2] += -c0*ddy;
        
    v_F[v2][1] += c0*ddx;
    v_F[v2][2] += c0*ddy;
    
    return (e_g[i])*lb;
}
