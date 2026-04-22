//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double area_facet_gradient_contribution(int j, int cellID, double _Cellarea, double x0, double y0, int coord){
    
    int vert_ref_id=basal_vertices[cellID][3];
    double vrefx = v[vert_ref_id][1];
    double vrefy = v[vert_ref_id][2];
    
    
    //vertices
    int v1, v2;
    if(j==2+basal_vertices[cellID][2]){
        v1=basal_vertices[cellID][j];
        v2=basal_vertices[cellID][3];
    }
    else {
        v1=basal_vertices[cellID][j];
        v2=basal_vertices[cellID][j+1];
    }
    
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
    
    //v2
    double v2x=v[v2][1];
    double v2y=v[v2][2];
    //x
    if(fabs(v2x-vrefx)>perioXYZ[0]/2.){
        if(v2x<vrefx) v2x+=perioXYZ[0];
        else if(v2x>vrefx) v2x-=perioXYZ[0];
    }
    //y
    if(fabs(v2y-vrefy)>perioXYZ[1]/2.){
        if(v2y<vrefy) v2y+=perioXYZ[1];
        else if(v2y>vrefy) v2y-=perioXYZ[1];
    }

    //v0
    double v0x=x0;
    double v0y=y0;
    //x
    if(fabs(v0x-vrefx)>perioXYZ[0]/2.){
        if(v0x<vrefx) v0x+=perioXYZ[0];
        else if(v0x>vrefx) v0x-=perioXYZ[0];
    }
    //y
    if(fabs(v0y-vrefy)>perioXYZ[1]/2.){
        if(v0y<vrefy) v0y+=perioXYZ[1];
        else if(v0y>vrefy) v0y-=perioXYZ[1];
    }
    
    
    double AA0=A01;
    if(c_type[cellID]==2) AA0=A02;
    double prefactor = 2.*c_kAb*(_Cellarea-AA0);
    
    double Lx = perioXYZ[0];
    double Ly = perioXYZ[1];

    double s1x = (v1x - v0x) / Lx;
    double s2x = (v2x - v0x) / Lx;

    double s1y = (v1y - v0y) / Ly;
    double s2y = (v2y - v0y) / Ly;

    if (coord == 1) return 2.*c_kAb*(_Cellarea-AA0)*(-(Ly*s1y*s2x) + Ly*s1x*s2y);
    if (coord == 2) return 2.*c_kAb*(_Cellarea-AA0)*(-(Lx*s1y*s2x) + Lx*s1x*s2y);
    else {cout << "no coord" << endl; return 0;}
    
}

//****************************************************************************
double area_grad(int cellID, int coord){
    double _Cellarea=CellArea_new(cellID);
    double x0 = center_coord(cellID, 1);
    double y0 = center_coord(cellID, 2);
    
    double out_sum = 0;
    for(int j=3; j<=2+basal_vertices[cellID][2]; j++ ){
        out_sum += area_facet_gradient_contribution(j,cellID,_Cellarea, x0, y0, coord);
    }

    return out_sum;
}
//****************************************************************************
double perimeter_grad(int i, int coord){

   //TENSIONS
   double p_current = CellPerimeter(i);

   double out_sum = 0;

   for (int j = 1; j <= basal_edges[i][2]; j ++)
   {
      double ddx,ddy;
      double *dxdydz = new double[2]; dxdydz[0]=0; dxdydz[1]=0;

      int iEdg = abs(basal_edges[i][j + 2]);

      int v1=e[iEdg][1], v2=e[iEdg][2];
      torus_dx_dy_dz(dxdydz,v1,v2);
      ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
      ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);

      double Lx = perioXYZ[0];
      double Ly = perioXYZ[1];

      double dsx = ddx/Lx;
      double dsy = ddy/Ly;

      double lb=sqrt(ddx*ddx+ddy*ddy);
      delete []dxdydz;

      if(lb<0.00000001) lb=0.00000001;

      if (coord == 1) out_sum += 2 * kPer * (p_current - P0) * Lx * dsx * dsx / lb; 
      if (coord == 2) out_sum += 2 * kPer * (p_current - P0) * Ly * dsy * dsy / lb; 
      
   }
   return out_sum;
}
//****************************************************************************
double tension_grad(int i, int coord){

   //TENSIONS
    double ddx,ddy;
    double *dxdydz = new double[2]; dxdydz[0]=0; dxdydz[1]=0;
    
    int v1=e[i][1], v2=e[i][2];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
    double lb=sqrt(ddx*ddx+ddy*ddy);
    
    delete []dxdydz;

    double Lx = perioXYZ[0];
    double Ly = perioXYZ[1];

    double dsx = ddx/Lx;
    double dsy = ddy/Ly;

    
    
    double c0;
    if(lb<0.00000001) c0=0;
    else c0 = (e_g[i])/lb;
    
    if (coord == 1) return e_g[i] * Lx * dsx * dsx / lb;
    if (coord == 2) return e_g[i] * Ly * dsy * dsy / lb;
    else {cout << "no coord" << endl; return 0;}
    
}
//****************************************************************************
double total_grad(int coord){
    double sum_grad = 0;
    for (int i = 1; i <= Nc; i++) if (basal_facets[i][1] > 0)
    {
        sum_grad += perimeter_grad(i, coord);
        sum_grad += area_grad(i, coord);
    }
    for (int i = 1; i <= Ne; i++) if (e[i][0] > 0)
    {
        sum_grad += tension_grad(i, coord);
    }
    return sum_grad;
}
//****************************************************************************
void compute_box_gradients(double& store_sig_xx, double& store_sig_yy){
    store_sig_xx = total_grad(1);
    store_sig_yy = total_grad(2);
}
//****************************************************************************
void rescale_box_by_gradients(double alpha, double sig_xx, double sig_yy){

    double old_Lx = perioXYZ[0];
    double old_Ly = perioXYZ[1];

    double new_Lx = old_Lx - alpha * h * sig_xx; 
    double new_Ly = old_Ly - alpha * h * sig_yy;

    double eta_xx = new_Lx / old_Lx; 
    double eta_yy = new_Ly / old_Ly; 

    for (int i = 1; i <= Nv; i++) if (v[i][0] > 0.5)
    {
        v[i][1] *= eta_xx;
        v[i][2] *= eta_yy;
    }

   perioXYZ[0] *= eta_xx;
   perioXYZ[1] *= eta_yy;

   for (int i = 1; i <= Nv; i++) if (v[i][0] > 0.5)
   {
      torus_vertex(i);
   }
   calc_central_vertices();


}
//****************************************************************************
