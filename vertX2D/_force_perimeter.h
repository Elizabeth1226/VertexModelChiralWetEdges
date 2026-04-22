//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double c_perimeter_force(int i){

   //TENSIONS
   double p_current = CellPerimeter(i);
   double eg = 2*kPer*(p_current - P0);

   for (int j = 1; j <= basal_edges[i][2]; j ++)
   {
      double ddx,ddy;
      double *dxdydz = new double[2]; dxdydz[0]=0; dxdydz[1]=0;

      int iEdg = abs(basal_edges[i][j + 2]);

      int v1=e[iEdg][1], v2=e[iEdg][2];
      torus_dx_dy_dz(dxdydz,v1,v2);
      ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
      ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);

      double lb=sqrt(ddx*ddx+ddy*ddy);
      delete []dxdydz;

      if(lb<0.00000001) lb=0.00000001;

      double c0 = -eg/lb;

      v_F[v1][1] += -c0*ddx;
      v_F[v1][2] += -c0*ddy;

      v_F[v2][1] += c0*ddx;
      v_F[v2][2] += c0*ddy;
   }
   return kPer * pow(p_current - P0, 2.);
}
//****************************************************************************
