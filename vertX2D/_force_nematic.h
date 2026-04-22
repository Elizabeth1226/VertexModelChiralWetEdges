//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
bool nonzero_vector(double vx, double vy)
{
   if (fabs(vx) > 1e-10 || fabs(vy) > 1e-10)
   {
      return true;
   }
   else return false;
}
//****************************************************************************
double center_coord(int i, int coord){
   int vert_id, vert_ref_id;
   double x_cent=0;
   double *dxdydz = new double[2];
   dxdydz[0]=0; dxdydz[1]=0;

   //vert_ref_id
   vert_ref_id=basal_vertices[i][3];
   x_cent+=v[vert_ref_id][coord];

   //vert_id
   for(int j=4; j <= 2+basal_edges[i][2]; j++){
      vert_id=basal_vertices[i][j];
      torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
      x_cent += v[vert_id][coord] + dxdydz[coord - 1];
   }

   //divided by number of edges
   x_cent/=(1.*basal_edges[i][2]);

   delete []dxdydz;

   return x_cent;
}
//****************************************************************************
double get_nemAB(int i, int coord_A, int coord_B, double cent_A, double cent_B)
{
   int vert_id, vert_ref_id;
   double x_cent=0, y_cent=0;
   double *dxdydz = new double[2];
   dxdydz[0]=0; dxdydz[1]=0;

   //vert_ref_id
   vert_ref_id=basal_vertices[i][3];

   double nemAB = (v[vert_ref_id][coord_A] - cent_A) * (v[vert_ref_id][coord_B] - cent_B);
   
   for(int j=4; j <= 2+basal_edges[i][2]; j++){
      vert_id=basal_vertices[i][j];
      torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
      nemAB += (v[vert_id][coord_A] + dxdydz[coord_A - 1] - cent_A) * (v[vert_id][coord_B] + dxdydz[coord_B - 1] - cent_B);
   }

   nemAB/=(1.*basal_edges[i][2]);

   delete []dxdydz;

   return nemAB;
}

//****************************************************************************
double compute_director(int i, double& nem_x, double& nem_y)
{
   double x_cent = center_coord(i, 1);
   double y_cent = center_coord(i, 2);

   double nemXX = get_nemAB(i, 1, 1, x_cent, x_cent); 
   double nemXY = get_nemAB(i, 1, 2, x_cent, y_cent); 
   double nemYY = get_nemAB(i, 2, 2, y_cent, y_cent); 

   double root_val = pow(nemXX,2.) + 4.*pow(nemXY,2.) - 2.*nemXX*nemYY + pow(nemYY,2.);
   if (root_val < 0.) root_val = 0.;

   double l_nem;

   //******************case 1: off diagonal term is small****************
   if (fabs(nemXY) < 1e-12)
   {
      if (fabs(nemXX-nemYY) < 1e-12)
      {
         nem_x = 0.;
         nem_y = 0.;
         l_nem = 0.;
      }
      else if (nemXX > nemYY)
      {
         nem_x = 1.;
         nem_y = 0.;
         l_nem = 1.;
      }
      else
      {
         nem_x = 0.;
         nem_y = 1.;
         l_nem = 1.;
      }
   }
   else
   {
      double eig1 = 0.5 * (nemXX + nemYY + sqrt(root_val));
      double eig2 = 0.5 * (nemXX + nemYY - sqrt(root_val));
      //******************case 2: the eigenvalues are the same****************
      if (fabs(eig1 - eig2) < 1e-12)
      {
         nem_x = 0.;
         nem_y = 0.;
         l_nem = 0.;
      }
      //******************case 3: general****************
      else
      {
         double temp_nem_x = -(-nemXX + nemYY - sqrt(root_val));
         double temp_nem_y = 2.*nemXY;

         l_nem = sqrt(temp_nem_x*temp_nem_x + temp_nem_y*temp_nem_y);
         nem_x = temp_nem_x/l_nem;
         nem_y = temp_nem_y/l_nem;
      }

   }

   return l_nem;
}
//****************************************************************************
int compute_Q_tensor(int i, double& Q_xx, double& Q_xy, double& Q_yy)
{
   double nem_x, nem_y;
   compute_director(i, nem_x, nem_y);
   int not_symmetric = 1;

   if (nonzero_vector(nem_x, nem_y))
   {
      Q_xx = 2. * (nem_x * nem_x - 0.5);
      Q_xy = 2. * (nem_x * nem_y);
      Q_yy = 2. * (nem_y * nem_y - 0.5);
   } else
   {
      Q_xx = 0.;
      Q_xy = 0.;
      Q_yy = 0.;
      not_symmetric = 1;
   }

   return not_symmetric;
}
bool opposite_sign(double val1, double val2)
{
   if (val1 * val2 < 0)
   {
      return true;
   }
   else return false;
}
//****************************************************************************
double compute_edge_direction(int v1, int v2, double& grad_phi_x, double& grad_phi_y)
{
   
   double ddx, ddy;
   double *dxdydz = new double[2];
   dxdydz[0]=0.; dxdydz[1]=0.;

   torus_dx_dy_dz(dxdydz,v1,v2);
   ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
   ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
   double lb=sqrt(ddx*ddx+ddy*ddy);

   delete []dxdydz;
      
   ddx/=lb;
   ddy/=lb;

   grad_phi_x = ddx;
   grad_phi_y = ddy;

   return lb;
}
//****************************************************************************
void get_Q(int i, double& Q_xx, double& Q_xy, double& Q_yy)
{

   double P_i = 0.;
   Q_xx = 0;
   Q_xy = 0;
   Q_yy = 0;

   for (int j = 1; j <= basal_edges[i][2]; j ++)
   {
      int iEdg = abs(basal_edges[i][j + 2]);
      int v1=e[iEdg][1], v2=e[iEdg][2];

      double t_x, t_y;
      double l_j = compute_edge_direction(v1, v2, t_x, t_y);
      //if (i == 100 && Time < h) cout << "{{" << v[v1][1] << "," << v[v1][2] << "},{" <<  v[v2][1] << "," << v[v2][2] << "}}," << endl;
      //if (i == 100 && Time < h) cout << l_j << endl;
      P_i += l_j;
      Q_xx += l_j * t_x * t_x;
      Q_xy += l_j * t_x * t_y;
      Q_yy += l_j * t_y * t_y;
      

   }

    //if (i == 100 && Time < h) cout <<"Qxx " << Q_xx << endl;

   if (P_i < 0.001) cout << "ERROR: P_i" << endl;
   
   Q_xx /= P_i;
   Q_xy /= P_i;
   Q_yy /= P_i;

   Q_xx -= 0.5;
   Q_yy -= 0.5;

   //if (i == 100 && Time < h) cout << Q_xx << " " << Q_xy << " " << Q_yy << endl;
}
//****************************************************************************
void compute_director_from_Q(int i, double& nem_x, double& nem_y)
{
   double Q_xx, Q_xy, Q_yy;
   get_Q(i, Q_xx, Q_xy, Q_yy);

   double qq = 2. * sqrt(Q_xx * Q_xx + Q_xy * Q_xy);

   double nem_theta = 0.5 * acos(2 * Q_xx / qq);

   if (opposite_sign(nem_theta, Q_xy)) nem_theta *= -1;
   nem_x = cos(nem_theta);
   nem_y = sin(nem_theta);
   
}

//****************************************************************************
double compute_edge_vector(int v1, int v2, double& grad_phi_x, double& grad_phi_y)
{
   
   double ddx, ddy;
   double *dxdydz = new double[2];
   dxdydz[0]=0.; dxdydz[1]=0.;

   torus_dx_dy_dz(dxdydz,v1,v2);
   ddx=v[v2][1]-(v[v1][1]+dxdydz[0]);
   ddy=v[v2][2]-(v[v1][2]+dxdydz[1]);
   double lb=sqrt(ddx*ddx+ddy*ddy);

   delete []dxdydz;

   grad_phi_x = ddx;
   grad_phi_y = ddy;

   return lb;
}
//****************************************************************************
double get_nematic_angle(double nem_x, double nem_y, double ee_x, double ee_y)
{
   double dot_prod = fabs(nem_x * ee_x + nem_y * ee_y);
   double abs_nem  = sqrt(nem_x * nem_x + nem_y * nem_y);
   double abs_ee   = sqrt(ee_x  * ee_x  + ee_y  * ee_y);


   double rescaled_prod = dot_prod / abs_nem / abs_ee;
   if (rescaled_prod > 1.)  {cout << "cos: " << rescaled_prod << endl; rescaled_prod = 1.;}
   if (rescaled_prod < -1.) {cout << "cos: " << rescaled_prod << endl; rescaled_prod = -1.;}
   return acos(rescaled_prod);
}
//****************************************************************************
void c_nematic_edge_tension(int i){

   double nem_x, nem_y;
   compute_director(i, nem_x, nem_y);

   if (nonzero_vector(nem_x, nem_y)){

      for (int j = 1; j <= basal_edges[i][2]; j ++)
      {
         int iEdg = abs(basal_edges[i][j + 2]);
         int v1=e[iEdg][1], v2=e[iEdg][2];

         double ee_x, ee_y;
         compute_edge_direction(v1, v2, ee_x, ee_y);

         if (nonzero_vector(ee_x, ee_y))
         {
            double edge_director_angle = get_nematic_angle(nem_x, nem_y, ee_x, ee_y);
            e_target[iEdg] += -1 * (zeta_inter / 2.) * cos(2 * edge_director_angle);
         }

      }
   }
}
//****************************************************************************
double cell_cell_angle(int i1, int i2, int& zero_count)
{
   zero_count = 0;

   double n1x, n1y, n2x, n2y;
   compute_director_from_Q(i1, n1x, n1y);
   compute_director_from_Q(i2, n2x, n2y);

   double dot_prod = fabs(n1x * n2x + n1y * n2y);
   double abs_n1   = sqrt(n1x * n1x + n1y * n1y);
   double abs_n2   = sqrt(n2x * n2x + n2y * n2y);

   if (abs_n1 < 1e-10) zero_count += 1;
   if (abs_n2 < 1e-10) zero_count += 1;


   double rescaled_prod = dot_prod / abs_n1 / abs_n2;
   if (rescaled_prod > 1.)  {cout << "cos: " << rescaled_prod << endl; rescaled_prod = 1.;}
   if (rescaled_prod < -1.) {cout << "cos: " << rescaled_prod << endl; rescaled_prod = -1.;}
   return acos(rescaled_prod);
}
//****************************************************************************
double cell_edge_angle(int i1, int i2)
{

   double n1x, n1y;
   compute_director_from_Q(i1, n1x, n1y);

   double n2x, n2y;
   compute_edge_direction(e[i2][1], e[i2][2], n2x, n2y);

   double dot_prod = fabs(n1x * n2x + n1y * n2y);
   double abs_n1   = sqrt(n1x * n1x + n1y * n1y);
   double abs_n2   = sqrt(n2x * n2x + n2y * n2y);

   double rescaled_prod = dot_prod / abs_n1 / abs_n2;
   if (rescaled_prod > 1.)  {cout << "cos: " << rescaled_prod << endl; rescaled_prod = 1.;}
   if (rescaled_prod < -1.) {cout << "cos: " << rescaled_prod << endl; rescaled_prod = -1.;}
   return acos(rescaled_prod);
}
//****************************************************************************
int other_cell(int cell_id, int edg_id)
{
   if (e_cell1[edg_id] == cell_id) return e_cell2[edg_id];
   if (e_cell2[edg_id] == cell_id) return e_cell1[edg_id];
   else {cout << "no cell error" << endl; return 0;}
}
//****************************************************************************
double nematic_order_parameter()
{
   double S_sum = 0;
   double S_count = 0;
   int isometric_cell = 0;
   for(int i=1; i<=Nc; i++){
      if(basal_edges[i][1]!=0){
         for (int j = 1; j<=basal_edges[i][2]; j++)
         {
             int other_cell_val = other_cell(i, abs(basal_edges[i][j + 2]));
             isometric_cell = 0;
             double angle_val = cell_cell_angle(i, other_cell_val, isometric_cell);
             if (isometric_cell == 0) {
                if (angle_val > PI/ 2.) angle_val = PI - angle_val;
                S_sum += cos(2 * angle_val);
                S_count += 1.;   
            }     
         }
      }
   }
   return S_sum/S_count;
}
//****************************************************************************
int edge_circle(int i, int n){
   if (i < 1.) {return i + n;}
   else if (i > n) {return i - n;}
   else return i;
}

//****************************************************************************
double cross_z(double ax, double ay, int coord){
   if (coord == 1) {return ay;}
   else if (coord == 2) {return -1 * ax;}
   else {cout << "ERROR: wrong coord!" << endl; return 0;}
}
//****************************************************************************
double A_dot_b(double Axx, double Axy, double Ayy, double bx, double by, int coord){
   if (coord == 1) {return Axx * bx + Axy * by;}
   else if (coord == 2) {return Axy * bx + Ayy * by;}
   else {cout << "ERROR: wrong coord!" << endl; return 0;}
}
//****************************************************************************
void c_nematic_stress_force(int i){

   double Q_xx, Q_xy, Q_yy;
   get_Q(i, Q_xx, Q_xy, Q_yy);

   double zeta_c = zeta_1;
   if (c_activity[i]==2) zeta_c = zeta_2;
   //if (i == 1) cout << Q_xx << endl;
   if (true){

      for (int j = 1; j <= basal_vertices[i][2]; j ++)
      {
         int vi      = basal_vertices[i][2 + j];

         int j_min1 = edge_circle(j - 1, basal_vertices[i][2]);
         int j_plus1 = edge_circle(j + 1, basal_vertices[i][2]);
         
         int v_min1  = basal_vertices[i][2 + j_min1];
         int v_plus1 = basal_vertices[i][2 + j_plus1];

         double ee_x, ee_y;
         compute_edge_vector(v_min1, v_plus1, ee_x, ee_y);

         double Ri_x = 0.5 * cross_z(ee_x, ee_y, 1);
         double Ri_y = 0.5 * cross_z(ee_x, ee_y, 2);
         //if (i == 100) cout << "R: " << Ri_x << " " << Ri_y << endl;


         double nem_Fx = zeta_c * A_dot_b(Q_xx, Q_xy, Q_yy, Ri_x, Ri_y, 1);
         double nem_Fy = zeta_c * A_dot_b(Q_xx, Q_xy, Q_yy, Ri_x, Ri_y, 2);

         //if (i == 100) cout << "F: " << nem_Fx << " " << nem_Fy << endl;

         v_F[vi][1] += nem_Fx;
         v_F[vi][2] += nem_Fy;

      }
   }
}



