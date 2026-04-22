//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

//****************************************************************************
//********************BASIC FUNCTIONS*****************************************
//****************************************************************************
double v_coord(int i, int coord){
   if (v[i][coord] < 0.) {return v[i][coord] + perioXYZ[coord - 1];}
   else if (v[i][coord] >= perioXYZ[coord - 1]) {return v[i][coord] - perioXYZ[coord - 1];}
   else return v[i][coord];
}
//****************************************************************************
int circle(int i, int n){
   if (i < 0.) {return i + n;}
   else if (i >= n) {return i - n;}
   else return i;
}
//****************************************************************************
void counterclockwise(int i0, int j0, int a, int& i1, int& j1){
   if (a == 0) {i1 = circle(i0    , N_grid); j1 = circle(j0 + 1, N_grid);}
   if (a == 1) {i1 = circle(i0 - 1, N_grid); j1 = circle(j0 + 1, N_grid);}
   if (a == 2) {i1 = circle(i0 - 1, N_grid); j1 = circle(j0    , N_grid);}
   if (a == 3) {i1 = circle(i0 - 1, N_grid); j1 = circle(j0 - 1, N_grid);}
   if (a == 4) {i1 = circle(i0    , N_grid); j1 = circle(j0 - 1, N_grid);}
   if (a == 5) {i1 = circle(i0 + 1, N_grid); j1 = circle(j0 - 1, N_grid);}
   if (a == 6) {i1 = circle(i0 + 1, N_grid); j1 = circle(j0    , N_grid);}
   if (a == 7) {i1 = circle(i0 + 1, N_grid); j1 = circle(j0 + 1, N_grid);}

   return;
}
//****************************************************************************
void get_vertex_indices(int i, int& x_coord, int& y_coord){
   double grid_step_x = perioXYZ[0] / ((double) N_grid);
   double grid_step_y = perioXYZ[1] / ((double) N_grid);

   x_coord = (int) (v_coord(i, 1) / grid_step_x);
   y_coord = (int) (v_coord(i, 2) / grid_step_y);

   return ;
}
//****************************************************************************
int get_sign(double x)
{
   if (x > 0) return 1;
   if (x < 0) return -1;
   return 0;
}
//****************************************************************************
//********************GRID FUNCTIONS******************************************
//****************************************************************************
void zero_grid(){
   for (int i = 0; i < N_grid; i ++ )
   {
      for (int j = 0; j < N_grid; j ++ )
      {
         g_Q_xx[i][j] = 0;
         g_Q_xy[i][j] = 0;
         g_beta[i][j] = 0;
         g_count[i][j] = 0;
         g_v_x[i][j] = 0;
         g_v_y[i][j] = 0;
         g_v_count[i][j] = 0;
      }
   }
   return;
}
//****************************************************************************
void sum_grid_Q_xx()
{
   int x_index, y_index;
   int not_symmetric;
   double Q_xx, Q_xy, Q_yy;

   //cout << perioXYZ[0] / ((double) N_grid) << endl;
   //cout << perioXYZ[1] / ((double) N_grid) << endl;

   vector <int> visited_grid_points = {};
      
   for (int i = 1; i <= Nc; i++) if (basal_facets[i][1] != 0)
   {
      not_symmetric = compute_Q_tensor(i, Q_xx, Q_xy, Q_yy);
      visited_grid_points = {};
      for (int j = 0; j < basal_facets[i][2]; j++)
      {
         int vert_in_grid = basal_vertices[i][3 + j];
         get_vertex_indices(vert_in_grid, x_index, y_index);
         //cout << v_coord(basal_vertices[i][3 + j], 1) << " " << v_coord(basal_vertices[i][3 + j], 2) << " " << x_index << " " << y_index << endl;

         int grid_linear = N_grid * x_index + y_index;

         //if (x_index == 0 && y_index == 0) cout << i << endl;



         if (std::count(visited_grid_points.begin(), visited_grid_points.end(), grid_linear) == 0) {
            visited_grid_points.push_back(grid_linear);
            g_Q_xx[x_index][y_index] += Q_xx;
            g_Q_xy[x_index][y_index] += Q_xy;
            g_count[x_index][y_index] += not_symmetric;
            //cout << "added to " << x_index << " " << y_index << " " << not_symmetric <<  endl;
         }

         g_v_x[x_index][y_index] += v_F[vert_in_grid][1];
         g_v_y[x_index][y_index] += v_F[vert_in_grid][2];
         g_v_count[x_index][y_index] += 1;

      }
   }

   //cout << "sum: " << g_Q_xx[0][0] << " " << g_count[0][0] << endl;

   return;
}
//****************************************************************************
void normalize_grid_Q_xx()
{
   for (int i = 0; i < N_grid; i++) 
   {
      for (int j = 0; j < N_grid; j++)
      {
         if (g_count[i][j] > 0.00000001) {
            g_Q_xx[i][j] = g_Q_xx[i][j]/g_count[i][j];
            g_Q_xy[i][j] = g_Q_xy[i][j]/g_count[i][j];
         } else 
         {
            g_Q_xx[i][j] = 0.;
            g_Q_xy[i][j] = 0.;
         }

         if (g_v_count[i][j] > 0.00000001) {
            g_v_x[i][j] /= g_v_count[i][j];
            g_v_y[i][j] /= g_v_count[i][j];
         } else 
         {
            g_v_x[i][j] = 0.;
            g_v_y[i][j] = 0.;
         }
      }
   }

   //cout << "sum: " << g_Q_xx[0][0] << endl;

   return;
}
//****************************************************************************
double beta_angle(int i, int j)
{
   double initial_beta = 0.5 * acos(2. * g_Q_xx[i][j]);
   if (get_sign(initial_beta) == get_sign(g_Q_xy[i][j])) {return initial_beta;} else {return -1 * initial_beta;}
}
//****************************************************************************
void compute_grid_beta()
{
   for (int i = 0; i < N_grid; i++) 
   {
      for (int j = 0; j < N_grid; j++)
      {
         g_beta[i][j] = beta_angle(i, j);
      }
   }

   //cout << "sum: " << g_beta[0][0] << endl;

   return;
}
//****************************************************************************
void fill_grid_arrays()
{
   zero_grid();
   sum_grid_Q_xx();
   normalize_grid_Q_xx();
   compute_grid_beta();

   return;
}

//****************************************************************************
//********************DEFECT DETECTION****************************************
//****************************************************************************
bool check_8_neighbor_directions(int i, int j)
{
   double out_count = 0.;
   int i1, j1;

   for (int a = 0; a < 8; a++)
   {
      counterclockwise(i, j, a, i1, j1);
      out_count += g_count[i][j];

   }

   if (out_count > 7.5) {return true;} else {return false;}
}
//****************************************************************************
bool check_4_neighbor_directions(int i, int j)
{
   double out_count = 0.;
   int i1, j1;

   out_count += g_count[circle(i    , N_grid)][circle(j    , N_grid)];
   out_count += g_count[circle(i    , N_grid)][circle(j - 1, N_grid)];
   out_count += g_count[circle(i + 1, N_grid)][circle(j - 1, N_grid)];
   out_count += g_count[circle(i + 1, N_grid)][circle(j    , N_grid)];

   if (out_count > 3.5) {return true;} else {return false;}
}
//****************************************************************************
double neighbor_angle(int i, int j, int a)
{
   int i1, j1;
   counterclockwise(i, j, a, i1, j1);

   return g_beta[i1][j1];
}
//****************************************************************************
double get_gamma(double beta_a, double beta_ap1)
{
   if (beta_ap1 - beta_a < -1 * PI / 2.) {return PI;}
   else if (beta_ap1 - beta_a > PI / 2.) {return -1 * PI;}
   else return 0;
}
//****************************************************************************
double winding_number(int i, int j)
{
   double beta_a, beta_ap1, gamma_a;
   double delta_beta = 0;

   for (int a = 0; a < 8; a++)
   {
      beta_a = neighbor_angle(i, j, a);
      beta_ap1 = neighbor_angle(i, j, circle(a + 1, 8));
      gamma_a = get_gamma(beta_a, beta_ap1);
      delta_beta += (beta_ap1 - beta_a + gamma_a);
   }

   return delta_beta / (2 * PI);
}
//****************************************************************************
double winding_number_grid_vertex(int i, int j, double& psi_direction)
{
   double beta_a, beta_ap1, gamma_a;
   double delta_beta    = 0;
   double theta_current = 0;
   double theta_sum    = 0;


      beta_a   = g_beta[circle(i    , N_grid)][circle(j    , N_grid)];
      beta_ap1 = g_beta[circle(i    , N_grid)][circle(j - 1, N_grid)];
      gamma_a  = get_gamma(beta_a, beta_ap1);
      delta_beta += (beta_ap1 - beta_a + gamma_a);

      theta_current = beta_a;
      theta_sum += theta_current;

      theta_current += delta_beta;
      theta_sum += theta_current;

      beta_a   = g_beta[circle(i    , N_grid)][circle(j - 1, N_grid)];
      beta_ap1 = g_beta[circle(i + 1, N_grid)][circle(j - 1, N_grid)];
      gamma_a  = get_gamma(beta_a, beta_ap1);
      delta_beta += (beta_ap1 - beta_a + gamma_a);

      theta_current += delta_beta;
      theta_sum += theta_current;

      beta_a   = g_beta[circle(i + 1, N_grid)][circle(j - 1, N_grid)];
      beta_ap1 = g_beta[circle(i + 1, N_grid)][circle(j    , N_grid)];
      gamma_a  = get_gamma(beta_a, beta_ap1);
      delta_beta += (beta_ap1 - beta_a + gamma_a);

      theta_current += delta_beta;
      theta_sum += theta_current;

      beta_a   = g_beta[circle(i + 1, N_grid)][circle(j    , N_grid)];
      beta_ap1 = g_beta[circle(i    , N_grid)][circle(j    , N_grid)];
      gamma_a  = get_gamma(beta_a, beta_ap1);
      delta_beta += (beta_ap1 - beta_a + gamma_a);

   psi_direction = 2 * theta_sum / 4. - PI;

   return delta_beta / (2 * PI);
}
//****************************************************************************
void grid_winding_number()
{  
   //N_grid = in_N_grid;
   if (N_grid > N_grid_max) cout << "grid too large" << endl;
   fill_grid_arrays();
   double defect_psi;

   for (int i = 0; i < N_grid; i++) 
   {
      for (int j = 0; j < N_grid; j++)
      {
         if (check_4_neighbor_directions(i, j)) {
            g_winding[i][j] = winding_number_grid_vertex(i, j, defect_psi);
            g_direction[i][j] = defect_psi;


         } else 
         {
            g_winding[i][j] = 0;
            g_direction[i][j] = 0;
         }
      }
   }

   return;
}
//****************************************************************************
void grid_winding_number_8()
{  
   //N_grid = in_N_grid;
   if (N_grid > N_grid_max) cout << "grid too large" << endl;
   fill_grid_arrays();
   double defect_psi;

   for (int i = 0; i < N_grid; i++) 
   {
      for (int j = 0; j < N_grid; j++)
      {
         if (check_8_neighbor_directions(i, j)) {
            g_winding[i][j] = winding_number(i, j);
            g_direction[i][j] = 0;//defect_psi;


         } else 
         {
            g_winding[i][j] = 0;
            g_direction[i][j] = 0;
         }
      }
   }

   return;
}
