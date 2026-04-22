//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

//****************************************************************************
//********************BASIC FUNCTIONS*****************************************
//****************************************************************************
void zero_vertex_neighbors(){
   for (int i = 1; i <= Nv; i++)
   {
      for (int j = 0; j<=4;j++) v_cell_neighbors[i][j] = 0;
   }
}
//****************************************************************************
void store_vertex_neighbors(){
   for (int i = 1; i <= Nc; i++) if (basal_edges[i][1] != 0)
   {
      for (int j = 1; j<=basal_edges[i][2];j++) 
      {
         int v1 = basal_vertices[i][j + 2];
         int current_found = v_cell_neighbors[v1][0];
         v_cell_neighbors[v1][current_found + 1] = i;
         v_cell_neighbors[v1][0] += 1;
      }
   }
}
//****************************************************************************
double dxy_vert_cell(int vert_id, int cell_id, int coord){
   double vert_x = v[vert_id][coord];
   double cell_x = center_coord(cell_id, coord);

   double dx = cell_x - vert_x;

    
   //x
   if(fabs(dx)>0.5*perioXYZ[coord-1]){
      if      (dx < 0.)   dx += perioXYZ[coord - 1];
      else if (dx > 0.)   dx -= perioXYZ[coord - 1];
   }

   return dx;
   
}
//****************************************************************************
double dxy_vert_vert(int v1_id, int v2_id, int coord){
   double v1_x = v[v1_id][coord];
   double v2_x = v[v2_id][coord];

   double dx = v2_x - v1_x;

    
   //x
   if(fabs(dx)>0.5*perioXYZ[coord-1]){
      if      (dx < 0.)   dx += perioXYZ[coord - 1];
      else if (dx > 0.)   dx -= perioXYZ[coord - 1];
   }

   return dx;
   
}
//****************************************************************************
double cell_angle(int cell_id, double& cell_nematic){
   double nem_x, nem_y;
   compute_director(cell_id, nem_x, nem_y);

   if (nonzero_vector(nem_x, nem_y)) 
   {
      cell_nematic = 1;
   }
   else
   {
      cell_nematic = 0;
   }

   return atan2(nem_y, nem_x);
}
//****************************************************************************
double get_vertex_defect_charge(int vert_id){
   vector<int> nearest = {};
   vector<int> next_nearest = {};
   vector<int> third_nearest = {};


   //**************FIND NEAREST NEIGHBORS*************************************
   for (int i = 1; i <= v_cell_neighbors[vert_id][0]; i++)
   {
      nearest.push_back(v_cell_neighbors[vert_id][i]);
   }

   //**************FIND NEXT NEAREST NEIGHBORS*************************************
   for (int i = 0; i < nearest.size(); i++ )
   {
      int cell_id = nearest[i];
      for (int j = 1; j <= basal_edges[cell_id][2]; j++)
      {
         int neighbor_cell_id = other_cell(cell_id, abs(basal_edges[cell_id][j + 2]));
         if (
            (count(nearest.begin(), nearest.end(), neighbor_cell_id) == 0) &&
            (count(next_nearest.begin(), next_nearest.end(), neighbor_cell_id) == 0)
            )
         {
            next_nearest.push_back(neighbor_cell_id);
         }
      }
   }

   //**************FIND THIRD NEAREST NEIGHBORS*************************************
   for (int i = 0; i < next_nearest.size(); i++ )
   {
      int cell_id = next_nearest[i];
      for (int j = 1; j <= basal_edges[cell_id][2]; j++)
      {
         int neighbor_cell_id = other_cell(cell_id, abs(basal_edges[cell_id][j + 2]));
         if (
            (count(nearest.begin(), nearest.end(), neighbor_cell_id) == 0) &&
            (count(next_nearest.begin(), next_nearest.end(), neighbor_cell_id) == 0)&&
            (count(third_nearest.begin(), third_nearest.end(), neighbor_cell_id) == 0)
            )
         {
            third_nearest.push_back(neighbor_cell_id);
         }
      }
   }

   //**************SORT CELLS*************************************
   vector<pair<double,int>> angle_id = {};

   double step_x;
   double step_y;
   double vert_cell_angle;

   for (int i = 0; i < third_nearest.size(); i++)
   {
      int cell_id = third_nearest[i];
      step_x = dxy_vert_cell(vert_id, cell_id, 1);
      step_y = dxy_vert_cell(vert_id, cell_id, 2);
      vert_cell_angle = atan2(step_y, step_x);
      angle_id.push_back(make_pair(vert_cell_angle, cell_id));
   }


   sort(angle_id.begin(), angle_id.end());


   vector<int> sorted_cells = {};
   for (int i=0;i < angle_id.size(); i++)
   {
      sorted_cells.push_back(angle_id[i].second);
   }
   //**************CALCULATE CHARGE*************************************
   int third_nearest_count = third_nearest.size();

   double beta_a, beta_ap1, gamma_a;
   double delta_beta    = 0;
   double theta_current = 0;
   double theta_sum    = 0;

   double sum_nematic_cells = 0;

   double cell_nematic_1 = 0;
   double cell_nematic_2 = 0;

   for (int i = 0; i < third_nearest_count; i++)
   {
      cell_nematic_1 = 0;
      cell_nematic_2 = 0;

      beta_a   = cell_angle(sorted_cells[circle(i,     third_nearest_count)], cell_nematic_1);
      beta_ap1 = cell_angle(sorted_cells[circle(i + 1, third_nearest_count)], cell_nematic_2);

      sum_nematic_cells += cell_nematic_1;

      gamma_a  = get_gamma(beta_a, beta_ap1);
      delta_beta += (beta_ap1 - beta_a + gamma_a);

      if (i == 0) {theta_current = beta_a;} else theta_current += delta_beta;
      theta_sum += theta_current;
   }

   if (fabs(sum_nematic_cells - (double) third_nearest_count) < 0.1) 
   {
      return delta_beta / (2 * PI);
   }
   else
   {
      return 0;
   }

}
//****************************************************************************
vector<int> get_third_order_vertices(int vert_id)
{
   vector<int> nearest = {};
   vector<int> next_nearest = {};
   vector<int> third_nearest = {};
   vector<int> vertices_out = {};


   //**************FIND NEAREST NEIGHBORS*************************************
   for (int i = 1; i <= v_cell_neighbors[vert_id][0]; i++)
   {
      nearest.push_back(v_cell_neighbors[vert_id][i]);
   }

   //**************FIND NEXT NEAREST NEIGHBORS*************************************
   for (int i = 0; i < nearest.size(); i++ )
   {
      int cell_id = nearest[i];
      for (int j = 1; j <= basal_edges[cell_id][2]; j++)
      {
         int neighbor_cell_id = other_cell(cell_id, abs(basal_edges[cell_id][j + 2]));
         if (
            (count(nearest.begin(), nearest.end(), neighbor_cell_id) == 0) &&
            (count(next_nearest.begin(), next_nearest.end(), neighbor_cell_id) == 0)
            )
         {
            next_nearest.push_back(neighbor_cell_id);
         }
      }
   }

   //**************FIND THIRD NEAREST NEIGHBORS*************************************
   for (int i = 0; i < next_nearest.size(); i++ )
   {
      int cell_id = next_nearest[i];
      for (int j = 1; j <= basal_edges[cell_id][2]; j++)
      {
         int neighbor_cell_id = other_cell(cell_id, abs(basal_edges[cell_id][j + 2]));
         if (
            (count(nearest.begin(), nearest.end(), neighbor_cell_id) == 0) &&
            (count(next_nearest.begin(), next_nearest.end(), neighbor_cell_id) == 0)&&
            (count(third_nearest.begin(), third_nearest.end(), neighbor_cell_id) == 0)
            )
         {
            third_nearest.push_back(neighbor_cell_id);
         }
      }
   }

   //**************add vertices*************************************
   //if (vert_id == 2) cout << "nearest: ";
   for (int i = 0; i < nearest.size(); i++)
   {
      int cell_id = nearest[i];
      //if (vert_id == 2) cout << cell_id << " (";
      for (int j = 1; j <= basal_edges[cell_id][2]; j++)
      {
         if (count(vertices_out.begin(), vertices_out.end(), basal_vertices[cell_id][j + 2]) == 0) 
            {
               vertices_out.push_back(basal_vertices[cell_id][j + 2]);
            }
      }
      
   }
   
   for (int i = 0; i < next_nearest.size(); i++)
   {
      int cell_id = next_nearest[i];
      for (int j = 1; j <= basal_edges[cell_id][2]; j++)
      {
         if (count(vertices_out.begin(), vertices_out.end(), basal_vertices[cell_id][j + 2]) == 0) 
            {
               vertices_out.push_back(basal_vertices[cell_id][j + 2]);
            }
      }
   }

   for (int i = 0; i < third_nearest.size(); i++)
   {
      int cell_id = third_nearest[i];
      for (int j = 1; j <= basal_edges[cell_id][2]; j++)
      {
         if (count(vertices_out.begin(), vertices_out.end(), basal_vertices[cell_id][j + 2]) == 0) 
            {
               vertices_out.push_back(basal_vertices[cell_id][j + 2]);
            }
      }
   }

   return vertices_out;

}
//****************************************************************************
void store_vertex_charges(){
   for (int i = 1; i <= Nv; i++) if (v[i][0] > 0.5)
   {
      v_charge[i] = get_vertex_defect_charge(i);
   }
}
//****************************************************************************
void zero_defect_arrays(){
   for (int i = 0; i < N_defects; i++) 
   {
      defect_charge[i] = 0;
      defect_x[i]      = 0;
      defect_y[i]      = 0;
   }
   N_defects = 0;
}
//****************************************************************************
void get_defects(){

   zero_vertex_neighbors();
   store_vertex_neighbors();
   store_vertex_charges();
   zero_defect_arrays();
   for (int i = 1; i <= Nv; i++) if (v[i][0] > 0.5)
   {
      if (fabs(v_charge[i]) > 0.001)
      {
         vector<int> neighboring_vertices = get_third_order_vertices(i);
         double store_charge = v_charge[i];
         double same_charge_count = 1;
         double sum_x = 0;
         double sum_y = 0;
         v_charge[i] = 0;
         for (int j = 0; j < neighboring_vertices.size(); j++)
         {
            int vert_now = neighboring_vertices[j];
            if (fabs(v_charge[vert_now] - store_charge) < 0.001)
            {
               sum_x += dxy_vert_vert(i, vert_now, 1);
               sum_y += dxy_vert_vert(i, vert_now, 2);
               same_charge_count += 1;
               v_charge[vert_now] = 0;
            }
         }

         defect_charge[N_defects] = store_charge;
         defect_x[N_defects] = v[i][1] + sum_x / same_charge_count;
         defect_y[N_defects] = v[i][2] + sum_y / same_charge_count;
         //cout << sum_x << " " << sum_y << " " << same_charge_count << endl;
         N_defects += 1;
      }
   }

   return;

}
//****************************************************************************
int count_defects(double target_charge){

   int count_val = 0;
   for (int i = 0; i < N_defects; i++) 
   {
      if (fabs(defect_charge[i] - target_charge) < 0.001) count_val += 1;
   }


   return count_val;

}
//****************************************************************************
int count_abs_larger_defects(double target_charge){

   int count_val = 0;
   for (int i = 0; i < N_defects; i++) 
   {
      if (fabs(defect_charge[i]) >  target_charge + 0.001) count_val += 1;
   }

   return count_val;

}

