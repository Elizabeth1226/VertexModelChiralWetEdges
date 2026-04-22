//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double compute_basal_center_x(int i){
    
    int vert_id, vert_ref_id, pass_id;
    double x_cent=0, y_cent=0;
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    //vert_ref_id
    vert_ref_id=basal_vertices[i][3];
    x_cent+=v[vert_ref_id][1];
    
    //vert_id
    for(int j=4; j <= 2+basal_edges[i][2]; j++){
        vert_id=basal_vertices[i][j];
        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        x_cent += v[vert_id][1] + dxdydz[0];
    }

    //divided by number of edges
    x_cent/=(1.*basal_edges[i][2]);
    
    
    
    delete []dxdydz;
    
    return x_cent;
}
//****************************************************************************
double compute_basal_center_y(int i){
    
    int vert_id, vert_ref_id, pass_id;
    double x_cent=0, y_cent=0;
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    //vert_ref_id
    vert_ref_id=basal_vertices[i][3];
    y_cent+=v[vert_ref_id][2];
    
    //vert_id
    for(int j=4; j <= 2+basal_edges[i][2]; j++){
        vert_id=basal_vertices[i][j];
        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        y_cent += v[vert_id][2] + dxdydz[1];
    }

    //divided by number of edges
    y_cent/=(1.*basal_edges[i][2]);
    
    
    
    delete []dxdydz;
    
    return y_cent;
}
//****************************************************************************
std::vector<int> find_adjacent_cells(int cell_i)
{
    std::vector<int> out_list = {};
    int edge_i;
    for (int i = 1; i <= basal_edges[cell_i][2]; i ++)
    {
        edge_i = abs(basal_edges[cell_i][2 + i]);
        //std::cout << edge_i <<std::endl;
        if (e_cell1[edge_i] != cell_i) {out_list.push_back(e_cell1[edge_i]);}
        else if (e_cell2[edge_i] != cell_i) {out_list.push_back(e_cell2[edge_i]);}
    }
    return out_list;
}
//****************************************************************************
double defect_density()
{
    double cell_count = 0.;
    double defect_count = 0.;
    for (int i = 1; i <= Nc; i ++ )
    {
        if (basal_vertices[i][1] != 0)
        {
            cell_count += 1.;
            if (basal_vertices[i][2] != 6) defect_count += 1.;
        }
    }
    return defect_count / cell_count;
}
//****************************************************************************
bool is_dislocation(int i)
{
    bool found_dislocation = true;
    if (basal_facets[i][2] == 5) //if the cell is 5-sided
    {
        std::vector<int> adjacent_cells = find_adjacent_cells(i);
        std::vector<int> adjacent_relevant = {};
        for (int j = 0; j < adjacent_cells.size(); j ++ )
        {
            if (basal_facets[adjacent_cells[j]][2] != 7 && basal_facets[adjacent_cells[j]][2] != 6) //if it has a neighbor that is neither 6- nor 7-sided
            {
                found_dislocation = false;
            }
            else if (basal_facets[adjacent_cells[j]][2] == 7) //if the cell is exactly 7-sided
            {
                adjacent_relevant.push_back(adjacent_cells[j]);
            }

            if (adjacent_relevant.size() == 1) //if we found exactly one 7-sided
            {
                int other_i = adjacent_relevant[0];
                std::vector<int> neighbor_adjacent = find_adjacent_cells(other_i);
                for (int j2 = 0; j2 < neighbor_adjacent.size(); j2 ++ )
                {
                    if (neighbor_adjacent[j2] != i) //if we are not looking at the original cell
                    {
                        if (basal_facets[neighbor_adjacent[j2]][2] != 6) found_dislocation = false;
                    }
                }
            }
            else
            {
                found_dislocation = false;
            }
        }                
    }
    else if (basal_facets[i][2] == 7) //if the cell is 5-sided
    {
        std::vector<int> adjacent_cells = find_adjacent_cells(i);
        std::vector<int> adjacent_relevant = {};
        for (int j = 0; j < adjacent_cells.size(); j ++ )
        {
            if (basal_facets[adjacent_cells[j]][2] != 5 && basal_facets[adjacent_cells[j]][2] != 6) //if it has a neighbor that is neither 6- nor 5-sided
            {
                found_dislocation = false;
            }
            else if (basal_facets[adjacent_cells[j]][2] == 5) //if the cell is exactly 5-sided
            {
                adjacent_relevant.push_back(adjacent_cells[j]);
            }

            if (adjacent_relevant.size() == 1) //if we found exactly one 5-sided
            {
                int other_i = adjacent_relevant[0];
                std::vector<int> neighbor_adjacent = find_adjacent_cells(other_i);
                for (int j2 = 0; j2 < neighbor_adjacent.size(); j2 ++ )
                {
                    if (neighbor_adjacent[j2] != i) //if we are not looking at the original cell
                    {
                        if (basal_facets[neighbor_adjacent[j2]][2] != 6) found_dislocation = false;
                    }
                }
            }
            else
            {
                found_dislocation = false;
            }
        }                
    }

    else
    {
        found_dislocation = false;
    }
    return found_dislocation;
}
//****************************************************************************
bool is_disclination(int i)
{
    bool found_disclination = true;
    if (basal_facets[i][2] == 5 || basal_facets[i][2] == 7)
    {
        std::vector<int> adjacent_cells = find_adjacent_cells(i);
        for (int j = 0; j < adjacent_cells.size(); j ++ )
        {
            if (basal_facets[adjacent_cells[j]][2] != 6) found_disclination = false;
        }                
    }
    else
    {
        found_disclination = false;
    }
    return found_disclination;
}
//****************************************************************************
double dislocation_density()
{
    double cell_count = 0.;
    double defect_count = 0.;
    for (int i = 1; i <= Nc; i ++ )
    {
        if (basal_vertices[i][1] != 0)
        {
            cell_count += 1.;
            if (is_dislocation(i) == true) defect_count += 1.;
        }
    }
    return defect_count / cell_count;
}
//****************************************************************************
double disclination_density()
{
    double cell_count = 0.;
    double defect_count = 0.;
    for (int i = 1; i <= Nc; i ++ )
    {
        if (basal_vertices[i][1] != 0)
        {
            cell_count += 1.;
            if (is_disclination(i) == true) defect_count += 1.;
        }
    }
    return defect_count / cell_count;
}
//****************************************************************************
double cell_cHexatic_re(int cellID){
    
    int adj_cell_id;
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    double xCenter = compute_basal_center_x(cellID); //center of the observed cell (x)
    double yCenter = compute_basal_center_y(cellID); //center of the observed cell (y)

    std::vector<int> adjacent_cells = find_adjacent_cells(cellID);

    double xCoord;
    double yCoord;
    double theta_j;

    double rePart = 0;
    double imPart = 0;

    //vert_id
    for(int j=0; j < adjacent_cells.size(); j++){
        adj_cell_id=adjacent_cells[j];
        xCoord = compute_basal_center_x(adj_cell_id);
        yCoord = compute_basal_center_y(adj_cell_id);
        
    
        dxdydz[0]=0;
        dxdydz[1]=0;

        //x
        if(fabs(xCoord - xCenter) > 0.5*perioXYZ[0]){
            if      (xCoord < xCenter)   dxdydz[0]=perioXYZ[0];
            else if (xCoord > xCenter)   dxdydz[0]=-perioXYZ[0];
        }
        else dxdydz[0]=0;
        //y
        if(fabs(yCoord - yCenter) > 0.5*perioXYZ[1]){
            if      (yCoord < yCenter)   dxdydz[1]=perioXYZ[1];
            else if (yCoord > yCenter)   dxdydz[1]=-perioXYZ[1];
        }
        else dxdydz[1]=0;

        xCoord += dxdydz[0];
        yCoord += dxdydz[1];

        theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
        rePart += cos(6 * theta_j);
        imPart += sin(6 * theta_j);
    }
    
    delete []dxdydz;
    return rePart / (double) basal_edges[cellID][2];
}
//****************************************************************************
double cell_cHexatic_im(int cellID){
    
    int adj_cell_id;
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    double xCenter = compute_basal_center_x(cellID); //center of the observed cell (x)
    double yCenter = compute_basal_center_y(cellID); //center of the observed cell (y)

    std::vector<int> adjacent_cells = find_adjacent_cells(cellID);

    double xCoord;
    double yCoord;
    double theta_j;

    double rePart = 0;
    double imPart = 0;

    //vert_id
    for(int j=0; j < adjacent_cells.size(); j++){
        adj_cell_id=adjacent_cells[j];
        xCoord = compute_basal_center_x(adj_cell_id);
        yCoord = compute_basal_center_y(adj_cell_id);
        
    
        dxdydz[0]=0;
        dxdydz[1]=0;

        //x
        if(fabs(xCoord - xCenter) > 0.5*perioXYZ[0]){
            if      (xCoord < xCenter)   dxdydz[0]=perioXYZ[0];
            else if (xCoord > xCenter)   dxdydz[0]=-perioXYZ[0];
        }
        else dxdydz[0]=0;
        //y
        if(fabs(yCoord - yCenter) > 0.5*perioXYZ[1]){
            if      (yCoord < yCenter)   dxdydz[1]=perioXYZ[1];
            else if (yCoord > yCenter)   dxdydz[1]=-perioXYZ[1];
        }
        else dxdydz[1]=0;

        xCoord += dxdydz[0];
        yCoord += dxdydz[1];
    
        theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
        rePart += cos(6 * theta_j);
        imPart += sin(6 * theta_j);
    }
    
    delete []dxdydz;
    return imPart / (double) basal_edges[cellID][2];
}
//****************************************************************************
double tissue_cHexatic()
{
    double t_hexatic_sum_re = 0;
    double t_hexatic_sum_im = 0;
    double t_hexatic_count = 0;

    double t_current_hexatic_re;
    double t_current_hexatic_im;
    for (int i = 1; i <= Nc; i ++)
    {
        if (basal_edges[i][1] > 0 )
        {
            t_current_hexatic_re = cell_cHexatic_re(i);
            t_current_hexatic_im = cell_cHexatic_im(i);
            t_hexatic_sum_re += t_current_hexatic_re;
            t_hexatic_sum_im += t_current_hexatic_im;
            t_hexatic_count += 1;
        }
    }
    return sqrt(t_hexatic_sum_re * t_hexatic_sum_re + t_hexatic_sum_im * t_hexatic_sum_im) / t_hexatic_count;
}

//****************************************************************************
double tissue_cHexatic_susceptibility()
{
    double t_hexatic_sum_re = 0;
    double t_hexatic_sum_im = 0;
    double t_hexatic_square = 0;
    double t_hexatic_count = 0;

    double t_current_hexatic_re;
    double t_current_hexatic_im;
    for (int i = 1; i <= Nc; i ++)
    {
        if (basal_edges[i][1] > 0 )
        {
            t_current_hexatic_re = cell_cHexatic_re(i);
            t_current_hexatic_im = cell_cHexatic_im(i);
            t_hexatic_sum_re += t_current_hexatic_re;
            t_hexatic_sum_im += t_current_hexatic_im;
            t_hexatic_square += t_current_hexatic_re * t_current_hexatic_re + t_current_hexatic_im * t_current_hexatic_im;
            t_hexatic_count += 1;
        }
    }
    return (t_hexatic_square / t_hexatic_count) - pow(sqrt(t_hexatic_sum_re * t_hexatic_sum_re + t_hexatic_sum_im * t_hexatic_sum_im) / t_hexatic_count , 2.);
}
