//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double cell_hexatic(int cellID){
    
    int vert_ref_id=basal_vertices[cellID][3];
    int vert_id;
    double vrefx = v[vert_ref_id][1];
    double vrefy = v[vert_ref_id][2];
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    double xCenter = compute_basal_center_x(cellID);
    double yCenter = compute_basal_center_y(cellID);

    double xCoord = vrefx;
    double yCoord = vrefy;


    double rePart = 0;
    double imPart = 0;
    
    double theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
    rePart += cos(6 * theta_j);
    imPart += sin(6 * theta_j);

    //vert_id
    for(int j=4; j <= 2+basal_edges[cellID][2]; j++){
        vert_id=basal_vertices[cellID][j];
        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        xCoord = v[vert_id][1] + dxdydz[0];
        yCoord = v[vert_id][2] + dxdydz[1];
        theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
        rePart += cos(6 * theta_j);
        imPart += sin(6 * theta_j);
    }
    
    delete []dxdydz;
    return sqrt(rePart * rePart + imPart * imPart) / (double) basal_edges[cellID][2];
}
//****************************************************************************
double cell_hexatic_re(int cellID){
    
    int vert_ref_id=basal_vertices[cellID][3];
    int vert_id;
    double vrefx = v[vert_ref_id][1];
    double vrefy = v[vert_ref_id][2];
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    double xCenter = compute_basal_center_x(cellID);
    double yCenter = compute_basal_center_y(cellID);

    double xCoord = vrefx;
    double yCoord = vrefy;


    double rePart = 0;
    
    double theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
    rePart += cos(6 * theta_j);

    //vert_id
    for(int j=4; j <= 2+basal_edges[cellID][2]; j++){
        vert_id=basal_vertices[cellID][j];
        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        xCoord = v[vert_id][1] + dxdydz[0];
        yCoord = v[vert_id][2] + dxdydz[1];
        theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
        rePart += cos(6 * theta_j);
    }
    
    delete []dxdydz;
    return rePart / (double) basal_edges[cellID][2];
}
//****************************************************************************
double cell_hexatic_im(int cellID){
    
    int vert_ref_id=basal_vertices[cellID][3];
    int vert_id;
    double vrefx = v[vert_ref_id][1];
    double vrefy = v[vert_ref_id][2];
    
    double *dxdydz = new double[2];
    dxdydz[0]=0; dxdydz[1]=0;
    
    double xCenter = compute_basal_center_x(cellID);
    double yCenter = compute_basal_center_y(cellID);

    double xCoord = vrefx;
    double yCoord = vrefy;


    double imPart = 0;
    
    double theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
    imPart += sin(6 * theta_j);

    //vert_id
    for(int j=4; j <= 2+basal_edges[cellID][2]; j++){
        vert_id=basal_vertices[cellID][j];
        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        xCoord = v[vert_id][1] + dxdydz[0];
        yCoord = v[vert_id][2] + dxdydz[1];
        theta_j = atan2(yCoord - yCenter, xCoord - xCenter);
        imPart += sin(6 * theta_j);
    }
    
    delete []dxdydz;
    return imPart / (double) basal_edges[cellID][2];
}
//****************************************************************************
double pair_hexatic_correlation_re(int cellID_1, int cellID_2)
{
	double re_1 = cell_hexatic_re(cellID_1);
	double im_1 = cell_hexatic_im(cellID_1);
	
	double re_2 = cell_hexatic_re(cellID_2);
	double im_2 = cell_hexatic_im(cellID_2);

	return (re_1*re_2 + im_1*im_2);
}
//****************************************************************************
double pair_hexatic_correlation_im(int cellID_1, int cellID_2)
{
	double re_1 = cell_hexatic_re(cellID_1);
	double im_1 = cell_hexatic_im(cellID_1);
	
	double re_2 = cell_hexatic_re(cellID_2);
	double im_2 = cell_hexatic_im(cellID_2);
	return (re_1*im_2 - im_1*re_2);
}
//****************************************************************************
double pair_distance(int cellID_1, int cellID_2)
{
	double x_1 = compute_basal_center_x(cellID_1);
	double y_1 = compute_basal_center_y(cellID_1);
	
	double x_2 = compute_basal_center_x(cellID_2);
	double y_2 = compute_basal_center_y(cellID_2);
	return sqrt((x_1 - x_2) * (x_1 - x_2) + (y_1 - y_2) * (y_1 - y_2));
}
//****************************************************************************
void compute_spatial_correlations()
{
	double r_distance;
	int r_index;
	for (int i = 0; i < 80; i ++ )
	{
		d_pair_sum_re[i] = 0.;
        d_pair_sum_im[i] = 0.;
        d_pair_sum_count[i] = 0.;
	}
	for (int i = 1; i <= Nc; i++ )
	{
		for (int j = i + 1; j <= Nc; j ++)
		{
			r_distance = pair_distance(i, j);
			r_index = (int) (2 * r_distance);
            if (r_index < 80) {
			d_pair_sum_re[r_index] += pair_hexatic_correlation_re(i, j);
			d_pair_sum_im[r_index] += pair_hexatic_correlation_im(i, j);
			d_pair_sum_count[r_index] += 1;
            }
		}
	}
}

//****************************************************************************
double tissue_hexatic()
{
    double hexatic_sum = 0;
    double hexatic_count = 0;

    double current_hexatic;
    for (int i = 1; i <= Nc; i ++)
    {
        if (basal_edges[i][1] > 0 )
        {
            current_hexatic = cell_hexatic(i);
            hexatic_sum += current_hexatic;
            hexatic_count += 1;
        }
    }
    return hexatic_sum / hexatic_count;
}//****************************************************************************
double tissue_hexatic_susceptibility()
{
    double hexatic_sum = 0;
    double hexatic_square = 0;
    double hexatic_count = 0;

    double current_hexatic;
    for (int i = 1; i <= Nc; i ++)
    {
        if (basal_edges[i][1] > 0 )
        {
            current_hexatic = cell_hexatic(i);
            hexatic_sum += current_hexatic;
            hexatic_square += pow(current_hexatic, 2.);
            hexatic_count += 1;
        }
    }
    return (hexatic_square / hexatic_count) - pow(hexatic_sum / hexatic_count , 2.);
}
