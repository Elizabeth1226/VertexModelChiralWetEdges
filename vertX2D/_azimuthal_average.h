//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
#include <cmath>

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

void get_chiral_cell_center(){
    double x_sum, y_sum;
    int chiral_cell_count = 0;

    x_sum = 0;
    y_sum = 0;

    for (int i=1; i<=Nc; i++){
        if (c_activity[i]==1){
            x_sum += center_coord(i,1);
            y_sum += center_coord(i,2);
            chiral_cell_count++;
        }
    }

    chiral_cell_center[0] = x_sum/chiral_cell_count;

    chiral_cell_center[1]= y_sum/chiral_cell_count;

}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void make_bins(){

    double grid_size;
    
    grid_size = 0.5*perioXYZ[1]/sqrt(Nc);

    int bin_num =0;
    for (double pos = 0; pos <= perioXYZ[1]; pos += grid_size){
        bin_positions[bin_num]=pos;
        bin_num++;
    }
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void make_fine_bins(){

    double grid_size;
    
    grid_size = 0.1*perioXYZ[1]/sqrt(Nc);

    int bin_num =0;
    for (double pos = 0; pos <= perioXYZ[1]; pos += grid_size){
        fine_bin_positions[bin_num]=pos;
        bin_num++;
    }
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void cell_position_check(){
//    if cellID is within grid_size of radius r, return 1; else return 0.
    double xc_val, yc_val, grid_size, rc_val;
    
    grid_size = 0.5*perioXYZ[1]/sqrt(Nc);
    for (int cellID=1; cellID<=Nc; cellID++){
        for (int bin_num=0; bin_num<=2.*sqrt(Nc); bin_num++){
            xc_val = center_coord(cellID,1)-chiral_cell_center[0];
            yc_val = center_coord(cellID,2)-chiral_cell_center[1];
            rc_val = sqrt(xc_val*xc_val+yc_val*yc_val);

            if ( abs(rc_val-bin_positions[bin_num]) <= grid_size ){
                c_bin_categories[cellID]=bin_num;
            }

        }
    }
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void edge_position_check(){
//    if cellID is within grid_size of radius r, return 1; else return 0.
    double xc_val, yc_val, grid_size, rc_val, central_coord_edge_x, central_coord_edge_y;
    int v1, v2;
    
    
    grid_size = 0.1*perioXYZ[1]/sqrt(Nc);
    for (int edgeID=1; edgeID<=Ne; edgeID++){
        v1=e[edgeID][1];
        v2=e[edgeID][2];
        central_coord_edge_x = (v[v2][1]+v[v1][1])/2;
        central_coord_edge_y = (v[v2][2]+v[v1][2])/2;
        xc_val = central_coord_edge_x-chiral_cell_center[0];
        yc_val = central_coord_edge_y-chiral_cell_center[1];
        rc_val = sqrt(xc_val*xc_val+yc_val*yc_val);
        for (int fine_bin_num=0; fine_bin_num<=10.*sqrt(Nc); fine_bin_num++){
            if ( abs(rc_val-fine_bin_positions[fine_bin_num]) <= grid_size ){
                e_bin_categories[edgeID]=fine_bin_num;
                azimuthal_edge_length [fine_bin_num] += EdgeLen(edgeID);
                azimuthal_edge_count [fine_bin_num] ++;
            }

        }
    }
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

double get_normal_rx(int cellID){
//  return the x component of unit radial vector at center of cell, radius r
    double xc_val, yc_val, rc_val, nx;

    xc_val = center_coord(cellID,1)-chiral_cell_center[0];
    yc_val = center_coord(cellID,2)-chiral_cell_center[1];
    rc_val = sqrt(xc_val*xc_val+yc_val*yc_val);
    nx = xc_val/rc_val;
    return nx;

}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

double get_normal_ry(int cellID){
//  return the y component of unit radial vector at center of cell, radius r
    double xc_val, yc_val, rc_val, ny;

    xc_val = center_coord(cellID,1)-chiral_cell_center[0];
    yc_val = center_coord(cellID,2)-chiral_cell_center[1];
    rc_val = sqrt(xc_val*xc_val+yc_val*yc_val);
    ny = yc_val/rc_val;
    return ny;
}


//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

double get_nematic_angle (int cellID){
// return value of phi of nematic tensor of cellID, phi has value [-0.5pi, 0.5pi]
    double Qxx, Qyy, Qxy, Q;

    get_Q(cellID, Qxx, Qxy, Qyy);
    Q = sqrt(-Qxx*Qyy+Qxy*Qxy);
    Qxx = Qxx/Q;
    Qxy = Qxy/Q;

    double nem_theta = 0.5 * acos(Qxx);

   if (opposite_sign(nem_theta, Qxy)) nem_theta *= -1;

    // if (Qxx<0){
    //     Qxx = -Qxx;
    //     Qxy = -Qxy;
    //     phi = 0.5*asin(Qxy);
    // } else {
    //     phi = 0.5*asin(Qxy);
    // }

    return nem_theta;
}


//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

double get_rotation_angle (int cellID){
// return the phi angle in cylindrical coordinate of cellID, has value [0, 2pi]
    double phi, rx, ry, pi;

    pi = 4*atan(1);
    rx = get_normal_rx(cellID);
    ry = get_normal_ry(cellID);

    if ((rx<0)&&(ry<0) ){
        phi = pi+acos(-rx);
    }else if ((rx>0)&&(ry<0)){
        phi = 2*pi+asin(ry);
    } else{
        phi = acos(rx);
    }

    return phi;
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

double rotated_nematic_tensor_Qxx (int cellID){
// rotate the nematic tensor by angle phi_rot
    double Qxx_rot, phi, phi_rot;
    double Q, Qxx, Qyy, Qxy;

    phi = get_nematic_angle(cellID);
    phi_rot = get_rotation_angle(cellID);
    phi -= phi_rot;
    get_Q(cellID, Qxx, Qxy, Qyy);
    Q = sqrt(-Qxx*Qyy+Qxy*Qxy);
    Qxx_rot = Q*cos(2*phi);

    return Qxx_rot;
}



//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double rotated_nematic_tensor_Qxy (int cellID){
// rotate the nematic tensor by angle phi_rot
    double Qxx_rot, phi, phi_rot, Qxy_rot;
    double Q, Qxx, Qyy, Qxy;

    phi = get_nematic_angle(cellID);
    phi_rot = get_rotation_angle(cellID);
    phi -= phi_rot;
    get_Q(cellID, Qxx, Qxy, Qyy);
    Q = sqrt(-Qxx*Qyy+Qxy*Qxy);
    Qxy_rot = Q*sin(2*phi);

    return Qxy_rot;
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double rotated_nematic_tensor_Qyy (int cellID){
// rotate the nematic tensor by angle phi_rot
    double Qxx_rot, phi, phi_rot, Qyy_rot;
    double Q, Qxx, Qyy, Qxy;

    phi = get_nematic_angle(cellID);
    phi_rot = get_rotation_angle(cellID);
    phi -= phi_rot;
    get_Q(cellID, Qxx, Qxy, Qyy);
    Q = sqrt(-Qxx*Qyy+Qxy*Qxy);
    Qyy_rot = -Q*cos(2*phi);

    return Qyy_rot;
}




//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double vx_cell_average (int cellID){
    double vx_sum;
    int num_vertex;
    int vertex_id;
    
    vx_sum = 0;
    num_vertex = basal_vertices[cellID][2];
    
    for (int i = 3; i <= 2+num_vertex; i++){
        vertex_id = basal_vertices[cellID][i];
        vx_sum += v_F[vertex_id][1];
    }
    return vx_sum/num_vertex;
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************


double vy_cell_average (int cellID){
    double vy_sum;
    int num_vertex;
    int vertex_id;
    
    num_vertex = basal_vertices[cellID][2];
    vy_sum = 0;
    
    for (int i = 3; i <= 2+num_vertex; i++){
        vertex_id= basal_vertices[cellID][i];
        vy_sum += v_F[vertex_id][2];
    }
    return vy_sum/num_vertex;
}

void compute_cell_force(int i, double& vel_x, double& vel_y);
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_averaged_v(){
    
    double vx_current, vy_current, v_current, v_sum;
    int cell_count;
    
    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        v_sum = 0;
        for (int j = 1; j <= Nc; j++){
            if (c_bin_categories[j] == bin_num){
                cell_count ++;
                compute_cell_force(j, vx_current, vy_current);
                v_current = sqrt(vx_current*vx_current+vy_current*vy_current);
                v_sum += v_current;
            }
        }
        if (cell_count == 0){
            v_azimuthal_ave[bin_num] = -1.;
        }else {
            v_azimuthal_ave[bin_num] = v_sum/cell_count;
        }
    }
    
    
    
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_averaged_vr(){
    
    double vr_sum, vr_current, vy_current, vx_current, rx, ry;
    int cell_count;
    
    
    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        vr_sum = 0;
        for (int j = 1; j <= Nc; j++){
            if (c_bin_categories[j] == bin_num){
                cell_count ++;
                compute_cell_force(j, vx_current, vy_current);
                rx = get_normal_rx(j);
                ry = get_normal_ry(j);
                vr_current = vx_current*rx+vy_current*ry;
                vr_sum += vr_current;
            }
        } if (cell_count == 0){
            vr_azimuthal_ave[bin_num] = 99999.;
        } else{
            vr_azimuthal_ave[bin_num]= vr_sum/cell_count;
        }
    }
    
}


//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_averaged_vphi(){
    
    double vphi_sum, vx_current, vy_current, vphi_current, rx, ry;
    int cell_count;
    
    
    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        vphi_sum = 0;
        for (int j = 1; j <= Nc; j++){
            if (c_bin_categories[j] == bin_num){
                cell_count ++;
                compute_cell_force(j, vx_current, vy_current);
                rx = get_normal_rx(j);
                ry = get_normal_ry(j);
                vphi_current = -ry*vx_current+rx*vy_current;
                vphi_sum += vphi_current;
            }
        } if (cell_count == 0){
            vphi_azimuthal_ave[bin_num] = 99999.;
        } else{
            vphi_azimuthal_ave[bin_num] = vphi_sum/cell_count;
        }
    }
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
    
// void azimuthal_averaged_Q(){

//     double Q, Q_sum, Q_current;
//     double Qxx, Qxy, Qyy;
//     int cell_count;

//     for (int bin_num =0; bin_num<=sqrt(Nc); bin_num++){
//         cell_count = 0;
//         Q_sum = 0;
//         for (int j = 1; j <= Nc; j++){
//             result = cell_position_check(j,i);
//             if (result == 1){
//                 cell_count ++;
//                 get_Q(j, Qxx, Qxy, Qyy);
//                 Q_current = sqrt(-Qxx*Qyy+Qxy*Qxy);
//                 Q_sum += Q_current;
//             }
//         }
//         if (cell_count == 0){
//             return 99999.;
//         } else{
//             Q = Q_sum/cell_count;
//             return Q;
//         }
//     }
    
// }

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
    
void azimuthal_averaged_Q_type1(){

    double Q, Q_sum, Q_current;
    double Qxx, Qxy, Qyy;
    int cell_count;

    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        Q_sum = 0;
        for (int j = 1; j <= Nc; j++){
            if ((c_bin_categories[j] == bin_num)&&(c_activity[j] == 1)){
                cell_count ++;
                get_Q(j, Qxx, Qxy, Qyy);
                Q_current = sqrt(-Qxx*Qyy+Qxy*Qxy);
                Q_sum += Q_current;
            }
        }
        if (cell_count == 0){
            Q_type1_azimuthal_ave[bin_num] = -1.;
        } else{
            Q_type1_azimuthal_ave[bin_num] = Q_sum/cell_count;
        }
    }
    
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
    
void azimuthal_averaged_Q_type2(){

    double Q, Q_sum, Q_current;
    double Qxx, Qxy, Qyy;
    int cell_count;

    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        Q_sum = 0;
        for (int j = 1; j <= Nc; j++){
            if ((c_bin_categories[j] == bin_num)&&(c_activity[j] == 2)){
                cell_count ++;
                get_Q(j, Qxx, Qxy, Qyy);
                Q_current = sqrt(-Qxx*Qyy+Qxy*Qxy);
                Q_sum += Q_current;
            }
        }
        if (cell_count == 0){
            Q_type2_azimuthal_ave[bin_num] = -1.;
        } else{
            Q_type2_azimuthal_ave[bin_num] = Q_sum/cell_count;
        }

    }
    
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************

// double azimuthal_averaged_Qphi(double i){
//     double Q, Qphi, grid_size, N_grid, Qxx_sum, Qxy_sum, Qyy_sum;
//     double Qxx, Qxy, Qyy;
//     int result, cell_count;


//         cell_count = 0;
//         Qxx_sum = 0;
//         Qxy_sum = 0;
//         Qyy_sum = 0;
//         for (int j = 1; j <= Nc; j++){
//             result = cell_position_check(j,i);
//             if (result == 1){
//                 cell_count ++;
//                 Qxx_sum += rotated_nematic_tensor_Qxx(j);
//                 Qxy_sum += rotated_nematic_tensor_Qxy(j);
//                 Qyy_sum += rotated_nematic_tensor_Qyy(j);
//             }
//         }
//         if (cell_count ==0){
//             return 99999.;
//         } else{
//             Q = sqrt(-Qxx_sum*Qyy_sum+Qxy_sum*Qxy_sum)/cell_count;
//             Qxx = Qxx_sum/cell_count/Q;
//             Qxy = Qxy_sum/cell_count/Q;

//             if (Qxx<0){
//                 Qxy = -Qxy;
//                 Qphi = 0.5*asin(Qxy);
//             } else {
//                 Qphi = 0.5*asin(Qxy);
//             }

//             return Qphi; 
//         }  
    
// }

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_averaged_Qphi_type1(){
    double Q, Qphi, Qxx_sum, Qxy_sum, Qyy_sum;
    double Qxx, Qxy, Qyy;
    int cell_count;
    
    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        Qxx_sum = 0;
        Qxy_sum = 0;
        Qyy_sum = 0;
        for (int j = 1; j <= Nc; j++){
            if ((c_bin_categories[j] == bin_num)&&(c_activity[j] == 1)){
                cell_count ++;
                Qxx_sum += rotated_nematic_tensor_Qxx(j);
                Qxy_sum += rotated_nematic_tensor_Qxy(j);
                Qyy_sum += rotated_nematic_tensor_Qyy(j);
            }
        }
        if (cell_count ==0){
            Qphi_type1_azimuthal_ave[bin_num] = 99999.;
        } else{
            Q = sqrt(-Qxx_sum*Qyy_sum+Qxy_sum*Qxy_sum)/cell_count;
            Qxx = Qxx_sum/cell_count/Q;
            Qxy = Qxy_sum/cell_count/Q;

            if (Qxx<0){
                Qxy = -Qxy;
                Qphi = 0.5*asin(Qxy);
            } else {
                Qphi = 0.5*asin(Qxy);
            }

            Qphi_type1_azimuthal_ave[bin_num] = Qphi; 
        }  
    }
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_averaged_Qphi_type2(){
    double Q, Qphi, Qxx_sum, Qxy_sum, Qyy_sum;
    double Qxx, Qxy, Qyy;
    int cell_count;
    
    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        Qxx_sum = 0;
        Qxy_sum = 0;
        Qyy_sum = 0;
        for (int j = 1; j <= Nc; j++){
            if ((c_bin_categories[j] == bin_num)&&(c_activity[j] == 2)){
                cell_count ++;
                Qxx_sum += rotated_nematic_tensor_Qxx(j);
                Qxy_sum += rotated_nematic_tensor_Qxy(j);
                Qyy_sum += rotated_nematic_tensor_Qyy(j);
            }
        }
        if (cell_count ==0){
            Qphi_type2_azimuthal_ave[bin_num] = 99999.;
        } else{
            Q = sqrt(-Qxx_sum*Qyy_sum+Qxy_sum*Qxy_sum)/cell_count;
            Qxx = Qxx_sum/cell_count/Q;
            Qxy = Qxy_sum/cell_count/Q;

            if (Qxx<0){
                Qxy = -Qxy;
                Qphi = 0.5*asin(Qxy);
            } else {
                Qphi = 0.5*asin(Qxy);
            }

            Qphi_type2_azimuthal_ave[bin_num] = Qphi; 
        }  
    }
}
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_averaged_cell_shape_index_type1(){
    double P_current, A_current, Kp_current, Kp_sum;
    int cell_count;

    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        Kp_sum = 0;
        for (int j=1; j<=Nc; j++){
            if ((c_bin_categories[j] == bin_num)&&(c_activity[j] == 1)){
                cell_count ++;
                A_current = CellArea_new(j);
                P_current = CellPerimeter(j);
                Kp_current = abs(P_current/sqrt(A_current));
                Kp_sum += Kp_current;
            }
        }
        if (cell_count == 0){
            Kp_type1_azimuthal_ave[bin_num] = -1.;
        }else{
            Kp_type1_azimuthal_ave[bin_num] = Kp_sum/cell_count;
        }
    }
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_averaged_cell_shape_index_type2(){
    double P_current, A_current, Kp_current, Kp_sum;
    int cell_count;

    for (int bin_num =0; bin_num<=2.*sqrt(Nc); bin_num++){
        cell_count = 0;
        Kp_sum = 0;
        for (int j=1; j<=Nc; j++){
            if ((c_bin_categories[j] == bin_num)&&(c_activity[j] == 2)){
                cell_count ++;
                A_current = CellArea_new(j);
                P_current = CellPerimeter(j);
                Kp_current = P_current/sqrt(A_current);
                Kp_sum += Kp_current;
            }
        }
        if (cell_count == 0){
            Kp_type2_azimuthal_ave[bin_num] = -1.;
        }else{
            Kp_type2_azimuthal_ave[bin_num] = Kp_sum/cell_count;
        }
    }
}

//********************************************************************************
void azimuthal_averaged_edge_length(){

for (int fine_bin_num=0; fine_bin_num<=10.*sqrt(Nc); fine_bin_num++){
    azimuthal_edge_length[fine_bin_num] = azimuthal_edge_length[fine_bin_num]/azimuthal_edge_count[fine_bin_num];
}

}

// //****************************************************************************
// void compute_director_from_Q_test(int i, double& nem_x, double& nem_y)
// {
//    double Q_xx, Q_xy, Q_yy;
//    get_Q(i, Q_xx, Q_xy, Q_yy);

//    double qq = 2. * sqrt(Q_xx * Q_xx + Q_xy * Q_xy);

//    double nem_theta = get_nematic_angle(i);

//    nem_x = cos(nem_theta);
//    nem_y = sin(nem_theta);
   
// }

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void azimuthal_ave_setup (){
    get_chiral_cell_center();
    make_bins();
    make_fine_bins();
    cell_position_check();
    edge_position_check();
    azimuthal_averaged_edge_length();   

}
