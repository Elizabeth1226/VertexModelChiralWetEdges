//****************************************************************************
//****************************************************************************
//*******************************TORUS****************************************
//****************************************************************************
//****************************************************************************
void torus_dx_dy_dz(double *dxdydz, int vert_id, int vert_ref_id){
    
    dxdydz[0]=0;
    dxdydz[1]=0;
    
    //x
    if(fabs(v[vert_id][1]-v[vert_ref_id][1])>0.5*perioXYZ[0]){
        if      (v[vert_id][1]<v[vert_ref_id][1])   dxdydz[0]=perioXYZ[0];
        else if (v[vert_id][1]>v[vert_ref_id][1])   dxdydz[0]=-perioXYZ[0];
    }
    else dxdydz[0]=0;
    //y
    if(fabs(v[vert_id][2]-v[vert_ref_id][2])>0.5*perioXYZ[1]){
        if      (v[vert_id][2]<v[vert_ref_id][2])   dxdydz[1]=perioXYZ[1];
        else if (v[vert_id][2]>v[vert_ref_id][2])   dxdydz[1]=-perioXYZ[1];
    }
    else dxdydz[1]=0;
}
//****************************************************************************
void torus_pass_dx_dy_dz(double *dxdydz, int vert_id, int vert_ref_id){
    
    dxdydz[0]=0;
    dxdydz[1]=0;
    
    //x
    if(fabs(v_pass[vert_id][1]-v_pass[vert_ref_id][1])>0.5*perioXYZ[0]){
        if      (v_pass[vert_id][1]<v_pass[vert_ref_id][1])   dxdydz[0]=perioXYZ[0];
        else if (v_pass[vert_id][1]>v_pass[vert_ref_id][1])   dxdydz[0]=-perioXYZ[0];
    }
    else dxdydz[0]=0;
    //y
    if(fabs(v_pass[vert_id][2]-v_pass[vert_ref_id][2])>0.5*perioXYZ[1]){
        if      (v_pass[vert_id][2]<v_pass[vert_ref_id][2])   dxdydz[1]=perioXYZ[1];
        else if (v_pass[vert_id][2]>v_pass[vert_ref_id][2])   dxdydz[1]=-perioXYZ[1];
    }
    else dxdydz[1]=0;
}
//****************************************************************************
void torus_passvert_dx_dy_dz(double *dxdydz, int vert_id, int vert_ref_id){
    
    dxdydz[0]=0;
    dxdydz[1]=0;
    
    //x
    if(fabs(v[vert_id][1]-v_pass[vert_ref_id][1])>0.5*perioXYZ[0]){
    	if      (v[vert_id][1]<v_pass[vert_ref_id][1])   dxdydz[0]=perioXYZ[0];
        else if (v[vert_id][1]>v_pass[vert_ref_id][1])   dxdydz[0]=-perioXYZ[0];
    }
    else dxdydz[0]=0;
    
    //y
    if(fabs(v[vert_id][2]-v_pass[vert_ref_id][2])>0.5*perioXYZ[1]){
        if      (v[vert_id][2]<v_pass[vert_ref_id][2])   dxdydz[1]=perioXYZ[1];
        else if (v[vert_id][2]>v_pass[vert_ref_id][2])   dxdydz[1]=-perioXYZ[1];
    }
    else dxdydz[1]=0;
}
//****************************************************************************
void torus_vertex(int i){
    //x
    if      (v[i][1]<0)             v[i][1] += perioXYZ[0];
    else if (v[i][1]>perioXYZ[0])   v[i][1] -= perioXYZ[0];
    //y
    if      (v[i][2]<0)             v[i][2] += perioXYZ[1];
    else if (v[i][2]>perioXYZ[1])   v[i][2] -= perioXYZ[1];
}
//****************************************************************************
void torus_vertex_pass(int i){
    
    //x
    if(v_pass[i][1]<0)              v_pass[i][1]+=perioXYZ[0];
    if(v_pass[i][1]>perioXYZ[0])    v_pass[i][1]-=perioXYZ[0];
    //y
    if(v_pass[i][2]<0)              v_pass[i][2]+=perioXYZ[1];
    if(v_pass[i][2]>perioXYZ[1])    v_pass[i][2]-=perioXYZ[1];
    
}
//****************************************************************************
void expand_box(double dpx, double dpy){
    
    //vertices
    for(size_t i = 1; i<=Nv; ++i){
        if(v[i][0]>0.5){
            v[i][1]*=(1.+dpx);
            v[i][2]*=(1.+dpy);
        }
    }
    
    //passive vertices
    for(size_t i = 1; i<=Nv_pass; ++i){
        if(v_pass[i][0]>0.5){
            v_pass[i][1]*=(1.+dpx);
            v_pass[i][2]*=(1.+dpy);
        }
    }
    
    perioXYZ[0]*=(1.+dpx);
    perioXYZ[1]*=(1.+dpy);
    
    printf("new px: %g  new py: %g\n", perioXYZ[0], perioXYZ[1]);
    
}
//****************************************************************************
