//****************************************************************************
//****************************************************************************
//***********************VERT, EDG, FAC***************************************
//****************************************************************************
//****************************************************************************
int make_vertex(double x, double y){
    
    //int v_id=v_free_id;
    int v_id=v_freeId.get(Nv+1);
    //if(Time>=0) printf("vertex %d created\n", v_id);
    v[v_id][0]=v_id;
    v[v_id][1]=x;
    v[v_id][2]=y;
    torus_vertex(v_id);
    
    v_edges[v_id][1]=v_id;
    v_cells[v_id][1]=v_id;
    v_type[v_id] = 0;
    
    if(v_id>Nv) Nv++;
    
    //v_free_id=v_freeId.get(Nv+1);
    
    return v_id;
}
//****************************************************************************
int make_vertex_pass(double x, double y){
    
    //int v_id=v_pass_free_id;
    int v_id=v_pass_freeId.get(Nv_pass+1);
    //if(Time>=0) printf("passive vertex %d created\n", v_id);
    v_pass[v_id][0]=v_id;
    v_pass[v_id][1]=x;
    v_pass[v_id][2]=y;
    
    if(v_id>Nv_pass) Nv_pass++;
    
    //v_pass_free_id=v_pass_freeId.get(Nv_pass+1);
    
    return v_id;
}
//****************************************************************************
int make_edge(int v1, int v2){
    
    //int e_id=e_free_id;
    int e_id=e_freeId.get(Ne+1);
    //if(Time>=0) printf("edge %d created\n", e_id);
    e[e_id][0]=e_id;
    e[e_id][1]=v1;
    e[e_id][2]=v2;
    
    //v_edges
    //v1
    v_edges[v1][2]++;
    v_edges[v1][2+v_edges[v1][2]]=e_id;
    //v2
    v_edges[v2][2]++;
    v_edges[v2][2+v_edges[v2][2]]=e_id;
    
    //attributes
    e_length[e_id]=0;//edge_length(e_id);
    e_dldt[e_id]=0;
    e_g[e_id]=gamma0;
    
    if(e_id>Ne) Ne++;
    
    //e_free_id=e_freeId.get(Ne+1);
    
    return e_id;
}
//****************************************************************************
int make_facet(int v1, int v2, int v3){
    
    //int f_id=f_free_id;
    int f_id=f_freeId.get(Nf+1);
    //if(Time>=0) printf("facet %d created\n", f_id);
    f[f_id][0]=f_id;
    f[f_id][1]=v1;
    f[f_id][2]=v2;
    f[f_id][3]=v3;
    
    if(f_id>Nf) Nf++;
    
    //f_free_id=f_freeId.get(Nf+1);
    
    return f_id;
}
//****************************************************************************
