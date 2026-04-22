//****************************************************************************
//****************************************************************************
//***************************BASAL NETWORK************************************
//****************************************************************************
//****************************************************************************
int make_cell(int poly_class, int e1, int e2, int e3, int e4, int e5, int e6, int e7, int e8, int e9, int e10, int e11, int e12){
    
    //int c_id=c_free_id;
    int c_id=c_freeId.get(Nc+1);
    //if(Time>=0) printf("cell %d created\n", c_id);
    basal_edges[c_id][1]=c_id;
    basal_edges[c_id][2]=poly_class;
    basal_edges[c_id][3]=e1;
    basal_edges[c_id][4]=e2;
    basal_edges[c_id][5]=e3;
    basal_edges[c_id][6]=e4;
    basal_edges[c_id][7]=e5;
    basal_edges[c_id][8]=e6;
    basal_edges[c_id][9]=e7;
    basal_edges[c_id][10]=e8;
    basal_edges[c_id][11]=e9;
    basal_edges[c_id][12]=e10;
    basal_edges[c_id][13]=e11;
    basal_edges[c_id][14]=e12;
    
    //ATTRIBUTES
    c_A0[c_id]=A0;
    c_type[c_id]=1;
    c_fixed[c_id] = 0;
    A0tot+=c_A0[c_id];
    
    if(c_id>Nc) Nc++;
    
    //c_free_id=c_freeId.get(Nc+1);
    
    return c_id;
}
//****************************************************************************
