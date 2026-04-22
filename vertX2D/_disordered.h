//****************************************************************************
void calc_forces_tension(){
    
    //AREA TERM
    wA=0;
    wP=0;
    for(int i = 1; i <= Nc; i++) if(basal_edges[i][1]!=0) 
    {
        wA += c_AreaCompressibility_force_New(i);
        wP += c_perimeter_force(i);
        //c_nematic_stress_force(i);
    }
    
    //PERIMETER TERM & TENSION DYNAMICS
    wl=0;
    for(int i = 1; i <= Ne; i++) if(e[i][0]!=0){
        
        //FORCES
        wl += e_BasalLength_force(i);
        
    }
}
//****************************************************************************
double eqOfMotion_tension(){
    
    //CALCULATES FORCES
    reset_forces();
    calc_forces_tension();
    
    //EQUATION OF MOTION
    for(int i = 1; i <= Nv; i++) if(v[i][0]>0.5){
            
            //MOVES VERTICES
            double pm_val = 0;

            if (v_type[i] != 1) for(int j = 1; j<=2; j++) 
            {
                v[i][j] += h*v_F[i][j];
                pm_val += pow(v_F[i][j], 2.);
            }
            pm_val = h*sqrt(pm_val);
            if (pm_val > max_move) max_move = pm_val;
            
            //PERIODIC BOUNDARY CONDITIONS
            torus_vertex(i);
    }
    
    //FIXES PASSIVE VERTICES
    calc_central_vertices();
    
    //OUTPUT
    //printf("%g \t\t wA=%.20g \t\t wP=%.20g \t\t w=%.20g \t\t %g\n", Time, wA, wP, wA+wP, max_move);
    
    //t=t+dt
    Time+=h;
    
    return 0;
}
//****************************************************************************
void make_disordered(double Sig_val, int tiling_seed)
{
    
    cout << "disordered tissue start" << endl;

    if (tiling_seed != 0)
    {
        //******************************************
        //**************RANDOM SEED*****************
        //******************************************
        ifstream seed_file;
        string filename = "seeds//seed_" + to_string(tiling_seed) + ".dat";
        seed_file.open (filename);
        seed_file >> gen_sigma;
        seed_file.close();
    }

    double tINITIAL = 1000;
    double P0_store = P0;
    double kPer_store = kPer;

    P0 = 3.6;
    kPer = 0.01;
    Sig = Sig_val;

    while(Time < tINITIAL - h /2.){  
        //MOTION
        line_tension_fluctuations();
        eqOfMotion_tension();//calculates forces and moves vertices, then increments time by h

        //TOPOLOGY
        update_topology(0.011, 0); //T1 transitions
    }

    reset_line_tensions();
    Sig = 0;
    Time = 0;
    wl = 0;
    P0 = P0_store;
    kPer = kPer_store;
    cout << "disordered tissue generated" << endl;
}
