//****************************************************************************
double A_dot_b_full(double Axx, double Axy, double Ayx, double Ayy, double bx, double by, int coord)
{
    if (coord == 1)
    {
        return Axx * bx + Axy * by;
    }
    else if (coord == 2)
    {
        return Ayx * bx + Ayy * by;
    }
    else
    {
        cout << "ERROR: wrong coord!" << endl;
        return 0;
    }
}
//****************************************************************************
void c_chiral_force(int i)
{

    double Q_xx = 0.;
    double Q_xy = 1.;
    double Q_yx = -1.;
    double Q_yy = 0.;

    double chiral_alpha_c = chiral_alpha_1 * c_activity[i];

    for (int j = 1; j <= basal_vertices[i][2]; j++)
    {
        int vi = basal_vertices[i][2 + j];

        int j_min1 = edge_circle(j - 1, basal_vertices[i][2]);
        int j_plus1 = edge_circle(j + 1, basal_vertices[i][2]);

        int v_min1 = basal_vertices[i][2 + j_min1];
        int v_plus1 = basal_vertices[i][2 + j_plus1];

        double ee_x, ee_y;
        compute_edge_vector(v_min1, v_plus1, ee_x, ee_y);

        double Ri_x = 0.5 * cross_z(ee_x, ee_y, 1);
        double Ri_y = 0.5 * cross_z(ee_x, ee_y, 2);
        // if (i == 100) cout << "R: " << Ri_x << " " << Ri_y << endl;

        double nem_Fx = -1. * chiral_alpha_c * A_dot_b_full(Q_xx, Q_xy, Q_yx, Q_yy, Ri_x, Ri_y, 1);
        double nem_Fy = -1. * chiral_alpha_c * A_dot_b_full(Q_xx, Q_xy, Q_yx, Q_yy, Ri_x, Ri_y, 2);

        // if (i == 100) cout << "F: " << nem_Fx << " " << nem_Fy << endl;

        v_F[vi][1] += nem_Fx;
        v_F[vi][2] += nem_Fy;
    }
}

//****************************************************************************
//***********************INITIAL CONDITIONS***********************************
//****************************************************************************
//****************************************************************************
void fix_cell_vertices(int i)
{   
    for (int j = 1; j <= basal_vertices[i][2]; j++)
    {
        v_type[basal_vertices[i][2 + j]] = 1;
    }
    c_fixed[i] = 1;

    return;
}

//****************************************************************************
void fix_disk()
{
    //*******************************************************
    //********************SET FIXED BOUNDARY*****************
    //*******************************************************

    double disk_center_x = perioXYZ[0] / 2;
    double disk_center_y = perioXYZ[1] / 2;
    cout << disk_center_x << " " << disk_center_y << endl;

    for (int i = 1; i <= Nc; i++)
    {
        double cell_x = center_coord(i, 1);
        double cell_y = center_coord(i, 2);
        double dist_val_x = cell_x - disk_center_x;
        double dist_val_y = cell_y - disk_center_y;
        double dr = sqrt(dist_val_x * dist_val_x + dist_val_y * dist_val_y);
        if (dr > 0.48 * perioXYZ[1])
            fix_cell_vertices(i);
    }
}

void node_initial(double node_radius, double node_center_y)
{
    for (int i = 1; i <= Nc; i++)
    {
        double cell_x = center_coord(i, 1);
        double cell_y = center_coord(i, 2);
        double dist_val_x = cell_x - perioXYZ[0] / 2.;
        double dist_val_y = cell_y - node_center_y;
        double dr = sqrt(dist_val_x * dist_val_x + dist_val_y * dist_val_y);
        if (dr < node_radius)
            c_activity[i] = 1.;
        else
            c_activity[i] = 0.;
    }
}

void fix_streak(double streak_length, double streak_width, double streak_center_y)
{
    for (int i = 1; i <= Nc; i++)
    {
        double cell_x = center_coord(i, 1);
        double cell_y = center_coord(i, 2);
        double dist_val_x = abs(cell_x - perioXYZ[0] / 2.);
        double dist_val_y = abs(cell_y - streak_center_y);
        if (dist_val_y < streak_length / 2. && dist_val_x < streak_width / 2.)
            fix_cell_vertices(i);
    }
}
//****************************************************************************
void chick_initial(double streak_length, double streak_width, double streak_center_y, double node_radius, double node_center_y)
{
    fix_streak(streak_length, streak_width, streak_center_y);
    fix_disk();

    node_initial(node_radius, node_center_y);

    return;
}

//****************************************************************************
void disk_initial(double radius_fraction, int radius_seed, int Nx, int transition_type)
{
    double xc_val, yc_val;

    double ddx = 1.0745699318;
    double ddy = 1.2408064788;

    double xmin = 0.08;
    double ymin = 1.0106048591;

    double x_mean = xmin + ddx * (((double)Nx) / 2. - 0.5);
    double y_mean = ymin + 0.5 * ddy + 3. / 2. * ddy * (((double)Nx) / 4. - 1.);

    cout << x_mean << " " << y_mean << endl;

    double radius_threshold = 0.5 * radius_fraction * perioXYZ[1];

    double radius_val;

    for (int i = 1; i <= Nc; i++)
    {
        xc_val = center_coord(i, 1);
        yc_val = center_coord(i, 2);

        radius_val = sqrt((xc_val - x_mean) * (xc_val - x_mean) + (yc_val - y_mean) * (yc_val - y_mean));

        if (transition_type == 1)
        {
            // smooth transition (cosine)
            if (radius_val < 0.5 * radius_threshold)
            {
                c_activity[i] = 1.0;
            }
            else if (radius_val < radius_threshold)
            {
                c_activity[i] = 0.5 * (cos((radius_val - 0.5 * radius_threshold) / radius_threshold * 2. * M_PI)) + 0.5;
            }
            else
            {
                c_activity[i] = 0.0;
            }
        }
        else if (transition_type == 0)
        {
            // sharp transition (transition_type == 0)
            if (radius_val < radius_threshold)
            {
                c_activity[i] = 1.;
            }
            else
            {
                c_activity[i] = 0.;
            }
        }
        else
        {
            cout << "ERROR: wrong transition type!" << endl;
        }
    }

    return;
}
//****************************************************************************
void mixed_initial(double population_fraction, int population_seed, int Nx)
{
    std::mt19937 gen(population_seed);
    std::uniform_int_distribution<> dis(1, (int)Nc);

    for (int i = 1; i <= Nc; i++)
    {
        c_activity[i] = 2;
    }

    double selected_cells = 0;
    int rnd_cell_id;

    while (selected_cells < population_fraction * ((double)Nc) - 0.00001)
    {
        rnd_cell_id = dis(gen);
        if (c_activity[rnd_cell_id] == 2)
        {
            c_activity[rnd_cell_id] = 1;
            selected_cells += 1.;
        }
    }

    return;
}
//****************************************************************************
void split_initial_x(double radius_fraction, int radius_seed, int Nx)
{
    double xc_val, yc_val;

    double x_threshold = radius_fraction * perioXYZ[0];

    for (int i = 1; i <= Nc; i++)
    {
        xc_val = center_coord(i, 1);

        if (xc_val < x_threshold)
        {
            c_activity[i] = 1;
        }
        else
        {
            c_activity[i] = 2;
        }
    }

    return;
}
//****************************************************************************
void split_initial_y(double radius_fraction, int radius_seed, int Nx)
{
    double xc_val, yc_val;

    double y_threshold = radius_fraction * perioXYZ[1];

    for (int i = 1; i <= Nc; i++)
    {
        yc_val = center_coord(i, 2);

        if (yc_val < y_threshold)
        {
            c_activity[i] = 1;
        }
        else
        {
            c_activity[i] = 2;
        }
    }

    return;
}
