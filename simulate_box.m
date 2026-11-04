function [tlist, Vlist, box_params] = simulate_box()
    %define system parameters
    LW = 10; LH = 1; LG = 3;
    k = 2; l0 = 1.5*LG;

    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;

    P_box = [Pbl_box,Pbl_box, ...
             Pbr_box,Pbr_box, ...
             Ptl_box, Ptl_box, ...
             Ptr_box, Ptr_box];

    Pbl1_world = Pbl_box + [-LG;-LG];
    Pbl2_world = Pbl_box + [LG;-LG];
    Pbr1_world = Pbr_box + [0;-l0];
    Pbr2_world = Pbr_box + [l0;0];
    Ptl1_world = Ptl_box + [-LG;LG];
    Ptl2_world = Ptl_box + [-LG;0];
    Ptr1_world = Ptr_box + [LG;LG];
    Ptr2_world = Ptr_box + [0;LG];
    
    P_world = [Pbl1_world,Pbl2_world, ...
               Pbr1_world,Pbr2_world, ...
               Ptl1_world,Ptl2_world, ...
               Ptr1_world, Ptr2_world];

    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];

    box_params = struct();
    box_params.m = 1;
    box_params.I = (1/12)*(LH^2+LW^2);
    box_params.g = 9.8;
    box_params.k_list = [.5*k,.5*k,2*k,3*k,3*k,.5*k,2*k,3*k];
    box_params.l0_list = [.5*LG,.5*LG,2*LG,5*LG,5*LG,.5*LG,2*LG,5*LG];
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;

    %load the system parameters into the rate function
    %via an anonymous function
    my_rate_func = @(t_in,V_in) box_rate_func(t_in,V_in,box_params);
    x0 = 0;
    y0 = 0;
    theta0 = 0;
    vx0 = 5;
    vy0 = 5;
    vtheta0 = pi/4;
    V0 = [x0;y0;theta0;vx0;vy0;vtheta0];
    tspan = [0 20];
    %run the integration
    [tlist,Vlist] = ode45(my_rate_func,tspan,V0);
end
