function modal_animation
    
    %define system parameters
    LW = 10; LH = 1; LG = 3;
    k = 20; l0 = 1.5*LG;

    Pbl_box = [-LW;-LH]/2;
    Pbr_box = [LW;-LH]/2;
    Ptl_box = [-LW;LH]/2;
    Ptr_box = [LW;LH]/2;

    % P_box = [Pbl_box,Pbl_box, ...
    %          Pbr_box,Pbr_box, ...
    %          Ptl_box, Ptl_box, ...
    %          Ptr_box, Ptr_box];

    P_box = [Pbl_box,Pbl_box,Pbr_box,Pbr_box];

    Pbl1_world = Pbl_box + [-LG;-LG];
    Pbl2_world = Pbl_box + [LG;-LG];
    Pbr1_world = Pbr_box + [0;-l0];
    Pbr2_world = Pbr_box + [l0;0];
    
    % P_world = [Pbl1_world,Pbl2_world, ...
    %            Pbr1_world,Pbr2_world, ...
    %            Ptl1_world,Ptl2_world, ...
    %            Ptr1_world, Ptr2_world];

    P_world = [Pbl1_world,Pbl2_world,Pbr1_world,Pbr2_world];

    boundary_pts = [Pbl_box,Pbr_box,Ptr_box,Ptl_box,Pbl_box];

    box_params = struct();
    box_params.m = 1;
    box_params.I = (1/12)*(LH^2+LW^2);
    box_params.g = 1;
    box_params.k_list = [.5*k,.5*k,2*k,5*k];   
    box_params.l0_list = l0*ones(1,4);                                 
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;

    my_rate_func = @(t,V) box_rate_func(t, V, box_params);

    %%Fehlburg
    Fehlberg = struct();
    Fehlberg.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
    Fehlberg.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
    25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    Fehlberg.A = [0,0,0,0,0,0;...
    1/4, 0,0,0,0,0;...
    3/32, 9/32, 0,0,0,0;...
    1932/2197, -7200/2197, 7296/2197, 0,0,0;...
    439/216, -8, 3680/513, -845/4104, 0,0;...
    -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
     BT_struct = Fehlberg;
    
    p=5;
    error_desired = 1e-8;
    h_ref = 0.01;
    
    %%
    rate_eq = @(V) my_rate_func(0, V);
    V_guess = zeros(6,1);

    solver_params = struct();
    solver_params.ftol = 1e-10;
    %solver_params.dxmin = 1e-12;
    solver_params.max_iter = 200;
    solver_params.numerical_diff = 1; 

    V_eq = multi_newton_solver(rate_eq, V_guess, solver_params);
    disp(norm(rate_eq(V_eq)));

    J_approx = approximate_jacobian(rate_eq, V_eq);
    %my_linear_rate = @(t_in,V_in) J_approx*(V_in - V_eq);

    dx0 = 1;
    dy0 = 0;
    dtheta0 = 0;
    vx0 = 0;
    vy0 = 0;
    vtheta0 = 0;
    dV = [dx0;dy0;dtheta0;vx0;vy0;vtheta0];

    %mode
    
    Q = -J_approx(4:6,1:3);
    Q = (Q+Q.')/2;
    [Umodes, D] = eig(Q);
    [lambda, idx] = sort(diag(D));
    Umodes = Umodes(:, idx);
    omega = sqrt(lambda);
    
%% animation
    mode = 1;
    epsilon_anim = 0.5;
    n_cycl = 4;
    
    T = 2*pi/omega(mode);
    V0 = V_eq + epsilon_anim*[Umodes(:,mode);0;0;0];
    [tA, VA] = explicit_RK_variable_step_integration(my_rate_func,[0, n_cycl*T], V0, h_ref, BT_struct, p, error_desired);
    t_uniform = linspace(0, n_cycl*T, n_cycl*60);
    V_uniform = interp1(tA, VA, t_uniform);
    %
    box_radius = max(vecnorm(boundary_pts));  
    x_all = [P_world(1,:), V_uniform(:,1)']; 
    y_all = [P_world(2,:), V_uniform(:,2)'];
    
    xlims = [min(x_all) - box_radius, max(x_all) + box_radius];
    ylims = [min(y_all) - box_radius, max(y_all) + box_radius];
    
    fig1 = figure();
    set(fig1,'units','pixels','position',[100 100 1440 1440]);
    set(fig1,'resize','off');

    hold on; 
    axis equal;
    xlim(xlims); 
    ylim(ylims);
    title("Mode Animation", 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('$x$ position (-)', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$y$ position (-)', 'Interpreter', 'latex','FontSize', 15);
    
    box_plot = plot(0, 0, 'k', 'LineWidth', 3);
    n = length(box_params.k_list);
    for i = 1:n
        spring_plot_structs(i) = initialize_spring_plot(5,.5);
    end
   
    vid = VideoWriter('box_mode_animation');
    open(vid);
    
    for k = 1:numel(t_uniform)
        x = V_uniform(k,1);
        y = V_uniform(k,2);
        theta = V_uniform(k,3);
        bp = compute_rbt(x, y, theta, boundary_pts);
        PB = compute_rbt(x, y, theta, P_box);
    
        set(box_plot,'xdata',bp(1,:),'ydata', bp(2,:));
        for i = 1:n
            update_spring_plot(spring_plot_structs(i), P_world(:,i), PB(:,i));
        end
    
        drawnow;
        writeVideo(vid, getframe(fig1));
    end
    close(vid);

end 