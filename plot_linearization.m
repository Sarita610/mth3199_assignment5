function plot_linearization()
    
    [tlist, Vlist, box_params] = simulate_box(); 
    
    x_guess = [0;0;0;0;0;0];
    my_rate_func = @box_rate_func;
    equil_func = @(V) my_rate_func(0,V);
    V_eq = multi_newton_solver(equil_func,x_guess, box_params);
    J_approx = approximate_jacobian(equil_func, V_eq);
    my_linear_rate = @(t_in,V_in) J_approx*(V_in-V_eq);


    % close ICs
    dx0 = 1;
    dy0 = 1;
    dtheta0 = 1;
    vx0 = 0;
    vy0 = 0;
    vtheta0 =0;
    

    %small number used to scale initial perturbation
    epsilon = 0.01;
    V0 = Veq + epsilon*[dx0;dy0;dtheta0;vx0;vy0;vtheta0];
    tspan = [0 20];
    %run the integration of nonlinear system
    [tlist_nonlinear,Vlist_nonlinear] = ode45(my_rate_func,tspan,V0);
    %run the integration of linear system
    [tlist_linear,Vlist_linear] = ode45(my_linear_rate,tspan,V0);

    t_uniform = linspace(tspan(1),tspan(end),500);
    V_uniform_nonlinear = interp1(tlist_nonlinear,Vlist_nonlinear,t_uniform_linear);
    V_uniform_linear = interp1(tlist_linear,Vlist_linear,t_uniform_linear);

    x_nonlinear = V_uniform_nonlinear(:,1); 
    y_nonlinear = V_uniform_nonlinear(:,2); 
    theta_nonlinear = V_uniform_nonlinear(:,3);        
    x_linear = V_uniform_linear(:,1); 
    y_linear = V_uniform_linear(:,2); 
    theta_linear = V_uniform_linear(:,3);        


    figure();
    hold on;
    subplot(3,1,1);
    plot(t_uniform, x_nonlinear)
    plot(t_uniform, x_linear)
    xlabel('Time', 'Interpreter','latex','FontSize',15);
    ylabel('$x$-position (-)', 'Interpreter','latex','FontSize',15);

    subplot(3,1,2);
    plot(t_uniform,y_nonlinear)
    plot(t_uniform,y_linear)
    xlabel('Time', 'Interpreter','latex','FontSize',15);
    ylabel('$y$-position (-)', 'Interpreter','latex','FontSize',15);

    subplot(3,1,3);
    plot(t_uniform,theta_nonlinear)
    plot(t_uniform,theta_linear)
    xlabel('Time', 'Interpreter','latex','FontSize',15);
    ylabel('Box Orientation (-)', 'Interpreter','latex','FontSize',15);

    Lgnd = legend('Nonlinear', 'Linearized', ...
        'Interpreter', 'latex', 'FontSize', 15);
    Lgnd.Position(1) = 0.01;
    Lgnd.Position(2) = 0.01;
end