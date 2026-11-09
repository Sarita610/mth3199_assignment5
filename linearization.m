[~, ~, box_params] = simulate_box();
my_rate_func = @(t_in, V_in) box_rate_func(t_in, V_in, box_params);
rate_eq = @(V) my_rate_func(0, V);

V_guess = zeros(6,1);
solver_params = struct();
solver_params.ftol = 1e-10;
solver_params.max_iter = 200;
solver_params.numerical_diff = 1; 
V_eq = multi_newton_solver(rate_eq, V_guess, solver_params);
disp(norm(rate_eq(V_eq)))

J_approx = approximate_jacobian(rate_eq, V_eq);
my_linear_rate = @(t_in,V_in) J_approx*(V_in-V_eq);

dx0 = 1;
dy0 = 0;
dtheta0 = 0;
vx0 = 0;
vy0 = 0;
vtheta0 = 0;
dV = [dx0;dy0;dtheta0;vx0;vy0;vtheta0];


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

 %small number used to scale initial perturbation
 epsilons = [0.01, 0.1, 0.5, 0.8];

 omega = abs(imag(eig(J_approx)));
 omega =omega(omega>1e-6);
 tspan = [0, 3*2*pi/min(omega)];

 for j = 1:length(epsilons)
    V0 = V_eq + epsilons(j)*dV;

    %run the integration of nonlinear system
    [tlist_nonlinear,Vlist_nonlinear] =...
    explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,BT_struct, p, error_desired);
    %run the integration of linear system
    [tlist_linear,Vlist_linear] =...
    explicit_RK_variable_step_integration(my_linear_rate,tspan,V0,h_ref,BT_struct, p, error_desired);

 end 