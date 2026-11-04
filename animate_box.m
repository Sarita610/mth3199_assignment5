function animate_box(tlist, Vlist, box_params)
    n = length(box_params.k_list);
    figure;
    hold on;
    axis equal;
    axis([-15, 15, -15, 15]);
    %plot for each spring
    for i = 1:n
        spring_plot_structs(i) = initialize_spring_plot(5,.5);
    end
    box_plot = plot(0, 0, 'k', 'linewidth', 2);

    t_uniform = linspace(tlist(1),tlist(end),500);
    V_uniform = interp1(tlist,Vlist,t_uniform);

    for j = 1:length(t_uniform)
        x = V_uniform(j,1); y = V_uniform(j,2); theta = V_uniform(j,3);

        % box outline
        bp = compute_rbt(x,y,theta,box_params.boundary_pts);
        set(box_plot,'xdata',bp(1,:),'ydata',bp(2,:));

        % springs
        for i = 1:n
            PA = box_params.P_world(:,i);
            PB = compute_rbt(x,y,theta,box_params.P_box(:,i));
            update_spring_plot(spring_plot_structs(i),PA,PB);
        end
        drawnow;
    end
end

