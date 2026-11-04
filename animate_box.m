function animate_box(tlist, Vlist, box_params)
    n = lenght(box_params.k_list);
    figure;
    hold on;
    axis equal;
    axis([-15, 15, -15, 15]);
    %plot for each spring
    for i = 1:n
        spring_plot_structs(i) = initialize_spring_plot(5,.5);
    end

    