function animate_box(tlist, Vlist, box_params)
    mypath1 = 'C:\Users\mepstein\OneDrive - Olin College of Engineering\Desktop\Applied Math\mth3199_assignment5\.git';
    fname='box_video.mp4';
    input_fname = [mypath1,fname];
    
    %create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);
    
    %must call open before writing any frames
    open(writerObj);

    %initialize the current figure and save as object
    fig1 = figure(1);
    xlabel('x-axis', 'Interpreter', 'latex','FontSize',15);
    ylabel('y-axis', 'Interpreter','latex','FontSize',15)
    title('Vibrating Box Animation (Fully Nonlinear)', 'Interpreter', 'latex','FontSize',20)

    %this is the critical line that forces the video to be high quality
    %adjust the numbers to adjust the position/size of the plotting window
    set(fig1,'units','pixels','position',[100 100 1440 1440]);
    set(fig1,'resize','off');

    n = length(box_params.k_list);
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

        %write the frame to the video
        pause(0.01);
        axis equal; 
        axis square
        axis([-15, 15, -15, 15])

        hold on;
        drawnow;

        %capture a frame (what is currently plotted)
        current_frame = getframe(gcf);
        writeVideo(writerObj,current_frame)
        drawnow;
    end
    close(writerObj)
end

