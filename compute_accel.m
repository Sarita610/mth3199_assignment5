%Computesthelinearandangularaccelerationofthebox
%givenitscurrentpositionandorientation
%INPUTS:
%x:currentxpositionofthebox
%y:currentypositionofthebox
%theta:currentorientationofthebox
%box_params:astructcontainingtheparametersthatdescribethesystem
%Fields:
%box_params.m:massofthebox
%box_params.I:momentofinertiaw/respecttocentroid
%box_params.g:accelerationduetogravity
%box_params.k_list:listofspringstiffnesses
%box_params.l0_list:listofspringnaturallengths
%box_params.P_world:2xnlistofstaticmounting
% pointsforthespring(intheworldframe)
%box_params.P_box:2xnlistofmountingpoints
% forthespring(intheboxframe)
%
%OUTPUTS
%ax:xaccelerationofthebox
%ay:yaccelerationofthebox
%atheta:angularaccelerationofthebox

function[ax,ay,atheta]=compute_accel(x,y,theta,box_params)
    % define function parameters
    m = box_params.m;
    I = box_params.I;
    g = box_params.g;
    k_list = box_params.k_list;
    l0_list = box_params.l0_list;
    P_box = box_params.P_box;
    P_world = box_params.P_world;

    total_F = [0;0];
    total_torque = 0;
    
    for i = 1:length(k_list)
        % compute P_lists
        PA = compute_rbt(x,y,theta,P_box(:,i));
        PB = compute_rbt(x,y,theta,P_world(:,i));

        % compute spring force
        F = compute_spring_force(k_list(i),l0_list(i),PA,PB);
        Fx = F(1);
        Fy = F(2);
        total_F = total_F+F;
    
        % compute accelerations
        % accels = [0; -g] + (1/m)*F;
        % ax(i) = accels(1);
        % ay(i) = accels(2);

        rA = PA-[x;y];
        rx = rA(1);
        ry = rA(2);
        torque = rx*Fy-ry*Fx;
        total_torque = total_torque+torque;
    end

    ax = (1/m)*total_F(1);
    ay = -g+(1/m)*total_F(2);
    atheta = total_torque/I;

end