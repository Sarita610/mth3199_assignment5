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
    k = box_params.k;
    l0 = box_params.l0;
    P_box = box_params.P_box;

    % compute P_lists
    PA = compute_rbt(x,y,theta,P_box(:,1));
    PB = compute_rbt(x,y,theta,P_box(:,2));
    
    % compute spring force
    F = compute_spring_force(k,l0,PA,PB);
    Fx = F(:,1);
    Fy = F(:,2);
    
    % compute accelerations
    accels = [0; -g] + (1/m)*sum(F, 2);
    ax = accels(:,1);
    ay = accels(:,2);
    
    atheta = zeros(size(ax));
    for i = 1:length(accels)
        atheta(i) = (1/I)*(PA(i)-x)*Fy-(PB(i)-y)*Fx;
    end

end