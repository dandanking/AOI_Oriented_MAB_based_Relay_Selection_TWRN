function draw_network_3D(wlans)
% DrawNetwork3D - Plots a 3D of the network 
%   INPUT: 
%       * wlan - contains information of each WLAN in the map. For instance,
%       wlan(1) corresponds to the first one, so that it has unique
%       parameters (x,y,z,BW,CCA,etc.)

    load('constants.mat')
    
    num_wlans = size(wlans, 2);
    
    for j = 1 : num_wlans
        x(j) = wlans(j).xn;
        y(j) = wlans(j).yn;
        z(j) = wlans(j).zn;
    end
    
    figure
    axes;
    set(gca,'fontsize',16);
    labels = num2str((1:size(y' ))','%d');    
    for i = 1 : num_wlans
        scatter3(wlans(i).x, wlans(i).y, wlans(i).z, 30, [0 0 0], 'filled');
        hold on;   
        scatter3(wlans(i).xn, wlans(i).yn, wlans(i).zn, 50, [0 0 1], 'filled');
        scatter3(wlans(i).xd, wlans(i).yd, wlans(i).zd, 30, [0 0 0], 'filled');
        line([wlans(i).x, wlans(i).xn], ...
            [wlans(i).y, wlans(i).yn], ...
            [wlans(i).z, wlans(i).zn], 'Color', [0.4, 0.4, 1.0], 'LineStyle', ':');
        line([wlans(i).xd, wlans(i).xn], ...
            [wlans(i).yd, wlans(i).yn], ...
            [wlans(i).zd, wlans(i).zn], 'Color', [0.4, 0.4, 1.0], 'LineStyle', ':');
    end
    text(x,y,z,labels,'horizontal','left','vertical','bottom') 
    xlabel('x [meters]','fontsize',14);
    ylabel('y [meters]','fontsize',14);
    zlabel('z [meters]','fontsize',14);
    axis([0 MaxX 0 MaxY 0 MaxZ])
end