%% plot the outline for the device
% -inputs:
% gfile                            structure,  gfile structure
% nw_refine/nh_refine              int, the refined grids
% position_coil                    matrix, position of the PF coils (the columns are:R, Z, Rwidth, Zwidth, flag for coil_grid_r_z.m, turns in R direction (used just for the number of grids in coil_grid_r_z.m,which can be arbitrary value), turns in Z direction, total turns.)
% vv                               vector, the dimensions of the VV, [defaut for EXL: 0.158 1.655 -1.405 1.405]
%rlim/zlimter

% -outputs:
% plot of the outline of the device

% Edited by Bin Chen in 2019/06/13
% Contact: chenbino@enn.cn
% ENN Group 1989-2019, all rights reserved.

% show the limiter just at the inner boundary
rlimter=[0.17 0.17 1.525 1.525 0.17 0.17];
zlimter=[0.000 1.245 1.245 -1.245 -1.245 0.000];
% rlimter(length(rlimter)+1)=rlimter(1);
% zlimter(length(zlimter)+1)=zlimter(1);

hold off;
plot(rlimter,zlimter,'color',[0.5 0.5 0.5],'linewidth',2);
set(gca,'fontweight','bold','fontsize',16,'fontname','times','XTick',[0:0.5:2.0],'YTick',[-2:0.5:2.0],'XTickLabel',{'0','0.5','1.0','1.5','2.0'}, ... 
        'YTickLabel',{'-2.0','-1.5','-1.0','-0.5','0.0','0.5','1.0','1.5','2.0'},'Color', 'None');
axis equal;
axis([0 2.2 -2.0 2.0]);
% axis([0 3.0 -3.0 3.0]);
xlabel('R (m)','FontSize',18);
ylabel('Z (m)','FontSize',18);

position_coil=[
    0.929	1.913	0.143	0.215	0 4 6 24 738
    0.929	-1.913	0.143	0.215	0 4 6 24 738
    2.108	1.335	0.143	0.215	0 4 6 24 738
    2.108	-1.335	0.143	0.215	0 4 6 24 738
    2.108	0.445	0.143	0.215	0 4 6 24 738
    2.108	-0.445	0.143	0.215	0 4 6 24 738
    ]; %the columns are:R, Z, Rwidth, Zwidth, flag for coil_grid_r_z.m, turns in R direction (used just for the number of grids in coil_grid_r_z.m,which can be arbitrary value), turns in Z direction, total turns.

hold on;
for i=1:length(position_coil(:,1))
    plot([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)-position_coil(i,4)/2],'b','linewidth',1.5);
    hold on;
    plot([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)+position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],'b','linewidth',1.5);
    plot([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)-position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],'b','linewidth',1.5);
    plot([position_coil(i,1)+position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],'b','linewidth',1.5);
    plot([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],'b','linewidth',1.5);
    plot([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)+position_coil(i,4)/2 position_coil(i,2)-position_coil(i,4)/2],'b','linewidth',1.5);
end

vv=[0.158 1.655 -1.405 1.405];  %%%the dimensions of the VV
plot([vv(1) vv(2)],[vv(4) vv(4)],'color',[0.667 0.667 1],'linewidth',2.5);
plot([vv(2) vv(2)],[vv(4) vv(3)],'color',[0.667 0.667 1],'linewidth',2.5);
plot([vv(2) vv(1)],[vv(3) vv(3)],'color',[0.667 0.667 1],'linewidth',2.5);
plot([vv(1) vv(1)],[vv(3) vv(4)],'color',[0.667 0.667 1],'linewidth',2.5);