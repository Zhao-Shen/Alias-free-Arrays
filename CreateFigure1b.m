%% Code to construct Figure 3
% Code to construct figure for array example:
% no fundamental parallelegram contains ball
% 2022
% Authors: Lee C. Potter


close all
clear all
clc

%% Create sensor positions for planar array
N=4;                                                                        %number of sensors
V=[2, 1; 0, 2];
Z=[0,0;1,3;2,2;3,1];                                                        %integers
B=2*Z*inv(V);
theta=atan(-3/2);
U=[cos(theta), sin(theta); -sin(theta), cos(theta)];
B = B*U';%rotate
V = U*V;
positions=B;
positions=positions - repmat(positions(3,:),4,1);
positions(1,:)=-positions(1,:);
% positions = B-repmat(mean(B,1),N,1);
% positions = positions-repmat([0, positions(2,2)],N,1);
% positions(:,2)=-positions(:,2);
positions = 0.95*positions;
V=V/0.95;
positions = [   0.26       2.11
                -1.71       0
                0           0
                1.71        0];


%% Figure for array geometry
figure;
plot(positions(:,1),positions(:,2),'kx',...
    'MarkerSize',10,'Linewidth',3);
grid;axis square;axis([-1,1,-1,1]*3)
title('Array Geometry','fontsize',16)
ax = gca;
set(ax,'XTick',-2:1:2)
set(ax,'YTick',-2:1:2)
ax.FontSize = 16;
pause(0.25)


%% Construct pairwise differences 
%(omit the factor of pi)
Nrows = N-1;
P = zeros(Nrows,2);
for k=2:N
        P(k-1,:)=(positions(k,:)-positions(1,:));
end


%% Theoretical value rather than solve via rational approximation
%V=diag([2, 2*sqrt(3)/3] )/(sidelength/2);
LatticeBasis=V;
display(LatticeBasis)
fprintf('Area of parallelogram: %g\n',abs(det(LatticeBasis)))
PVby2 = P*V/2;
display(PVby2)
display(positions)

%% Plot lattice and ball
figure;
ninepoints=[-1, -1; -1, 0; -1, 1; ...
             0, -1;  0, 0;  0, 1; ...
             1, -1;  1, 0;  1, 1];
x=ninepoints*[V(1,1);V(1,2)];
y=ninepoints*[V(2,1);V(2,2)];
p1=plot(x,y,'ro','MarkerSize',10,'Linewidth',3);hold on
viscircles([x, y], ones(size(x)),'Color','b');
p2=plot(-1000,1000,'b-','Linewidth',3);                                     %legend entry for viscircles
p3=quiver([0,0],[0,0],V(1,:),V(2,:),'off',...
    'Color','k','LineWidth',2);
% parallelogram
tmp=zeros(2,3);
tmp(:,1)=V(:,1);
tmp(:,3)=V(:,2);
tmp = tmp -(V(:,1)+V(:,2))/2;
p4=plot(tmp(1,:),tmp(2,:),'Color','g','LineWidth',2);
% remainder of parallelogram (<1 part)
tmp(:,1)=V(:,1);
tmp(:,2)=V(:,1)+V(:,2);
tmp(:,3)=V(:,2);
tmp = tmp -(V(:,1)+V(:,2))/2;
p5=plot(tmp(1,:),tmp(2,:),'g--','LineWidth',2);
% finish labeling graph
legend([p1,p2,p3,p4],...
    {'Lattice','Boundary of $\mathcal{B}$','Basis','Parallelogram'},...
    'Location','southwest','Interpreter','latex')
axis square;axis([-1,1,-1,1]*2.5);                                          %max(max(V))*1.5)
%title('Lattice','fontsize',16)
ax = gca;
tmp=floor(max(max(V))*1.5);
set(ax,'XTick',-3:1:3)
set(ax,'YTick',-3:1:3)
ax.FontSize = 16;

ax = gca;
% Requires R2020a or later
exportgraphics(ax,'Figure1b.eps');
