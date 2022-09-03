%% Construct figure for alias-free array example
%  This example is alias-free for the subset of DoAs with an azimuth angle |theta|<60 degrees
% 2022
% Authors: Lee C. Potter

close all
clear all
clc

%% Create sensor positions for planar array
L=4; % Number of sensors
theta_max = deg2rad(60);

v1 = 1.4*[0.40,sqrt(3)/2].'; v2 = [2.25,0].';
V = [v1, v2];
Z = [1,3;    2,2;    3,1];
P=2*Z*inv(V); %Construct pairwise differences (omit the factor of pi)
P=round(P,2);
V=2*inv(P'*P)*P'*Z;
positions = zeros([L,2]); positions(2:end,:)=P;                             %Element positions (in half-wavelengths)

% Parallogram vertices
vertices=zeros(4,2);
vertices(1,:)=V(:,1)-(V(:,1)+V(:,2))/2;
vertices(2,:)=-(V(:,1)+V(:,2))/2;
vertices(3,:)=V(:,2)-(V(:,1)+V(:,2))/2;
vertices(4,:)=V(:,1)+V(:,2)- (V(:,1)+V(:,2))/2;

ij = nchoosek(1:4,2);
minnorm_vertices = min( sqrt(sum((vertices(ij(:,2),:) - vertices(ij(:,1),:)).^2,2)));
disp("Minimum separation between parallogram vertices: "+minnorm_vertices)
minnorm_positions = min( sqrt(sum((positions(ij(:,2),:)-positions(ij(:,1),:)).^2,2)));
disp("Minimum separation between sensor elements (in half-wavelengths): "+minnorm_positions);

%% Figure for array geometry
figure(1);clf;
plot(positions(:,1),positions(:,2),'kx',...
    'MarkerSize',10,'Linewidth',3);
grid;axis square;axis([-1,1,-1,1]*3)
title('Array Geometry','fontsize',16)
ax = gca;
set(ax,'XTick',-2:1:2)
set(ax,'YTick',-2:1:2)
ax.FontSize = 16;

%% Theoretical value rather than solve via rational approximation
%V=diag([2, 2*sqrt(3)/3] )/(sidelength/2);
LatticeBasis=V;
display(LatticeBasis)
fprintf('Area of parallelogram: %g\n',abs(det(LatticeBasis)))
PVby2 = P*V/2;
display(PVby2)
disp("Array positions (units of half-wavelength)")
display(positions)

%% Plot lattice and ball
figure(2);clf;hold on;
ninepoints=[-1, -1; -1, 0; -1, 1; ...
             0, -1;  0, 0;  0, 1; ...
             1, -1;  1, 0;  1, 1];
x=ninepoints*[V(1,1);V(1,2)];
y=ninepoints*[V(2,1);V(2,2)];
r=1;

for k=1:size(ninepoints,1)
    th=linspace(-theta_max,theta_max,100).';
    % Construct hourglass centered on (x(k),y(k))
    seg = [x(k),y(k); r*cos(th)+x(k), r*sin(th)+y(k); x(k),y(k); r*cos(th+pi)+x(k), r*sin(th+pi)+y(k); x(k),y(k)];
    
    %p2=patch('Faces',1:size(seg,1),'Vertices',seg,'EdgeColor','none','FaceColor','b','FaceAlpha',.5);
    p2=plot(seg(:,1),seg(:,2),'Color','b','LineWidth',2);
end
p1=plot(x,y,'ro','MarkerSize',10,'Linewidth',2);
p3=quiver([0,0],[0,0],V(1,:),V(2,:),1,'k','Linewidth',2);

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
h_legend = legend([p1,p2,p3,p4],...
    {'Lattice','Boundary of $\mathcal{B}$','Basis','Parallelogram'},...
    'Location','southwest','Interpreter','latex')
axis square;axis([-1,1,-1,1]*2.5);

ax = gca;
set(ax,'Box','on')
tmp=floor(max(max(V))*1.5);
set(ax,'XTick',-3:1:3)
set(ax,'YTick',-3:1:3)
ax.FontSize = 16;

% Requires R2020a or later
exportgraphics(ax,'Figure1c.eps');
