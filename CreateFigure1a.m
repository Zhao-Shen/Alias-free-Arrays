%% Code to construct Figure 1
% Code to construct figure for uniform circular array
% 2022
% Authors: Lee C. Potter

close all
clear all
clc


%% Create sensor positions for circular array
sidelength=2;                                                               %units of 1/2 wavelength
%sidelength=2*sqrt(3)/3;                                                    %for ball diameter = 2
N=6; %number of sensors
anglesDeg = 360/N*(0:(N-1))';
positions=sidelength*[ cosd(anglesDeg) sind(anglesDeg)];


%% Moment of inertia
positions = positions - mean(positions,1);
MomentInertia = positions'*positions/N;
display(MomentInertia)
disp('(scaled identity for isotropic array)')
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
V=diag([2, 2*sqrt(3)/3] )/(sidelength/2);
LatticeBasis=V;
display(LatticeBasis)
fprintf('Area of parallelogram: %g\n',abs(det(LatticeBasis)))

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
exportgraphics(ax,'Figure1a.eps')
