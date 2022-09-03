%% 3D example
% create figure 2
% 2022
% Authors: Lee C. Potter

close all
clear all
clc

%% Create sensor positions for circular array
scale=2*diag([1,1/sqrt(3),1/sqrt(3)/2]);
N=5; %number of sensors
%tetrahedron
positions = [   0, 0 , sqrt(3);
                0, 0 , -sqrt(3);
                cosd(0)  ,sind(0)  ,0;
                cosd(120),sind(120),0;
                cosd(240),sind(240),0];
positions=positions*scale;


%% Moment of inertia
positions = positions - mean(positions,1);
MomentInertia = positions'*positions/N;
display(MomentInertia)
disp('(scaled identity for isotropic array)')
%% Figure for array geometry
figure;
f=plot3(positions(:,1),positions(:,2),positions(:,3),'ro',...
    'MarkerSize',12,'Linewidth',4);
grid;axis equal;hold on;
VizEdges1=[   ...
    0, 0 , sqrt(3);
    1, 0 ,  0;
    0, 0 , -sqrt(3);
    cosd(240),sind(240),0];
VizEdges2=[   ...
    cosd(240),sind(240),0;
    cosd(120),sind(120),0; 
    0, 0 , sqrt(3);
    cosd(240),sind(240),0;
    1, 0 ,  0];
HidEdges1=[   ...
    cosd(120),sind(120),0;
    0, 0 ,  -sqrt(3)];
HidEdges2=[   ...
    cosd(120),sind(120),0;
    1, 0 , 0];
VizEdges1=VizEdges1*scale;VizEdges2=VizEdges2*scale;
HidEdges1=HidEdges1*scale;HidEdges2=HidEdges2*scale;
plot3(VizEdges1(:,1),VizEdges1(:,2),VizEdges1(:,3),'k-','LineWidth',2)
plot3(VizEdges2(:,1),VizEdges2(:,2),VizEdges2(:,3),'k-','LineWidth',2)
plot3(HidEdges1(:,1),HidEdges1(:,2),HidEdges1(:,3),'k--','LineWidth',2)
plot3(HidEdges2(:,1),HidEdges2(:,2),HidEdges2(:,3),'k--','LineWidth',2)
view(15,22) ;
xticks([-1,0,2])
yticks([-1,0,1])
zticks([-1,0,1])
set(gca,'fontsize',18)
xlabel('x','fontsize',21)
ylabel('y','fontsize',21)
zlabel('z','fontsize',21)
%title('Array Geometry','fontsize',21)
ax = gca;
% Requires R2020a or later
exportgraphics(ax,'TetraGeom5.eps')


%% Construct pairwise differences 
%(omit the factor of pi)
Nrows = N-1;
P = zeros(Nrows,3);
for k=2:N
        P(k-1,:)=(positions(k,:)-positions(1,:));
end

%% Theoretical value rather than solve via rational approximation
V =2*eye(3);
LatticeBasis=V;
display(LatticeBasis)
fprintf('Volume of parallelogram: %g\n',abs(det(LatticeBasis)))
PtimesV=P*V;
display(PtimesV)
disp('Should be 0 mod 2')


%% Plot lattice and ball
figure;
plotcube([2 2 2],[0 0 0],.2,[0 1 0]); 
hold on;

quiver3(zeros(3,1),zeros(3,1),zeros(3,1),...
    V(:,1),V(:,2),V(:,3),'off',...
    'Color','k','LineWidth',4);

% Make base sphere data
radius=1;
[xb, yb, zb] = sphere(90);
surf(radius*xb+1, radius*yb+1, radius*zb+1, 'facecolor', 'b', 'edgealpha', 0);
% smooth and shaded
light;
lighting gouraud;
view(30,15) ;
axis off
text(2.15,0,0,'[2,0,0]','FontSize',21)
text(-0.25,0,2.15,'[0,0,2]','FontSize',21)
ax = gca;
% Requires R2020a or later
exportgraphics(ax,'TetraLattice5.eps')