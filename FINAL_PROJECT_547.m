clearvars
N=4;
%% Input from the user
prompt="Please input desired position as x y z ";
inp=input(prompt,'s');
st=split(inp,' ');
dub=str2double(st);
% clearAllMemoizedCaches

%% Dh parameters initialized
d=[0 0 0.15005 0.4318 ];
alpha=[pi/2 0 -pi/2 pi/2];
a=[0 0.4318 0.0203 0];
% for i=1:N 
L1=Link('revolute','d', 0,'a', 0,'alpha', pi/2);
L2=Link('revolute','d', 0,'a', 0.4318,'alpha',0);
L3=Link('revolute','d', 0.15005,'a', 0.0203,'alpha', -pi/2);
L4=Link('revolute','d', 0.4318,'a',0,'alpha', pi/2);
L=[L1 L2 L3 L4];
% end
R = SerialLink(L,'name','4 Link Arm');
goal_pos=[dub(1),dub(2),dub(3)];
% present the workspace for the user
workspaceforproj(goal_pos)
% hold on
flag=false;
prompt="Check if point is in the workspace of Robot. Do you want to continue with this ?(y/n) ";

inp=input(prompt,'s');
if prompt=='n'
  fprintf('Thank you! Please restart the software')
end

q_initial=[0 0 0 0];
start_pos=forward_kinematics(q_initial(1), q_initial(2), q_initial(3), q_initial(4) ); 
start_pos=start_pos(1:3)';
% rrt star planning with an obstacle is done
path=rrt_star_plan_obs(start_pos, goal_pos);
% spline interpolation done
way_points=cscvn(path(:,[1:end]));
num_of_points=400;

syms t
x1=[];
y1=[];
z1=[];
vx=[];
vy=[];
vz=[];
%Finding the equations for trajectory from polynomial coefficients
for iter=1:length(path)-1
    t1=linspace(way_points.breaks(iter), way_points.breaks(iter+1), num_of_points);
    eqn_x=way_points.coefs((iter-1)*3+1,4)+(way_points.coefs((iter-1)*3+1, 3)*(t-way_points.breaks(iter)))+(way_points.coefs((iter-1)*3+1, 2)*((t-way_points.breaks(iter))^2))+(way_points.coefs((iter-1)*3+1, 1)*(t-way_points.breaks(iter))^3);
    eqn_y=way_points.coefs((iter-1)*3+2,4)+(way_points.coefs((iter-1)*3+2, 3)*(t-way_points.breaks(iter)))+(way_points.coefs((iter-1)*3+2, 2)*((t-way_points.breaks(iter))^2))+(way_points.coefs((iter-1)*3+2, 1)*(t-way_points.breaks(iter))^3);   
    eqn_z=way_points.coefs((iter-1)*3+3,4)+(way_points.coefs((iter-1)*3+3, 3)*(t-way_points.breaks(iter)))+(way_points.coefs((iter-1)*3+3, 2)*((t-way_points.breaks(iter))^2))+(way_points.coefs((iter-1)*3+3, 1)*(t-way_points.breaks(iter))^3);   
    vel_x=diff(eqn_x,t);
    vel_y=diff(eqn_y,t);
    vel_z=diff(eqn_z,t);
    xxx=subs(eqn_x,t,t1);
    yyy=subs(eqn_y,t,t1);
    zzz=subs(eqn_z,t,t1);
    x1=[x1 xxx];
    y1 =[y1 yyy];
    z1= [z1 zzz];
     xxx=subs(vel_x,t,t1);
    yyy=subs(vel_y,t,t1);
    zzz=subs(vel_z,t,t1);
    x1=double(x1)
    y1=double(y1)
    z1=double(z1)
    vx=double([vx xxx]);
    vy =double([vy yyy]);
    vz= double([vz zzz]);
end
clear t
% Using time parameter to come up with discretized trajectory.
tot_points=num_of_points*(length(path)-1);
t=linspace(0, way_points.breaks(length(path)), tot_points);
dt=way_points.breaks(length(path))/tot_points;
%% Code for finding Transformation matrix and Jacobian
% T=eye(4);
% phi=0;
% t=sym('theta', N);
% 
% for iter=1:N
% 
%     fprintf('hello')
%     T_temp=dh_trans(t(iter),d(iter),alpha(iter), a(iter));
%     phi=phi+t(iter);
% 
%     T=T*T_temp;
%     
% end
% 
% pe=[T(1:3,4); phi];
% J_a=[];
% for iter=1:N
%     j=diff(pe, t(iter));
%     J_a=[J_a j];
% end

%% Do inverse kinematics on trajectory calculated
%now using existing trajectory
% dt=0.01;
% t=linspace(0,5,5/dt);
l1=length(x1);

q_final=[];
pd=[x1',y1',z1',(t*pi/10)'];%+5*pi/12
pd(1:length(t),:)=pd(length(t):-1:1,:)
% pd_dot=[-0.04*pi*sin(t*pi/10)', 0.04*pi*cos(pi*t/10)',0.1*cos(t)', (ones(1,length(t))*pi/10)'];
pd_dot=[vx',vy',vz',(ones(1,length(t))*pi/10)']
pd_dot(1:length(t),:)=pd_dot(length(t):-1:1,:);
I = eye( 3 );
% I = eye( 4 );
k=30*I; %Proportional gain

ki=0.05*I; %Integral  gain
q=[0, 0, 0, 0]';
actual_pos_x=zeros(length(t), 3);
e1=zeros(1,4);
t_iter_sum=0;
e_sum=0;
e_norm_fin=[]
for i= 1:length(t)
%     break
    pd1=pd(i,1:3)';
    xd_dot=pd_dot(i,1:3)';
    

    
    for iter=1:1000
        
        pe1=forward_kinematics(q(1),q(2),q(3),q(4)); 
        
        pe1=pe1(1:3,1);
    %     e=pd1-pe1(1:3);
        e=pd1-pe1;
        e_sum=(e_sum+e);
        e_norm=norm(e);
        if e_norm<0.1 %check if the norm of error is less than 0.1
            break
        end
        J1=analytical_jacobian_val(q(1),q(2), q(3), q(4));
        
        %the new state is determined based on PI terms and iterative
        %process
    
        dq_dt=pinv(J1)*(k*e + ki*e_sum);%+kd*e_dot +ki*e_sum);%+ xd_dot)%kd*e_dot);%+k*e +ki*e_sum);
    %     q
        q=q+dq_dt*dt
    end
    e_norm_fin=[e_norm_fin e_norm];
    t_iter_sum=t_iter_sum+dt;
    q_final=[q_final,q];
    
    actual_pos_x(i,:)=pe1';        
    fprintf(' %g th point reached!',i)
    
end
subplot(2,2,1)
plot(t,pd(:,1),t,actual_pos_x(:,1));%plotting error in x
xlabel('time step')
legend('Desired position', 'actual position ')
ylabel('position in x')
xlim([0 2])
ylim([-0.4 2])
subplot(2,2,2)
plot(t,pd(:,2),t,actual_pos_x(:,2)); %plotting error in y
xlabel('time step')
legend('Desired position', 'actual position ')
ylabel('position in y')
xlim([0 2])
ylim([-0.4 2])
subplot(2,2,3)
plot(t,pd(:,3),t,actual_pos_x(:,3));%plotting error in z
xlabel('time step')
legend('Desired position', 'actual position ')
ylabel('position in z')
xlim([0 2])
ylim([-0.4 2])
subplot(2,2,4)% plot(t,pd(:,4),t,actual_pos_x(:,4));
plot(t,e_norm_fin); %plotting norm of e

figure()
% f=subs(T,{theta1,theta2,d3,theta3},{pi/3,pi/2, 0.2,pi/3})
%for animation
% for i = 1:10
plot3(x1,y1,z1,'o','color','k')
% plot3(pd(:,1)',pd(:,2)',pd(:,3)','o','color','k')
plot(R,q_final','workspace',[-1 1 -1 1 -1 1], 'delay',0.01);
% hold on
pause(1)
theta1=q(1)
theta2=q(2)
v=sym('v',2);
ve=v(:,1);
L1=1;
L2=1;
J_a=double(analytical_jacobian_val(theta1, theta2, 0,0))
%% Finding velocity ellipsoid for final joint position.

A1=inv(J_a*J_a')
[x,y,z]=Ellipse_plot(A1,[0,0,0]');

hold on
surf(x,y,z)
% end
function J= analytical_jacobian_val(theta1, theta2, d3, theta3)
    x1=[- 0.5000*sin(theta1) - 0.3000*cos(theta1)*sin(theta2) - 0.3000*cos(theta2)*sin(theta1), - 0.3000*cos(theta1)*sin(theta2) - 0.3000*cos(theta2)*sin(theta1), 0, 0];
    x2=[0.5000*cos(theta1) + 0.3000*cos(theta1)*cos(theta2) - 0.3000*sin(theta1)*sin(theta2),   0.3000*cos(theta1)*cos(theta2) - 0.3000*sin(theta1)*sin(theta2), 0, 0];
    x3=[0, 0, 1, 0];
%     x4=[1,1, 0, 1];
    J=[x1;x2;x3];
end
function xyz_path=rrt_star_plan_obs(start_pos, goal_pos ) 
close all
x_max = 1;
y_max = 1;
z_max = 1;
EPS = 20;
numNodes = 500;        
q_start.coord = start_pos;
q_start.cost = 0;
q_start.parent = 0;
q_goal.coord = goal_pos;
q_goal.cost = 0;
nodes(1) = q_start;
figure(1)
xlabel('x'),ylabel('y'),zlabel('z')

%hold on
obstacle1 = [0.5,0.2,0.2] ;   % origin point 
obstacle2 = [0.5,0.5,0.5] ;  % your cube dimensions 
%O = P-L/2 ;% Get the origin of cube so that P is at center 
plotcube(obstacle1,obstacle2,1,[0 0 1]);% use function plotcube 
hold on
for i = 1:1:numNodes
    q_rand = [rand(1)*x_max rand(1)*y_max rand(1)*z_max];
    %plot3(q_rand(1), q_rand(2), q_rand(3), 'x', 'Color',  [0 0.4470 0.7410])
    
    % Break if goal node is already reached
    for j = 1:1:length(nodes)
        if nodes(j).coord == q_goal.coord
            break
        end
    end
    
    % Pick the closest node from existing list to branch out from
    ndist = [];
    
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = dist_3d(n.coord, q_rand);
        ndist = [ndist tmp];
    end
    [val, idx] = min(ndist);
    q_near = nodes(idx);
    
    q_new.coord = steer3d(q_rand, q_near.coord, val, EPS);
  if nooCollision(q_rand, q_near.coord)%, obstacle1,obstacle2) %&& noCollision(q_rand, q_near.coord,obstacle2) %&& noCollision(q_rand, q_near.coord, obstacle3) && noCollision(q_rand, q_near.coord, obstacle4)
    line([q_near.coord(1), q_new.coord(1)], [q_near.coord(2), q_new.coord(2)], [q_near.coord(3), q_new.coord(3)], 'Color', 'k', 'LineWidth', 2);
    drawnow
    hold on
    q_new.cost = dist_3d(q_new.coord, q_near.coord) + q_near.cost;
    
    % Within a radius of r, find all existing nodes
    q_nearest = [];
    r = 90;
    neighbor_count = 1;
    for j = 1:1:length(nodes)
        if nooCollision(nodes(j).coord, q_new.coord) && (dist_3d(nodes(j).coord, q_new.coord)) <= r
            %noCollision(nodes(j).coord, q_new.coord, obstacle2) && noCollision(nodes(j).coord, q_new.coord, obstacle3) && noCollision(nodes(j).coord, q_new.coord, obstacle4) && (dist_3d(nodes(j).coord, q_new.coord)) <= r
            q_nearest(neighbor_count).coord = nodes(j).coord;
            q_nearest(neighbor_count).cost = nodes(j).cost;
            neighbor_count = neighbor_count+1;
        end
    end
 
    
    
    % Initialize cost to currently known value
    q_min = q_near;
    C_min = q_new.cost;
    
    % Iterate through all nearest neighbors to find alternate lower
    % cost paths
    
    for k = 1:1:length(q_nearest)
        if  nooCollision(q_nearest(k).coord, q_new.coord) && q_nearest(k).cost + dist_3d(q_nearest(k).coord, q_new.coord) < C_min
            %&& noCollision(q_nearest(k).coord, q_new.coord, obstacle2) && noCollision(q_nearest(k).coord, q_new.coord, obstacle3) && noCollision(q_nearest(k).coord, q_new.coord, obstacle4) && q_nearest(k).cost + dist_3d(q_nearest(k).coord, q_new.coord) < C_min
            q_min = q_nearest(k);
            C_min = q_nearest(k).cost + dist_3d(q_nearest(k).coord, q_new.coord);
            %line([q_min.coord(1), q_new.coord(1)], [q_min.coord(2), q_new.coord(2)], [q_min.coord(3), q_new.coord(3)], 'Color', 'g');            
            hold on
        end
    end
    
    % Update parent to least cost-from node
    for j = 1:1:length(nodes)
        if nodes(j).coord == q_min.coord
            q_new.parent = j;
        end
    end
    
    % Append to nodes
    nodes = [nodes q_new];
  end
end
D = [];
for j = 1:1:length(nodes)
    tmpdist = dist_3d(nodes(j).coord, q_goal.coord);
    D = [D tmpdist];
end
% Search backwards from goal to start to find the optimal least cost path
[val, idx] = min(D);
q_final = nodes(idx);
q_goal.parent = idx;
q_end = q_goal;
nodes = [nodes q_goal];
xyz_path=[q_end.coord'];
while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)], [q_end.coord(3), nodes(start).coord(3)], 'Color', 'r', 'LineWidth', 4);
    hold on
    xyz_path=[xyz_path nodes(start).coord'];
%     q_end.coord
    q_end = nodes(start);
end
end

function workspaceforproj(pos)
syms t1;
syms t2
syms t3
a1=0.3;
a2=0.5;
a3=0.2;
L1=Link('revolute','d', 0,'a', 0,'alpha', pi/2);
L2=Link('revolute','d', 0,'a', 0.4318,'alpha',0);
L3=Link('revolute','d', 0.15005,'a', 0.0203,'alpha', -pi/2);
L4=Link('revolute','d', 0.4318,'a',0,'alpha', pi/2);
L=[L1 L2 L3 L4];
% end
R = SerialLink(L,'name','4 Link Arm');
figure
plot(R,[1,1,1,1]);
t_1=linspace(-120,120,15)*pi/180;
t_2=linspace(-120,120,15)*pi/180;
t_3=linspace(-120,120,15)*pi/180;
t_4=linspace(-120,120,15)*pi/180;
[theta1,theta2,theta3,theta4]=ndgrid(t_1,t_2,t_3, t_4);  % This will create matrices of 100x100x100 for each variable
xM = 0.1501*sin(theta1) + 0.4318*cos(theta1).*cos(theta2) + 0.0203*cos(theta3).*(cos(theta1).*cos(theta2) - 6.1232e-17*sin(theta1).*sin(theta2)) - 0.4318*cos(theta3).*(cos(theta1).*sin(theta2) + 6.1232e-17*cos(theta2).*sin(theta1)) - 0.4318*sin(theta3).*(cos(theta1).*cos(theta2) - 6.1232e-17*sin(theta1).*sin(theta2)) - 0.0203*sin(theta3).*(cos(theta1).*sin(theta2) + 6.1232e-17*cos(theta2).*sin(theta1)) - 2.6440e-17*sin(theta1).*sin(theta2);
yM = 0.0203*cos(theta3).*(6.1232e-17*cos(theta1).*sin(theta2)) + cos(theta2).*sin(theta1) - 0.1501*cos(theta1) + 0.4318*cos(theta3).*(6.1232e-17*cos(theta1).*cos(theta2) - sin(theta1).*sin(theta2)) + 2.6440e-17*cos(theta1).*sin(theta2) + 0.4318*cos(theta2).*sin(theta1) - 0.4318*sin(theta3).*(6.1232e-17*cos(theta1).*sin(theta2) + cos(theta2).*sin(theta1)) + 0.0203*sin(theta3).*(6.1232e-17*cos(theta1).*cos(theta2) - sin(theta1).*sin(theta2));
zM = 0.4318*sin(theta2) + 0.4318*cos(theta2).*cos(theta3) + 0.0203*cos(theta2).*sin(theta3) + 0.0203*cos(theta3).*sin(theta2) - 0.4318*sin(theta2).*sin(theta3) + 9.1879e-18;
hold on
plot3(xM(:),yM(:),zM(:),'g') % This is the plot type you should be using.
plot3(pos(1),pos(2),pos(3),'r*','MarkerSize', 50,'LineWidth',2)
hold off
end

function P= forward_kinematics(theta1, theta2, theta3, theta4)
    x1=0.1501*sin(theta1) + 0.4318*cos(theta1)*cos(theta2) + 0.0203*cos(theta3)*(cos(theta1)*cos(theta2) - 6.1232e-17*sin(theta1)*sin(theta2)) - 0.4318*cos(theta3)*(cos(theta1)*sin(theta2) + 6.1232e-17*cos(theta2)*sin(theta1)) - 0.4318*sin(theta3)*(cos(theta1)*cos(theta2) - 6.1232e-17*sin(theta1)*sin(theta2)) - 0.0203*sin(theta3)*(cos(theta1)*sin(theta2) + 6.1232e-17*cos(theta2)*sin(theta1)) - 2.6440e-17*sin(theta1)*sin(theta2);
    x2=0.0203*cos(theta3)*(6.1232e-17*cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) - 0.1501*cos(theta1) + 0.4318*cos(theta3)*(6.1232e-17*cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + 2.6440e-17*cos(theta1)*sin(theta2) + 0.4318*cos(theta2)*sin(theta1) - 0.4318*sin(theta3)*(6.1232e-17*cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + 0.0203*sin(theta3)*(6.1232e-17*cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2));
    x3=0.4318*sin(theta2) + 0.4318*cos(theta2)*cos(theta3) + 0.0203*cos(theta2)*sin(theta3) + 0.0203*cos(theta3)*sin(theta2) - 0.4318*sin(theta2)*sin(theta3) + 9.1879e-18;
    x4=theta1 + theta2 + theta1 + theta4;
    P=[x1;x2;x3;x4];
end

function [XX,YY,ZZ]=Ellipse_plot(A, C, varargin)

N = 20; % Default value for grid

% See if the user wants a different value for N.
%----------------------------------------------
if nargin > 2
 	N = varargin{1};
end


% check the dimension of the inputs: 2D or 3D
%--------------------------------------------
if length(C) == 3,
    Type = '3D';
elseif length(C) == 2,
    Type = '2D';
else
    display('Cannot plot an ellipse with more than 3 dimensions!!');
    return
end

% "singular value decomposition" to extract the orientation and the
% axes of the ellipsoid
[U D V] = svd(A);

if strcmp(Type, '2D'),
    % get the major and minor axes
    %------------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));

    theta = [0:1/N:2*pi+1/N];

    % Parametric equation of the ellipse
    %----------------------------------------
    state(1,:) = a*cos(theta); 
    state(2,:) = b*sin(theta);

    % Coordinate transform 
    %----------------------------------------
    X = V * state;
    X(1,:) = X(1,:) + C(1);
    X(2,:) = X(2,:) + C(2);
    
elseif strcmp(Type,'3D'),
    % generate the ellipsoid at (0,0,0)
    %----------------------------------
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));
    c = 1/sqrt(D(3,3));
    [X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);
    
    %  rotate and center the ellipsoid to the actual center point
    %------------------------------------------------------------
    XX = zeros(N+1,N+1);
    YY = zeros(N+1,N+1);
    ZZ = zeros(N+1,N+1);
    for k = 1:length(X),
        for j = 1:length(X),
            point = [X(k,j) Y(k,j) Z(k,j)]';
            P = V * point;
            XX(k,j) = P(1)+C(1);
            YY(k,j) = P(2)+C(2);
            ZZ(k,j) = P(3)+C(3);
        end
    end
end
end


