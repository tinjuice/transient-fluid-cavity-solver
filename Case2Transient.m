clear
close all

% User Inputs
param.nnx			=	6;			%Defines grid points in x : 2^nnx+1
param.nny			=	12;			%Defines grid points in y   WONT THE ELEMENTS BE TOO THIN?
param.Ly			=	0.4;		%(m) Domain Height
param.Lx			=	0.02;		%(m) Domain Width
param.dt			=	0.01;		%(s) Timestep 
param.tmax			=	50;			%(s) End time
param.Th			=	30;			%(C) Temperature on Hot Wall
param.Ti            =   20;         %(C) Initial ambient temperature

param.gravity		=	9.81;		%(m/s^2) Gravity
param.beta			=	2.59E-4;	%(1/K) Fluid Prperty, thermal expansion Coeff
param.visc			=	0.893E-6;	%(m^2/s) Fluid Prperty, kinematic viscosity
param.alpha 		= 	1.43E-7;    %(m^2/s) Fluid Prperty, thermal diffusivity
param.k				=	0.607;      %(W/m/K) thermal conductivity
param.gangle 		= 	-85;        %(deg) gravity angle from -y axis, +ve tilts g toward +x
 
gx_hat =  sind(param.gangle);      % unit gravity components (cavity frame)
gy_hat = -cosd(param.gangle);

Tmon 				= 	[0.005,param.Ly];	%monitor point (example)

%**** NS Solver Controls ***
param.Stream_tol	= 	1e-9; 		%Stream Function solver tolerance
param.Utime_tol 	= 	1e-6; 		%Velocity loop solver tolerance
param.wits_min		=	8;			%minimum number of iterations in vorticity solver
param.stream_its_min=	5;			%minimum number of iterations in stream function solver
param.mom_its_min 	= 	3;			%minimum number of iterations through momentum solver

gr	=	param.gravity*param.beta;

nx = 2^param.nnx + 1;
dx = param.Lx / (nx-1);
ny = 2^param.nny + 1;
dy = param.Ly / (ny-1);


mx(1:nx) = 0:dx:param.Lx;
my(1:ny) = 0:dy:param.Ly;
for ii = 1:nx
    for jj = 1:ny
	XX(ii,jj) = mx(ii);
	YY(ii,jj) = my(jj);
    end 
end

MonInd = [find(int16(mx./dx)==int16(Tmon(1)./dx)),find(int16(my./dy)==int16(Tmon(2)./dy))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise
stream 	= zeros(nx, ny);
w 		= zeros(nx, ny);
w_new 	= zeros(nx, ny);
wo		= zeros(nx, ny);
u 		= zeros(nx, ny);
u_new 	= zeros(nx, ny);
v 		= zeros(nx, ny);
v_new 	= zeros(nx, ny);
T 		= zeros(nx, ny);
To      = zeros(nx, ny);
T(:,:) 	= param.Ti;
To(:,:) = param.Ti;
buoy    = zeros(nx, ny);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary Conditions
T(1,1:ny)	=	param.Th;			%HOT WALL 		
% T(nx,1:ny)	= 	param.Tc;			%COLD WALL 
w_wall=0; 					%streamfunction value at the wall

To(1,1:ny) = param.Th;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xdv = param.visc*param.dt/dx^2;
ydv = param.visc*param.dt/dy^2; 
Mii = (1+ 2*param.visc*param.dt/dx.^2 + 2*param.visc*param.dt/dy.^2);
idtdx=param.dt./dx;
idtdy=param.dt./dy;
ajj=(2/dx^2 +2/dy^2);

timestep = 0;
momit=0;
t=0;

tic

while t < param.tmax
	t=t+param.dt;
    timestep = timestep + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Call temperature solver tn+1
%	Step 2.1: Your Job ....
d_x = param.alpha*param.dt/dx^2;
d_y = param.alpha*param.dt/dy^2;
for i=2:nx-1
    for j=2:ny-1
        CFL_x = u(i,j)*param.dt/dx;
        CFL_y = v(i,j)*param.dt/dy;

        % Gauss-Seidel
        T(i,j) = To(i,j) + d_x*(To(i+1,j) - 2*To(i,j) + To(i-1,j)) + d_y*(To(i,j+1) - 2*To(i,j) + To(i,j-1)) - 0.5*CFL_x*(To(i+1,j)-To(i-1,j)) - 0.5*CFL_y*(To(i,j+1)-To(i,j-1));

	    buoy(i,j) = gr*( -gy_hat*(T(i+1,j)-T(i-1,j))/(2*dx) ...
                   + gx_hat*(T(i,j+1)-T(i,j-1))/(2*dy) ); 	% T here is the newest temp value at Tn+1
    end
end
To = T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	Iterate on Momentum Solver
%	Step 2.2: You dont have to touch this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	du_it=0;
	delta_U = param.Utime_tol*10;
	while delta_U > param.Utime_tol |   du_it < param.mom_its_min
	momit = momit + 1;
	du_it = du_it+1;

	%	Enforce Boundary conditions for Vorticity
	%	Step 2.2.1
    w(:,1) 	=  		2 *(w_wall - stream(:, 2)) / (dy^2);
    w(1,:) 	=  		2 *(w_wall - stream(2, :) )/ (dx^2);
    w(nx,:) = 		2 *(w_wall - stream(nx-1, :)) / (dx^2);
    w(:,ny) = 		2 *(w_wall - stream(:, ny-1)) / (dy^2) ;

	%	Update Vorticity at n+1
	%	Step 2.2.2
	w_new	=	w;
	for wits = 1:param.wits_min 
		for ii = 2:nx-1
		for jj = 2:ny-1
			w_new(ii,jj)=	( w(ii+1,jj)+w(ii-1,jj) )*xdv/Mii + ...
							( w(ii,jj+1)+w(ii,jj-1) )*ydv/Mii + wo(ii,jj)/Mii  + buoy(ii,jj)*param.dt/Mii ...
							- 0.5*idtdx*u(ii,jj)*(w(ii+1,jj)-w(ii-1,jj))/Mii  + ...
							- 0.5*idtdy*v(ii,jj)*(w(ii,jj+1)-w(ii,jj-1))/Mii 	; %changed dt to param.dt
		end
		end
		w=w_new;
	end

	% 	Update Streamfunction to enforce continuity
	%	Step 2.2.3
    delta_stream = 1;
    sits = 0;
    stream_n = stream;
    while delta_stream > param.Stream_tol | sits < param.stream_its_min
        sits = sits + 1;
        stream_n(i,j) = ( w_new(i,j) + (stream(i+1, j) + stream(i-1, j))/dx^2 + (stream(i, j+1) + stream(i,j-1))/dy^2 )/ajj;

        delta_stream = sum(sum(abs(stream_n - stream)))/(nx*ny);

        stream = stream_n;
    end

    i = 2:nx-1;
    j = 2:ny-1;

    %   Update Velocity at n+1
	%	Step 2.2.4
    u_new(i, j) = (stream(i, j+1) - stream(i, j-1)) / (2 * dy);
    v_new(i, j) = (-stream(i+1, j) + stream(i-1, j)) / (2 * dx);

	%	Test for velocity change
	%	Step 2.2.5
	delta_U = sum(sum(abs(u_new - u)))/(nx*ny);
	u=u_new;
	v=v_new;

   	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Swap old and new values of Vorticity
    w = w_new;
    wo = w_new;

    contourf(XX(1:nx,1:ny),YY(1:nx,1:ny),T(1:nx,1:ny));
    axis equal,clim([0 param.Th]),colorbar
	set(gca,'fontsize',18,'linewidth',2);
    title(sprintf('t=%.4g s, Pe = %.4f',t, sqrt(max(u.^2+v.^2,[],'all'))*dx/param.alpha), 'FontSize', 15)
	drawnow			
	% if mod(timestep,400)==0
    % 
	% 	fprintf('Time Step %d; time %7.3e s; Max v %7.3e m/s; \n', timestep,t, max(v(:)));		
    % 
	% 	hFig = figure(1);
    % 
	% 	subplot(1,2,1);
	% 	contourf(XX(1:nx,1:ny),YY(1:nx,1:ny),T(1:nx,1:ny));
	% 	axis equal,clim([param.Tc param.Th]),colorbar
	% 	set(gca,'fontsize',18,'linewidth',2)
    % 
	% 	subplot(1,2,2);
	% 	contourf(XX(1:nx,1:ny),YY(1:nx,1:ny),v(1:nx,1:ny));
	% 	axis equal;
    %     colorbar
	% 	set(gca,'fontsize',18,'linewidth',2)
	% 	drawnow;	
	% end
end

elapsedTime = toc;
fprintf('Elapsed time: %.4f seconds\n', elapsedTime);

save('Case2_Cavity.mat','u','v','T','stream','w','param')
 
Tmonitor = T(MonInd(1),MonInd(2));
fileID = fopen('Dataout.txt','w');
fprintf(fileID,'%16.12e %16.12e %16.12e %16.12e e\n',[Tmonitor; dx ; dy; param.dt;]);
fclose(fileID); 

%%
hFig = figure(1);

subplot(1,2,1);
contourf(XX(1:nx,1:ny),YY(1:nx,1:ny),T(1:nx,1:ny));
axis equal,clim([0 param.Th]),colorbar
set(gca,'fontsize',18,'linewidth',2)
title('Temp field')

subplot(1,2,2);
contourf(XX(1:nx,1:ny),YY(1:nx,1:ny),v(1:nx,1:ny));
axis equal;
colorbar
set(gca,'fontsize',18,'linewidth',2)
title('v field')
drawnow;	

CFL = sqrt(max(u.^2+v.^2,[],'all'))*param.dt/dx;
Peclet = sqrt(max(u.^2+v.^2,[],'all'))*dx/param.alpha;
fprintf('Final Pe = %.4f\nFinal CFL = %.4f\n', Peclet, CFL)
