%method of characteristics solution, following Johnston et al "Use of Pipeline Wave
%Propagation Model for Measuring Unsteady Flow Rate" 2013

%% Liquid properties

nu = 100e-6; %(m^2/s) kinematic viscosity
rho = 870; %(kg/m^3) density
K = 1.5e9; %(Pa) bulk modulus

%% Pipeline Properties
% Elastic Modulus (Young's)
% Material: Steel
E =190e9; %Pa
nu_p = 0.3; %Poisson's ratio
% Axial effects due to Poisson's ratio (3 models) (Assumes inertial and acceleration effects
% of pipe negligable)
% 1 = Anchored upstream end only
% 2 = Anchored throughout to prevent axial movement
% 3 = Pipe with expansion joints throughout
% From: Ghidaoui et al. (2005)
axial_effect = 2; % Assume pipe does not move axially
if axial_effect == 1
    alpha = 1 - (nu_p/2);
elseif axial_effect == 2
    alpha = 1 - nu_p^2;
else
    alpha = 1;
end

%% Pipeline Dimensions

L=1000;%(m) pipe length
OD=8*25.4e-3;%(m) pipe outer diameter
e1=1/8*25.4e-3;%(m) pipe wall thickness
e2=e1*0.1;%(m) pipe wall thickness

r1=OD/2-e1;%(m) inner radius
r2=OD/2-e2;%(m) inner radius



%% MOC params
N_cycles=50;%number of cycles to calculate
N_x=100;%number of x grid points
N_t=N_x*N_cycles*2;%number of time points


p_IC=0;%(Pa) initial pressure throughout
q_IC=0;%(m^3/s) initial flow throughout

p_BC=[1e6 0];%(Pa) pressure boundary conditions (nan if flow or RL BC)
q_BC=[nan nan];%(m^3/s) flow boundary conditions (nan if pressure or RL BC)
%RL_BC=[nan nan];%(Pa/(m^3/s) resistive boundary condition (set below)

r=@(x) r1+(r2-r1)/L*x;%radius function
e=@(x) e1+(e2-e1)/L*x;%pipe wall thickness function
c=@(x) sqrt(K/rho./(1+alpha*2*K/E*r(x)./e(x)));%(m/s) wave speed function

%% friction
%steady only
% n=0;
% m=0;

%Trikha
% n=[26.4 200 8000];
% m=[1.0 8.1 40.0];

%Johnston 2006
beta_f=2;
m1=1.4064;
m2=2.5200;
n1=33.104;
r_approx=(r(0)+r(L))/2;
dt_approx=1.4e-3;
k=ceil((2*log(r_approx)-log(n1*nu*dt_approx))/(2*log(beta_f)));

m=nan(k,1);
n=nan(k,1);
m(1)=m1;
m(2)=m2;
n(1)=n1;

for i=3:k
    m(i)=beta_f*m(i-1);
end

for i=2:k
    n(i)=beta_f^2*n(i-1);
end


%% solve MOC solution
%initialize, including finding x grid spacing (uneven if c is not constant)
[ x,t,Zc,c_bar ] = MOCinit( N_x,N_t, L, c, rho, r  );

RL=0.5*Zc(L);
RL_BC=[nan RL];%set load resistance at exit, relative to Zc

tic
%solve
[ p, q, y ] =  MOCsolverR(x, t, p_IC, q_IC, p_BC, q_BC,RL_BC, Zc, r, nu, n, m  );
dt=toc;

%% steady resistance

R=@(x) 8*nu*rho/(pi*r(x)^4);%(Pa/(m^3/s)/m) static laminar resistance per unit length
R1=R(0)*L;%resistance if all pipe at r1
R2=R(L)*L;%resistance if all pipe at r2
q_ss1=(p(end,1)-p(end,end))/R1;%flow for R1
q_ss2=(p(end,1)-p(end,end))/R2;%flow for R2






fprintf('dt=%f s\n',dt)



figure(1)
pcolor(x/L,t/(2*L/c_bar),p*1e-6)
shading interp
h=colorbar;
ylabel(h,'P (MPa)')
xlabel('x/L')
ylabel('t/(2*L/c)')



figure(2)
plot(t/(2*L/c_bar),[q(:,1) q(:,end)]*60000)
hold all
plot(xlim,q_ss1*[1 1]*60000,'--')
plot(xlim,q_ss2*[1 1]*60000,'--')

hold off
xlabel('t/(2*L/c)')
ylabel('q (L/min)')
legend({'inlet','outlet'},'location','best')

figure(3)
idx=round(N_x/2);
plot(t/(2*L/c_bar),p(:,idx)*1e-6);
xlabel('t/(2*L/c)')
ylabel('p (MPa) midpoint')
