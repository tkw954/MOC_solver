%method of characteristics solution, following Johnston et al "Use of Pipeline Wave
%Propagation Model for Measuring Unsteady Flow Rate" 2013

%% Liquid properties

nu = 100e-6; %(m^2/s) kinematic viscosity
rho = 870; %(kg/m^3) density
K = 1.5e9; %(Pa) bulk modulus
c=sqrt(K/rho); % (m/s) wave speed

%% Pipeline Dimensions

L=1000;%(m) pipe length
r=4*25.4e-3/2;%(m) pipe inner radius



%% MOC params
N_cycles=20;%number of cycles to calculate
N_x=100;%number of x grid points
N_t=N_x*N_cycles*2;%number of time points


p_IC=0;%(Pa) initial pressure throughout
q_IC=0;%(m^3/s) initial flow throughout

p_BC=[1e6 nan];%(Pa) pressure boundary conditions (nan if flow or RL BC)
q_BC=[nan 0];%(m^3/s) flow boundary conditions (nan if pressure or RL BC)
R_BC=[nan nan];%(Pa/(m^3/s) resistive bounary condition (nan if P or Q BC)


%% friction
%steady only
 n=0;
 m=0;


%% solve MOC solution
%initialize, including finding x grid spacing (uneven if c is not constant)
[ x,t,Zc,c_bar ] = MOCinit( N_x,N_t, L, c, rho, r  );


%solve
[ p, q, y ] =  MOCsolverR(x, t, p_IC, q_IC, p_BC, q_BC, R_BC, Zc, r, nu, n, m  );



figure(1)
pcolor(x/L,t/(2*L/c_bar),p*1e-6)
shading interp
h=colorbar;
ylabel(h,'P (MPa)')
xlabel('x/L')
ylabel('t/(2*L/c)')



figure(2)
plot(t/(2*L/c_bar),[q(:,1) q(:,end)]*60000)
xlabel('t/(2*L/c)')
ylabel('q (L/min)')
legend({'inlet','outlet'},'location','best')

figure(3)
idx=round(N_x/2);
plot(t/(2*L/c_bar),p(:,idx)*1e-6);
xlabel('t/(2*L/c)')
ylabel('p (MPa) midpoint')
