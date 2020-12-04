%method of characteristics solution for a pipe with arbitrarily changing
%wall thickness

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
%e1=1/8*25.4e-3;%(m) pipe wall thickness
%e2=e1*0.1;%(m) pipe wall thickness

e_x=1/8*25.4e-3*[1 1 0.1 0.1 1 1];%(m)lookup table for wall thickness
x_e=L*[0 0.3 0.45 0.55 0.7 1];%(m) positions for lookup table




%% MOC params
N_cycles=5;%number of cycles to calculate
N_x=100;%number of x grid points
N_t=N_x*N_cycles*2;%number of time points


p_IC=0;%(Pa) initial pressure throughout
q_IC=0;%(m^3/s) initial flow throughout

p_BC=[1e6 0];%(Pa) pressure boundary conditions (nan if flow BC)
q_BC=[nan nan];%(m^3/s) flow boundary conditions (nan if pressure BC)
R_BC=[nan nan];%(Pa/(m^3/s) resistive bounary condition (nan if P or Q BC)


e=@(x) interp1(x_e,e_x,x,'pchip');%interp is slow, is there a faster method?
r=@(x) OD/2-e(x);
c=@(x) sqrt(K./rho./(1+alpha*2*K/E*r(x)./e(x)));%(m/s) wave speed function


%% solve MOC solution
%initialize, including finding x grid spacing (uneven if c is not constant)
[ x,t,Zc,c_bar ] = MOCinit( N_x,N_t, L, c, rho, r  );


tic
%solve
[ p, q, y ] =  MOCsolverR(x, t, p_IC, q_IC, p_BC, q_BC, R_BC, Zc, r, nu );
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
hold all
plot(x/L,zeros(size(x)),'kx')
hold off
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

figure(4)
subplot(3,1,1)
plot(x/L,r(x),'x-');
ylim([0 inf])
xlabel('x/L')
ylabel('r (m)')

subplot(3,1,2)
plot(x/L,c(x),'x-');
ylim([0 inf])
xlabel('x/L')
ylabel('c (m/s)')

subplot(3,1,3)
plot(x/L,Zc(x),'x-');
ylim([0 inf])
xlabel('x/L')
ylabel('Zc (Ps/(m^3/s))')
