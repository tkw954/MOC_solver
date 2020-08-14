function [ p, q, y ] = MOCsolverR(x, t, p_IC, q_IC, p_BC, q_BC, R_BC, Zc, r, nu, n, m  )
%[ p, q ] = MOCsolverF(x, t, p_IC, q_IC, p_BC, q_BC, R_BC, Zc, r, nu, n, m  )
%   Method of characteristics pipeline dynamics solver, following Johnston,
%   2006 and 2014
%  This function allows for resistive boundary conditions
% x,t: spatial and temporal grid, output from MOCinit
% p_IC, q_IC: pressure and flow initial conditions. May be scalar or vector size(x)
% p_BC, q_BC, R_BC: pressure, flow or resistance boundary conditions. May be 1 x 2 or 
%  N_t x 2 matrices. Nan for the unknown BC
% Zc: characteristic impedance, f(x). Output from MOCinit
% r: iner radius, f(x)
% nu: kinematic viscosity
% n,m: friction parameters as defined by Johnston 2006. n=0 m=0 for steady
%  friction only
%
% Returns:
% p: pressure
% q: flow
%
% Reference: 
% Johnston, N., 2006, Efficient Methods for Numerical Modeling of Laminar 
%  Friction in Fluid Lines, Proc. IMechE Vol. 220 Part I: J. Systems 
%  and Control Engineering
% Johnston, N. 2014, Use of Pipeline Wave Propagation Model for
%  Measuring Unsteady Flow Rate, ASME Journal of Fluids Engineering



N_x=numel(x);
N_t=numel(t);

%% check inputs
if numel(p_IC)==1
    p_IC=p_IC*ones(size(x));
end

if numel(q_IC)==1
    q_IC=q_IC*ones(size(x));
end

if numel(p_IC)~=N_x
    error('x and p initial conditions must be the same size')
end

if numel(q_IC)~=N_x
    error('x and q initial conditions must be the same size')
end

if isempty(q_BC)
    q_BC=[nan nan];
end

if isempty(p_BC)
    p_BC=[nan nan];
end

if isempty(R_BC)
    R_BC=[nan nan];
end

if numel(q_BC)==2
    q_BC=repmat(q_BC(:)',N_t,1);
end

if numel(p_BC)==2
    p_BC=repmat(p_BC(:)',N_t,1);
end

if numel(R_BC)==2
    R_BC=repmat(R_BC(:)',N_t,1);
end

if ~isa(r, 'function_handle')
    r=@(x) r*ones(size(x));
end

%% Allocate matrices

p=nan(N_t,N_x);
q=nan(N_t,N_x);
k=numel(n);
y=nan(N_t,N_x,k);

%% initial conditions
p(1,:)=p_IC(:);
q(1,:)=q_IC(:);
y(1,:,:)=0;

%% boundary conditions
idx=isfinite(p_BC(:,1));
p(idx,1)=p_BC(idx,1);

idx=isfinite(p_BC(:,2));
p(idx,end)=p_BC(idx,2);

idx=isfinite(q_BC(:,1));
q(idx,1)=q_BC(idx,1);

idx=isfinite(q_BC(:,2));
q(idx,end)=q_BC(idx,2);



%grid at midpoint
x_midgrid=(x(1:(end-1))+x(2:end))/2;


% Get Zc at element midpoint
if isa(Zc, 'function_handle')
    Zc=Zc(x_midgrid);
elseif numel(Zc)==1
    Zc=Zc*ones(1,N_x-1);
end


tic
for k=2:N_t
    dt=t(k)-t(k-1);%time step, should be constant
    %fprintf('%d/%d\n',k,N_t)
    i=1;
    
    xB=x_midgrid(i);
    ZB=Zc(i);
    pB=p(k-1,i+1);
    qB=q(k-1,i+1);
    yB=y(k-1,i+1,:);
    
    if k>2
        qBold=q(k-2,i+1);
    else
        qBold=qB;
    end
    
    [fB, yB]=unsteady_update(qB,qBold, yB,n,m,nu,r(xB),dt);
    
    if isfinite(p(k,i))
        pP=p(k,i);
        q(k,i)=qB+1/ZB*(pP-pB)-dt*fB;%eq 6
    elseif isfinite(q(k,i))
        qP=q(k,i);
        p(k,i)=pB+ZB*(qP-qB+dt*fB);
    elseif isfinite(R_BC(k,1))
        
        q(k,i)=(qB-1/ZB*pB-dt*fB)/(1+R_BC(k,1)/ZB);%eq 6 with pP=-R_BC*qP
        p(k,i)=-R_BC(k,1)*q(k,i);
    else
        error('Missing boundary condition')
    end
    y(k,i+1,:)=yB;
    
    for i=2:(N_x-1)
        
        xA=x_midgrid(i-1);
        xB=x_midgrid(i);
        ZA=Zc(i-1);
        ZB=Zc(i);
        pA=p(k-1,i-1);
        pB=p(k-1,i+1);
        qA=q(k-1,i-1);
        qB=q(k-1,i+1);
        yA=y(k-1,i-1,:);
        yB=y(k-1,i+1,:);
        if k>2
            qAold=q(k-2,i-1);
            qBold=q(k-2,i+1);
        else
            qAold=qA;
            qBold=qB;
        end
        
        psi=1/(1+ZA/ZB);
        
        [fA, yA]=unsteady_update(qA,qAold, yA,n,m,nu,r(xA),dt);
        [fB, yB]=unsteady_update(qB,qBold, yB,n,m,nu,r(xB),dt);
        
        p(k,i)=psi*(pA+ZA/ZB*pB+ZA*(qA-qB)+dt*ZA*(fB-fA));%eq 12
        
        q(k,i)=(1-psi)*qA+psi*qB+1/ZA*(1-psi)*pA-1/ZB*psi*pB-dt*((1-psi)*fA+psi*fB);%eq 13
        
        y(k,i-1,:)=yA;
        y(k,i+1,:)=yB;
    end
    
    i=N_x;
    xA=x_midgrid(i-1);
    % xB=x_midgrid(i);
    ZA=Zc(i-1);
    %ZB=Zc(i);
    pA=p(k-1,i-1);
    %pB=p(k-1,i+1);
    qA=q(k-1,i-1);
    %qB=q(k-1,i+1);
    yA=y(k-1,i-1,:);
    
    if k>2
        qAold=q(k-2,i-1);
    else
        qAold=qA;
    end
    
    [fA, yA]=unsteady_update(qA,qAold, yA,n,m,nu,r(xA),dt);
    
    if isfinite(p(k,i))
        
        pP=p(k,i);
        q(k,i)=qA-1/ZA*(pP-pA)-dt*fA;%eq 5
    elseif isfinite(q(k,i))
        qP=q(k,i);
        p(k,i)=pA-ZA*(qP-qA+dt*fA);
    elseif isfinite(R_BC(k,2))
        
        q(k,i)=(qA+1/ZA*(pA)-dt*fA)/(1+R_BC(k,2)/ZA);%eq 5 with Pp=R*qp
        p(k,i)=R_BC(k,2)*q(k,i);
    else
        error('Missing boundary condition')
    end
    
    y(k,i-1,:)=yA;
end



end


function [f, y_new, f_unsteady]=unsteady_update(q_new,q_old, y_old,n,m,nu,r,dt)
k=numel(n);

y_new=nan(size(y_old));
for i=1:k
    y_new(i)=y_old(i).*exp(-n(i)*nu*dt/r^2)+m(i)*(q_new-q_old)*exp(-n(i)*nu*dt/(2*r^2));
end

f_unsteady=4*nu/r^2*sum(y_new);

f=8*nu/r^2*q_new+f_unsteady;
end
