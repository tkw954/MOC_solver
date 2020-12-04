function [ p, q, y ] = MOCsolverR_vectorize(x, t, p_IC, q_IC, p_BC, q_BC, R_BC, Zc, r, nu, n, m  )
%[ p, q ] = MOCsolverF(x, t, p_IC, q_IC, p_BC, q_BC, R_BC, Zc, r, nu, n, m  )
%   Method of characteristics pipeline dynamics solver, following Johnston,
%   2006 and 2014
%  This function allows for resistive boundary conditions
% x,t: spatial and temporal grid, output from MOCinit
% p_IC, q_IC: pressure and flow initial conditions. May be scalar or vector size(x)
% p_BC, q_BC, R_BC: pressure, flow or resistance boundary conditions. May be 1 x 2 or
%  N_t x 2 matrices. Nan for the unknown BC
% Zc: characteristic impedance, f(x). Output from MOCinit
% r: inner radius, f(x)
% nu: kinematic viscosity
% n,m: friction parameters as defined by Johnston 2006. n=0 m=0 for steady
%  friction only. Leave blank for frequency-dependent friction.
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
if nargin<11||isempty(n)||isempty(m)
    %default to Johnston 2006 method
        beta_f=2;
        m1=1.4064;
        m2=2.5200;
        n1=33.104;
        r_max=max(r(x));
        dt=t(2)-t(1);
        k=ceil((2*log(r_max)-log(n1*nu*dt))/(2*log(beta_f)));
        k=max(k,3);
        
        m=nan(k,1);
        n=nan(k,1);
        m(1)=m1;
        m(2)=m2;
        n(1)=n1;
        
        for i_f=3:k
            m(i_f)=beta_f*m(i_f-1);
        end
        
        for i_f=2:k
            n(i_f)=beta_f^2*n(i_f-1);
        end
end

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

% if ~isa(r, 'function_handle') %calculate at grid rather than function
%     r=@(x) r*ones(size(x));
% end



%% Allocate matrices

p=nan(N_t,N_x);
q=nan(N_t,N_x);
k=numel(n);
y=nan(N_x,k,N_t);

%% initial conditions
p(1,:)=p_IC(:);
q(1,:)=q_IC(:);
y(:,:,1)=0;

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

if isa(r, 'function_handle')
    r=r(x);
elseif numel(r)==1
    r=r*ones(1,N_x);
end

tic
for k=2:N_t
    dt=t(k)-t(k-1);%time step, should be constant
    
    % friction
    if k>2
        q_old=q(k-2,:);
    else
        q_old=q(k-1,:);
    end
    
    [f, y_new]=unsteady_update(q(k-1,:),q_old, y(:,:,k-1),n,m,nu,r,dt);
    y(:,:,k)=y_new;
    
    %% inlet BC
    i=1;
    
    ZB=Zc(i);
    pB=p(k-1,i+1);
    qB=q(k-1,i+1);
    fB=f(i+1);
    
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
    
    %% interior points
    ZA=Zc(1:N_x-2);
    ZB=Zc(2:N_x-1);
    pA=p(k-1,1:N_x-2);
    pB=p(k-1,3:N_x);
    qA=q(k-1,1:N_x-2);
    qB=q(k-1,3:N_x);
    
    fA=f(1:N_x-2);
    fB=f(3:N_x);
    
    psi=1./(1+ZA./ZB);
    
    p(k,2:(N_x-1))=psi.*(pA+ZA./ZB.*pB+ZA.*(qA-qB)+dt.*ZA.*(fB-fA));%eq 12
    q(k,2:(N_x-1))=(1-psi).*qA+psi.*qB+1./ZA.*(1-psi).*pA-1./ZB.*psi.*pB-dt.*((1-psi).*fA+psi.*fB);%eq 13
    
    
    
    %% exit boundary
    i=N_x;
    
    ZA=Zc(i-1);
    pA=p(k-1,i-1);
    qA=q(k-1,i-1);
    fA=f(i-1);
    
    
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
    
end


end


function [f, y_new, f_unsteady]=unsteady_update(q_new,q_old, y_old,n,m,nu,r,dt)
k=numel(n);

N_x=numel(r);

y_new=nan(N_x,k);
for i=1:k %could vecotrize this but probably not much speedup possible as k is small
    y_new(:,i)=y_old(:,i).*exp(-n(i)*nu*dt./r(:).^2)+m(i)*(q_new(:)-q_old(:)).*exp(-n(i)*nu*dt./(2*r(:).^2));
end

f_unsteady=4*nu./r.^2.*sum(y_new,2)';

f=8.*nu./r.^2.*q_new+f_unsteady;

end
