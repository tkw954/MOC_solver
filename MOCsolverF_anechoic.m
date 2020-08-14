function [ p, q, y ] = MOCsolverF_anechoic(x, t, p_IC, q_IC, p_BC, q_BC, Zc, r, nu, n, m  )
%% solver for anechoic BC



N_x=numel(x);
N_t=numel(t);

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

if numel(q_BC)==2
    q_BC=repmat(q_BC(:)',N_t,1);
end

if numel(p_BC)==2
    p_BC=repmat(p_BC(:)',N_t,1);
end

if ~isa(r, 'function_handle')
    r=@(x) r*ones(size(x));
end

% if (isfinite(p_BC(1)) && isfinite(q_BC(1))) || ...
%     (isfinite(p_BC(2)) && isfinite(q_BC(2)))
%     error('Boundary conditions must be p OR q')
% end



p=nan(N_t,N_x);
q=nan(N_t,N_x);
k=numel(n);
y=nan(N_t,N_x,k);

%initial conditions
p(1,:)=p_IC(:);
q(1,:)=q_IC(:);
y(1,:,:)=0;

%boundary conditions
idx=isfinite(p_BC(:,1));
p(idx,1)=p_BC(idx,1);

idx=isfinite(p_BC(:,2));
p(idx,end)=p_BC(idx,2);

idx=isfinite(q_BC(:,1));
q(idx,1)=q_BC(idx,1);

idx=isfinite(q_BC(:,2));
q(idx,end)=q_BC(idx,2);

%grid
x_midgrid=(x(1:(end-1))+x(2:end))/2;

%friction
% if nargin(f)==1
%     f=@(q,x) f(q);
% end
%f=@(q,x) 8*nu/r(x)^2*q;%steady friction model (eq 7)



% Zc
if isa(Zc, 'function_handle')
    Zc=Zc(x_midgrid);
elseif numel(Zc)==1
    Zc=Zc*ones(1,N_x-1);
end


tic
for k=2:N_t
    dt=t(k)-t(k-1);%should be constant
    
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
    
%     if isfinite(p(k,i))
%         
%         pP=p(k,i);
%         q(k,i)=qA-1/ZA*(pP-pA)-dt*fA;%eq 5
%     elseif isfinite(q(k,i))
%         qP=q(k,i);
%         p(k,i)=pA-ZA*(qP-qA+dt*fA);
%     else
%         error('Missing boundary condition')
%     end
    R=1*ZA;
    q(k,i)=1/(1+R/ZA)*(qA+1/ZA*pA-dt*fA);
    p(k,i)=q(k,i)*R;
    
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
