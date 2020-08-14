function [ x,t,Zc,c_bar ] = MOCinit( N_x,N_t, x_lim, c, rho, r  )
%[ x,t,Zc,c_bar ] = MOCinit( N_x,N_t, x_lim, c, rho, r  )
%   Initializes Method of Characteristics (MOC) solver. This selects grid
%   points for varying wave speed.
% N_x: number of spatial grid points
% N_t: number of temporal grid points
% x_lim: scalar length or two element vector with start and end position
% c: speed of sound. May be scalar or function c(x)
% rho: fluid density. May be scalar or function rho(x)
% r: inner radius. May be scalar or function r(x)
%
% Returns:
% x: grid positions
% t: time
% Zc: characteristic impedance. f(x)
% c_bar: average sonic speed


if ~isa(c, 'function_handle')
    c=@(x) c*ones(size(x));
end

if ~isa(r, 'function_handle')
    r=@(x) r*ones(size(x));
end

if ~isa(rho, 'function_handle')
    rho=@(x) rho*ones(size(x));
end

if numel(x_lim)==1
    x_lim=[0 x_lim];
end



%calculate wave progression
N_int=N_x*10;%number of points to integrate over
x_int=linspace(x_lim(1),x_lim(2),N_int);

t_wave=cumtrapz(x_int,1./c(x_int));%time to reach x_int points

c_bar=(x_lim(2)-x_lim(1))/t_wave(end);%average wve speed

t_grid=linspace(0,t_wave(end),N_x);%evenly spaced times
dt=t_grid(2)-t_grid(1);%time step
t=(0:(N_t-1))*dt;%time

x=interp1(t_wave,x_int,t_grid);%x grid for MOC, with grid points equally spaced in wave time

if nargout>2
    
    
    Zc=@(x) rho(x).*c(x)./(pi*r(x).^2);
    
end

end

