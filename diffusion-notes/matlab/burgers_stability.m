function [Acrit,iem,ie,iep] = burgers_stability(selection,period)
% GAJ 26/02/2015
% Burgers' eq. with unit viscosity, u_t=u_xx-uu_x.
% Solve on a [0,2pi] periodic domain.
% Use finite difference approximations for u_x and u_xx.
%----------------------------------------------------------------
% Set up domain:
init_domain(period);
global T Dt
T=10;
Dt=which_dudt(selection);
% Run search:
Acrit=fzero(@search, [1 50]);
iem=0.5-search(Acrit-0.01);
ie=0.5-search(Acrit-0.01);
iep=0.5-search(Acrit+0.01);
end
%----------------------------------------------------------------
% Search function:
function v=search(A)
    global T Dt
    [~,~,ie] = integ(T,A,Dt);
    if isempty(ie), ie=0; end
    v=0.5-ie;
end
%----------------------------------------------------------------
% Temporal derivatives of u(x,t):
function u_t=dudt_std(t,u)
    global H D2 MD
    u_t=D2*u/H^2-(u.*(MD*u))/H;
end
function u_t=dudt_forn(t,u)
    global H D2 MD
    u_t=D2*u/H^2-1/3*(u.*(MD*u)+MD*(u.*u))/H;
end
function u_t=dudt_holi(t,u)
    global H D2 MD S
    u_t=S*(D2*u/H^2-1/3*(u.*(MD*u)+MD*(u.*u))/H);
end
function dudt=which_dudt(sel)
    switch sel
    case 'std'
        dudt = @dudt_std;
    case 'forn'
        dudt = @dudt_forn;
    case 'holi'
        dudt = @dudt_holi;
    end
end
%----------------------------------------------------------------
% Solve PDE from t=0 to t=T with initial wave of amplitude A:
% Stop early if an extreme event occurs.
function [t,u,ie] = integ(T, A, dudt)
    opts=odeset('Events',@events);
    [t,u,te,ye,ie]=ode15s(dudt,[0 T],u0(A),opts);
end
%----------------------------------------------------------------
% Initialisation of u(x,t) at t=0:
function u=u0(A)
    global L x
    u=A*sin(x(2:L+1));
end
%----------------------------------------------------------------
% Integration events:
function [val,isfin,dirn]=events(t,y)
    isfin=[1, 1];
    % Detect any non-monotonic modes.
    v=sum(abs(diff(sign(diff(y)))))/4;
    % Detect extreme wave heights.
    h=max(abs(y));
    % Event occurs if either quantity below is 0 or changes sign.
    val=[1.01-v, 1e3-h];
    dirn=[0, 0];
end