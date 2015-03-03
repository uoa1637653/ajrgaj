function [Acrit,Tcrit,ie] = burgers_stability(selection,period)
% GAJ 26/02/2015
% Burgers' eq. with unit viscosity, u_t=u_xx-uu_x.
% Solve on a [0,2pi] periodic domain.
% Use finite difference approximations for u_x and u_xx.
%----------------------------------------------------------------
    % Set up domain:
    init_domain(period);
    global T Dt
    T = 10;
    Dt = which_dudt(selection);
    % Run search:
    Acrit = fzero(@search, [1 50]);
    [t, ~, ie] = integ(T, Acrit, Dt);
    if ie == 0
        Acrit = Acrit + 1e-3;
        [t, ~, ie] = integ(T, Acrit, Dt);
    end
    Tcrit = t(end);
end
%----------------------------------------------------------------
% Search function:
function v=search(A)
    global T Dt
    [~,~,ie] = integ(T,A,Dt);
    v=0.5-ie;
end
%----------------------------------------------------------------
% Temporal derivatives of u(x,t):
function u_t=dudt_adv(t,u)
    global H D2 MD
    u_t=D2*u/H^2-(u.*(MD*u))/H;
end
function u_t=dudt_cons(t,u)
    global H D2 MD
    u_t=D2*u/H^2-MD*(u.*u/2)/H;
end
function u_t=dudt_mix(t,u)
    global H D2 MD
    u_t=D2*u/H^2-1/3*(u.*(MD*u)+MD*(u.*u))/H;
end
function u_t=dudt_hol(t,u)
    global H D2 MD S
    u_t=S*(D2*u/H^2-1/3*(u.*(MD*u)+MD*(u.*u))/H);
end
function u_t=dudt_hol2(t,u)
    global H MD S SD2 II
    u_t=SD2*(II+(7*II-2*S)*SD2/60)*u/H^2-1/3*S*(u.*(MD*u)+MD*(u.*u))/H;
end
function dudt=which_dudt(sel)
    switch sel
    case 'adv'
        dudt = @dudt_adv;
    case 'cons'
        dudt = @dudt_cons;
    case 'mix'
        dudt = @dudt_mix;
    case 'hol'
        dudt = @dudt_hol;
    case 'hol2'
        dudt = @dudt_hol2;
    end
end
%----------------------------------------------------------------
% Solve PDE from t=0 to t=T with initial wave of amplitude A:
% Stop early if an extreme event occurs.
function [t,u,ie] = integ(T, A, dudt)
    opts=odeset('Events',@events);
    [t,u,~,~,ie]=ode15s(dudt,[0 T],u0(A),opts);
    if isempty(ie), ie=0; end
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