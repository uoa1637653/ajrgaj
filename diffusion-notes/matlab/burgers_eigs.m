function [Acrit] = burgers_eigs(selection,period)
% GAJ 17/03/2015
% Burgers' eq. with unit viscosity, u_t=u_xx-uu_x.
% Solve on a [0,2pi] periodic domain.
% Use finite difference approximations for u_x and u_xx.
%----------------------------------------------------------------
    % Set up domain:
    init_domain(period);
    global Dt
    Dt = which_dudt(selection);
    % Run search:
    Acrit = fzero(@search, [1 50]);
end
%----------------------------------------------------------------
% Search function:
function v=search(A)
    global Dt
    [~,lam]=calc_eigs(u0(A),Dt);
    if sum(real(lam) > 1e-5)
        v = 1;
    else
        v = -1;
    end
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
    global H D2 MD Sinv
    u_t=Sinv\(D2*u/H^2-1/3*(u.*(MD*u)+MD*(u.*u))/H);
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
% Initialisation of u(x,t) at t=0:
function u=u0(A)
    global L x
    u=A*sin(x(2:L+1));
end