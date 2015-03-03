function burgers_stab_plot(Ls)
%----------------------------------------------------------------
% GAJ 02/03/2015
%----------------------------------------------------------------
    figure
    hold on
    plot_it('adv', Ls);
    plot_it('cons', Ls);
    plot_it('mix', Ls);
    plot_it('hol', Ls);
    plot_it('hol2', Ls);
    xlabel('#intervals, L')
    ylabel('critical amplitude, A')
    legend('advective','conservative','mixture','holistic','holistic2')
end

function plot_it(selection,Ls)
    N = length(Ls);
    As = zeros(1, N);
    Ts = zeros(1, N);
    ies = zeros(1, N);
    for i = 1:N
        [Acrit, Tcrit, ie] = burgers_stability(selection, Ls(i));
        As(i) = Acrit;
        Ts(i) = Tcrit;
        ies(i) = ie;
    end
    if sum(ies == 1) == N
        sym = '+-';
    elseif sum(ies == 2) == N
        sym = 'x-';
    else
        sym = 'o-';
    end
    plot(Ls, As, sym);
end
