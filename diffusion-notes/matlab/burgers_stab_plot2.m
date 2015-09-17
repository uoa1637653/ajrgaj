function burgers_stab_plot2(M)
%----------------------------------------------------------------
% GAJ 15/09/2015
%----------------------------------------------------------------
    figure
    hold on
    plot_it('adv', [3,5:M], [1 0 0]*0.7); 
    plot_it('cons', 6:M, [1 1 0]*0.7);
    plot_it('mix', 5:M, [0 1 0]*0.7);
    plot_it('hol', 5:M, [0 1 1]*0.7);
    xlabel('#intervals, L')
    ylabel('critical amplitude, A')
    legend('advective','conservative','mixture','holistic')
    hold off
end

function plot_it(selection,Ls,rgb)
    N = length(Ls);
    As = zeros(1, N);
    Ts = zeros(1, N);
    ies = zeros(1, N);
    for i = 1:N
        [Acrit, Tcrit, ie] = burgers_stability(selection, Ls(i));
        As(i) = Acrit;
        Ts(i) = Tcrit;
        ies(i) = ie;
    if ie == 1
        sym = '+-';
    elseif ie == 2
        sym = 'x-';
    else
        sym = 'o-';
    end
        plot(Ls(i), Acrit, sym, 'Color', rgb);
    end
    plot(Ls, As, '-', 'Color', rgb);
end
