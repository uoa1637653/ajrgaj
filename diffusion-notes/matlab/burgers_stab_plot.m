function burgers_stab_plot(Ls)
%----------------------------------------------------------------
% GAJ 02/03/2015
%----------------------------------------------------------------
    figure
    hold on
    plot_it('adv', Ls,[1 0 0]*0.7); 
    plot_it('cons', Ls,[1 1 0]*0.7);
    plot_it('mix', Ls,[0 1 0]*0.7);
    plot_it('hol', Ls,[0 1 1]*0.7);
%    plot_it('hol2', Ls,[0 0 1]*0.7);
    xlabel('#intervals, L')
    ylabel('critical amplitude, A')
%    legend('advective','conservative','mixture','holistic','holistic2')
    legend('advective','conservative','mixture','holistic')
    plot_it('adv', Ls+1,[1 0 0]*0.7); 
    plot_it('cons', Ls+1,[1 1 0]*0.7);
    plot_it('mix', Ls+1,[0 1 0]*0.7);
    plot_it('hol', Ls+1,[0 1 1]*0.7);
%    plot_it('hol2', Ls+1,[0 0 1]*0.7);
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
    end
    if sum(ies == 1) == N
        sym = '+-';
    elseif sum(ies == 2) == N
        sym = 'x-';
    else
        sym = 'o-';
    end
    plot(Ls, As, sym,'Color',rgb);
end
