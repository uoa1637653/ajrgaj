function burgers_stab_plot2(Ls)
%----------------------------------------------------------------
% GAJ 17/03/2015
%----------------------------------------------------------------
    figure
    hold on
    plot_it('adv', Ls,[1 0 0]*0.7); 
    %plot_it('cons', Ls,[1 1 0]*0.7);
    %plot_it('mix', Ls,[0 1 0]*0.7);
    %plot_it('hol', Ls,[0 1 1]*0.7);
    xlabel('#intervals, L')
    ylabel('critical amplitude, A')
    %legend('advective','conservative','mixture','holistic')
    legend('advective')
    hold off
end

function plot_it(selection,Ls,rgb)
    N = length(Ls);
    As = zeros(1, N);
    for i = 1:N
        [Acrit] = burgers_eigs(selection, Ls(i));
        As(i) = Acrit;
    end
    plot(Ls, As, 'o-', 'Color', rgb);
end
