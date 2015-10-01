function burgers_stab_plot2(M)
%----------------------------------------------------------------
% GAJ 01/10/2015
%----------------------------------------------------------------
    figure
    hold on
    p1 = plot_it('adv', [3,5:M], [1 0.1 0.1]*0.7); 
    p2 = plot_it('cons', [3,5:M], [1 1 0.1]*0.7);
    p3 = plot_it('mix', 5:M, [0.1 1 0.1]*0.7);
    p4 = plot_it('hol', 5:M, [0.1 1 1]*0.7);
    xlabel('number of intervals, N')
    ylabel('critical amplitude, A')
    legend([p1,p2,p3,p4],'advective','conservative','mixture','holistic')
    hold off
end

function pid = plot_it(selection,Ls,rgb)
    N = length(Ls);
    As = zeros(1, N);
    for i = 1:N
        try
            [Ap, ~, iep] = burgers_stability(selection, Ls(i));
        catch
            Ap = +Inf;
            iep = 0;
        end
        try
            [Am, ~, iem] = burgers_stability_neg(selection, Ls(i));
        catch
            Am = -Inf;
            iem = 0;
        end
        Acrit = min(abs(Am), Ap);
        As(i) = Acrit;
        if Acrit == Ap, ie = iep; else ie = iem; end
        if ie == 1
            sym = '+-';
        elseif ie == 2
            sym = 'x-';
        else
            sym = 'o-';
        end
        plot(Ls(i), Acrit, sym, 'Color', rgb);
    end
    pid = plot(Ls, As, '-', 'Color', rgb);
end
