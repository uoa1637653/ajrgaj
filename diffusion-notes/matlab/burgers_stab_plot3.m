function burgers_stab_plot3(M)
%----------------------------------------------------------------
% GAJ 01/10/2015
%----------------------------------------------------------------
    figure
    hold on
    p1 = plot_it('adv', 3:M, [1 0.1 0.1]*0.7); 
    p2 = plot_it('cons', 3:M, [1 1 0.1]*0.7);
    p3 = plot_it('mix', 3:M, [0.1 1 0.1]*0.7);
    p4 = plot_it('hol', 3:M, [0.1 1 1]*0.7);
    xlabel('number of intervals, N')
    ylabel('critical amplitude, A')
    legend([p1,p2,p3,p4],'advective','conservative','mixture','holistic')
    hold off
end

function pid = plot_it(selection,Ls,rgb)
    N = length(Ls);
    Ls2 = [];
    As = [];
    for i = 1:N
        try
            [Ap, ~, iep] = burgers_stability(selection, Ls(i));
            if iep > 0
            	if iep == 1; sym = '+-'; else sym = 'x-'; end
            	plot(Ls(i), Ap, sym, 'Color', rgb);
                if isempty(As)
                    As = Ap;
                    Ls2 = Ls(i);
                end
            end
        catch
        end
        try
            [Am, ~, iem] = burgers_stability_neg(selection, Ls(i));
            if iem > 0
            	if iem == 1; sym = '+-'; else sym = 'x-'; end
            	plot(Ls(i), Am, sym, 'Color', rgb);
                if isempty(As)
                    As = Am;
                    Ls2 = Ls(i);
                end
            end
        catch
        end
    end
    pid = plot(Ls2, As, '-', 'Color', rgb);
end
