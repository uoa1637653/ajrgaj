function burgers_stab_plot2(M)
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
		Acrit = Inf;
		ie = -1;
        try
            [Ap, ~, iep] = burgers_stability(selection, Ls(i));
            if iep > 0
                Acrit = Ap;
                ie = iep;
            end
        catch
        end
        try
            [Am, ~, iem] = burgers_stability_neg(selection, Ls(i));
            Am = abs(Am);
            if iem > 0 && Am < Acrit
            	Acrit = Am;
            	ie = iem;
           	end
        catch
        end
        if ie > 0
	        As = [As; Acrit];
	        Ls2 = [Ls2; Ls(i)];
        	if ie == 1
            	sym = '+-';
        	elseif ie == 2
            	sym = 'x-';
        	else
            	sym = 'o-';
        	end
        	plot(Ls(i), Acrit, sym, 'Color', rgb);
        end
    end
    pid = plot(Ls2, As, '-', 'Color', rgb);
end
