function burgers_stab_plot4(M)
%----------------------------------------------------------------
% GAJ 09/02/2016
%----------------------------------------------------------------
    figure
    hold on
    handles = [];
    labels = {};
    [handles, labels] = plot_it(handles, labels, 'advective', 3:M, [1 0.1 0.1], ['s','x']); 
    [handles, labels] = plot_it(handles, labels, 'conservative', 3:M, [0.1 1 0.1], ['d','+']);
    [handles, labels] = plot_it(handles, labels, 'mixture', 3:M, [0.1 0.1 1], ['o','*']);
    [handles, labels] = plot_it(handles, labels, 'holistic', 3:M, [1 1 0.1]*0.7, ['v','.']);
    xlabel('number of intervals, N')
    ylabel('critical amplitude, A')
    legend(handles, labels, 'Location', 'East');
    hold off
end

function [handles, labels] = plot_it(handles, labels, selection, Ls, rgb, syms)
    [L1,A1,L2,A2] = compute(selection,Ls);
    if ~isempty(L1)
        handles(1+length(handles)) = plot(L1, A1, syms(1), 'Color', rgb);
        labels{1+length(labels)} = strcat(selection,', non-monotonic');
    end
    if ~isempty(L2)
        handles(1+length(handles)) = plot(L2, A2, syms(2), 'Color', rgb);
        labels{1+length(labels)} = strcat(selection,', instability');
    end
end

function [L1,A1,L2,A2] = compute(selection,Ls)
    N = length(Ls);
    L1 = [];
    A1 = [];
    L2 = [];
    A2 = [];
    for i = 1:N
        try
            [Ap, ~, iep] = burgers_stability(selection, Ls(i));
            if iep > 0
                if iep == 1 % non-monotonic
                    L1 = [L1, Ls(i)];
                    A1 = [A1, Ap];
                else        % instability
                    L2 = [L2, Ls(i)];
                    A2 = [A2, Ap];
                end
            end
        catch
        end
        try
            [Am, ~, iem] = burgers_stability_neg(selection, Ls(i));
            if iem > 0
                if iem == 1 % non-monotonic
                    L1 = [L1, Ls(i)];
                    A1 = [A1, Am];
                else        % instability
                    L2 = [L2, Ls(i)];
                    A2 = [A2, Am];
                end
            end
        catch
        end
    end
end
