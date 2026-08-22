% function export_run(logsout, fname)
% %EXPORT_RUN  Write logged signals to a CSV on a common time grid.
% names = {'delta_f','P_BESS','P_Wind','SoH','Shortfall'};
% t = (0:0.05:500)';                    % 0.05 s grid -> ~10k rows
% T = table(t,'VariableNames',{'time'});
% for k = 1:numel(names)
%     ts = logsout.get(names{k}).Values;
%     T.(names{k}) = interp1(ts.Time, squeeze(ts.Data), t, 'linear', 'extrap');
% end
% writetable(T, fname);
% fprintf('wrote %s (%d rows)\n', fname, height(T));
% end

% function export_run(logsout, fname)
% idx   = [3 1 9 2 8];                     % delta_f, P_BESS, P_Wind, SoH, Shortfall
% names = {'delta_f','P_BESS','P_Wind','SoH','Shortfall'};
% t = unique([ (0:0.1:30)' ; (30:1:500)' ]);
% T = table(t,'VariableNames',{'time'});
% for k = 1:numel(idx)
%     ts = logsout{idx(k)}.Values;
%     T.(names{k}) = interp1(ts.Time, squeeze(ts.Data), t, 'linear', 'extrap');
% end
% writetable(T, fname);
% fprintf('wrote %s (%d rows)\n', fname, height(T));
% end

function export_run(logsout, fname)
%EXPORT_RUN  Write logged DVPP signals to CSV on a common time grid.

names = {'delta_f','delta_f_COI','P_BESS','P_Wind','SoH','Shortfall'};

% fine over the transient, coarse after
t = unique([ (0:0.1:40)' ; (40:1:500)' ]);

T = table(t, 'VariableNames', {'time'});
missing = {};

for k = 1:numel(names)
    try
        ts = logsout.get(names{k}).Values;
        T.(names{k}) = interp1(ts.Time, squeeze(ts.Data), t, 'linear', 'extrap');
    catch
        missing{end+1} = names{k};  %#ok<AGROW>
    end
end

if ~isempty(missing)
    warning('export_run:missing', 'not logged: %s', strjoin(missing, ', '));
end

writetable(T, fname);
fprintf('wrote %s  (%d rows, %d signals)\n', fname, height(T), width(T)-1);
end