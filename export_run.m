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

function export_run(logsout, fname)
idx   = [3 1 9 2 8];                     % delta_f, P_BESS, P_Wind, SoH, Shortfall
names = {'delta_f','P_BESS','P_Wind','SoH','Shortfall'};
t = unique([ (0:0.1:30)' ; (30:1:500)' ]);
T = table(t,'VariableNames',{'time'});
for k = 1:numel(idx)
    ts = logsout{idx(k)}.Values;
    T.(names{k}) = interp1(ts.Time, squeeze(ts.Data), t, 'linear', 'extrap');
end
writetable(T, fname);
fprintf('wrote %s (%d rows)\n', fname, height(T));
end