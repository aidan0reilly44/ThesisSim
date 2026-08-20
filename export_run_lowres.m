function export_runlowres(logsout, fname)
names = {'delta_f','P_BESS','P_Wind','SoH','Shortfall'};
t = unique([ (0:0.1:30)' ; (30:1:500)' ]);   % fine over the transient, coarse after
T = table(t,'VariableNames',{'time'});
for k = 1:numel(names)
    ts = logsout.get(names{k}).Values;
    T.(names{k}) = interp1(ts.Time, squeeze(ts.Data), t, 'linear', 'extrap');
end
writetable(T, fname);
fprintf('wrote %s (%d rows)\n', fname, height(T));
end