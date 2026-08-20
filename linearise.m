mdl = 'IEEE9BusSystemDVPP';
old = get_param(mdl,'Solver');
set_param(mdl,'Solver','ode23t');

io = getlinio(mdl);
io(1).Type = 'looptransfer';      % was 'openoutput'
L = linearize(mdl, io, 5);

set_param(mdl,'Solver',old);
size(L)
[Gm, Pm, wcg, wcp] = margin(L)
figure; bode(L); grid on