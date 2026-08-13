This repository contains my currently incomplete DVPP honours thesis work.
NEM_OPEN_DVPP contains my open loop testing and IEEE9BusSystemDVPP contains the DVPP scaled up to a 9 bus system
Thesis Current contains my current thesis writing, although it is not up to date.

To change between modes in the simulation type this command into the command window, the number after DVPP_MODE can be either 1 or 2.
mw = get_param('IEEE9BusSystemDVPP','ModelWorkspace');
mw.assignin('DVPP_MODE', 2);
