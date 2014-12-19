%%%%%%%%%%
write "Generating matlab function file holpde.m";
linelength 60$
jj:=(stenwidth-1)/2$
load_package scope;
out "holred";
off nat;
write "function uudot=holpde(t,uu)$
global h"$
if jj=1 then write "% periodic BCs
uu=[uu(end);uu;uu(1)]"
else write "% periodic BCs
uu=[uu(end-",jj-1,":end);uu;uu(1:",jj,")]"$
write "m=length(uu)$
j=(",jj+1,":m-",jj,")'"$
optimize uudot:=:gj iname o$
on nat;
shut "holred";
showtime;
write "finished generating matlab file";

%echo convert file to matlab with
%sed -f red2mat holred > ../Sites/RWX/holpde.m
