function [header,myline1,myline2] = Header
% Set stuff
myline1 = [repmat('=',1,32),'\n'];
myline2 = [repmat('-',1,32),'\n'];
header  = [' iter  |   Obj   |   alpha  | time (s) |\n'];