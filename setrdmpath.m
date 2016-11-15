% setup the path for rdm folders 

rdmhome=pwd; 

addpath(rdmhome);
addpath([rdmhome,'/Integrators']);
addpath([rdmhome,'/ModelPool']);
addpath([rdmhome,'/TestPool']);
addpath([rdmhome,'/ModifyState']);
addpath([rdmhome,'/Utilities']);
addpath([rdmhome,'/Bifurcation']);
addpath([rdmhome,'/Continuation']);
addpath([rdmhome,'/Cluster']);

clear rdmhome;
