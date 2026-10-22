function [root,base]=setup
root=fileparts(fileparts(mfilename('fullpath')));
base=fullfile(fileparts(root),'6.位置反馈方案');
assert(exist(fullfile(base,'+posfb','config.m'),'file')==2);
addpath(root,base,fullfile(base,'lib'));maxNumCompThreads(1);
end
