function result=run_upa_ue(varargin)
% Run the UPA UE-side co-phase feedback branch using the shared v2 protocol.
result=run_upa_cpu(varargin{:},'PanelBranch','UE');
end
