

function result = get_ABIC_for_min(mu)

%% Set parameters.
parm_num = 2; % For sigma and omega.
un = 0; % Rank of B*B, i.e., number of edges inside.
nkm = 0; % Number of conditions + number of 0-singular value of A + number
% of 0-eigenvalue of B, i.e., number of cells.
sSa = [];
sSb0 = [];

%% Get ABIC.
result = ABIC(mu,parm_num,un,nkm,sSa,sSb0);

end
