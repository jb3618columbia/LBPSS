function [ marginalsLBP, marginalsJT ] = get_LBP_marginals( a , b )
% Computes loopy belief propagation marginal approximations

% Likelihoods are of the form \prod_{i=1}^d \exp{a_i s_i} \prod_{i < j} \exp{b_{i,j} s_i s_j}
% where s_i \in \{+1,-1\}
% require b to be symmetric, i.e. b_{i,j}=b_{j,i}, if not will make b=b+b'
% note that the edges b_{i,j} are NOT double counted since the product is
% over {i < j}
% Returns marginal vector with i^{th} coordinate = Pr(s_i=1)

nStates = 2;

% Make dimensions for vector a correct
[a_rows, a_col] = size(a);
if a_col > a_rows
    a = a';
end

% Make b symmetric if not already
if ~isequal(b,b') % If not symmetric
    disp('b matrix not symmetric, making it symmetric')
    b = b+b'; % Make symmetric
end

% Create structures
nodePot = [exp(a),exp(-a)]; % Node potentials
adj = b~=0; % Adjacency matrix
edgeStruct = UGM_makeEdgeStruct(adj,nStates); % Generate edgeStruct
maxState = max(edgeStruct.nStates); % Number of edges
edgePot = zeros(maxState,maxState,edgeStruct.nEdges); % Initialize edge potentials
for e = 1:edgeStruct.nEdges
    indices = edgeStruct.edgeEnds(e,:);
   edgePot(:,:,e) = [exp(b(indices(1),indices(2))) exp(-b(indices(1),indices(2))) ; exp(-b(indices(1),indices(2))) exp(b(indices(1),indices(2)))]; 
end

% Find the marginals
% Using LBP
[nodeBelLBP,~,~] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);
marginalsLBP = nodeBelLBP(:,1);

% Using Junction Tree (exact)
[nodeBelJT,~,~] = UGM_Infer_Junction(nodePot,edgePot,edgeStruct); 
marginalsJT = nodeBelJT(:,1);

end

