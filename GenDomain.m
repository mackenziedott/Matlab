%%%
% Function to generate the a N dimensionsl domain - based on nodal 
% distribution along the various dimensions. It executes in polynomial
% time. (C) Mrinal Kumar Feb 5, 2006.
%%%
function [glo_domain Nodespecs] = GenDomain(numnodes, Domlimits)
global N

%%%
% Build the node distribution vector if uniform distribution is asked for
%%%
for j = 1:N %good loop! - over the dimensions
    Nodespecs(j).ns = (linspace(Domlimits(j,1), Domlimits(j,2), numnodes(j)))';
end

%%%
% Main loop - for building the index-matrix - this matrix basically
% contains all the possible combinations of indices from 1 to numnodes
% along the respective directions. e.g. for 2 nodes each in 2D, the result
% would be [1 1; 2 1; 1 2; 2 2]
%%%
canon0 = (1:numnodes(1))';  %short for canonical_0 - first dimension's nodes: this will be loooped through the dimensions
for ct = 2:N    %good loop! - over the dimensions
    repel = canon0; %REPetition-ELement
    repsize = length(canon0(:,1));  %REPetition SIZE
    repwith = ones(repsize,1);  %REPeat WITH this structure
    for rs = 2:numnodes(ct)
        repwith = [repwith; ones(repsize,1)*rs];
    end
    canon0 = [repmat(repel,numnodes(ct),1), repwith];
end
%%%
% Build the global domain out of the index matrix
%%%
glo_domain = [];
for j = 1:N %good loop! - run through the number of dimensions!
    glo_domain = [glo_domain, Nodespecs(j).ns(canon0(:,j))];
end