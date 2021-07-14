function subm = submeshFromV(mesh, vidx)
    % get a submesh by specifying vertex indices
    %
    % subm = submeshFromV(mesh, vidx)
    %
    % mesh : struct with fields v, f, vn, (optionally u)
    % subm : same
    % vidx : logical index

    if size(vidx,2) > 1
        vidx = vidx';
    end
    
    assert(islogical(vidx), 'vidx should be logical');
    
    % determine the corresponding faces in the submesh
    fidx = ismember(mesh.f, find(vidx));
    submFidx = all(fidx, 2);

    % now create a proper submesh
    newv = mesh.v(vidx,:);
    newvn = mesh.vn(vidx,:);
    
    if isfield(mesh, 'u')
        newu = {};
        for i = 1:length(mesh.u)
            newu{i} = mesh.u{i}(vidx,:);
        end
    else
        newu = [];
    end

    old2new = zeros(size(vidx));
    old2new(vidx) = 1:sum(vidx);

    newf = mesh.f(submFidx,:);
    for j = 1:3
        newf(:,j) = old2new(newf(:,j));
    end

    % submesh boundaries
    b = compute_boundaries(newf);

    subm = struct('v', newv, 'f', newf, 'vn', newvn, 'u', {newu}, 'b', {b});
end
            