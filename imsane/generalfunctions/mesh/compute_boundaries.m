function boundaries=compute_boundaries(faces)
% compute boundaries of a mesh
%
% boundaries=compute_boundary(faces);
%
% faces:        triangulation
% boundaries:   cell array containing ordered list of vertices for each
%               boundary
% 
% based on function by Gabriel Peyre, modified to deal with multiply
% connected meshes

nvert=max(max(faces));
nface=size(faces,1);

% each position in the matrix A is a directed edge connecting two vertices
% the value is increased if it is part of a face
% if it part of only one face, it is on the boundary
A=sparse(nvert,nvert);
for i=1:nface
    f=faces(i,:);
    A(f(1),f(2))=A(f(1),f(2))+1;
    A(f(1),f(3))=A(f(1),f(3))+1;
    A(f(3),f(2))=A(f(3),f(2))+1;
end
A=A+A';

% idenitify vertices on the boundary
bdryVertIdx = false([nvert 1]);
for i=1:nface
    
    f=faces(i,:);

    if A(f(3),f(2)) == 1
        bdryVertIdx(f(2)) = true;
        bdryVertIdx(f(3)) = true;
    end
    if A(f(1),f(2)) == 1
        bdryVertIdx(f(1)) = true;
        bdryVertIdx(f(2)) = true;
    end
    if A(f(1),f(3)) == 1
        bdryVertIdx(f(1)) = true;
        bdryVertIdx(f(3)) = true;
    end
end

bdryVertIdx = find(bdryVertIdx);
nBdries = 0;
boundaries = {};
while ~isempty(bdryVertIdx)
    
    nBdries = nBdries + 1;
    
    % starting point
    i=bdryVertIdx(1);
    u=find(A(i,:)==1);
    boundary=[i u(1)];
    
    % the rest
    s=boundary(2);
    i=2;
    while(i<=nvert)
        u=find(A(s,:)==1);
        % if a boundary vertex is connected to more than 2 other boundary
        % vertices, something is wrong (i.e. a triangle sharing only one vertex
        % with the rest of the mesh)
        if length(u)~=2
            warning('problem in boundary');
        end
        if u(1)==boundary(i-1)
            s=u(2);
        else
            s=u(1);
        end
        if s~=boundary(1)
            boundary=[boundary s];
        else
            break;
        end
        i=i+1;
    end
    
    bdryVertIdx = setdiff(bdryVertIdx, boundary);
    boundaries{nBdries} = boundary;
end