function [newface, newfaceIdx] = clip_mesh(face)
% remove dangling triangles and ears from a mesh
%
% [newface, newfaceIdx] = clip_mesh(face)
%
% face:         triangulation
% newface:      clipped triangulation 
% newfaceIdx:   indices of new faces in old triangulation
%
% partially based on compute_boundary by Gabriel Peyre

if size(face,1)<size(face,2)
    face=face';
end

nvert=max(max(face));
nface=size(face,1);

% each position in the matrix A is a directed edge connecting two vertices
% the value is increased if it is part of a face
% if it part of only one face, it is on the boundary
A=sparse(nvert,nvert);
for i=1:nface
    f=face(i,:);
    A(f(1),f(2))=A(f(1),f(2))+1;
    A(f(1),f(3))=A(f(1),f(3))+1;
    A(f(3),f(2))=A(f(3),f(2))+1;
end
A=A+A';

% find ears (bdry triangles sharing one edge with mesh) and triangles
% sharing only a vertex with the mesh
% the former have S == 4, the latter S==3, the good ones S == 6
newfaceIdx = false([nface 1]);
for i=1:nface
    
    f=face(i,:);
    S = A(f(1),f(2)) + A(f(1),f(3)) + A(f(3),f(2));
    newfaceIdx(i) = S>=5;
end

% make new face list
newface = face(newfaceIdx,:);


    
