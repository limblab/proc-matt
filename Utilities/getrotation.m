     function [R,N]=getrotation(U,v)
     %U - a matrix of size nx(n-2) whose columns are an orthogonal basis 
     %    for a subspace of dimension n-2
     %
     %v - a given vector
     %
     %R - final rotation matrix
     %N - axis of rotation
        v=v(:)/norm(v);
        N=null([U,v].'); %axis of rotation
        q=U*(U.'*v);
          q=q/norm(q);
        A=null([N,v].');
        B=null([N,q].');
        AA=[N,v,A];
        BB=[N,q,B].';
         AA(:,3)=AA(:,3)*sign(det(BB*AA));
        R=AA*BB; %final rotation