function L=linegraph(A)
% Finds linegraph for a graph G represented by its adjacency matrix A
% L=linegraph(A) returns adjacency matrix L of the linegraph G
n=size(A);
if (n(1)~=n(2))
    display('Matrix is no t square');
   
else
    n=n(1);
    for i=1:n
        A(i)=0;
    end
    
    
    tar=[];
    nei=[];
    for i=1:n
        for j=i+1:n
            if(A(i,j)==1)
                tar=[tar,i];
                nei=[nei,j];
            end
        end
    end
    m=length(tar);
    L=zeros(m);
    
        for j=m
            for k=1:m
                if(tar(j)==tar(k) || nei(j)==nei(k))
                    L(j,k)=1;
                end
            end
            L(j,j)=0;
        end
end
end
    