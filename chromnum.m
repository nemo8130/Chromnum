function [crn,colour,order] = chromnum(adj,flag,num)
% Chromnum computes the chromatic number of an "UNDIRECTED" graph
% Chromatic Number of a graph is the minimum number of colors by whcih
% All nodes of the graph could be exhaustively mapped so that no two
% adjacent nodes share the same color.
% Chromnum implements the "trailpathSA" algorithm
% The function takes the 1-0 adjacency matrix stored in a text file as its sole input argument
% And returns the chromatic number of the corresponding graph
%
% Usage: [crn,colour,order] = chromnum('adj.inp')
% Usage: [crn,colour,order] = chromnum('adj.inp',flag) 
% for graphix display flag=1, flag=0 will not display graphix
% Usage: [crn,colour,order] = chromnum('adj.inp',flag,num) num is the minimun no of iteration
%
% It also generates one color map for the graph and tabulates the same
% Please note that this color map could be degenerate but the chromatic
% number should be identical
% Furthermore it also returns the order at which the nodes have been
% colored. Again this may have degenerate solutions.
%
% Example adjacency matrix file :
% (Each row should contain the adjacencies of one node in the graph)
% (The matrix entries should be seperated by a single white-space)
%
% 0 1 1 0 1
% 1 0 0 1 1
% 1 0 0 0 1
% 0 1 0 0 1
% 1 1 1 1 0
%
% OUTPUT:
%--------
% Number of Nodes in the given undirected graph : 5
% Number of edges : 7 
% Chromatic Number : 3 

% Proposed Colormap:
% ------------------
% Node-1 : color-3
% Node-2 : color-2
% Node-3 : color-2
% Node-4 : color-3
% Node-5 : color-1
%
%
%
% Authors: Sankar Basu 1* and Abhirup Bandyopadhyay 2
% 
% 1. Department of Biochemistry, University of Calcutta
% Circular Rd, Ballygunge, Kolkata, 
% West Bengal 700019
%
% Current affiliation: Department of Chemistry
% University of Delhi (North Campus)
% Delhi 110007 
%
% 2. Department of Mathematics, National Institute of Technology, Durgapur
% Mahatma Gandhi Avenue, Durgapur 713209, West Bengal, India
%
% Email: Sankar Basu (nemo8130@gmail.com, scbasu@chemistry.du.ac.in)
% Abhirup Bandyopadhyay (ab.13math1110@phd.nitdgp.ac.in, abhirupnit@gmail.com)
% 
% (C) All copyrights reserved to the authors; any violations may be treated as illegal  
%
%

tic
if(nargin>3)
    sprintf('Enter the adjacency matrix as ONLY ONE iput file- .txt etc......');
else
    A0=load(adj,'s');
    no=size(A0);
    
    A=A0;
    nn=size(A);
    n=nn(1);

    % ZERO PADDING OF THE DIAGONAL ELEMENTS
    
    for i = 1:n
        A(i,i) = 0;        
    end
    
    A0=A;
    
    % CHECK FOR SYMMETRY (i.e., WHETHER THE GRAPH IS UNDIRECTED)
    
    csym = 0;
    matent = n*n;    
    nelmND = (n*n)-n;
    nones = 0;
    
    for i = 1:n
        for j = 1:n
            if (A0(i,j) == A0(j,i))
                csym = csym + 1;
            end
            if (A0(i,j) == 1)
                nones = nones + 1;
            end
        end
    end
    
    if (csym == matent)
    fprintf('The Adjacency Matrix is symmetric:\nThe graph is undirected;\nProgram will proceed\n\n');
    else
        fprintf('The Adjacency Matrix is not symmetric:\nThe graph is NOT undirected !!!\nProgram will Cease.\n\n');
        return;
    end
    
    % QUICK CHECK AND RETURN FOR COMPLETE GRAPHS
    
    if (nones == nelmND)
        fprintf ('%s %d %s \n%s %d\n','Its a complete graph of ',n,' nodes','The chromatic number therefore will trivially be ',n);
        crn = n;
        colour = [1:n]';
        order = [1:n];
        return;
    end
    
    if(length(no)~=2)
        sprintf('Enter a TWO DIMENSONAL adjacency matrix as an iput file- .txt etc......');
    else
        if(no(1)~=no(2))
            sprintf('Enter a 2-dimensonal SQUARE adjacency matrix as an iput file- .txt etc......');
        else
            n=no(1);
            size(A0);
            
            if(nargin<3)
                num=1;
            end
            %==============================================================            
            %==============================================================
            [crn,colour,order]=chromatic_no_regcor(A0,num);
            %==============================================================
            %==============================================================
                        
            targ = [];
            nei = [];
            for i = 1:n
                for j = (i+1):n
                    if (A0(i,j)==1)
                        targ = [targ;i];
                        nei = [nei;j];
                    end
                end
            end
            
            GNCF = zeros(n,3);
            rtrip = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0];
            nc=size(rtrip);
            nc=nc(1);
            for i=1:n
                if(crn<i^3)
                    no=i;
                    break;
                end
            end
            diff=(1/no)-(1/no)^10;
            
            for i = 1:n
                if(colour(i)<=nc)
                    GNCF(i,:) = rtrip(colour(i),:);
                    
                else
                    
                    test=zeros(1,n);
                    tst=0;
                    while(tst<diff)
                        rr=rand(1,3);
                        for j=1:n
                            test(j)=sqrt(((GNCF(j,1)-rr(1))^2)+((GNCF(j,2)-rr(2))^2)+((GNCF(j,3)-rr(3))^2));
                        end
                        tst=min(test);
                    end
                    GNCF(i,:)=rr;
                    
                end
                
            end
            GEC = 'black';
            
            y=zeros(1,n);
            le=sum(sum(A0))/2;
            GEW = zeros(le,1)';
            length(GEW);
            ld = le/(n*(n-1)/2);

            fprintf('\nNumber of Nodes in the given undirected graph : %d',n);
            fprintf('\nNumber of edges : %d ',le);
            fprintf('\nLink Density: %f ',ld);
            fprintf('\nChromatic Number : %d \n',crn);
            
            fprintf ('\nProposed Colormap:');
            fprintf ('\n------------------\n');
            
            for i = 1:n            
            fprintf('Node-%d : color-%d\n',i,colour(i));
            end
            
            fprintf ('\n');
        end
    end
                        

       
end

if (nargin>1 && flag==1)
                figure
            hold on
            
            for i=1:crn
                b=num2str(i);
                c='colour-';
                name=[c,b];
                for j=1:n
                    if(colour(j)==i)
                        y(i)=plot(0,0,'color',GNCF(j,:),'DisplayName',name);
                        break
                    end
                end
            end
            
            G = graph(targ,nei,GEW,'OmitSelfLoops');
            
            title('Chromnum: Colour hierarchy is as follows as in the legend')
            plot(G,'Layout','circle','NodeColor',GNCF,'EdgeColor',GEC,'MarkerSize',25,'LineWidth',2.0);%,'EdgeLabel',GEW)%    Layout: circle
            legend('show');
            title('Chromnum: Colour hierarchy is as follows as in the legend')
            set (gca,'FontSize',20)
            axis equal
 
            figure
            hold on
            
            for i=1:crn
                b=num2str(i);
                c='colour-';
                name=[c,b];
                for j=1:n
                    if(colour(j)==i)
                        y(i)=plot(0,0,'color',GNCF(j,:),'DisplayName',name);
                        break
                    end
                end
            end
                                   
            plot(G,'Layout','force','NodeColor',GNCF,'EdgeColor',GEC,'MarkerSize',25,'LineWidth',2.0);%,'EdgeLabel',GEW)%    Layout: circle
            legend('show');
            title('Chromnum: Colour hierarchy is as follows as in the legend')
            set (gca,'FontSize',20)
            axis equal

end
toc
end


%======================================================================================================

function [chnum,colour,sequence] = chromatic_no_regcor(A0,num)
%tstart=tic;
dim=size(A0);
dim=dim(2);

indicator=0;

%==============================================================================
% SETTING UP PARAMETERS FOR NUMBER OF ITERATIONS:  nmax >= ni MUST BE SATISFIED
%==============================================================================

if(num==1)
    ni=1; 
    nmax=1; 
else
    nmax=round(sqrt(10*num))+1;
    if(nmax==0)
        nmax=1;
    end
    ni=round(nmax/10)+1;
    if(ni==0)
        ni=1;
    end
end

fprintf('Total number of iteration is %d',ni*nmax);


ntest=25;
testi=10;
counter=1;

while(indicator==0)
    nr=ni;
    
    if(nr>nmax)
        str=num2str(counter);
        %sprintf('\nError Message: step %s : error-00 Convergence not attained.\nConvergence could not be attend through the program.\nFor random simulation upto %d times state %d',str,nmax,counter)
        indicator=1;
    else
        str=num2str(counter);
        number=zeros(1,ntest);
        for loop=1:ntest
            val=dim+1;
%	fprintf('%s %d\n','val:',val)
            for index=1:nr
                n=dim;
                A=A0;
                for i = 1:n
                    A(i,i) = 0;
                end
                
                cn=n;
                ba=zeros(n,cn);
                ca=zeros(n,1);
                bfa=zeros(n,1);
                order=zeros(1,n);
                od=zeros(1,n);
                cntr=0;
                while(min(ca)==0)
                    m=max(sum(A));
                    cnt=0;
                    
                    while (m>0)
                        %%                  FINDING THE INDEX NODE ID WHICH SHOULD BE COLOURED IN EACH STEP
                        dseq=sum(A);
                        m=max(dseq);
                        if (cnt==0)
                            [m,id]=max(dseq);
%================================================ FROM CPOL ==========================
                            tpm3=zeros(1,n);
                            inx3=[];
                        for i=1:n
                            if (dseq(i)==m)
                                tpm3(i)=dseq(i);
                                inx3=[inx3,i];
                            end
                        end
                        nnx=length(inx3);
                        idd=randi([1,nnx],1);
                        inx3(idd);
                        id=inx3(idd);
                        inx3=[];
%====================================================================================
                        else
                            [ttp,ii]=max(bfa);
                            
                            tpm1=zeros(1,n);
                            inx1=zeros(1,n);
                            
                            tpm1(1)=dseq(ii);
                            inx1(1)=ii;
                            for i=1:n
                                if (bfa(i)==ttp)
                                    tpm1(i)=dseq(i);
                                    inx1(i)=i;
                                end
                            end
                            [mn,ij]=max(tpm1);
                            if(mn==0)
                                break
                            else
                                ii=inx1(ij);
                                
                                inx=zeros(1,n);
                                cnt=0;
                                for i=1:n
                                    if (tpm1(i)==mn)
                                        inx(i)=inx1(i);
                                        cnt=cnt+1;
                                    end
                                end
                                
                                if(cnt>1)
                                    r=randi(cnt,1);
                                    j=0;
                                    for i=1:n
                                        if(inx(i)>0)
                                            j=j+1;
                                        end
                                        if(j==r)
                                            break
                                        end
                                    end
                                    id=inx(i);
                                else
                                    id=ii;
                                end
                            end
                        end
                        
                        %%                  FINDING THE SET OF NEIGHBOURS OF INDEX NODE ID
                        
                        idx=0;
                        for i=1:n
                            if (A(id,i)==1)
                                idx=idx+1;
                            end
                        end
                        nbl=zeros(1,idx);
                        
                        nn=idx;
                        idx=0;
                        for i=1:n
                            if (A(id,i)==1)
                                idx=idx+1;
                                nbl(idx)=i;
                            end
                        end
                        
                        
                        while (isempty(nbl)==0)
                            %%                         COLOURINR THE INDEX NODE ID
                            
                            for i=1:cn
                                if (ba(id,i)==0)
                                    ca(id)=i;
                                    break
                                end
                            end
                            cntr=cntr+1;
                            order(cntr)=id;
                            A(id,:)=zeros(1,n);
                            A(:,id)=zeros(n,1);
                            for i=1:nn
                                ba(nbl(i),ca(id))=1;
                            end
                            for i=1:nn
                                if (dseq(nbl(i))==1)
                                    for k=1:cn
                                        if (ba(nbl(i),k)==0)
                                            ca(nbl(i))=k;
                                            break
                                        end
                                    end
                                    cntr=cntr+1;
                                    order(cntr)=nbl(i);
                                    A(nbl(i),:)=zeros(1,n);
                                    A(:,nbl(i))=zeros(n,1);
                                end
                            end
                            
                            %%                   BFA NO. OF THE COLOURS THAT ARE ACCUPIED BY THE NEIGHBOURS
                            for i=1:n
                                bfa(i)=sum(ba(i,:));
                            end
                            
                            %%                        FINDING THE NEW INDEX NODE ID WHICH SHOULD BE COLOURED IN EACH STEP
                            tmpd=zeros(1,nn);
                            tmpb=zeros(1,nn);
                            for i=1:nn
                                tmpd(i)=dseq(nbl(i));
                                tmpb(i)=bfa(nbl(i));
                            end
                            [mm,ii]=max(tmpb);
                            
                            tpm2=zeros(1,n);
                            tpm2(1)=tmpd(ii);
                            inx2=zeros(1,n);
                            inx2(1)=ii;
                            for i=1:nn
                                if (tmpb(i)==mm && i~=ii)
                                    tpm2(i)=tmpd(i);
                                    inx2(i)=i;
                                end
                            end
                            [mn,ij]=max(tpm2);
                            ii=inx2(ij);
                            
                            inx=zeros(1,n);
                            inx(1)=ii;
                            cnt=0;
                            for i=1:n
                                if (tpm2(i)==mn)
                                    inx(i)=inx2(i);
                                    cnt=cnt+1;
                                end
                            end
                            if(cnt>1)
                                r=randi(cnt,1);
                                j=0;
                                
                                for i=1:n
                                    if(inx(i)>0)
                                        j=j+1;
                                    end
                                    if(j==r)
                                        break
                                    end
                                end
                                id=nbl(inx(i));
                            else
                                id=ii;
                            end
                            
                            dseq=sum(A);
                            m=max(dseq);
                            
                            
                            idx=0;
                            for i=1:n
                                if (A(id,i)==1)
                                    idx=idx+1;
                                end
                            end
                            nbl=zeros(1,idx);
                            
                            nn=idx;
                            idx=0;
                            for i=1:n
                                if (A(id,i)==1)
                                    idx=idx+1;
                                    nbl(idx)=i;
                                end
                            end
                            
                            if(nn==0)
                                for i=1:cn
                                    if (ba(id,i)==0)
                                        ca(id)=i;
                                        break
                                    end
                                end
                                cntr=cntr+1;
                                order(cntr)=id;
                                A(id,:)=zeros(1,n);
                                A(:,id)=zeros(n,1);
                                for i=1:nn
                                    ba(nbl(i),ca(id))=1;
                                end
                            end
                        end
                        
                        cnt=cnt+1;
                        
                    end
                    
                    if(m==0)
                        for j=1:n
                            if(ca(j)==0)
                                for i=1:cn
                                    if (ba(j,i)==0)
                                        ca(j)=i;
                                        break
                                    end
                                end
                                cntr=cntr+1;
                                order(cntr)=j;
                            end
                        end
                    end
                end
                idx=1;
                for i=1:length(order)
                    flag=0;
                    if(i==1)
                        od(i)=order(i);
                    else
                        for j=1:i-1
                            if(order(i)==order(j))
                                flag=1;
                                break
                            end
                        end
                        if(flag==0)
                            idx=idx+1;
                            od(idx)=order(i);
                        end
                    end
                end
                
                
                crn=max(ca);
                if (crn<val)
                    val=crn;
                    colour=ca;
                    sequence=od;
                end
                
            end
            number(loop)=val;
        end
        if(min(number)==max(number))
            chnum=min(number);
%	fprintf('%s %d\n','val:',val)
            indicator=1;
        else
            ntest=ntest+counter*testi;
            counter=counter+1;
            ni=ni*counter;
            indicator=0;
            chnum=min(number);
            colour=ca;
            sequence=od;
        end
%	fprintf('%s %d\n','indicator:',indicator)
    end
end
%fprintf('%s %d\n','ntest:',length(number))
telp=toc;
fprintf('Elapsed time for the sub-routine: %f sec.\n',telp);
end
