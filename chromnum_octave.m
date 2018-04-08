function [crn,colour,order] = find_chromatic_number(adj)
% Chromnum computes the chromatic number of an "UNDIRECTED" graph
% Chromatic Number of a graph is the minimum number of colors by whcih
% All nodes of the graph could be exhaustively mapped so that no two
% adjacent nodes share the same color.
% Chromnum implements the "trailpathSA" algorithm
% The function takes the 1-0 adjacency matrix stored in a text file as its sole input argument
% And returns the chromatic number of the corresponding graph
% Usage: crn = chromnum_octave('adj.inp')
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
% Authors: Abhirup Bandyopadhyay and Sankar Basu 
% 
% Email: Sankar Basu (nemo8130@gmail.com)
% Abhirup Bandyopadhyay (abhirupnit@gmail.com)
% 
% (c) Indian Copyright filed, Diary Number. 4955/2018-CO/SW, Filed. April, 2018
%
%

tic
if(nargin~=1)
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
            %==============================================================            
            %==============================================================
            [crn,colour,order]=chromatic_no_regcor(A0);
            %==============================================================
            %==============================================================

            le=sum(sum(A0))/2;
            ld = le/(n*(n-1)/2);
            
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
            
            fprintf('\nNumber of Nodes in the given undirected graph : %d',n);
            fprintf('\nNumber of edges : %d ',le);
	    fprintf('\nLink Density: %f ',ld);
            fprintf('\nChromatic Number : %d \n',crn);
            
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

            % Colormap is stored in a file and called to execute drawgraph 
            
            fprintf ('\nProposed Colormap:');
            fprintf ('\n------------------\n');
            delete colormap.inp;
            fidcm = fopen('colormap.inp','w');
            for i = 1:n
            fprintf(fidcm,'%f %f %f\n',GNCF(i,:));
            fprintf('Node-%d : color-%d\n',i,colour(i));
            end
            fclose(fidcm);
 
            %================ LEGEND ========================
            
          
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

            
            %=============== DRAWGRAPH ======================
            
            legend('show');
            title('Chromnum: Colour hierarchy is as follows as in the legend')                     
            drawgraph(adj,'colormap.inp');
            set (gca,'FontSize',15)
            axis equal 
        end
    end
    
end
toc
end

%=======================================================================================================

% Draw a graph in a circular layout given an adjacency matrix
% Color the nodes as specified in a colormap
% provide the matrix file_name as the input argument
% Usage: drawgraph('adj.inp','colormap.inp')
%
function drawgraph(adjfile,colormap)

adj = load (adjfile,'s');
colmap = load (colormap,'s');

sa = size(adj);
sc = size(colmap);
N = sa(1);
N2 = sc(1);

if (N ~= N2)
    fprintf('Dimension Missmatch between the adjacency matrix and the colormap\n');
    disp(N)
    disp(N2)
    return
end

deg2rad = pi/180.0;
rbig = 15.0;

XX = [];
YY = [];

for theta = 0:1:360
    theta_rad = theta*deg2rad;
    xcir = rbig*cos(theta_rad);
    ycir = rbig*sin(theta_rad);
    XX = [XX;xcir];
    YY = [YY;ycir];
end

%figure
hold on
%plot (XX,YY,'.')
axis equal

% Draw Nodes in a circular layout

deltheta = 360.0/N;

xdts = [];
ydts = [];

for i = 1:N
    thN = deltheta*i;
    thNrad = thN*deg2rad;
    xdt = rbig*cos(thNrad);
    ydt = rbig*sin(thNrad);
    plot(xdt,ydt,'o')
    cm = colmap(i,:);
    D=filled_circle(xdt,ydt,cm);
    xdts = [xdts;xdt];
    ydts = [ydts;ydt];
    text((xdt+0.5),(ydt+0.5),num2str(i),'FontSize',15)
end

% Draw Edges: Connect Nodes 

Nlinks = 0;

for i = 1:N
    for j = (i+1):N
        xp = [];
        yp = [];
        if (adj(i,j)==1)
            xp = [xp;xdts(i);xdts(j)];
            yp = [yp;ydts(i);ydts(j)];
            plot(xp,yp,'k-','LineWidth',2.5)
            Nlinks = Nlinks + 1;
        end
    end
end
        
%disp(Nlinks)        
       
end

function Y = filled_circle(xo,yo,cm)

    XXf = [];
    YYf = [];
    r = 0.5;

    deg2rad = pi/180.0;
    
    for th = 0:1:360
        th_rad = th*deg2rad;
        for rs = 0:0.01:r
            xc = rs*cos(th_rad);
            yc = rs*sin(th_rad);
            XXf = [XXf;(xc+xo)];
            YYf = [YYf;(yc+yo)];
        end
    end
    
    plot(XXf,YYf,'Color',cm,'Marker','.')

    Y=1;
    
end


%=======================================================================================================


function [chnum,colour,sequence] = chromatic_no_regcor(A0)
tstart=tic;
dim=size(A0);
dim=dim(2);

indicator=0;
ni=50;
ntest=25;
testi=10;
nmax=1000;
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


