clear all;

%TM

A = 2e-2; x0=0; y0=0;
gd = [1 ;
      x0 ;
      y0 ;
      A ];
 
d1 = decsg(gd);

[p,e,t] = initmesh(d1);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
figure;
%pdeplot(p,e,t); axis equal; axis tight;



Nn = size(p,2); % Number of nodes
Ne = size(t,2); % Number of elements (triangles)
Nd = size(e,2); % Number of (boundary) edges



node_id = ones(Nn,1); % Inialization of node flag (1 is unknown, 0 is Dirichlet node)
X0 = zeros(Nn,1);

for id = 1:Nn
    x1 = p(1,id) ; y1= p(2,id);
    
    if(sqrt(x1^2 + y1^2) >= 2e-2 - 1e-6)
        node_id(id)=0; 
        
    end

end


% % Define unknown numbering
index = zeros(Nn,1); % Define index vector with unknown's numbering for each node

ic = 0; %Define counter to count unknowns
for in=1:Nn 
    if (node_id(in) == 1) % Node is unknown
        ic = ic + 1; 
        index(in) = ic;
    end
end
Nf = ic; % Total number of unknowns


for in=1:Nn text(p(1,in),p(2,in),num2str(index(in))); end;


S = spalloc(Nf,Nf,7*Nf);  T = spalloc(Nf,Nf,7*Nf);

for ie = 1:Ne % Scan all elements
    n(1:3) = t(1:3,ie); % n(1), n(2), n(3) are the three nodes of element ie
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3)); % Nodes' coordinates
    D = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    b(1) = (y(2)-y(3))/D; b(2) = (y(3)-y(1))/D; b(3) = (y(1)-y(2))/D;
    c(1) = (x(3)-x(2))/D; c(2) = (x(1)-x(3))/D; c(3) = (x(2)-x(1))/D;
    Ae = abs(D)/2; % Element area
    for i=1:3
        for j=1:3
            Se(i,j) = (b(i)*b(j) + c(i)*c(j))*Ae;
                if(i==j)
                    Te(i,j) = Ae/6;
                else
                    Te(i,j) = Ae/12;
                end
        
            if (node_id(n(i))~= 0 && node_id(n(j)) ~=0)
                
                    S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
                    T(index(n(i)), index(n(j))) =  T(index(n(i)), index(n(j))) + Te(i,j);
            end
                
        end          
    end
end

 [V,D] = eigs(S,T,6,0);

for k = 1:6
    X1 = V(:,k);
    for in = 1:Nn
        if(node_id(in) ~= 0)
            X0(in)=X1(index(in));
        end
    end
    
    [Ex,Ey] =  pdegrad(p,t,-X0);
    figure; pdeplot(p,e,t,'xydata',X0,'contour','on','FlowData',[Ex;Ey]); axis equal; axis tight; hold on; colormap jet;

    %figure; pdeplot(p,e,t,'xydata',X0, 'contour', 'on', 'mesh', 'off'); axis equal; axis tight; hold on; colormap jet;
end


 
 
