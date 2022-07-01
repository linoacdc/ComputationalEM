clear all;
x0 = 0; y0=0; bb=1.75e-3; aa=0.76e-3; e0 = 8.85418782e-12;
A = 2e-2; 
gd = [1 1;
      x0 x0;
      y0 y0;
      bb aa];
ns = [82 82; 49 50]; sf = 'R1-R2'; 
d1 = decsg(gd,sf,ns);
%d1 = decsg(gd);
[p,e,t] = initmesh(d1);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);

pdeplot(p,e,t); axis equal; axis tight;



Nn = size(p,2); % Number of nodes
Ne = size(t,2); % Number of elements (triangles)
Nd = size(e,2); % Number of (boundary) edges

node_id = ones(Nn,1); % Inialization of node flag (1 is unknown, 0 is Dirichlet node)
X0 = zeros(Nn,1); V0 = 100;

for id = 1:Nd
    n(1:2) = e(1:2,id); x(1:2) = p(1,n(1:2)); y(1:2) = p(2,n(1:2));  
    rg1 = e(6,id); rg2 = e(7,id);
    if (rg1==0 || rg2==0)        
        node_id(n(1))=0; node_id(n(2))=0;
        radius(1:2) = sqrt(x(1:2).^2 + y(1:2).^2);
        if (radius(1) > (aa+bb)/2) X0(n(1)) = 0; else X0(n(1)) = 100; end;
        if (radius(2) > (aa+bb)/2) X0(n(2)) = 0; else X0(n(2)) = 100; end;   
    end;
end

% Define unknown numbering
ic = 0; %Define counter to count unknowns
index = zeros(Nn,1); % Define index vector with unknown's numbering for each node
for in=1:Nn 
    if (node_id(in) == 1) % Node is unknown
        ic = ic + 1; 
        index(in) = ic;
    end;
end;
Nf = ic % Total number of unknowns

%for in=1:Nn text(p(1,in),p(2,in),num2str(index(in))); end;

S = spalloc(Nf,Nf,7*Nf); B = zeros(Nf,1);

for ie = 1:Ne % Scan all elements
    n(1:3) = t(1:3,ie); % n(1), n(2), n(3) are the three nodes of element ie
    rg = t(4,ie); % Element region  
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3)); % Nodes' coordinates
    D = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    b(1) = (y(2)-y(3))/D; b(2) = (y(3)-y(1))/D; b(3) = (y(1)-y(2))/D;
    c(1) = (x(3)-x(2))/D; c(2) = (x(1)-x(3))/D; c(3) = (x(2)-x(1))/D;
    Ae = abs(D)/2; % Element area
    for i=1:3
        for j=1:3
            Se(i,j) = (b(i)*b(j) + c(i)*c(j))*Ae;
            if (node_id(n(i)) == 1)
                if (node_id(n(j)) == 1)
                    S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
                else
                    B(index(n(i))) = B(index(n(i))) - Se(i,j)*X0(n(j));
                end
            end
        end
    end
end

tic;
X = S\B; % System solution via Gauss-Jordan elimination
%X=gmres(S,B,100);
%X=bicg(S,B,1e-6,100);
%ms = (toc * 1000)
toc

for in=1:Nn
    if (node_id(in) == 1)
        X0(in) = X(index(in));
    end
end

% pdeplot(p,e,t,'xydata',X0,'contour','off','mesh','off'); axis equal; axis tight; hold on; colormap jet;

X0_exact = zeros(Nn,1); 
c1 = V0/log(aa/bb); c2 = -c1*log(bb);
for in=1:Nn
    x = p(1,in); y = p(2,in);
    rho = sqrt(x^2+y^2);
    X0_exact(in) = c1*log(rho) + c2;
end


figure;     
pdeplot(p,e,t,'xydata',X0_exact,'contour','off','mesh','off'); axis equal; axis tight; hold on; colormap jet;

Canalytic = 2*pi*e0/log(bb/aa)      %Calculate C analytically

[Ex,Ey] =  pdegrad(p,t,-X0);        %Get Ex, Ey





Welement = zeros(Ne,1);


for ie = 1:Ne
    n(1:3) = t(1:3,ie);
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3));
    D = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    a(1) = (x(2)*y(3) - x(3)*y(2))/D; a(2) = (x(1)*y(3) - x(3)*y(1))/D; a(3) = (x(1)*y(2) - x(2)*y(1))/D;
    b(1) = (y(2)-y(3))/D; b(2) = (y(3)-y(1))/D; b(3) = (y(1)-y(2))/D;
    c(1) = (x(3)-x(2))/D; c(2) = (x(1)-x(3))/D; c(3) = (x(2)-x(1))/D;
    
    for k = 1:3
        z(k) = a(k) + b(k)* x(k) + c(k)*y(k);
    end
    Ae = abs(D)/2; % Element area

    
    %Calculate W for this element
    for i = 1:3
        for j = 1:3
            Se(i,j) = (b(i)*b(j) + c(i)*c(j))*Ae;       %Local S
            Welement(ie) = Welement(ie) + X0(n(i))* Se(i,j) * X0(n(j)); %Add contribution from nodes i, j
        
        end 
    end
    Welement(ie) = 1/2 * e0 * Welement(ie);     %Make W the energy of the whole element

end

Wtotal = sum(Welement);     %Sum the energy of all elements
C = 2 * Wtotal/ (V0^2)      %Calculate C

error = abs((Canalytic - C) / Canalytic)        %Calculate error %

figure(3);
pdeplot(p,e,t,'xydata',X0,'FlowData',[Ex;Ey]);      %Plot E



 
 
