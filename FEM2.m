clear all;
x0 = 0; y0 = 0; w=4e-2; h=2e-3;d=1e-2; V0=100; A=9*w; B=A; e0 = 8.85418782e-12;
er = [2.2,1];

gd = [3 3 3 3;
      4 4 4 4;
      -A/2 w/2 w/2 w/2;
      A/2 w/2 w/2 w/2;
      A/2 -w/2 -w/2 -w/2;
      -A/2 -w/2 -w/2 -w/2;
      -B/2 d/2 -d/2 -d/2;
      -B/2 d/2+h -d/2-h d/2;
      B/2 d/2+h -d/2-h d/2;
      B/2 d/2 -d/2 -d/2];
ns = [82 82 82 82; 49 50 51 52]; sf = 'R1-R2-R3+R4';    %Remove parallel plates from mesh
d1 = decsg(gd,sf,ns);
%d1 = decsg(gd);
[p,e,t] = initmesh(d1);
[p,e,t] = refinemesh(d1,p,e,t);
[p,e,t] = refinemesh(d1,p,e,t);
%[p,e,t] = refinemesh(d1,p,e,t);

figure(4); pdeplot(p,e,t); axis equal; axis tight;



Nnodes = size(p,2); % Number of nodes
Nelements = size(t,2); % Number of elements (triangles)
Nedges = size(e,2); % Number of (boundary) edges

node_id = ones(Nnodes,1); % Inialization of node flag (1 is unknown, 0 is Dirichlet node)
X0 = zeros(Nnodes,1); 

for id = 1:Nedges
    n1 = e(1,id); n2 = e(2,id); x1 = p(1,n1); y1 = p(2,n1);
    x2 = p(1,n2); y2 = p(2,n2);
    r1 = e(6,id); r2 = e(7,id);
    if (r1==0 || r2==0) && (y1 > 0 && y1 < B/4 && abs(x1)<1.2*w/2 && y2 > 0 && y2 < B/4 && abs(x2)<1.2*w/2)
        node_id(n1)=0; node_id(n2)=0; X0(n1)=V0/2; X0(n2)=V0/2;
    end
    if (r1==0 || r2==0) && (y1 < 0  && y1 > -B/4 && abs(x1)<1.2*w/2 && y2 < 0 && y2 > - B/4 && abs(x2)<1.2*w/2)
        node_id(n1)=0; node_id(n2)=0; X0(n1)=-V0/2; X0(n2)=-V0/2;
    end
end





ic = 0; %Define counter to count unknowns
index = zeros(Nnodes,1); % Define index vector with unknown's numbering for each node
for in=1:Nnodes 
    if (node_id(in) == 1) % Node is unknown
        ic = ic + 1; 
        index(in) = ic;
    end;
end;
Nf = ic % Total number of unknowns


% for in = 1:Nnodes
%     x = p(1,in); y = p(2,in);
%     text(x,y,num2str(index(in)))
% end

% for ie = 1:Nelements
%     n(1:3) = t(1:3,ie); x(1:3) = p(1,n(1:3)); y(1:3) = p(2, n(1:3));
%     
%     region = t(4,ie);
%     xg=sum(x)/3;
%     yg=sum(y)/3;
%     text(xg,yg,num2str(region));
% end


S = spalloc(Nf,Nf,7*Nf); B = zeros(Nf,1);


for ie = 1:Nelements % Scan all elements
    n(1:3) = t(1:3,ie); % n(1), n(2), n(3) are the three nodes of element ie
    rg = t(4,ie); % Element region  
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3)); % Nodes' coordinates
    D = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    b(1) = (y(2)-y(3))/D; b(2) = (y(3)-y(1))/D; b(3) = (y(1)-y(2))/D;
    c(1) = (x(3)-x(2))/D; c(2) = (x(1)-x(3))/D; c(3) = (x(2)-x(1))/D;
    Ae = abs(D)/2; % Element area
    er1 = er(rg);
    for i=1:3
        for j=1:3
            Se(i,j) = er1*(b(i)*b(j) + c(i)*c(j))*Ae;
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

X = S\B;

for in=1:Nnodes
    if (node_id(in) == 1)
        X0(in) = X(index(in));
    end
end

Welement = zeros(Nelements,1);

for ie = 1:Nelements
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

    
    
    for i = 1:3
        for j = 1:3
            Se(i,j) = (b(i)*b(j) + c(i)*c(j))*Ae;
            Welement(ie) = Welement(ie) + X0(n(i))* Se(i,j) * X0(n(j)); 
        
        end 
    end
    Welement(ie) = 1/2 * e0 * Welement(ie);

end

Wtotal = sum(Welement);
C = 2 * Wtotal/ (V0^2)

Canalytic = e(1) * e0 * w / d
error = abs((Canalytic - C) / Canalytic)

figure ;pdeplot(p,e,t,'xydata',X0,'contour','on','mesh','off'); axis equal; axis tight; hold on; colormap jet;

[Ex,Ey] =  pdegrad(p,t,-X0);

figure(3);
pdeplot(p,e,t,'xydata',X0,'FlowData',[Ex;Ey]); axis equal; axis tight; hold on; colormap jet;

