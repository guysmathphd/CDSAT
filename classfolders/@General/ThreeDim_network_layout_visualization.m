function hax = ThreeDim_network_layout_visualization(d,v,n,fname,hax)
%% Input of the function 
% d is the degree vector of the network
% v is the quantity on the z axis
% n is the number of nodes (Please try to give n maximum 1000)


%% You may change the color codes from here
c1 = [102, 205, 170]/255; c2 = [0, 102, 51]/255;
%c1=[250,128,114]/255;c2=[153,0,0]/255; %hub-loc

[dval, ~] = max(d);
rng(1);

%% Node coordinate creation
xi = zeros(n,1);
yi = zeros(n,1);
 for i=1:n
   if dval== d(i)
     r = 0; % Keeping hub node to the center
   else
     r = (1/d(i))+1./(exp(-1./(d(i)))+exp(-1./log10(d(i)))); %(1/(4+exp(-1/(d(i))))); %(1/(1+exp(-1/log10(d(i)))));  %(1/(1+exp(-d(i)^(1/2)))); %exp(-d(i)^(1/2)); 
   end  
   a = randi([0, 360],1);
   xi(i) = r*cos((pi*a)/180);
   yi(i) = r*sin((pi*a)/180);
 end

% After generating the coordinate make more spaceing between the nodes 
m1 = xi.*20; 
m2 = yi.*20; 

%% Scaling the node size in the coordinate axes
low_d = 1;
high_d = 5;
md = low_d+(high_d-low_d).*(d./max(d));

clear xi yi  

%% Creating color coding based on the node degree but we can change based on our requirement
depth = size(unique(d),1);
grad = [linspace(c1(1),c2(1),depth)', linspace(c1(2),c2(2),depth)', linspace(c1(3),c2(3),depth)'];

%% iterate over all the nodes
for i=1:n   
    mu = [m1(i);m2(i)];
    if d(i) == 2   
        %r = d(i);
        r = md(i);  % scaling the degree of a node
        sigma = 0.8;  % sqrt(r/2)*log2(r); % how to decide sigma as a function of r still unsolve and I am trying
        fcolor = grad(1,:); % assigning color 
        x1 = mu(1)-r:0.01:r+mu(1);  % deciding the 2-D surface plot area
        x2 = mu(2)-r:0.01:r+mu(2);  
    elseif d(i)==3  
        %r = d(i);
        r = md(i);
        sigma = 0.8;  % sqrt(r/2)*log2(r);
        fcolor = grad(2,:);
        x1 = mu(1)-r:0.1:r+mu(1);
        x2 = mu(2)-r:0.1:r+mu(2);
    elseif d(i) == 4  
        %r = d(i)/2;
        r = md(i);
        sigma = 0.8;  % sqrt(r/2)*log2(r);
        fcolor = grad(3,:);
        x1 = mu(1)-r:0.1:r+mu(1);
        x2 = mu(2)-r:0.1:r+mu(2);
    elseif d(i)==5  
        %r = d(i)/2;
        r = md(i);
        sigma = 0.8;  % sqrt(r/2)*log2(r);
        fcolor = grad(4,:); 
        x1 = mu(1)-r:0.1:r+mu(1);
        x2 = mu(2)-r:0.1:r+mu(2);
    elseif d(i)== 6  
        %r = d(i)/3;
        r = md(i);
        sigma = 0.8;  % sqrt(r/2)*log2(r);
        fcolor = grad(5,:);
        x1 = mu(1)-r:0.1:r+mu(1);
        x2 = mu(2)-r:0.1:r+mu(2);
    elseif d(i)== 7  
        %r = d(i)/3;
        r = md(i);
        sigma = 0.8;  % sqrt(r/2)*log2(r);
        fcolor = grad(6,:);
        x1 = mu(1)-r:0.1:r+mu(1);
        x2 = mu(2)-r:0.1:r+mu(2);
    elseif d(i) == 8  
        %r = d(i)/3;
        r = md(i);
        sigma = 0.8;  % sqrt(r/2)*log2(r);
        fcolor = grad(7,:);
        x1 = mu(1)-r:0.1:r+mu(1);
        x2 = mu(2)-r:0.1:r+mu(2);
        
    elseif d(i) == 9  
        %r = d(i)/3;
        r = md(i);
        sigma = 0.8;  % sqrt(r/2)*log2(r);
        fcolor = grad(8,:); 
        x1 = mu(1)-r:0.1:r+mu(1);
        x2 = mu(2)-r:0.1:r+mu(2);
    elseif d(i) == 10  
        %r = d(i)/3;
        r = md(i);
        sigma = 0.8; % sqrt(r/2)*log2(r);
        fcolor = grad(9,:); 
        x1 = mu(1)-r:0.5:r+mu(1);
        x2 = mu(2)-r:0.5:r+mu(2);
       
    elseif  10 < d(i) && d(i)< 20
       %r = d(i)/4;
       r = md(i);
       sigma = sqrt(r/3)*log2(r);
       fcolor = grad(10,:); 
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);    
       
    elseif d(i)== 20
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(11,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);
       
   elseif d(i)== 21
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(12,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);
 
    elseif d(i)== 22
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(13,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);
       
   elseif d(i)== 23
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(14,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);    
    
    elseif d(i) == 24  
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(15,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);   
    
    elseif d(i) == 25  
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(16,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);   
    
    elseif d(i) == 26 
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(17,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);   
    
    elseif d(i) == 27 
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(18,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);  
    
    elseif d(i) == 28 
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(19,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);    
       
    elseif d(i) == 29 
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(20,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2); 
    
    elseif  29 < d(i) && d(i)< 41
       %r = d(i)/4;
       r = md(i);
       sigma = 8; %sqrt(r/3)*log2(r);
       fcolor = grad(21,:); 
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);      
    
    elseif  40 < d(i) && d(i)< 51
       %r = d(i)/4;
       r = md(i);
       sigma = 8; %sqrt(r/3)*log2(r);
       fcolor = grad(22,:); 
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);         
       
    elseif d(i) == 60  
       %r = d(i)/3;
       r = md(i);
       sigma = 8; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(23,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);  
    else 
%        fprintf('special: %d\n',d(i));
       %r = d(i)/14;
       r = md(i);
       sigma =3; % sqrt(r/2)*log2(r)+2;
       fcolor = grad(24,:);   
       x1 = mu(1)-r:0.1:r+mu(1);
       x2 = mu(2)-r:0.1:r+mu(2);        
    end
    
    [X,Y] = meshgrid(x1,x2);
    clear x1 x2;
    temp = [X(:)-mu(1) Y(:)-mu(2)];
    temp1 = temp(:,1).*temp(:,1)+temp(:,2).*temp(:,2);
    clear temp
    pdf = v(i)*exp(-(1/sigma)*temp1);
    clear temp1
    %pdf = v(i)*exp(-(1/sigma)*diag(temp*temp'));
    
    %clear temp
    pdf = reshape(pdf,size(X));
    
    n1 = size(pdf,1);
    for i1=1:n1
      for j1=1:n1
        dist = sqrt((X(i1,j1)-mu(1))^2+(Y(i1,j1)-mu(2))^2);
          if dist>r
            pdf(i1,j1) = NaN;
          end
      end
    end
%     f = figure(20);f.Name = fname;f.NumberTitle = 'off';
    
    s = surf(hax,X, Y, pdf);
    alpha 0.1;
    s.FaceColor = fcolor;
    s.EdgeColor = 'None';
    grid off
    box on
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    hold on
      
    view(55,15)

    clear X Y pdf
end
set(gca,'FontSize',25);
set(gca,'linewidth',1.5)
zlim([-0.01 1.0])
clear  m1 m2

end