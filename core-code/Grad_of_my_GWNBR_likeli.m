function [result,g] = Grad_of_my_GWNBR_likeli(x,u,N,w,G,y,i)
% Anotation:
% u y N G w i are the global coefficients,x is the parameters to be solved.
% w is weighted matrix,u is scRNA features,y is spatial matrix, 
% N is the spots num,G is the genes num, i is target index,
% n is the nums of neighbors of the destination spot, j is gene index.

if nargout > 1 % with jacob
    xindex=size(u,1);
    result=0;
    theta_result=0;
    beta_result=zeros(xindex,1);
    for n=1:size(N,1)
        % neighbors or not
        if w(i,n)==0
            continue;
        end
        %if true,start computing likelihood
        neihelog=0;
        theta_neihelog=0;
        beta_neihelog=zeros(xindex,1);
        for j=1:G
            nsum=dot(x(1:xindex),u(:,j));
            first=0;
            first_3=0;
            if y(n,j)==0
                first=0;
                first_3=0;
            else
                temp=0:(y(n,j)-1);
                first=sum(log(temp+x(xindex+1)));
                first_3=sum(1./(temp+x(xindex+1)));
               
            end
            neihelog=neihelog+first+...
                x(xindex+1)*(log(x(xindex+1))-log(x(xindex+1)+N(n)*nsum))+...
                y(n,j)*(log(N(n)*nsum)-log(N(n)*nsum+x(xindex+1)));
            theta_neihelog=theta_neihelog+first_3+log(x(xindex+1))-log(x(xindex+1)+N(n)*nsum)+1-(x(xindex+1)+y(n,j))/(x(xindex+1)+N(n)*nsum);
            for k=1:xindex
                beta_neihelog(k)=beta_neihelog(k)+(-1)*x(xindex+1)*N(n)*u(k,j)/(x(xindex+1)+N(n)*nsum)+y(n,j)*(u(k,j)/nsum-N(n)*u(k,j)/(N(n)*nsum+x(xindex+1)));
            end
        end
        %weighted
        result=result+neihelog*w(i,n);
        theta_result=theta_result+theta_neihelog*w(i,n);
        beta_result=beta_result+beta_neihelog*w(i,n);
    end
    result=(-1)*result;
    theta_result=(-1)*theta_result;
    beta_result=(-1)*beta_result;
    g=[beta_result;theta_result];

end
end

