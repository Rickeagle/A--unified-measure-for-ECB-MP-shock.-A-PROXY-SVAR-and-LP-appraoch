function [Y X] = LagY(y,p)
         Y=y(p+1:end,:);
         T=size(Y,1);
         X=ones(T,1);
        for j=1:p
            X=[X y(p+1-j:end-j,:)]; 
        end  
        
