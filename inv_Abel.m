function Data_inv = inv_Abel(Data)

I_pro = Data(:,(0.5*size(Data,2)+1):end);

DY = diff(I_pro,1,2);

DY(:,end+1) = zeros;

x = linspace(0,0.5*size(Data,1),0.5*size(Data,1));

Data_inv = zeros(size(I_pro));

for i = 1:(0.5*size(Data,1))
    
    int = -(1/pi)*(heaviside(x-i)).*DY./(sqrt(x.^2 - i.^2));
    
    Data_inv(:,i) = sum(int,2);
    
end

Data_inv(isinf(Data_inv)==1) = 0;
Data_inv(isnan(Data_inv)==1) = 0;

Data_inv = [fliplr(Data_inv), Data_inv];

end

