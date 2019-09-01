function [A]=constuct_dis_Matrix(x,y)
totalNum=length(x);
Edm=[];
for i=1:totalNum % for each cell
    for j=1:totalNum
        Edm(i,j)=(sqrt((x(i)-x(j))^2+(y(i)-y(j))^2));
    end
end
A=Edm;%>20;