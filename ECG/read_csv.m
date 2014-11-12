function [X] = read_csv (ruta,N,L,offset)
data = csvread(ruta,2,0);
X = data(1+offset:offset+N,2:L+1);
end