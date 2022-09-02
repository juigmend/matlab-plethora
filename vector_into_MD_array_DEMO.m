% place a vector into a multidimensional array DEMO

clc
test_data = 1:6;   % <--- Vector to put in any 2&UP dimension of a multidimensional array. The vector will be the first dimension.
MD2UP_dim = [4,2]; % <--- 2&UP dimensions of the multidimensional array. The first dimension is the length of the vector.
data_indx = [1,2]; % <--- INDEX of 2&UP dimensions of the multidimensional array where to put the vector.

data_length = length(test_data);
MD_dim = [data_length,MD2UP_dim];
MD_array = zeros(MD_dim); 
if ( length(MD2UP_dim) ~= length(data_indx) )
    disp('INDEX and 2&UP dimensions should have the same amount of elements')
    return
end
if sum(MD2UP_dim  <  data_indx) ~= 0
    disp('INDEX values have to be less or equal that 2&UP dimensions')
    return
end
beginning = 1;
for i = 1:length(data_indx)
    beginning = beginning + ( (data_indx(i)-1) * prod(MD_dim(1:i)) );
end
ending = data_length + beginning - 1;

MD_array(beginning:ending) = test_data % result

% ending = data_length;
% 
% % beginning = 1;
% 
% for i = 1:length(data_indx)
% %     beginning = beginning + data_indx(i)*sum(MD_dim(i));
% ending = ending + (data_indx(i)-1)*sum(MD2UP_dim(i));
% end
% 
% % ending = data_length + beginning - 1
% 
% beginning = ending - data_length + 1
% 
% MD_array(beginning:ending) = test_data % result

% beginning = 1;
% for i = 1:length(MD_dim)
%     beginning = beginning + (sum(MD_dim_minu(1:i)))*data_length - data_length;
% end
% 
% beginning
% ending = data_length + beginning - 1
