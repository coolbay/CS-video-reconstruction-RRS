
%


function y = BCS_Encoder(current_image, Phi, block_size)

x = im2col(current_image, [block_size block_size], 'distinct');

y = Phi * x;
