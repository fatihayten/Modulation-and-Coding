function [bits] = random_bit_generator(number_of_bits)
% This functions generates random bits, at amount of input
bits = randi([0 1],number_of_bits,1);
end