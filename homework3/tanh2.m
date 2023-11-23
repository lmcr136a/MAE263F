function A = tanh2(Z)
% Evelyn Kim, 706180341
% 
% tanh2: a function that applies the tanh activation function on the input to get a nonlinear output.
%    Inputs:
%         Z: A M x N matrix representing the output of the neurons, 
%       which serves as the input of the activation function. 
%       M is the number of neurons, and N is the number of examples
%   Outputs:
%         A: a M x N matrix representing the output after the tanh activation function

    A = 2./(1+exp(-2*Z))-1;
end