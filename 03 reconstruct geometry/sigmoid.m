% Compute sigmoid function
function sig = sigmoid(varargin)

if nargin == 1
    z = varargin{1};
    a = 1;
elseif nargin == 2
    z = varargin{1};
    a = varargin{2};
end

sig = 1./(1+exp(-a*z));

end