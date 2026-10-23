function y = cochlear_map(x, direction)
%COCHLEAR_MAP Eq. (1), and its algebraic inverse. Angles are in degrees.
if nargin < 2, direction = 'frequency'; end
validateattributes(x, {'numeric'}, {'real','finite'});
switch lower(direction)
    case 'frequency'
        y = 165.4 .* (3251^2.1 .* (x+177.3).^(-2.1*1.149)-0.88);
    case 'angle'
        if any(x(:)<=-165.4*0.88), error('CT:Map','Invalid frequency.'); end
        y = (((x./165.4)+0.88)./3251^2.1).^(-1/(2.1*1.149))-177.3;
    otherwise
        error('CT:Map','Direction must be frequency or angle.');
end
end
