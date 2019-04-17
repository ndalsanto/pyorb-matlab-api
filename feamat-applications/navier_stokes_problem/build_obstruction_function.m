function [cc] = build_obstruction_function( )

    cc = @(x, param) 1000 * ( ( 10*(x(1)-1.5).^2 + (x(2)-0.46).^2) <  param(2)^2 ) ...
                   + 10^-10 * ( ( 10*(x(1)-1.5).^2 + (x(2)-0.46).^2) <  0.4^2 );

end