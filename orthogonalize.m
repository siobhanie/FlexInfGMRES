function [h,y] = orthogonalize(Q,w)
    h=Q'*w;
    y = w - Q*h;
    g = Q'*y;
    y = y - Q*g;
    h = h + g;
end
