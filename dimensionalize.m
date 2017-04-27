function [cc] = dimensionalize(c)
global qdp h Tinf k L D dx dy row col nodetype dxnd dynd
    cc = c*Tinf;
end