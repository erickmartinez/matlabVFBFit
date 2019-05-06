function [res] = VFBSH_Residue(x)
    global Vsh; % The flatband shift
    global weights;
    Vsh_model = getVFBSH_Matlab(10.^x);
    res = (Vsh_model - Vsh).*weights;
end