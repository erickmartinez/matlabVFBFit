function [VFB_SH] = getVFBSH_Matlab_callback(x,timeArray)
    global thickness;
    global tempC;
    global VBias;
    global simulationTime;
    global time_exp;
    global area_m;
    global CN;
%     global C0;
    D   = 10.^x(1);
    C0   = 10.^x(2);
    
    [C_t,VFB,time_fd,depth_um] = FDNP_SiNxDevice1Dirichlet(D,C0,thickness,tempC,VBias,simulationTime,...
    CN,area_m);
    
    %Interpolate simulated time and flatband array from
    %experimental time
   
    VFB_SH = (interp1(time_fd, VFB, time_exp));
    if isrow(VFB_SH)
        VFB_SH = VFB_SH';
    end
end
    