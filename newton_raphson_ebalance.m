
function Xs = newton_raphson_ebalance(X0,Err,imax,Lin_coarse,...
    Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,emissivity,Cp,xLs,...
    rho_air,De_h,B1,B2,Kd,mode,tau,cc,opt_input)
% matlab coder doesn't support anonymous function hence the repeated code
Xs=NaN;
Xi = NaN(1,1,'single');

switch mode
    case 'normal'
        Xest=X0;
        Xs=X0;
        for i = 1:imax
            Xi(1) = single(Xest - ebalance(Xest,Lin_coarse,...
            Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,emissivity,Cp,xLs,...
            rho_air,De_h,B1,B2,Kd,mode,tau,cc,[])/...
                d_ebalance(Xest,T_fine,Vf,pres_fine,ea,sig,emissivity,...
                Cp,xLs,rho_air,De_h,B1,B2,tau,cc));
            if abs((Xi - Xest)/Xest) < Err
                Xs = Xi;
                break
            end
            Xest = Xi;
        end
    case 'debris'
        d=opt_input;
        Xest=X0;
        Xs=X0;
        for i = 1:imax
            Xi(1) = Xest - ebalance(Xest,Lin_coarse,...
                Zdiff,T_fine,Vf,albedo,pres_fine,Sin,ea,sig,...
                emissivity,Cp,xLs,...
                rho_air,De_h,B1,B2,Kd,mode,tau,cc,d)/...
                d_ebalance_debris(Xest,T_fine,Vf,sig,emissivity.ems,...
                Cp,rho_air,De_h,B1,B2,Kd,d);
            if abs((Xi - Xest)/Xest) < Err
                Xs = Xi;
                break
            end
            Xest = Xi;
        end
end
