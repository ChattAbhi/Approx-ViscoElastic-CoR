function thc=theta_gtilde(alp,gtil)
% theta_gtilde(alp,gtil) computes the value of $\theta$ for a corresponding
% value $\tilde{g}$, given a $\alpha$ value. The inputs alp and gtil
% represent the parameters $\alpha$ and $\tilde{g}$, respectively.
max_iter=100;
iter=0;

tol=1e-8; 

tha=eps; thb=1; thc=0.5;
while iter<=max_iter
    gtila=gtilde_theta(alp,tha);
    gtilb=gtilde_theta(alp,thb);
    gtilc=gtilde_theta(alp,thc);
    
    if abs(gtilc-gtil)<=tol
        break;
    else
        if gtilc>gtil
            if gtila>gtilc && gtilb<gtilc
                tha=thc;
            elseif gtilb>gtilc && gtila<gtilc
                thb=thc;
            else
                error('Something Wrong!')
            end
        elseif gtilc<gtil
            if gtila>gtilc && gtilb<gtilc
                thb=thc;
            elseif gtilb>gtilc && gtila<gtilc
                tha=thc;
            else
                error('Something Wrong!')
            end
        end
    end
    thc=(tha+thb)/2;
end

end