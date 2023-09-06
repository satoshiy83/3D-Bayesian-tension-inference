% Functoin to get a parameter mu for the minimum ABIC.
% Written by Satoshi Yamashita.

function result = get_mu_for_min_ABIC(param_num,un,nkm,sSa,sSb0)
% Functoin to return mu = sigma^2/omega^2 for the minimum ABIC.
% result = get_mu_for_min_ABIC(parm_num,un,nkm,sSa,sSb0)
% Argument param_num is a number of hyper-parameter, 2 for sigma and omega.
% Argument un is a number of edges inside.
% Argument nkm is number of conditions + number of 0-singular value of
% balance matrix A + number of cells.
% Argument sSa = [A, 0].
% Argument sSb = [B, g], where B is an m times m matrix whose diagonal
% elements corresponding to edges are 1 and the rest are 0, g is an
% m-dimensional vector whose elements corresponding to edges are 1 and the
% rest are 0.
% Return value is mu for the minimum ABIC.

result = fminsearch(@ABIC,1);

    function sesult = ABIC(mu)
        if mu <= 0
            sesult = 1.0e10;
            return
        end

        tau = sqrt(mu);
        sSb = tau*sSb0;
        sS = cat(1,sSa,sSb);
        
        [~,R] = qr(sS);
        
        odh = diag(R);

        F = odh(end) * odh(end); % F is the minimum value of ||Ap - b||^2 + mu||Bp - g||^2.
        odhh = odh(1:end - 1);
        dh = abs(odhh(odhh ~= 0));
        detlA = 2 * sum(log(dh)); % detlA is 2log(det[A^TA + mu B^TB]).
        detlB = un * log(mu); % detlB is log(det[mu B^TB]).
        
        sesult = nkm + nkm * log(2.0 * pi() * F / nkm) + detlA - detlB + 2 * param_num;
    end
end
