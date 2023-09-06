

function result = ABIC(mu,parm_num,un,nkm,sSa,sSb0)

if mu <= 0
    result = 1.0e10;
    return
end

tau = sqrt(mu);
sSb = tau*sSb0;
sS = cat(1,sSa,sSb);

[~,R,E] = qr(sS);

dhs = diag(R);
odh = dhs' * E;

F = odh(end) * odh(end); % F is the minimum value of ||Ap - b||^2 + mu||Bp - g||^2.
odhh = odh(1:end - 1);
dh = abs(odhh(odhh ~= 0));
detlA = 2 * sum(log(dh)); % detlA is 2log(det[A^T A + mu B^T B]).
detlB = un * log(mu);

ABIC = nkm + nkm * log(2.0 * pi() * F / nkm) + detlA - detlB + 2 * parm_num;

result = ABIC;
end
