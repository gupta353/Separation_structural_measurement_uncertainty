% Log of Beta pdf (a faster compuattion than that of MATALB's in-built function)
% inputs: smp = point in support at which log of density is to be computed
%         alpha  = first shape parameter
%         betav = second shape parameter
% outputs: log_pdf = log of beta density at smp

function log_pdf = betapdfAG(smp,alpha,betav)
    if smp<=1 && smp>=0
        log_pdf = gammaln(alpha+betav) - gammaln(alpha) - gammaln(betav) + (alpha-1)*log(smp) + (betav-1)*log(1-smp);
    else
        log_pdf = -inf;
    end
end