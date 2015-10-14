module Bernstein

export berntuples, flattuples, stroud_eval1d, stroud_eval2d, evalmat

iter2list(a) = [a[i] for i=1:length(a)]


function berntuples(sdim, deg)
    @assert sdim >= 1
    @assert deg >= 0
    if sdim == 1
        result = [(i, deg-i) for i=0:deg]
    else
        result = []
        for i=0:deg
            bts = berntuples(sdim-1, deg-i)
            for bt in bts
                item = ntuple((j)->vcat([i], iter2list(bt))[j],
                              1+length(bt))
                push!(result, item)
            end
        end
    end
    return result
end

function flattuples(sdim, num)
    @assert sdim >= 1 && num > 0
    if sdim == 1
        result = [(i,) for i=1:num]
    else
        result = []
        for i=1:num
            fts = flattuples(sdim-1, num)
            for ft in fts
                push!(result,
                      ntuple((j)->vcat([i], iter2list(ft))[j],
                             1+length(ft)))
            end
        end
    end
    return result        
end

function berntuple_indices(sdim, deg)
    bts = berntuples(sdim, deg)
    D = Dict()
    for (i, bt) in enumerate(bts)
        D[bt] = i
    end
    return D
end

multi_binomial(a, b) = prod([binomial(ai, bi) for (ai, bi) in zip(a,b)])
polydim(sdim, deg) = binomial(sdim + deg, sdim)

bern_dumb_1d(deg, i, x) = binomial(deg, i) * x^i * (1-x)^(deg-i)

# this is transposed from Python!
function evalmat(sdim, deg, qpts)
    nqp1d = size(qpts, 1)
    mat = zeros(nqp1d^sdim, polydim(sdim, deg))
    alphas = berntuples(sdim, deg)
    Is = flattuples(sdim, nqp1d)
    for j=1:length(alphas)
        alpha = alphas[j]
        offset = cumsum(vcat([0],iter2list(alpha[1:length(alpha)-1])))
        for i=1:length(Is)
            mat[i,j] = prod([bern_dumb_1d(deg-offset[k],
                                          alpha[k],
                                          qpts[Is[i][k], k])
                             for k=1:sdim])
        end
    end

    return mat                
end

# Note: can't work with 0 as a qpt.
function stroud_eval1d(deg, coeffs, result, qpts)
    nqp1d = size(qpts, 1)
    for i=1:nqp1d
        result[i] = 0.0
        xi = qpts[i, 1]
        s = 1-xi
        r = xi / s
        w = s^deg
        for alpha=0:deg
            result[i] += w * coeffs[alpha+1]
            w *= r * (deg-alpha) / (alpha+1.0)
        end
    end
    return
end

function stroud_eval2d(deg, coeffs, result, qpts)
    nqp1d = size(qpts, 1)
    result[:] = 0.0
    tmp = zeros(nqp1d, deg+1)
    for i2=1:nqp1d
        xi = qpts[i2, 2]
        s = 1 - xi
        r = xi / s
        for alpha1=0:deg
            w = s^(deg-alpha1)
            aind = polydim(2,deg) - polydim(2, deg-alpha1) + 1
            for alpha2=0:deg-alpha1
                tmp[i2, alpha1+1] += w * coeffs[aind+alpha2]
                w = w * r * (deg-alpha1-alpha2) / (1.0 + alpha2)
            end
        end
    end
    for i1=1:nqp1d
        xi = qpts[i1, 1]
        s = 1-xi
        r = xi / s
        w = s^deg
        for alpha1=0:deg
            for i2=1:nqp1d
                result[(i1-1)*nqp1d + i2] += w*tmp[i2, alpha1+1]
            end
            w = w * r * (deg-alpha1) / (1. + alpha1)
        end
    end
    
    return
end

function stroud_eval3d(deg, coeffs, result, qpts)
    nqp1d = size(qpts, 1)
    tmp = zeros(nqp1d, deg+1, deg+1, 2)
    for i=1:length(result)
        result[i] = 0.0
    end
    for i3=1:nqp1d
        xi=qpts[i3, 3]
        s = 1-xi
        r = xi/s
        for alpha1=0:deg
            for alpha2=0:(deg-alpha1)
                w = s^(deg-alpha1-alpha2)
                aind = polydim(3, deg) - polydim(3, deg-alpha1) + polydim(2, deg-alpha1) - polydim(2, deg-alpha1 - alpha2) + 1
                for alpha3=0:(deg-alpha1-alpha2)
                    tmp[i3, alpha2+1, alpha1+1, 1] += w * coeffs[aind+alpha3]
                    w = w * r * (deg-alpha1-alpha2-alpha3)/(1.+alpha3)
                end
            end
        end
    end

    for i2=1:nqp1d
        xi = qpts[i2, 2]
        s = 1-xi
        r = xi/s
        for alpha1=0:deg
            w=s^(deg-alpha1)
            for alpha2=0:(deg-alpha1)
                for i3=1:nqp1d
                    tmp[i3, i2, alpha1+1, 2] += w * tmp[i3, alpha2+1, alpha1+1, 1]
                end
                w = w * r * (deg-alpha1-alpha2) / (1. + alpha2)
            end
        end
    end

    for i1=1:nqp1d
        xi=qpts[i1, 1]
        s = 1-xi
        r = xi/s
        w=s^deg
        for alpha1=0:deg
            for i2=1:nqp1d
                for i3=1:nqp1d
                    result[(i1-1)*nqp1d^2 + (i2-1)*nqp1d + i3] += w * tmp[i3,i2,alpha1+1,2]
                end
            end
            w = w * r * (deg-alpha1)/(1.+alpha1)
        end
    end

end

function stroud_integrate1d(deg, vals, result, qpts, qwts)
    nqp1d = size(qwts, 1)
    result[:] = 0.0
    for i=1:nqp1d:
        xi = qpts[i,1]
        s = 1 - xi
        r = xi / s
        w = s^deg
        for alpha=0:deg
            result[alpha+1] += qwts[i,1] * w * vals[i]
            w = w * r * (deg - alpha) / (1. + alpha)
        end
    end

    return
end

function stroud_integrate2d(deg, vals, result, qpts, qwts)
    nqp1d = size(qwts, 1)
    tmp = zeros(nqp1d, deg+1)
    result[:] = 0.0
    
    for i1=1:nqp1d
        xi = qpts[i1, 1]
        s = 1.0 - xi
        r = xi / s
        w = (s^deg) * qwts[i1, 1]
        for a1=0:deg
            for i2=1:nqp1d
                tmp[i2, a1+1] += w * vals[(i1-1)*nqp1d+i2]
            end
            w = w * r * (deg-a1) / (1. + a1)
        end
    end

    for i2=1:nqp1d
        xi = qpts[i2, 2]
        s = 1.0 - xi
        r = xi / s
        for a1=0:deg
            w = s^(deg-a1)*qwts[i2, 2]
            aind = polydim(2, deg) - polydim(2, deg-a1) + 1
            for a2=0:deg-a1
                result[aind+a2] += w * tmp[a1+1, i2]
                w = w * r * (deg-a1-a2) / (1. + a2)
            end
        end

    return
end

function stroud_eval3d(deg, vals, result, qpts, qwts)
    nqp1d = shape(qwts, 1)
    tmp = zeros(nqp1d, deg+1, deg+1, 2)
    result[:] = 0.0
    
    for i1=1:nqp1d
        xi = qpts[i1, 1]
        s = 1.0 - xi
        r = xi / s
        w = qwts[i1, 1] * s^deg
        for a1=0:deg
            for i2=1:nqp1d
                for i3=1:nqp1d
                    tmp[i3, i2, a1+1, 1] += w * vals[(i1-1)*nqp1d^2+(i2-1)*nqp1d+i3]
                    w = w * r * (deg-a1) / (1.+a1)
                end
            end
        end
    end

    for i2=1:nqp1d
        xi = qpts[i2, 2]
        s = 1.0 - xi
        r = xi / s
        for a1=0:deg
            w = qwts[i2, 1] * s^(deg-a1)
            for a2=0:deg-a1
                for i3=1:nqp1d
                    tmp[i3,a2+1,a1+1,2] += w*tmp[i3,i2,a1+1,1]
                end
                w = w * r * (deg-a1-a2) / (1. + a2)
            end
        end
    end

    for i3=1:nqp1d
        xi = qpts[i3, 3]
        s = 1.0 - xi
        r = xi / s
        for a1=0:deg
            for a2=0:deg-a1
                w=qwts[i3, 3] * s^(deg-a1-a2)
                aind = polydim(3, deg) - polydim(3, deg-a1) + polydim(2, deg-a1) - polydim(2, deg-a1-a2) + 1
                for a3=0:deg-a1-a2
                    result[aind+a3] += w * tmp[i3, a2+1, a1+1, 2]
                    w = w * r * (deg-a1-a2-a3) / (1.+a3)
                end
            end
        end
    end

                    

    return
end

end
