module OpenChannel

using Roots

abstract Section

immutable Trapezoidal <: Section
    b::Real
    kl::Real
    kr::Real
end

immutable Triangular <: Section
    kl::Real
    kr::Real
end

function Fr2(section::Section, Q, g, y)
    Q2 = Q*Q
    a3 = area(section, y) * area(section, y) * area(section, y)
    return Q2*tw(section, y) / (g * a3)
end

function tw(section::Trapezoidal, y::Real)
    return section.b + y*(section.kl + section.kr)
end

function area(section::Trapezoidal, y::Real)
    return y * ( section.b + (section.kl + section.kr)*y/2 )
end

function yc(section::Triangular, Q, g)
    Q2 = Q*Q
    k2 = (section.kl + section.kr)*(section.kl + section.kr)
    return ( 8*Q2 / (k2*g) )^(1.0/5.0)
end

function yc(section::Trapezoidal, Q, g)
    tri = Triangular(section.kl, section.kr)
    guess = yc(tri, Q, g)
    #Review the bracket
    return fzero(x -> Fr2(section, Q, g, x) - 1, 0, guess)
end

end
