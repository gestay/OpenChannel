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

function Fr2(section::Section, Q::Real, g::Real, y::Real)
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

function wp(section::Trapezoidal, y::Real)
    temp_a = sqrt( 1 + section.kl*section.kl )
    temp_b = sqrt( 1 + section.kr*section.kr )
    return section.b + y*( temp_a + temp_b )
end

#Critical depth
function yc(section::Triangular, Q::Real, g::Real)
    Q2 = Q*Q
    k2 = (section.kl + section.kr)*(section.kl + section.kr)
    return ( 8*Q2 / (k2*g) )^(1/5)
end

function yc(section::Trapezoidal, Q::Real, g::Real)
    tri = Triangular(section.kl, section.kr)
    guess = yc(tri, Q, g)
    #Review the bracket
    return fzero(x -> Fr2(section, Q, g, x) - 1, 0, guess)
end

#Normal depth
function yn(section::Triangular, Q::Real, s::Real, n::Real)
    temp_a = ( Q*n/sqrt(s) )^(3/8)
    temp_b = ( sqrt(1 + section.kl*section.kl)
               + sqrt(1 + section.kr*section.kr) )^(1/4)
    temp_c = ( (section.kl + section.kr)/2 )^(5/8)

    return temp_a*temp_b/temp_c
end

end
