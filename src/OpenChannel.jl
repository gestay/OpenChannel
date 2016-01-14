module OpenChannel

using Roots

abstract Section

immutable Trapezoidal <: Section
    b::Real
    kl::Real
    kr::Real
end

function tw(section::Trapezoidal, y::Real)
    return section.b + y*(section.kl + section.kr)
end

function area(section::Trapezoidal, y::Real)
    return y * ( section.b + (section.kl + section.kr)*y/2 )
end

end
