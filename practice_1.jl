# № 1 ----------------------------------
function gcd_(a::T, b::T) where T<:Integer
    while !iszero(b)
        a, b = b, a % b
    end
    return abs(a)
end

# № 2 ----------------------------------
function gcdx_(a::T, b::T) where T<:Integer 
    u, v = one(T), zero(T) 
    u1, v1 = zero(T), one(T)
    while !iszero(b)
        k = div(a, b)
        a, b = b, a - k * b
        u, v = u1, v1
        u1, v1 = u - k * u1, v - k * v1
    end
    if isnegarive(a) 
        a, u, v = -a, -u, -v
    end
    return a, u, v 
end

isnegative(a::Integer) = (a < 0)

# № 4 ----------------------------------
function diophant(x::T, y::T, z::T) where T<:Integer
    a, b, c = gcdx_(x,y)
    return (b * z)÷a , (c * z)÷a
end

# № 7 ----------------------------------
struct Polynom{T}
    a::AbstractVector{T}
    Polynom{T}(a::AbstractVector{T}) where {T<:Integer} = new(a)
end

function Base. +(a::Polynom{T}, b::Polynom{T}) where {T<:Integer} 
    m = Polynom{T}(zeros(T, max(length(a.a), length(b.a))))
    for i in range(1, max(length(a.a), length(b.a)))
        m.a[i] = 0
        if i <= length(a.a)
            m.a[i] += a.a[i]
        end
        if i <= length(b.a)
            m.a[i] += b.a[i]
        end
    end
    return m
end
Base. -(a::Polynom{T}) where {T<:Integer} = Polynom{T}(-a.a)
Base. -(a::Polynom{T}, b::Polynom{T}) where {T<:Integer} = a + (-b)

function Base. *(a::Polynom{T}, b::Polynom{T}) where {T<:Integer}
    m = Polynom{T}(zeros(T, (length(a.a)-1) + (length(b.a)-1) + 1))
    for i in range(1, length(a.a))
        for j in range(1, length(b.a))
            m.a[(i-1)+(j-1)+1] += a.a[i] * b.a[j]
            println((i-1)+(j-1)+1)
        end
    end
    return m
end

function divmod(a::Polynom{T}, b::Polynom{T}) where {T<:Integer} # Made using Bezout's sheme
    # assuming length(a) >= length(b)
    a_a = copy(a.a)
    b_a = zeros(T, length(a.a) - length(b.a))
    m = Polynom{T}(zeros(T, length(a.a) - length(b.a) + 1))
    append!(b_a, b.a)
    for i in length(a.a):-1:2
        n = div(a_a[i], b_a[i])
        a_a -= (b_a * n)
        popfirst!(b_a)
        push!(b_a, 0)
        if i-length(b.a)+1 != 0
            m.a[i-length(b.a)+1] = n
        end
    end
    while a_a[length(a_a)] == 0
        pop!(a_a)
        if length(a_a) == 0
            a_a = [0]
            break
        end
    end
    return (m, Polynom{T}(a_a))
end

Base. div(a::Polynom{T}, b::Polynom{T}) where {T<:Integer} = divmod(a, b)[1]
Base. mod(a::Polynom{T}, b::Polynom{T}) where {T<:Integer} = divmod(a, b)[2]
Base. rem(a::Polynom{T}, b::Polynom{T}) where {T<:Integer} = divmod(a, b)[2]

function Base. display(a::Polynom{T}) where {T<:Integer}
    println(string(a.a))
end
function Base. print(a::Polynom{T}) where {T<:Integer}
    print(string(a.a))
end
function Base. println(a::Polynom{T}) where {T<:Integer}
    println(string(a.a))
end

# № 5,6 ----------------------------------
struct Z{T,N}
    a::T
    Z{T,N}(a::T) where {T<:Integer, N} = new(mod(a, N))
    Z{T,N}(a::T) where {T<:Polynom{Int}, N} = new(mod(a, Polynom{Int}(Vector(collect(N)))))
end


# № 3 ----------------------------------
function Base. *(a::Z{T,N}) where {T<:Integer, N}
    if gcd(a.a, N) != 1 
        return -1
    else
        f, s, d = gcdx_(a.a, N)
        return Z{T,N}(s)
    end 
end


Base. +(a::Z{T,N}, b::Z{T,N}) where {T<:Any, N} = Z{T,N}(a.a + b.a)
Base. -(a::Z{T,N}, b::Z{T,N}) where {T<:Any, N} = Z{T,N}(a.a - b.a)
Base. *(a::Z{T,N}, b::Z{T,N}) where {T<:Any, N} = Z{T,N}(a.a * b.a)
Base. -(a::Z{T,N}) where {T<:Any, N} = Z{T,N}(-a.a)

F = Z{Int, 5}(7)
Q = Z{Int, 5}(9)

function Base. display(a::Z{T,N}) where {T<:Any, N}
    println(string(a.a))
end
function Base. print(a::Z{T,N}) where {T<:Any, N}
    print(string(a.a))
end
function Base. println(a::Z{T,N}) where {T<:Any, N}
    println(string(a.a))
end

# Кортеж коэффицентов

K = Z{Polynom{Int}, (1, 1)}(Polynom{Int}([1, 1, 1]))