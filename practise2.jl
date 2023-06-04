# № 1 ----------------------------------
function fast_pow(a::T, n::Int) where T<:Any
    k=n
    p=a
    t=1
    #ИНВАРИАНТ: p^k*t=a^n
    while k>0
        if iseven(Integer(k))
            k /=2
            p *= p # - это преобразование следует из инварианта
        else
            k -= 1
            t *= p
        end
    end
    return t
end

# № 2 ----------------------------------
function fast_fib(a::Integer)
    matr = [1 1; 1 0]
    return (fast_pow(matr, a))[1,2]
end
# № 3 ----------------------------------
"""
z=x; t=1; y=0
#ИНВАРИАНТ:  x = z^t * a^y
while z < 1/a || z > a || t > ε 
    if z < 1/a
        z *= a # это перобразование направлено на достижения условия окончания цикла
        y -= t # тогда необходимрсть этого преобразования следует из инварианта цикла
    elseif z > a
        z /= a # это перобразование направлено на достижения условия окончания цикла
        y += t # тогда необходимрсть этого преобразования следует из инварианта цикла
    elseif t > ε
        t /= 2 # это перобразование направлено на достижения условия окончания цикла
        z *= z # тогда необходимрсть этого преобразования следует из инварианта цикла
    end
end
# y - искомое приближенное значение
"""
function log(a, x, e) # a > 1    z^t * a^y = x    
    z = x
    t = 1
    y = 0
    while z < 1/a || z > a || t > e  # написать инвариант
        if z < 1/a
            z *= a 
            y -= t 
        elseif z > a
            z /= a
            y += t
        elseif t > e
            t /= 2 
            z *= z 
        end
    end
    return y
end

# № 4 ----------------------------------
function bisection(f::Function, a, b, epsilon)
    @assert f(a)*f(b) < 0 
    @assert a < b

        f_a = f(a)
        #ИНВАРИАНТ: f_a*f(b) < 0
        while b-a > epsilon
            t = (a+b)/2
            f_t = f(t)
            if f_t == 0
                return t
            elseif f_a*f_t < 0
                b=t
            else
                a, f_a = t, f_t
            end
        end  
        return (a+b)/2

end

# № 5 ----------------------------------
println(bisection(xz -> cos(xz) - xz , 0, 1, 1e-8))

# № 6 ----------------------------------
function newton(r::Function, x; epsilon = 1e-8, num_max = 10)
    dx = r(x); x += dx; k=1
    while abs(dx) > epsilon && k < num_max
        dx = r(x); x += dx; k += 1
    end
    abs(dx) > epsilon && @warn("Требуемая точность не достигнута")
    return x
end

function newton(r::Function, x, epsilon, num_max = 10)
    dx = -r(x)
    k=0
    while abs(dx) > epsilon && k <= num_max
        x += dx
        dx = -r(x)
        k += 1
    end
    k > num_max && @warn("Требуемая точность не достигнута")
    return x
end

# № 7 ----------------------------------
f(x) = cos(x) - x

r(x) = -f(x)/(sin(x)+1)


# № 8 ----------------------------------
p(x) = 6*x^5 - 23*x^4 + 12*x^2 + 86

rp(x) = p(x) / (30*x^4 - 92*x^3 + 24*x)


ad = [86, 0, 12, 0, 23, 6]
function gorner(poly::AbstractVector{T}, x::T) where T<:Number
    n = length(poly)
    if n == 0
        return zero(T)
    end
    p = 0
    dp = 0
    for i in 1:n - 1
        dp = dp * x + p
        p = p * x + poly[i]
    end
    return p, dp
end

println(gorner(ad, 0))
println(newton(x->\(gorner(ad,x)...)))