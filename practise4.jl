using LinearAlgebra
using Plots
 
#1. Написать функцию, вычисляющую n-ю частичную сумму ряда Телора
# (Маклорена) функции для произвольно заданного значения аргумента x.
# Сложность алгоритма должна иметь оценку O(n).
function exp_taylor(x::T,n::Integer=10)::T where T
    value = oneunit(T)
    fact = oneunit(T)  
    x_base = x
    f = oneunit(T) 
 
    for i in 1:n
        value += x / fact
        x *= x_base
        fact *= f + oneunit(T)
        f += 1
    end
 
    return value
end

#2. Написать функцию, вычиляющую значение exp(x) с машинной точностью (с
# максимально возможной в арифметике с плавающей точкой).
function exp_ideal(x::T)::T where T
    return (Float64(ℯ) ^ Int(trunc(x))) * exp_taylor(x - trunc(x))
end

println(exp_taylor(5.3))
println(exp(5.3))
println(exp_ideal(5.3))

#3. Написать функцию, вычисляющую функцию Бесселя заданного целого
# неотрицательного порядка по ее ряду Тейлора с машинной точностью. Для
# этого сначала вывести соответствующую рекуррентную формулу,
# обеспечивающую возможность эффективного вычисления.
function besselj(M::Integer, x::Real)
    sqrx = x*x
    a = 1/factorial(M)
    m = 1
    s = 0 
    
    while s + a != s
        s += a
        a = -a * sqrx /(m*(M+m)*4)
        m += 1
    end
    
    return s*(x/2)^M
end
 
#3. Построить семейство графиков этих функций для нескольких порядков, начиная с нулевого порядка.
values = 0:0.1:20
myPlot = plot()
for m in 1:5
    plot!(myPlot, values, besselj.(m, values))
end
display(myPlot)
 
#4. Реализовать алгорим, реализующий обратный ход алгоритма Жордана-Гаусса

function reverse_gauss(A::AbstractMatrix{T}, b::AbstractVector{T}) where T
    x = similar(b) #вектор решений
    N = size(A, 1) #количчество строк
    for k in 0:N-1
        x[N-k] = (b[N-k] - sum(A[N-k,N-k+1:end] .* x[N-k+1:end])) / A[N-k,N-k] #классическое вычисление x_k
    end
    return x
end

#5. Реализовать алгоритм, осуществляющий приведение матрицы матрицы к ступенчатому виду
function transform_to_steps!(A::AbstractMatrix; epsilon::Real = 1e-7)
    @inbounds for k ∈ 1:size(A, 1)
        absval, Δk = findmax(abs, @view(A[k:end,k]))
        (absval <= epsilon) && throw("Вырожденая матрица")
        Δk > 1 && swap!(@view(A[k,k:end]), @view(A[k+Δk-1,k:end]))
        for i ∈ k+1:size(A,1)
            t = A[i,k]/A[k,k] 
            @. @views A[i,k:end] = A[i,k:end] - t * A[k,k:end]
        end
    end
    return A
end

#6. Реализовать алгоритм, реализующий метод Жордана-Гаусса решение СЛАУ для произвольной невырожденной матрицы (достаточно хорошо обусловленной).

function solve_sla(A::AbstractMatrix{T}, b::AbstractVector{T}) where T 
    Ab = [A b]
    transform_to_steps!(Ab; epsilon = 10*sqrt(eps(T))*maximum(abs,A))
    return reverse_gauss(Ab)
end
   
 
#7. Постараться обеспечить максимально возможную производительность алгорима решения СЛАУ;
# провести временные замеры с помощью макроса @time для систем большого размера (порядка 1000)

#более эфективынй гаусс
function Res_Gauss(A::AbstractMatrix{T}, b::AbstractVector{T}) where T
    x = similar(b)
    N = size(A, 1)
    for k in 0:N-1
        x[N-k] = (b[N-k] - sumprod(@view(A[N-k,N-k+1:end]), @view(x[N-k+1:end]))) / A[N-k,N-k]
    end
    return x
end

@inline function sumprod(A::AbstractVector{T}, B::AbstractVector{T}) where T
    s = T(0)
    @inbounds for i in eachindex(A)
        s = fma(A[i], B[i], s)
    end
    return s
end

for n in 1000:500:5000
    println("Матрица порядка ",n,"×",n,":")
    @time Res_Gauss(randn(n,n),randn(n))
    @time reverse_gauss(randn(n,n),randn(n))
    println("--------------------")
end

#8. Написать функцию, возвращающую ранг произвольной прямоугольной матрицы (реализуется на базе приведения матрицы к ступенчатому виду).
function GetRank(matrix::AbstractMatrix{T}) where T
    transform_to_steps!(Matrix)
    i = 1
    while matrix[i,i] != zero(T)
        i+=1
    end
    return i-1
end
 
#9. Написать функцию, возвращающую определитель произвольной квадратной матрицы (реализуется на основе приведения матрицы к ступенчатому виду).
function CalcDeterminant(matrix::AbstractMatrix{T}) where T
    transform_to_steps!(matrix)
    det = oneunit(T)
    i = 1
    while i <= size(matrix, 1)
        if matrix[i, i] == zero(T)
            det = 0
            break
        end
        det *= matrix[i, i]
        i += 1
    end
    return det
end
