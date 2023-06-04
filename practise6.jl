using Plots
using LinearAlgebra
 
# 1. Спроектировать типы Vector2D и Segment2D с соответсвующими функциями.
Vector2D{T <: Real} = NamedTuple{(:x, :y), Tuple{T,T}}
Base. +(a::Vector2D{T},b::Vector2D{T}) where T = Vector2D{T}(Tuple(a) .+ Tuple(b))
Base. -(a::Vector2D{T}, b::Vector2D{T}) where T = Vector2D{T}(Tuple(a) .- Tuple(b))
Base. *(α::T, a::Vector2D{T}) where T = Vector2D{T}(α.*Tuple(a))
LinearAlgebra.norm(a::Vector2D) = norm(Tuple(a))
LinearAlgebra.dot(a::Vector2D{T}, b::Vector2D{T}) where T = dot(Tuple(a), Tuple(b))
Base. cos(a::Vector2D{T}, b::Vector2D{T}) where T = dot(a,b)/(norm(a)*norm(b))
xdot(a::Vector2D{T}, b::Vector2D{T}) where T = a.x*b.y-a.y*b.x
Base.sin(a::Vector2D{T}, b::Vector2D{T}) where T = xdot(a,b)/(norm(a)*norm(b))
Base.angle(a::Vector2D{T}, b::Vector2D{T}) where T = atan(sin(a,b),cos(a,b))
Base.sign(a::Vector2D{T}, b::Vector2D{T}) where T = sign(sin(a,b))
Segment2D{T <: Real} = NamedTuple{(:A, :B), NTuple{2,Vector2D{T}}}
 
# 2. Написать функцию, проверяющую, лежат ли две заданные точки по одну сторону от заданной прямой (прямая задается некоторым содержащимся в ней отрезком).
function oneside(P::Vector2D{T}, Q::Vector2D{T}, s::Segment2D{T})::Bool where T
    # l - направляющий вектор прямой
    l = s.B - s.A
    # Тогда, точки , лежат по одну сторону от прямой <=> когда углы между векторами имеют один и тот же знак (отложены в одну и ту же сторону от прямой)
    return sin(l, P-s.A) * sin(l,Q-s.A) > 0
end

# 3. Написать функцию, проверяющую, лежат ли две заданные точки по одну сторону от заданной кривой (кривая задается уравнением вида F(x,y) = 0).

oneside(F::Function, P::Vector2D, Q::Vector2D)::Bool = ( F(P...) * F(Q...) > 0 )
 
# 4. Написать функцию, возвращающую точку пересечения (если она существует) двух заданных отрезков.

isinner(P::Vector2D, s::Segment2D) = (s.A.x <= P.x <= s.B.x || s.A.x >= P.x >= s.B.x)  &&  (s.A.y <= P.y <= s.B.y || s.A.y >= P.y >= s.B.y)

function intersection(s1::Segment2D{T},s2::Segment2D{T})::Union{Vector2D{T},Nothing} where T
    A = [s1.B[2]-s1.A[2] s1.A[1]-s1.B[1] 
        s2.B[2]-s2.A[2] s2.A[1]-s2.B[1]]

    #проверка матрицы A на вырожденность, если вырождена, то прямые паралельны, если нет, то всё хорошо
    if (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) == 0
        return nothing
    end

    b = [s1.A[2]*(s1.A[1]-s1.B[1]) + s1.A[1]*(s1.B[2]-s1.A[2])
        s2.A[2]*(s2.A[1]-s2.B[1]) + s2.A[1]*(s2.B[2]-s2.A[2])]
    
    x,y = A\b

    if isinner((;x, y), s1) == false || isinner((;x, y), s2) == false
        return nothing
    end

    return (;x, y) #Vector2D{T}((x,y))
end
 
println("Пересечение: ",intersection( (A=(x=-1.0,y=-1.0),B=(x=1.0,y=2.0)) , (A=(x=1.0,y=-1.0),B=(x=-1.0,y=3.0)) ))

# 5. Написать функцию, проверяющую лежит ли заданная точка внутри заданного многоугольника.
function isinside(point::Vector2D{T},polygon::AbstractArray{Vector2D{T}})::Bool where T
    @assert length(polygon) > 2
    sum = zero(Float64)
 
    # Вычислить алгебраическую (т.е. с учетом знака) сумму углов, между направлениями из заданной точки на каждые две сосоедние вершины многоугольника.
    # Далее воспользоваться тем, что, если полученная сумма окажется равнной нулю, то точка лежит вне многоугольника, а если она окажется равной 360 градусам, то - внутри.
    for i in firstindex(polygon):lastindex(polygon)
        sum += angle(polygon[i] - point, polygon[i % lastindex(polygon) + 1] - point)
    end

    if abs(sum) == 2 * pi
        return true
    end
    return false
end

println("Внутри: ",isinside( (x=0,y=0),[(x=0,y=1),(x=1,y=-1),(x=-1,y=-1)] ))
println("Внутри: ",isinside( (x=5,y=0),[(x=0,y=1),(x=1,y=-1),(x=-1,y=-1)] ))

# 6. Написать функцию, проверяющую, является ли заданный многоугольник выпуклым.
function isconvex(polygon::AbstractArray{Vector2D{T}})::Bool where T
    @assert length(polygon) > 2

    for i in firstindex(polygon):lastindex(polygon)
        # У выпуклого многоугольника все внутренние углы будут меньше 180 градусов.
        # А у не выпуклого многоугольника обязательно найдутся, как углы меньшие, так и большие 180 градусов
        k = angle( polygon[i > firstindex(polygon) ? i - 1 : lastindex(polygon)] - polygon[i] , polygon[i % lastindex(polygon) + 1] - polygon[i])
        if (k >= 0 && k >= pi) || (k < 0 && 2 * pi + k >= pi)
            return false
        end
    end
    return true
end
 
println("Выпуклый: ",isconvex([(x=0,y=1),(x=1,y=-1),(x=-1,y=-1)]))
println("Выпуклый: ",isconvex([(x=0,y=0),(x=2,y=2),(x=5,y=-2),(x=3,y=5)]))

# 7. Выпуклая оболочка по Джарвису
function orientation(p1::Vector2D{T}, p2::Vector2D{T}, p3::Vector2D{T}) where T
    d = (p3.y - p2.y) * (p2.x - p1.x) - (p2.y-p1.y) * (p3.x - p2.x)
    if d > 0
        return 1
    elseif d < 0
        return -1
    else
        return 0
    end
end
function jarvis(points::AbstractArray{Vector2D{T}})::AbstractArray{Vector2D{T}} where T
    on_hull = minimum(points)
    hull = []
    while true
        push!(hull, on_hull)
        next_point = points[1]
        for point in points
            o = orientation(on_hull, next_point, point)
            if next_point == on_hull || o == 1 || (o == 0 && norm(on_hull - point) > norm(on_hull - next_point))
                next_point = point
            end
        end
        on_hull = next_point
        if on_hull == hull[1]
            break
        end
    end
    return hull
end

test_hull = [(x=0.0,y=0.0),(x=5.0,y=1.0),(x=4.0,y=3.0),(x=2.0,y=2.0),(x=3.0,y=-2.0),(x=2.0,y=-2.0),(x=1.0,y=9.0),(x=-3.0,y=8.0),(x=-5.0,y=2.0),(x=-5.0,y=1.0),(x=-2.0,y=3.0)]
println("Алгоритм Джарвиса: ",jarvis(test_hull))
 
# 8. Написать функцию, реализующую алгоритм Грехома построения выпуклой оболочки заданных точек плоскости.
function graham_scan(points::AbstractArray{Vector2D{T}})::AbstractArray{Vector2D{T}} where T
    p0 = points[1] #нижняя левая точка
    for i in firstindex(points):lastindex(points) 
        if (points[i].y < p0.y) || (points[i].y == p0.y && points[i].x < p0.x) 
            p0 = points[i] 
        end 
    end

    sort!(@view(points[1:end]), by=(point -> angle((x=oneunit(T), y=zero(T)), point - p0))) #сортируем в зависимости от угла от нормали к p0

    hull = []
    for i in firstindex(points):lastindex(points) 
        while length(hull) >= 2 && orientation(hull[lastindex(hull) - 1], hull[lastindex(hull)], points[i]) != 1
            pop!(hull) #удаляем последний элемент из оболочки если ориентация не против часовой стрелки
        end
        push!(hull, points[i])
    end
    return hull
end
println("Алгоритм Грехома: ",graham_scan(test_hull))
# 9. Написать функцию вычисляющую площадь (ориентированную) заданного многоугольника методом трапеций.
function area_trapeze(poly::AbstractArray{Vector2D{T}})::T where T
    res = zero(T)
    # area = (yk + yk+1)(xk+1 − xk)/2
    for i in 1:length(poly)
        j = mod1(i + 1, length(poly))
        res += (poly[i].y + poly[j].y) * (poly[j].x - poly[i].x) / 2
    end
    return res
end
test_area = [(x=2.0,y=-1.0),(x=1.0,y=2.0),(x=-1.0,y=3.0),(x=-3.0,y=-1.0),]
println("Площадь (Трапеция): ",area_trapeze(test_area))
# 10. Написать функцию вычисляющую площадь (ориентированную) заданного многоугольника методом треугольников.
function area_triangle(poly::AbstractArray{Vector2D{T}})::T where T
    res = zero(T)
 
    # area = (yk + yk+1)(xk+1 − xk)/2
    for i in firstindex(poly)+1:lastindex(poly)-1
        res += xdot(poly[1] - poly[i], poly[i+1] - poly[i])/2
    end
 
    return res
end
println("Площадь (Треугольники): ",area_triangle(test_area))






stored_lims = [0,0,0,0]
 
function lims!(x1,y1,x2,y2)
    stored_lims[1] = min(x1-1,stored_lims[1])
    stored_lims[2] = min(y1-1,stored_lims[2])
    stored_lims[3] = max(x2+1,stored_lims[3])
    stored_lims[4] = max(y2+1,stored_lims[4])
 
    xlims!(stored_lims[1], stored_lims[3])
    ylims!(stored_lims[2], stored_lims[4])
end
 
lims!(x,y) = lims!(x,y,x,y)
 
function draw(vertices::AbstractArray{Vector2D{T}}) where T
    vertices = copy(vertices)
    push!(vertices,first(vertices))
 
    x = [v.x for v in vertices]
    y = [v.y for v in vertices]
 
    plot(x, y, color=:blue, legend=false)
 
    lims!( minimum(x) , minimum(y) , maximum(x) , maximum(y) )
end
 
function draw(point::Vector2D{T}) where T
    scatter!([point.x,point.x], [point.y,point.y], color=:red, markersize=5, legend=false)
 
    lims!( point.x , point.y )
end

function draw_a(points::AbstractArray{Vector2D{T}}) where T
    for point in points
        scatter!([point.x,point.x], [point.y,point.y], color=:red, markersize=5, legend=false)
        lims!( point.x , point.y )
    end
end
function clear()
    fill!(stored_lims,0)
    xlims!(0,1)
    ylims!(0,1)
    plot!()
end
 

draw(graham_scan(test_hull))
draw_a(test_hull)
savefig("D:\\grekhom.png")
clear()

draw(jarvis(test_hull))
draw_a(test_hull)
savefig("D:\\jarvis.png")
clear()
