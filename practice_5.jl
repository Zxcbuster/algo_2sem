using BenchmarkTools
# № 1 ----------------------------------

function sort_perm!(a)
    indexes = collect(firstindex(a):lastindex(a))
    n = length(a)
    for i in 1:n
        for j in i:n
            if a[i]>a[j]
                a[i], a[j] = a[j], a[i]
                indexes[i], indexes[j] = indexes[j], indexes[i]
            end
        end
    end
    return indexes
end
 
sort_perm(a)=sort_perm!(copy(a))

# № 2 ----------------------------------
function bubble_sort!(a)
    n = length(a)
    for k in 1:n-1
        istranspose = false
        for i in 1:n-k
            if a[i]>a[i+1]
                a[i], a[i+1] = a[i+1], a[i]
                istranspose = true
            end
        end
        if istranspose == false
            break
        end
    end
    return a
end

function comb_sort!(a; factor=1.2473309)
    step = length(a)
    while step >= 1
        for i in 1:length(a)-step
            if a[i] > a[i+step]
                a[i], a[i+step] = a[i+step], a[i]
            end
        end
        step = Int(floor(step/factor))
    end
    bubble_sort!(a)
end

f = randn(1000)
@time bubble_sort!(f)
f = randn(1000)
@time comb_sort!(f)
#как мы видим, уже на 1000 элементах сортировка расчёсыванием перегоняет пузырьковую.
#сложность пузырьковой O(n^2) во всех случаях, а у сортировки расчёсыванием в среднем же O(n^2/2^p)


# № 3 ----------------------------------

function insert_sort!(vector)
    n = 1
    while n < length(vector) 
        n += 1
        i = n
        while i > 1 && vector[i-1] > vector[i]
            vector[i], vector[i-1] = vector[i-1], vector[i]
            i -= 1
        end
    end
    return vector
end

function shell_sort!(a;  step_series = (length(a)÷2^i for i in 1:Int(floor(log2(length(a))))) )
    for step in step_series
        for i in firstindex(a):lastindex(a)-step
            j = i
            while j >= firstindex(a) && a[j] > a[j+step]
                a[j], a[j+step] = a[j+step], a[j]
                j -= step
            end
        end
    end
    return a
end

function shell_sort2(arr)
    n = length(arr)
    gap = div(n, 2)
    
    while gap > 0
        for i = gap+1:n
            temp = arr[i]
            j = i
            
            while j > gap && arr[j-gap] > temp
                arr[j] = arr[j-gap]
                j -= gap
            end
            
            arr[j] = temp
        end
        
        gap = div(gap, 2)
    end
    
    return arr
end

#=
f = rand(Int, 10000)
@time shell_sort!(f)
f = rand(Int, 10000)
@time insert_sort!(f)
=#
#на значении в 10000 элементов сортировка шела становится быстрее сортировки вставками. 
#сортировка шела немного лучше чем обычная сортировка вставками, работающая во всех случаях за O(n^2), 
#так как это вероятностный алгоритм, в лучшем случае, который может выдавать даже O(n * (logn)^2)

# № 4 ----------------------------------

@inline function Base.merge!(a1, a2, a3)::Nothing 
    i1, i2, i3 = 1, 1, 1
    @inbounds while i1 <= length(a1) && i2 <= length(a2) 
        if a1[i1] < a2[i2]
            a3[i3] = a1[i1]
            i1 += 1
        else
            a3[i3] = a2[i2]
            i2 += 1
        end
        i3 += 1
    end
    @inbounds if i1 > length(a1)
        a3[i3:end] .= @view(a2[i2:end]) 
    else
        a3[i3:end] .= @view(a1[i1:end])
    end
    nothing
end

function merge_sort!(a)
    b = similar(a) 
    N = length(a)
    n = 1 
    @inbounds while n < N
        K = div(N,2n) 
        for k in 0:K-1
            merge!(@view(a[(1:n).+k*2n]), @view(a[(n+1:2n).+k*2n]), @view(b[(1:2n).+k*2n]))
        end
        if N - K*2n > n
            merge!(@view(a[(1:n).+K*2n]), @view(a[K*2n+n+1:end]), @view(b[K*2n+1:end]))
        elseif 0 < N - K*2n <= n 
            b[K*2n+1:end] .= @view(a[K*2n+1:end])
        end
        a, b = b, a
        n *= 2
    end
    if isodd(log2(n)) 
        b .= a 
        a = b
    end
    return a
end
#=
f = rand(Int, 100000)
@time bubble_sort!(f)
f = rand(Int, 100000)
@time comb_sort!(f)
f = rand(Int, 100000)
@time shell_sort!(f)
f = rand(Int, 100000)
@time insert_sort!(f)
f = rand(Int, 100000)
@time merge_sort!(f)
=#
#видно, что на маленьких значениях сортировка слиянием показывает себя хуже, но на больших значениях она выигрывает, так как её сложность O(n*logn)

# № 5 ----------------------------------
function part_sort!(A, b)
    N = length(A)
    K=0
    L=0
    M=N
    while L < M 
        if A[L+1] == b
            L += 1
        elseif A[L+1] > b
            A[L+1], A[M] = A[M], A[L+1]
            M -= 1
        else
            L += 1; K += 1
            A[L], A[K] = A[K], A[L]
        end
    end
    return K, M+1 
end

function quick_sort!(A)
    if isempty(A)
        return A
    end
    N = length(A)
    K, M = part_sort!(A, A[rand(1:N)])
    quick_sort!(@view A[1:K])
    quick_sort!(@view A[M:N])
    return A
end

#=
f = rand(Int, 10000000)
@time comb_sort!(f)
f = rand(Int, 10000000)
@time shell_sort!(f)
f = rand(Int, 10000000)
@time merge_sort!(f)
f = rand(Int, 10000000)
@time quick_sort!(f)
=#
#в данной ситуации получилось, что сортировка хоара не настолько эффективна, как другие влоть до 1+e^6, но её эффективность раскрывается на больших значениях
#так как её сложность в среднем равна O(N*log(N)) хоть и с достаточно большим коэфициентом 

# № 6 ----------------------------------
function findMedian(a, n)
    if (n % 2 == 0)
        part_sort!(a, n÷2)
        part_sort!(a, (n - 1) ÷ 2)
        return (a[(n - 1) ÷ 2] + a[n ÷ 2]) ÷ 2
    else 
        part_sort!(a, n÷2)
        return a[n ÷ 2]
    end
end
f = rand(Int, 1000)
@time findMedian(f, length(f))


# № 7 ----------------------------------
function counting_sort!(arr::Vector{T}) where T <: Integer
    max_val = maximum(arr)
    count_arr = zeros(T, max_val + 1)
    
    for val in arr
        count_arr[val + 1] += 1
    end
    for i in 2:length(count_arr)
        count_arr[i] += count_arr[i-1]
    end
    sorted_arr = zeros(T, length(arr))
    for val in arr
        sorted_arr[count_arr[val + 1]] = val
        count_arr[val + 1] -= 1
    end
    
    return sorted_arr
end

x = rand(1:1000, 1000) 
@time bubble_sort!(copy(x))
@time comb_sort!(copy(x))
@time insert_sort!(copy(x))
@time shell_sort!(copy(x))
@time merge_sort!(copy(x))
@time quick_sort!(copy(x))
@time counting_sort!(copy(x))
println()
f = rand(1:1000, 1000) 
@time bubble_sort!(copy(f))
@time comb_sort!(copy(f))
@time insert_sort!(copy(f))
@time shell_sort!(copy(f))
@time merge_sort!(copy(f))
@time quick_sort!(copy(f))
@time counting_sort!(copy(f)) 
println()
m = rand(1:10000, 1000) 
@time bubble_sort!(copy(m))
@time comb_sort!(copy(m))
@time insert_sort!(copy(m))
@time shell_sort!(copy(m))
@time merge_sort!(copy(m))
@time quick_sort!(copy(m))
@time counting_sort!(copy(m)) 
println()
d = rand(1:10000, 100000) 
@time bubble_sort!(copy(d))
@time comb_sort!(copy(d))
@time insert_sort!(copy(d))
@time shell_sort!(copy(d))
@time merge_sort!(copy(d))
@time quick_sort!(copy(d))
@time counting_sort!(copy(d)) 