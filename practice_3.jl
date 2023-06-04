# № 1 ----------------------------------
function isprime(n::IntType) where IntType <: Integer # является ли заданное число простым
    for d in 2:IntType(ceil(sqrt(n)))
        if n % d == 0
            return false
        end
    end
    return true
end

# № 2 ----------------------------------
function eratosphenes_sieve(n::Integer)
    prime_indexes = ones(Bool, n)
    prime_indexes[begin] = false
    i = 2
    prime_indexes[i^2:i:n] .= false # - четные вычеркнуты
    i=3
    #ИНВАРИАНТ: i - простое нечетное
    while i <= n
        prime_indexes[i^2:2i:n] .= false
        # т.к. i^2 - нечетное, то шаг тут можно взять равным 2i, т.к. нечетное+нечетное=четное, а все четные уже вычеркнуты
        i+=1
        while i <= n && prime_indexes[i] == false
            i+=1
        end
        # i - очередное простое (первое не вычеркунутое)
    end
    return findall(prime_indexes)
end

# № 3 ----------------------------------
function factorize(n::IntType) where IntType <: Integer
    list = NamedTuple{(:div, :deg), Tuple{IntType, IntType}}[]
    for p in eratosphenes_sieve(Int(ceil(n/2)))
        k = degree(n, p) # кратность делителя
        if k > 0
            push!(list, (div=p, deg=k))
        end
    end
    return list
end

function degree(n, p) # кратность делителя `p` числа `n`
    k=0
    n, r = divrem(n,p)
    while n > 0 && r == 0
        k += 1
        n, r = divrem(n,p)
    end
    return k
end

# № 4 ----------------------------------
function meanstd(V :: Vector)
    m, s = 0, 0
    for i in V
        m += i
        s += i*i
    end
    return s / length(V) - (m/length(V))^2
end

# № 5 ----------------------------------
#5 Взаимные преобразования различных способов представления деревьев
struct Node
    index :: Int
    children :: Vector{Union{Nothing,Node}}
end
function convert!( arr :: Vector, tree :: Dict{Int,Vector})
    isempty(arr) && return

    list = []

    for subarr in arr[1:end-1]
        if isempty(subarr)
            push!(list,nothing)
            continue
        end
        if typeof(subarr) <: Int
            push!(list,subarr)
            continue
        end
        push!(list,subarr[end])
        convert!(subarr,tree)
    end

    tree[arr[end]] = list

    return tree
end

function convert!(tree :: Dict{Int,Vector}; root ::  Union{Int,Nothing}) :: Union{Vector,Int}
    arr = []
    isnothing(root) && return []
    !(root in keys(tree)) && return root
    for subroot in tree[root]
        push!(arr,convert!(tree; root = subroot))
    end
    push!(arr,root)
    return arr
end

function convert!( tree :: Dict{Int,Vector}, root :: Union{Int,Nothing}) ::Union{Node,Nothing}

    isnothing(root) && return nothing
    !(root in keys(tree)) && return Node(root,[])
    node = Node(root,[])

    for sub_root in tree[root]
        push!(node.children, convert!(tree, sub_root))
    end

    return node
end

function convert!( node :: Node) :: Union{Vector,Int}
    arr = []
    length(node.children)==0 && return node.index
    for child in node.children
        if isnothing(child)
            push!(arr, [])
            continue
        end
        push!(arr,convert!(child))
    end
    push!(arr,node.index)
    return arr

end
function convert!(node :: Node, tree :: Dict{Int, Vector}) :: Union{Dict{Int,Vector},Int}
    list = []
    for child in node.children
        if isnothing(child)
            push!(list, nothing)
            continue
        end
        push!(list,child.index)
        length(child.children) != 0 && convert!(child,tree)
    end
    tree[node.index] = list
    return tree
end

#----------------------------------------------------------------------------------------------------------#
arr = [[[[],[],6], [], 2], [[10,11,4], [[],[],5], 3],1]
tree = Dict{Int,Vector}();
tree = convert!(arr, tree)

display(tree)

_arr = convert!(tree; root = 1)
println(_arr)

node = convert!(tree, 1)
println(node)

_arr = convert!(node)
println(_arr)

tree = convert!(node,tree)
display(tree)

# № 6 Рекурсивные алгоритмы дерева

# Высота дерева
function depth( arr :: Union{Vector,Int}, i :: Int = 1)
    length(arr)==0 && return i-1
    typeof(arr) <: Int && return i
    return max([depth(subarr,i+1) for subarr in arr[1:end-1]]...)
end

function leafs( arr :: Union{Vector,Int})
    length(arr)==0  && return 0
    typeof(arr) <: Int && return 1
    all(length(subarr)==0 for subarr in arr[1:end-1]) && return 1 
    s=0
    for subarr in arr[1:end-1]
        s+=leafs(subarr)
    end
    return s
end

function nodes( arr :: Union{Vector,Int})
    typeof(arr) <: Int && return 1
    s = 0
    for subarr in arr
        s += nodes(subarr)
    end
    return s
end

function treepower( arr :: Union{Vector,Int}, root :: Bool = true)
    length(arr)==0  && return 1
    typeof(arr) <: Int && return 1
    return max(1,Int(!root)+length(arr[1:end-1]),[treepower(subarr, false) for subarr in arr[1:end-1]]...)
end

function meantrail(_arr :: Union{Vector,Int})
    function recursive(arr; depth :: Int = 0)
        isempty(arr) && return 0,0
        typeof(arr) <: Int && return depth,1
        s,n = depth, 1
        for subarr in arr[1:end-1]
            _s, _n = recursive(subarr; depth = depth+1)
            s+=_s
            n+=_n
        end
        return s,n
    end
    paths,verts = recursive(_arr)
    return paths/verts
end