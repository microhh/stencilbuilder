## Packages
using BenchmarkTools
using LoopVectorization
using Tullio


## Macros
function make_index(a, arrays, i, j, k)
    if a in arrays
        i += (a in [ Symbol("u"), Symbol("ut")]) ? 0.5 : 0
        j += (a in [ Symbol("v"), Symbol("vt")]) ? 0.5 : 0
        k += (a in [ Symbol("w"), Symbol("wt")]) ? 0.5 : 0

        if i < 0
            i_int = convert(Int, abs(i))
            ex_i = :( i-$i_int )
        elseif i > 0
            i_int = convert(Int, i)
            ex_i = :( i+$i_int )
        else
            ex_i = :( i )
        end

        if j < 0
            j_int = convert(Int, abs(j))
            ex_j = :( j-$j_int )
        elseif j > 0
            j_int = convert(Int, j)
            ex_j = :( j+$j_int )
        else
            ex_j = :( j )
        end

        if k < 0
            k_int = convert(Int, abs(k))
            ex_k = :( k-$k_int )
        elseif k > 0
            k_int = convert(Int, k)
            ex_k = :( k+$k_int )
        else
            ex_k = :( k )
        end

        return :( $a[ $ex_i, $ex_j, $ex_k] )
    else
        return :( $a )
    end
end

function process_expr(ex, arrays, i, j, k)
    n = 1

    if (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("gradx"))
        ex.args[1] = Symbol("gradx_")
        ex = :( $ex * dxi )
    elseif (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("grady"))
        ex.args[1] = Symbol("grady_")
        ex = :( $ex * dyi )
    elseif (isa(ex.args[1], Symbol) && ex.args[1] == Symbol("gradz"))
        ex.args[1] = Symbol("gradz_")
        ex = :( $ex * dzi )
    end

    args = ex.args
    while n <= length(args)
        if isa(args[n], Expr)
            args[n] = process_expr(args[n], arrays, i, j, k)
            n += 1
        elseif isa(args[n], Symbol)
            if args[n] == Symbol("gradx_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = process_expr(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = process_expr(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = process_expr(args[n+3], arrays, i+1.5, j, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = make_index(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = make_index(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = make_index(args[n+3], arrays, i+1.5, j, k)
                end
                args[n  ] = :( (  1//24) * $(args[n  ])  )
                args[n+1] = :( (-27//24) * $(args[n+1])  )
                args[n+2] = :( ( 27//24) * $(args[n+2])  )
                args[n+3] = :( ( -1//24) * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("grady_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = process_expr(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = process_expr(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = process_expr(args[n+3], arrays, i, j+1.5, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = make_index(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = make_index(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = make_index(args[n+3], arrays, i, j+1.5, k)
                end
                args[n  ] = :( (  1//24) * $(args[n  ])  )
                args[n+1] = :( (-27//24) * $(args[n+1])  )
                args[n+2] = :( ( 27//24) * $(args[n+2])  )
                args[n+3] = :( ( -1//24) * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("gradz_")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = process_expr(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = process_expr(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = process_expr(args[n+3], arrays, i, j, k+1.5)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = make_index(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = make_index(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = make_index(args[n+3], arrays, i, j, k+1.5)
                end
                args[n  ] = :( (  1//24) * $(args[n  ])  )
                args[n+1] = :( (-27//24) * $(args[n+1])  )
                args[n+2] = :( ( 27//24) * $(args[n+2])  )
                args[n+3] = :( ( -1//24) * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("interpx")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = process_expr(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = process_expr(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = process_expr(args[n+3], arrays, i+1.5, j, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i-1.5, j, k)
                    args[n+1] = make_index(args[n+1], arrays, i-0.5, j, k)
                    args[n+2] = make_index(args[n+2], arrays, i+0.5, j, k)
                    args[n+3] = make_index(args[n+3], arrays, i+1.5, j, k)
                end
                args[n  ] = :( (-1//16) * $(args[n  ])  )
                args[n+1] = :( ( 9//16) * $(args[n+1])  )
                args[n+2] = :( ( 9//16) * $(args[n+2])  )
                args[n+3] = :( (-1//16) * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("interpy")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = process_expr(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = process_expr(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = process_expr(args[n+3], arrays, i, j+1.5, k)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j-1.5, k)
                    args[n+1] = make_index(args[n+1], arrays, i, j-0.5, k)
                    args[n+2] = make_index(args[n+2], arrays, i, j+0.5, k)
                    args[n+3] = make_index(args[n+3], arrays, i, j+1.5, k)
                end
                args[n  ] = :( (-1//16) * $(args[n  ])  )
                args[n+1] = :( ( 9//16) * $(args[n+1])  )
                args[n+2] = :( ( 9//16) * $(args[n+2])  )
                args[n+3] = :( (-1//16) * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            elseif args[n] == Symbol("interpz")
                if isa(args[n+1], Expr)
                    args[n] = copy(args[n+1])
                    push!(args, copy(args[n+1]), copy(args[n+1]))
                    args[n  ] = process_expr(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = process_expr(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = process_expr(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = process_expr(args[n+3], arrays, i, j, k+1.5)
                elseif isa(args[n+1], Symbol)
                    args[n] = args[n+1]
                    push!(args, args[n+1], args[n+1])
                    args[n  ] = make_index(args[n  ], arrays, i, j, k-1.5)
                    args[n+1] = make_index(args[n+1], arrays, i, j, k-0.5)
                    args[n+2] = make_index(args[n+2], arrays, i, j, k+0.5)
                    args[n+3] = make_index(args[n+3], arrays, i, j, k+1.5)
                end
                args[n  ] = :( (-1//16) * $(args[n  ])  )
                args[n+1] = :( ( 9//16) * $(args[n+1])  )
                args[n+2] = :( ( 9//16) * $(args[n+2])  )
                args[n+3] = :( (-1//16) * $(args[n+3])  )
                insert!(args, n, Symbol("+"))
                n += 5
            else
                args[n] = make_index(args[n], arrays, i, j, k)
                n += 1
            end
        else
            n += 1
        end
    end
    return ex
end

macro fd(arrays, ex)
    i = (ex.args[1] in [ Symbol("u"), Symbol("ut")]) ? -0.5 : 0
    j = (ex.args[1] in [ Symbol("v"), Symbol("vt")]) ? -0.5 : 0
    k = (ex.args[1] in [ Symbol("w"), Symbol("wt")]) ? -0.5 : 0

    if isa(arrays, Symbol)
        arrays = :( [ $arrays ] )
    end

    ex = process_expr(ex, arrays.args, i, j, k)

    println("Generated stencil: ")
    println(ex)
    println("")

    return esc(ex)
end


## Advection, diffusion, time kernel.
function kernel!(
        ut, u, v, w,
        visc, dxi, dyi, dzi, dt,
        is, ie, js, je, ks, ke)

    @tturbo for k in ks:ke
        for j in js:je
            for i in is:ie
                @fd (ut, u, v, w) ut += (
                    - gradx(interpx(u) * interpx(u)) + visc * (gradx(gradx(u))) )
                @fd (ut, u, v, w) ut += (
                    - grady(interpx(v) * interpy(u)) + visc * (grady(grady(u))) )
                @fd (ut, u, v, w) ut += (
                    - gradz(interpx(w) * interpz(u)) + visc * (gradz(gradz(u))) )
            end
        end
    end

    @tturbo for k in ks:ke
        for j in js:je
            for i in is:ie
                @fd (ut, u) u += dt*ut
                @fd (ut, u) ut = 0
            end
        end
    end
end


## Initialize the grid.
itot = 512; jtot = 512; ktot = 512
igc = 4; jgc = 4; kgc = 4

dx = 1/itot; dy = 1/jtot; dz = 1/ktot
dxi = 1/dx; dyi = 1/dy; dzi = 1/dz
x = dx*collect(0:itot-1)
y = dy*collect(0:jtot-1)
z = dz*collect(0:ktot-1)


## Solve the problem in double precision.
visc = 1.5
dt = 1.e-3

u = zeros(Float64, (itot+2*igc, jtot+2*kgc, ktot+2*kgc))
v = zeros(Float64, (itot+2*igc, jtot+2*kgc, ktot+2*kgc))
w = zeros(Float64, (itot+2*igc, jtot+2*kgc, ktot+2*kgc))
ut = zeros(Float64, (itot+2*igc, jtot+2*kgc, ktot+2*kgc))


## Initialize with a sinus.
n_waves = 3
is = igc+1; ie = igc+ktot; js = jgc+1; je = jgc+jtot; ks = kgc+1; ke = kgc+ktot
uc = @view u[is:ie, js:je, ks:ke]
@tullio uc = sin(n_waves*2*pi*x[i]) + cos(n_waves*2*pi*y[j]) + sin(n_waves*2*pi*z[k])


## Run kernel.
is = igc+1; ie = igc+itot; js = jgc+1; je = jgc+jtot; ks = kgc+1; ke = kgc+ktot
@btime kernel!(
        $ut, $u, $v, $w,
        $visc, $dxi, $dyi, $dzi, $dt,
        $is, $ie, $js, $je, $ks, $ke)

