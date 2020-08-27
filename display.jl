using Plots

include("evo_analysis_functions.jl")

mutable struct Node
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
    fx::Float64
    fy::Float64
    activated::Bool
end

function force_sim(GRN)
    N = length(GRN[1,:])
    arr = Array{Any, 1}()
    print(N)
    print(arr)
    for i in 0:N-1
        x = sqrt(i)*cos(i)/sqrt(N)
        y = sqrt(i)*sin(i)/sqrt(N)
        push!(arr, Node(x,y,0,0,0,0,false))
    end
    for step in 1:5000
        for i in 1:N
            arr[i].fx = 0
            arr[i].fy = 0
        end
        for i in 1:N
            for j in i+1:N
                dx = arr[j].x - arr[i].x
                dy = arr[j].y - arr[i].y

                nx = dx/pyth(dx,dy)
                ny = dy/pyth(dx,dy)

                repulsion = 10/((pyth(dx,dy)*sqrt(N))^3)

                dfx = repulsion * nx
                dfy = repulsion * ny

                arr[i].fx -= dfx
                arr[j].fx += dfx
                arr[i].fy -= dfy
                arr[j].fy += dfy
                if GRN[i,j] != 0 || GRN[j,i] != 0
                    # println(pyth(dx,dy))
                    attraction = 1.5*sqrt(pyth(dx,dy))

                    dfx = attraction * nx
                    dfy = attraction * ny

                    arr[i].fx += dfx
                    arr[j].fx -= dfx
                    arr[i].fy += dfy
                    arr[j].fy -= dfy
                end
            end
            edgex = -5*(arr[i].x)^7
            edgey = -5*(arr[i].y)^7

            arr[i].fx += edgex
            arr[i].fy += edgey
        end
        delta = 0.01
        coeff= 0.7
        for i in 1:N
            node = arr[i]
            # friction
            nvx = node.vx/(pyth(node.vx,node.vy)+0.00001)
            nvy = node.vy/(pyth(node.vx,node.vy)+0.00001)

            node.fx -= nvx*coeff
            node.fy -= nvy*coeff
            # println(node)

            node.vx += node.fx * delta
            node.vy += node.fy * delta

            node.x += node.vx * delta
            node.y += node.vy * delta
        end
        if step % 100 == 0
            display_graph(GRN, arr)
        end
    end
end

function pyth(x,y)
    return sqrt(x^2 + y^2)
end

function display_graph(GRN, arr, title="test")

    on = "#ff0000"
    off = "#00ffff"
    activate = "#00ff00"
    repress = "#ff0000"
    buffer = 0.06

    N = length(arr)
    # p = plot(title=title,size=(800,800),lims=(-1,1))
    p = plot(title=title,size=(800,800))
    for i in 1:N
        c = arr[i].activated ?  on : off
        scatter!(p, [arr[i].x],[arr[i].y],label="",c=c,markersize=18)
        # quiver!(p, [arr[i].x],[arr[i].y],quiver=([arr[i].fx*0.1],[arr[i].fy*0.1]),c="#000000")
        # quiver!(p, [arr[i].x],[arr[i].y],quiver=([arr[i].vx*0.1],[arr[i].vy*0.1]),c="#ff0000")
        for j in 1:N
            if GRN[i,j] != 0 && i!=j
                dx = arr[j].x - arr[i].x
                dy = arr[j].y - arr[i].y

                nx = dx/pyth(dx,dy)
                ny = dy/pyth(dx,dy)

                dx = dx - 2*buffer*nx
                dy = dy - 2*buffer*ny

                c =  GRN[i,j] == 1 ? activate : repress

                quiver!(p, [arr[i].x + buffer*nx],
                           [arr[i].y + buffer*ny],
                           quiver=([dx], [dy]), c=c,
                           linestyle=:dashdot,
                           linewidth=3)

            end
        end
    end
    display(p)
end

GRN = load("results_new/Size50/replicate_101/adaptive_combined_evolved_GRN.dat")
force_sim(GRN)
