module ArbolGrafosTikz

using LightGraphs, TikzGraphs

# Tipos y estructuras, no se porque se usa el abstract primero pero bueno.

abstract type LabelledDiGraph
end

export LabelledTree, Nodo, Subarbol, Arbol,
       walk_tree, walk_tree!, tikz_representation

struct LabelledTree <: LabelledDiGraph
    g::DiGraph
    labels::Vector{Any}
end

mutable struct Nodo
    nombre::Any
end

mutable struct Subarbol
    raiz::Nodo
    hijos::Union{Nothing, Array{Nodo, 1}}
end

mutable struct Arbol
    nodostodos::Vector{Subarbol}
end


# Funciones. ----------------------------------------------------

add_numbered_vertex!(g) = (add_vertex!(g); top = nv(g)) 

function walk_tree!(g, labels, arbol)

    ll=length(arbol.nodostodos)
    raiz_vertex=0
    
    for l in 1:ll
        nodolabel=arbol.nodostodos[l].raiz.nombre
       
        if !(nodolabel in labels)
            raiz_vertex = add_numbered_vertex!(g)
            push!(labels,nodolabel)
        else
            raiz_vertex=findfirst(x->x==nodolabel, labels)
        end
        
        ramas=arbol.nodostodos[l].hijos
        for r in ramas
            nodolabel=r.nombre
            rnum=add_numbered_vertex!(g)
            push!(labels, nodolabel)
            add_edge!(g, raiz_vertex, rnum)
        end
    end
   
end

function walk_tree(arbol::Arbol)
    g = DiGraph()
    labels = Any[]

    walk_tree!(g, labels, arbol)

    return LabelledTree(g, labels)

end

function tikz_representation(tree::LabelledDiGraph)
    labels = String[string(x) for x in tree.labels]
    return TikzGraphs.plot(tree.g, labels)
end
    


end #module
