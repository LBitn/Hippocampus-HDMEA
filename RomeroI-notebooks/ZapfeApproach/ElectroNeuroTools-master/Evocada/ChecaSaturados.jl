using JLD,HDF5

listacompleta=readdir()

println("tuabuela")

for i in listacompleta
  #  println("listacompleta no esta vacia")
  #  println(i)
    if (i[end-2:end]=="jld" && i[1:3]=="LFP")
        try
            saturados=load(i)["CanalesSaturados"]
            buenos=load(i)["Canalesrespuesta"]
            nsatu=length(saturados)
            nbuen=length(buenos)
            println(i, " ", nsatu, " ",nbuen)
        catch
            println(i, " esto no tiene lo que buscas")
        end
    end
end
