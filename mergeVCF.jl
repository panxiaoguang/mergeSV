using CSV
using DataFrames
using Fire


readin(x::String, y::String, z::String) = DataFrame(CSV.File("$(z)/format.$(x).sv.$(y).txt", delim="\t"))

function mergeSV(sample::String, path::String, IO::IO)
    alls = foldl((x, y) -> vcat(x, y), readin.(Ref(sample), ["delly", "novobreak", "manta", "svaba"], Ref(path)))
    sort!(alls, [:svclass, :chrom1, :start1])
    alls.label = collect(1:nrow(alls))
    shunxu = Dict("Manta" => 1, "Svaba" => 2, "Delly" => 3, "NovoBreak" => 4)
    buyao = []
    HEADER = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod"
    println(IO, HEADER)
    for row in eachrow(alls)
        if row.label in buyao
            continue
        else
            panxuan = (alls.chrom1 .== row.chrom1) .& (alls.chrom2 .== row.chrom2) .& (row.start1 - 100 .< alls.start1 .< row.start1 + 100) .& (row.start2 - 100 .< alls.start2 .< row.start2 + 100)
            need = alls[panxuan, :]
            need = sort(need, order(:svmethod, by=(x) -> shunxu[x]))
            if length(unique(need.svmethod)) > 1
                need.mix = fill(join(unique(need.svmethod), ","), nrow(need))
                woyao = need[1, :]
                println(IO, woyao.chrom1, "\t", woyao.start1, "\t", woyao.end1, "\t", woyao.chrom2, "\t", woyao.start2, "\t", woyao.end2, "\t", woyao.sv_id, "\t", woyao.pe_support, "\t", woyao.strand1, "\t", woyao.strand2, "\t", woyao.svclass, "\t", woyao.mix)
                append!(buyao, need.label)
            end
        end
    end
end

@main function main(; path::String="none", samples::String="none")
    for line in eachline(samples)
        sm = line
        open("$(sm).merged.sv.txt", "w+") do io
            mergeSV(sm, path, io)
        end
    end
end

