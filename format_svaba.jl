
using Fire
function IDdct(fs::String)
    dct = Dict{String,String}()
    for line in eachline(fs)
        if startswith(line, "#")
            continue
        else
            _, _, ID, _, _, _, _, _, _, _, _ = split(line, "\t")
            ID = convert(String, ID)
            dct[ID] = line
        end
    end
    dct
end

function parse_class(alt::AbstractString, alt2::AbstractString)
    class = ""
    strand = ""
    if match(r"\D\].+\]", alt) != nothing || match(r"\D\].+\]", alt2) != nothing
        class = "INV"
        strand = "3to3"
    elseif match(r"\[.+\[\D", alt) != nothing || match(r"\[.+\[\D", alt2) != nothing
        class = "INV"
        strand = "5to5"
    elseif match(r"\D\[.+\[", alt) != nothing || match(r"\].+\]\D", alt2) != nothing
        class = "DEL"
        strand = "3to5"
    elseif match(r"\].+\]\D", alt) != nothing || match(r"\D\[.+\[", alt2) != nothing
        class = "DUP"
        strand = "5to3"
    end
    (class, strand)
end


function svclass(ln::AbstractString, dct::Dict{String,String})
    shunxu = Dict("chr1" => 1, "chr2" => 2, "chr3" => 3, "chr4" => 4, "chr5" => 5, "chr6" => 6, "chr7" => 7, "chr8" => 8, "chr9" => 9, "chr10" => 10, "chr11" => 11, "chr12" => 12, "chr13" => 13, "chr14" => 14, "chr15" => 15, "chr16" => 16, "chr17" => 17, "chr18" => 18, "chr19" => 19, "chr20" => 20, "chr21" => 21, "chr22" => 22, "chrX" => 23, "chrY" => 24, "chrM" => 25)
    CHROM, POS, ID, _, alt, _, _, INFO, _, _, _ = split(ln, "\t")
    mateid = string(replace(match(r"MATEID=[\d:]+", INFO).match, "MATEID=" => ""))
    CHROM2, POS2, _, _, alt2, _, _, _, _, _, _ = split(dct[mateid], "\t")
    pesupport = 0
    if CHROM == CHROM2   ### intrachromosomal
        if POS > POS2
            tmp_chrom = CHROM2
            tmp_pos = POS2
            tmp_alt = alt2
            CHROM2 = CHROM
            POS2 = POS
            alt2 = alt
            CHROM = tmp_chrom
            POS = tmp_pos
            alt = tmp_alt
        end
        svtype, strand = parse_class(alt, alt2)
    else
        if shunxu[CHROM] > shunxu[CHROM2]
            tmp_chrom = CHROM2
            tmp_pos = POS2
            tmp_alt = alt2
            CHROM2 = CHROM
            POS2 = POS
            alt2 = alt
            CHROM = tmp_chrom
            POS = tmp_pos
            alt = tmp_alt
        end
        _, strand = parse_class(alt, alt2)
        svtype = "TRA"
    end
    POS = parse(Int64, POS)
    POS2 = parse(Int64, POS2)
    (CHROM, POS, CHROM2, POS2, strand, pesupport, svtype, ID, mateid)
end

function parse_svaba(fs::String, IO)
    mateids = []
    BNDs = IDdct(fs)
    HEADER = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod"
    println(IO, HEADER)
    for line in eachline(fs)
        if startswith(line, "#")
            continue
        else
            CHROM, POS, CHROM2, POS2, strand, pesupport, svtype, ID, mateid = svclass(line, BNDs)
            if ID in mateids
                continue
            else
                push!(mateids, mateid)
                if svtype == "TRA"
                    if strand == "5to5"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", "TRA", "\t", "Svaba")
                    elseif strand == "3to3"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", "TRA", "\t", "Svaba")
                    elseif strand == "5to3"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "TRA", "\t", "Svaba")
                    elseif strand == "3to5"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", "TRA", "\t", "Svaba")
                    else
                        @warn "Something was wrong.  You should pay attention"
                    end
                elseif svtype == "INV"
                    if strand == "5to5"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", svtype, "\t", "Svaba")
                    elseif strand == "3to3"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", svtype, "\t", "Svaba")
                    else
                        @warn "Something was wrong.  You should pay attention"
                    end
                elseif svtype == "DEL"
                    println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", svtype, "\t", "Svaba")
                else
                    println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", svtype, "\t", "Svaba")
                end
            end
        end
    end
end

@main function main(; input::String="none")
    for line in eachline(input)
        prx, fs = split(line, "\t")
        open("format.$(prx).sv.svaba.txt", "w+") do io
            parse_svaba(fs, io)
        end
    end
end





