using Fire

function parse_line(ln::AbstractString)
    CHROM, POS, ID, _, _, _, _, INFO, _, _, ot... = split(ln, "\t")
    svclass = replace(match(r"SVTYPE=\w+", INFO).match, "SVTYPE=" => "")
    pesupport = ot[4]
    strand = replace(match(r"CT=\Sto\S", INFO).match, "CT=" => "")
    POS2 = replace(match(r"END=\d+", INFO).match, "END=" => "")
    CHROM2 = replace(match(r"CHR2=chr[\d\w_]+", INFO).match, "CHR2=" => "")
    POS2 = parse(Int64, POS2)
    POS = parse(Int64, POS)
    (CHROM, POS, CHROM2, POS2, strand, pesupport, svclass, ID)
end

function parse_novobreak(fs::String, IO)
    shunxu = Dict("chr1" => 1, "chr2" => 2, "chr3" => 3, "chr4" => 4, "chr5" => 5, "chr6" => 6, "chr7" => 7, "chr8" => 8, "chr9" => 9, "chr10" => 10, "chr11" => 11, "chr12" => 12, "chr13" => 13, "chr14" => 14, "chr15" => 15, "chr16" => 16, "chr17" => 17, "chr18" => 18, "chr19" => 19, "chr20" => 20, "chr21" => 21, "chr22" => 22, "chrX" => 23, "chrY" => 24, "chrM" => 25)
    HEADER = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod"
    println(IO, HEADER)
    for line in eachline(fs)
        if startswith(line, "#")
            continue
        else
            CHROM, POS, CHROM2, POS2, strand, pesupport, svclass, ID = parse_line(line)
            if CHROM ∉ keys(shunxu) || CHROM2 ∉ keys(shunxu)
                continue
            else
                if svclass == "TRA"
                    if strand == "5to5"
                        if shunxu[CHROM] > shunxu[CHROM2]
                            println(IO, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", "TRA", "\t", "NovoBreak")
                        else
                            println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", "TRA", "\t", "NovoBreak")
                        end
                    elseif strand == "3to3"
                        if shunxu[CHROM] > shunxu[CHROM2]
                            println(IO, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", "TRA", "\t", "NovoBreak")
                        else
                            println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", "TRA", "\t", "NovoBreak")
                        end
                    elseif strand == "5to3"
                        if shunxu[CHROM] > shunxu[CHROM2]
                            println(IO, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", "TRA", "\t", "NovoBreak")
                        else
                            println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "TRA", "\t", "NovoBreak")
                        end
                    elseif strand == "3to5"
                        if shunxu[CHROM] > shunxu[CHROM2]
                            println(IO, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "TRA", "\t", "NovoBreak")
                        else
                            println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", "TRA", "\t", "NovoBreak")
                        end
                    else
                        @warn "Something was wrong.  You should pay attention"
                    end
                elseif svclass == "INV"
                    if strand == "5to5"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", svclass, "\t", "NovoBreak")
                    elseif strand == "3to3"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", svclass, "\t", "NovoBreak")
                    else
                        @warn "Something was wrong.  You should pay attention"
                    end
                elseif svclass == "DEL"
                    println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", svclass, "\t", "NovoBreak")
                else
                    println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "DUP", "\t", "NovoBreak")
                end
            end
        end
    end
end


@main function main(; input::String="none")
    for line in eachline(input)
        prx, fs = split(line, "\t")
        open("format.$(prx).sv.novobreak.txt", "w+") do io
            parse_novobreak(fs, io)
        end
    end
end
