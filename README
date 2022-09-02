# mergeSV

1. you can use the merged result to call chromothripsis based ShatterSeek.

2. you can use the scripts to classify the BND from svaba into SVtype.


## sv callers we can use 

- Delly
- NovoBreak
- Manta
- Svaba

## Usage

#### formatting the four vcf to bed format


```bash
julia format_delly.jl --input input_delly_path
julia format_manta.jl --input input_manta_path
julia format_novobreak.jl --input input_novobreak_path
julia format_svaba.jl --input input_svaba_path
```

**Note:**  input file should be a tab delimated file with sample name and vcf path, eg: ```sample24 examples/24_delly.vcf```

#### merge sv from four callers

*please notice that we only keep the SVs exited more than one caller!*

```bash
julia mergeVCF.jl --path merged --samples example/samples
```
`path` should be the file path after formatting

`samples` should be a file that contains all names


*for manta, you should firstly use their scrpts to transfer the BND into INV*

