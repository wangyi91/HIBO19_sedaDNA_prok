using Glob
using CodecZlib
using StatsBase
using Plots
using StatsPlots
using CSV
using DataFrames
using FileIO
using Measures
using BioAlignments, XAM
import XAM.BAM

# Load metadata
metadata = CSV.read("./metadata/HIBO_library_metadata.csv", DataFrame)

# Extract XX from Label and flag controls
metadata[!, :XX] = replace.(string.(metadata.Label), r"HB_" => "")
metadata[!, :is_control] = ismissing.(metadata.yrBP) .|| string.(metadata.yrBP) .== "NA"

# Create a mapping from "XX" to Label and control/sample type
label_map = Dict{String, Tuple{String, Bool}}()
for row in eachrow(metadata)
    xx = row.XX
    is_ctrl = row.is_control
    label = is_ctrl ? "$(row.Label) (control)" : row.Label
    label_map[xx] = (label, is_ctrl)
end

# Function to get read lengths from a gzipped FASTQ file
function get_read_lengths_fastq(file_path::String)
    lengths = Int[]
    open(GzipDecompressorStream, file_path) do io
        line_count = 0
        for line in eachline(io)
            line_count += 1
            if line_count % 4 == 2
                push!(lengths, length(line))
            end
        end
    end
    return lengths
end

# Function to get read lengths from a BAM file
function get_read_lengths_bam(file_path::String)
    lengths = Int[]
    reader = open(XAM.BAM.Reader, file_path)
    for record in reader
        if !XAM.BAM.ismapped(record)
            continue
        end
        seq = XAM.BAM.sequence(record)
        push!(lengths, length(seq))
    end
    close(reader)
    return lengths
end
function get_read_lengths_txt(file_path::String)
    return parse.(Int, readlines(file_path))
end
# Map fastq files: XX => file path
fastq_files = sort(glob("../trim2ends/*.fq.gz"))
xx_fastq_map = Dict{String, String}()
for file in fastq_files
    m = match(r"AVXF-1-(\d+)", basename(file))
    if m !== nothing
        xx = m.captures[1]
        xx_fastq_map[xx] = file
    end
end

# Map bam files: XX => file path
#bam_files = sort(glob("../EnvironmenTracker/results/taxonomic-profiling/*.dedup.bam"))
txt_files = sort(glob("./LibrarySummary/input/*read_len.txt"))
xx_bam_map = Dict{String, String}()
for file in txt_files
    m = match(r"HB_(\d+)\.read_len\.txt", basename(file))
    if m !== nothing
        xx = m.captures[1]
        xx_bam_map[xx] = file
    end
end

# Generate subplots
plots_list = Plots.Plot[]
Threads.@threads for row in eachrow(metadata[20:21,:])
    xx = row.XX
    label, is_control = label_map[xx]

    # Skip if no fastq file
    if !haskey(xx_fastq_map, xx)
        @warn "No FASTQ file found for $xx"
        continue
    end
    fastq_file = xx_fastq_map[xx]
    fastq_lengths = get_read_lengths_fastq(fastq_file)

    color_fastq = is_control ? :red : :blue
    color_bam = :black  # BAM always black overlay

    # Create base plot with FASTQ
    p = density(fastq_lengths;
                title=label,
                linewidth=1.5,
                linecolor=color_fastq,
                grid=false,
                legend=false,
                xlabel="",
                ylabel="")

    # If BAM file exists, add BAM read length overlay
    if haskey(xx_bam_map, xx)
        txt_file = xx_bam_map[xx]
        bam_lengths = get_read_lengths_txt(txt_file)
        density!(p, bam_lengths;
              linewidth=1.5,
              linecolor=color_bam,
              linestyle=:dot)
    else
        @warn "No BAM read lengths txt file found for $xx"
    end

    push!(plots_list, p)
end

# Final plot grid
final_plot = plot(plots_list..., layout=(11, 4), left_margin=5mm, size=(1000, 1600), dpi=300)

# Save to PDF
savefig(final_plot, "./LibrarySummary/output/read_length_distributions_raw_vs_mapped.pdf")

