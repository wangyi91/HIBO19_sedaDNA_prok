using Glob
using StatsBase
using Plots
using StatsPlots
using CSV
using DataFrames
using FileIO
using Measures
using Base.Threads
# Load metadata
metadata = CSV.read("./metadata/HIBO_library_metadata.csv", DataFrame)

# Extract XX from Label and flag controls
metadata[!, :XX] = replace.(string.(metadata.Label), r"HB_" => "")
metadata[!, :is_control] = ismissing.(metadata.yrBP) .|| string.(metadata.yrBP) .== "NA"
metadata[!, :depth] = [s == "NA" ? missing : ceil(Int, parse(Float64, s)) for s in metadata.Middle_depth]

# Create a mapping from "XX" to Label and control/sample type
label_map = Dict{String, Tuple{String, Bool}}()
for row in eachrow(metadata)
    xx = row.XX
    is_ctrl = row.is_control
    label = is_ctrl ? "$(row.Label) (control)" : "$(row.Label) ($(row.depth)cm)"
    label_map[xx] = (label, is_ctrl)
end

function get_read_lengths_txt(file_path::String)
    return parse.(Int, readlines(file_path))
end

# Map fastq files: XX => file path
#fastq_files = sort(glob("../trim2ends/*.fq.gz"))
txt_files_all = sort(glob("./LibrarySummary/input/*.all_read_len.txt"))

xx_fastq_map = Dict{String, String}()
for file in txt_files_all
    m = match(r"HB_(\d+)\.all_read_len\.txt", basename(file))
    if m !== nothing
        xx = m.captures[1]
        xx_fastq_map[xx] = file
    end
end

# Map bam files: XX => file path
#bam_files = sort(glob("../EnvironmenTracker/results/taxonomic-profiling/*.dedup.bam"))
txt_files_mapped = sort(glob("./LibrarySummary/input/*.filtered_bam_read_len.txt"))

xx_bam_map = Dict{String, String}()
for file in txt_files_mapped
    m = match(r"HB_(\d+)\.filtered_bam_read_len\.txt", basename(file))
    if m !== nothing
        xx = m.captures[1]
        xx_bam_map[xx] = file
    end
end

# Generate subplots
plots_list = Plots.Plot[]
@threads for row in eachrow(metadata)
    xx = row.XX
    label, is_control = label_map[xx]

    # Skip if no fastq file
    if !haskey(xx_fastq_map, xx)
        @warn "No FASTQ file found for $xx"
        continue
    end
    all_file = xx_fastq_map[xx]
    all_lengths = get_read_lengths_txt(all_file)

    color_all = is_control ? :darkred : :darkblue
    color_mapped = is_control ? :red : :blue

    # Create base plot with FASTQ
    p = density(all_lengths; trim = true, bandwidth = 0.5,
                title=label,
                linewidth=1.5, linecolor=color_all, linestyle=:dot,
                grid=false,
                legend=false,
                xlabel="", ylabel="")

    # If BAM file exists, add BAM read length overlay
    if haskey(xx_bam_map, xx)
        mapped_file = xx_bam_map[xx]
        mapped_lengths = get_read_lengths_txt(mapped_file)
        density!(p, mapped_lengths;
                 trim = true,
                 bandwidth = 0.5,
                 linewidth=1.5,
              linecolor=color_mapped)
    else
        @warn "No BAM read lengths txt file found for $xx"
    end

    push!(plots_list, p)
end

# Final plot grid
final_plot = plot(plots_list..., layout=(11, 4), left_margin=5mm, size=(1000, 1600), dpi=300)

# Save to PDF
savefig(final_plot, "./LibrarySummary/output/read_length_distributions_raw_vs_mapped.pdf")

