using Glob
using StatsBase
#using Plots
using Base.Threads
using StatsPlots
default(fontfamily = "DejaVu Sans")
using CSV
using DataFrames
using FileIO
using Measures

# Load metadata
metadata = CSV.read("./metadata/HIBO_library_metadata.csv", DataFrame);

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

function get_read_lengths_txt(file_path::String)
    return parse.(Int, readlines(file_path))
end

# Map fastq files: XX => file path
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
txt_files_mapped = sort(glob("./LibrarySummary/input/*.mapped_read_len.txt"));
txt_files_mapped = sort(glob("./LibrarySummary/input/*.filtered_bam_read_len.txt"));

xx_bam_map = Dict{String, String}()
for file in txt_files_mapped
    m = match(r"HB_(\d+)\.filtered_bam_read_len\.txt", basename(file))
    if m !== nothing
        xx = m.captures[1]
        xx_bam_map[xx] = file
    end
end

# Generate summary plots
# --- Merge all "sample" data into one combined plot ---

# Collect all read lengths from sample entries
nthd = nthreads()
local_all = [Vector{Int}() for _ in 1:nthd];
local_mapped = [Vector{Int}() for _ in 1:nthd];

# Use index-based iteration (safe) rather than eachrow()
N = size(metadata, 1)
@threads for idx in 1:N
    tid = threadid()
    xx = metadata[idx, :XX]           # or metadata.XX[idx]
    _, is_control = label_map[xx]

    if is_control
        continue
    end

    if haskey(xx_fastq_map, xx)
        append!(local_all[tid], get_read_lengths_txt(xx_fastq_map[xx]))
    else
        @warn "Missing FASTQ for sample $xx"
    end

    if haskey(xx_bam_map, xx)
        append!(local_mapped[tid], get_read_lengths_txt(xx_bam_map[xx]))
    else
        @warn "Missing BAM for sample $xx"
    end
end

# combine results once on the main thread
sample_all_lengths = vcat(local_all...);
sample_mapped_lengths = vcat(local_mapped...);


# Create density plot for combined sample data
combined_plot = plot(
    title="",
    xlabel="Read Length",
    ylabel="Density",
    legend=:topright,
    grid=false,
    xtickfont=font(8),
    ytickfont=font(8),
    left_margin=5mm, bottom_margin=5mm,
    size=(400, 300)
)

density!(combined_plot, sample_all_lengths; 
         bandwidth = 0.5, linewidth=2, linecolor=:black, linestyle=:dot, 
         label="All Reads"
        )

density!(combined_plot, sample_mapped_lengths; 
         bandwidth = 0.5, linewidth=2, linecolor=:blue, 
             label="Mapped Reads, Passed Quality Filtration")

# Save combined plot
savefig(combined_plot, "./LibrarySummary/output/read_length_distribution_combined_samples_only.pdf")

