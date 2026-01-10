using Glob
using CodecZlib
using StatsBase
using Plots
using StatsPlots
using CSV
using DataFrames
using FileIO
using Measures

# Load metadata
metadata = CSV.read("./metadata/HIBO_library_metadata.csv", DataFrame)

# Extract XX from Label and flag controls
metadata[!, :XX] = replace.(string.(metadata.Label), r"HB_" => "")
metadata[!, :is_control] =ismissing.(metadata.yrBP) .|| string.(metadata.yrBP) .== "NA"

# Create a mapping from "XX" to Label and control/sample type
label_map = Dict{String, Tuple{String, Bool}}()  # "XX" => (Label, is_control)
for row in eachrow(metadata)
    xx = row.XX
    is_ctrl = row.is_control
    label = is_ctrl ? "$(row.Label) (control)" : row.Label
    label_map[xx] = (label, is_ctrl)
end

# Function to get read lengths from a gzipped FASTQ file
function get_read_lengths(file_path::String)
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

# Build a Dict: XX => file path
fastq_files = sort(glob("../trim2ends/*.fq.gz"))

xx_file_map = Dict{String, String}()
for file in fastq_files
    m = match(r"AVXF-1-(\d+)", basename(file))
    if m !== nothing
        xx = m.captures[1]
        xx_file_map[xx] = file
    end
end

# Generate subplots in metadata order
plots_list = Plots.Plot[]
Threads.@threads for row in eachrow(metadata)
    xx = row.XX

    file = xx_file_map[xx]
    label, is_control = label_map[xx]
    color = is_control ? :red : :blue

    lengths = get_read_lengths(file)

    p = density(lengths,
                title=label,
                linewidth=1.5,
                linecolor=color,
                grid=false,
                legend=false,
                xlabel="",
                ylabel="")
    push!(plots_list, p)
end


# Plot in 10 rows x 5 columns
final_plot = plot(plots_list..., layout=(11, 4), left_margin = 5mm, size=(1000, 1600), dpi=300)

# Save to PDF
savefig(final_plot, "./LibrarySummary/output/read_length_distributions.pdf")

