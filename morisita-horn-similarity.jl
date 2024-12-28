# using DataFrames

# data = groupby(side_data, :sample_age_ma)

function morisitahorn(data)
    """
    Calculates the Morisita-Horn index for a given diversity dataset and returns similarity matrix.

    Parameters:
    - data (GroupedDataFrame):
        A `GroupedDataFrame` object obtained by grouping a `DataFrame` by samples (e.g., ages). Each group should include the following fields:
        - 'taxon_name': species name
        - 'value': Occurrence of each species or taxon in the sample.

    - `Matrix{Float64}`: A square matrix of size `dim × dim`, where `dim` is the number of samples (groups) in `data`.
        Each element `[i, j]` represents the Morisita-Horn similarity index between samples `i` and `j`, ranging from:
        - `0`: No similarity.
        - `1`: Complete similarity.
    """
    dim = length(data)
    matrix = zeros(Float64, dim, dim)

    for i in 1:dim, j in 1:dim

        common = intersect(data[i].taxon_name, data[j].taxon_name)

        Na = sum(data[i].value)
        Nb = sum(data[j].value)

        ga = filter(row -> row.taxon_name in common, data[i])
        gb = filter(row -> row.taxon_name in common, data[j])

        # Ensure ga and gb are sorted based on taxon_name

        sort!(ga, :taxon_name)
        sort!(gb, :taxon_name)

        da = sum((ga.value) .^ 2) / Na^2
        db = sum((gb.value) .^ 2) / Nb^2

        g = (ga.value .* gb.value)

        a = (2 * sum(g)) / ((da + db) * (Na * Nb))

        matrix[i, j] = a
    end

    return matrix
end
