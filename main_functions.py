# import pandas as pd
# import numpy as np

def speceies_continuity_extended(sitedata, ci):
    """
    Calculate species estimated stratigraphic ranges based on sampling gaps and species occurrence gaps.

    Parameters:
    - sitedata (pd.DataFrame):
        A DataFrame containing diversity data for the targeted site.
        Required columns include:
        - 'taxon_name': The species name.
        - 'sample_age_ma': The age of each sample in million years (Ma).
    - ci (float):
        A confidence interval value (between 0 and 1) used to extend the stratigraphic range
        by quantiles of the sampling gaps.

    Returns:
    - continuity_df (pd.DataFrame), its each column name explained below:
        A DataFrame summarizing the stratigraphic range and continuity of each species, with the following columns:
        - taxon_name: The name of the species.
        - continuity_percentage: The percentage of occurrence of the species within its observed stratigraphic range.
        - observed_datums: A list of ages where the species was observed.
        - LAD: The observed last appearance datum of the species (minimum age observed).
        - FAD: The observed first appearance datum of the species (maximum age observed).
        - LAD_extended: The estimated last appearance datum, extended based on the confidence interval.
        - FAD_extended: The estimated first appearance datum, extended based on the confidence interval.

    Notes:
    - The function assumes the input DataFrame is sorted by 'sample_age_ma'.
        If not, it will sort the data in-place.
    - Species with a single occurrence (where FAD equals LAD) are excluded from the results.
    - Gaps between consecutive observations are calculated to derive the extensions using the specified confidence interval.
    """
    sitedata.sort_values('sample_age_ma', inplace = True)
    grouped_data = sitedata.groupby('taxon_name')

    results = []
    for species_name, group in grouped_data:
        species_ages = group['sample_age_ma'].values
        FAD = max(species_ages)
        LAD = min(species_ages)

        if FAD == LAD:
           continue

        sample_ages_in_range = sitedata[(sitedata['sample_age_ma'] <= FAD) & (sitedata['sample_age_ma'] >= LAD)]['sample_age_ma'].unique()
        overlap = set(sample_ages_in_range).intersection(species_ages)
        overlap_percentage = (len(overlap) / len(sample_ages_in_range)) * 100

        LAD_new = LAD - np.quantile(np.diff(species_ages), ci)
        FAD_new = FAD + np.quantile(np.diff(species_ages), ci)

        results.append({'taxon_name': species_name, 'continuity_percentage' : overlap_percentage, 'observed_datums' : species_ages, 'LAD' : LAD, 'FAD' : FAD, 'LAD_extended' : LAD_new, 'FAD_extended' : FAD_new})

    continuity_df = pd.DataFrame(results)
    return continuity_df

def create_bins(start_ma, end_ma, bin_width):
    """
    Generate a series of bins (time intervals) based on the specified start and end times and bin width.

    Parameters:
    - start_ma (float):
        The start of the time period.
    - end_ma (float):
        The end of the time period. Naturally, it must be greater than or equal to `start_ma`.
    - bin_width (float):
        The width of each bin. Guess what? Must be positive.

    Returns:
    - bins (list of lists):
        A list of time intervals (bins), where each bin is represented as a two-element list [start, end].
        The bins are constructed from `start_ma` to `end_ma`, inclusive of the start and ensuring the last bin does not exceed `end_ma`.

    Notes:
    - The function rounds bin boundaries to four decimal places to prevent floating-point inaccuracies.
    - The final bin's end time is adjusted to exactly match `end_ma` if it would otherwise exceed the specified range.
    - A small offset (`0.0001 Ma`) is applied to handle edge cases and prevent overlapping or skipped bins due to rounding errors.

    Example:
    >>> create_bins(0, 10, 3)
    [[0, 3], [3.0001, 6], [6.0001, 9], [9.0001, 10]]
    """
    bins = []
    current_bin_start = start_ma

    # Create bins until the current bin start is less than the end_ma
    while current_bin_start < end_ma:
        # Calculate the end of the current bin
        bin_end = round((current_bin_start -0.0001) + bin_width, 4)

        # If bin_end is less than end_ma, adjust it to end_ma
        if bin_end > end_ma:
            bin_end = end_ma  # Ensure the last bin does not exceed the end range

        # Append the bin as [start, end]
        bins.append([current_bin_start, bin_end])

        # Update the start of the next bin
        current_bin_start = round(bin_end + 0.0001, 4)

    return bins


def foote_updated(data, start_ma, end_ma, bin_width):
    """
    Calculate Foote's origination and extinction rates using First Appearance Datum (FAD) and Last Appearance Datum (LAD) estimates.

    Parameters:
    - data (pd.DataFrame):
        A DataFrame containing the following required columns:
        - 'FAD_extended': Estimated First Appearance Datum (FAD).
        - 'LAD_extended': Estimated Last Appearance Datum (LAD).
    - start_ma (float):
        The start of the time period.
    - end_ma (float):
        The end of the time period.
    - bin_width (float):
        The width of each time bin in million years.

    Returns:
    - counters (pd.DataFrame):
        A DataFrame with the following columns:
        - `mids`: Midpoints of each bin.
        - `Bins`: Time bins represented as [start, end] sublists.
        - `tOrix`: Count of origination events in each bin.
        - `tExtx`: Count of extinction events in each bin.
        - `tThroughx`: Count of taxa present throughout each bin.
        - `extinction_rate`: Raw extinction rate for each bin (calculated using the Foote method).
        - `origination_rate`: Raw origination rate for each bin (calculated using the Foote method).
        - `foote_extinction_rate_normalized`: Extinction rate normalized by the bin interval length.
        - `foote_origination_rate_normalized`: Origination rate normalized by the bin interval length.

    Notes:
        - The function uses the `create_bins` utility to generate time bins and their midpoints.
        - Extinction and origination rates are calculated followin Foote, 2000.
    """
    # Create bins as (start, end) tuples
    bins = create_bins(start_ma, end_ma, bin_width)

    # Initialize counters DataFrame with bins and their midpoints
    counters = pd.DataFrame({
        'mids': [(b[0] + b[1]) / 2 for b in bins],  # Calculate midpoint of each bin
        'Bins': bins,
        'tOrix': 0,  # Initialize origination count
        'tExtx': 0,  # Initialize extinction count
        'tThroughx': 0  # Initialize taxa-through count
    })

    # Iterate through each taxon in the dataset
    for index, row in data.iterrows():
        # Find the bin index for FAD and LAD
        bin_index_fad = next((i for i, b in enumerate(bins) if b[0] <= row['FAD_extended'] <= b[1]), None)
        bin_index_lad = next((i for i, b in enumerate(bins) if b[0] <= row['LAD_extended'] <= b[1]), None)

        # Check if both FAD and LAD fall within the bins
        if bin_index_fad is not None and bin_index_lad is not None:
            # Increment origination count at the FAD bin
            counters.at[bin_index_fad, 'tOrix'] += 1
            # Increment extinction count at the LAD bin
            counters.at[bin_index_lad, 'tExtx'] += 1

            # Increment taxa-through count for all bins between LAD and FAD
            if bin_index_lad + 1 < bin_index_fad:  # Ensure at least one bin between LAD and FAD
                counters.loc[bin_index_lad + 1 : bin_index_fad - 1, 'tThroughx'] += 1
            elif bin_index_lad == bin_index_fad:
                counters.loc[bin_index_lad, 'tThroughx'] += 1

    # Calculate extinction and origination rates
    counters['extinction_rate'] = -np.log(counters['tThroughx'] / (counters['tExtx'] + counters['tThroughx']))
    counters['origination_rate'] = -np.log(counters['tThroughx'] / (counters['tOrix'] + counters['tThroughx']))

    # Normalize the rates by the time interval length
    counters['foote_extinction_rate_normalized'] = counters['extinction_rate'] / bin_width
    counters['foote_origination_rate_normalized'] = counters['origination_rate'] / bin_width

    # Replace any NaN values resulting from log(0) or division by zero
    counters.fillna(0, inplace=True)

    return counters

def chao1(data):
    """
    Calculate the Chao1 diversity estimate for a sample, including lower and upper confidence intervals.

    Parameters:
    - data (pd.DataFrame):
        A DataFrame containing sample diversity data with at least one column:
        - 'value': 'number of times observed' of each species in the sample.

    Returns:
    - Schao1 (float):
        The Chao1 diversity estimate, which accounts for unseen species based on rare species in the sample.
    - LCI (float):
        The lower 95% confidence interval for the Chao1 estimate.
    - UCI (float):
        The upper 95% confidence interval for the Chao1 estimate.

    Notes:
    - Chao1 estimates are based on the observed species richness and the number of rare species with occurrences of 1 (`n1`: singleton) and 2 (`n2`: doubleton).
    - If `n2` (doubleton) is 0, an alternative variance formula is applied.
    - Confidence intervals (LCI and UCI) are calculated using a log-normal distribution approach.
    """
    n1 = np.sum(data['value'] == 1)
    n2 = np.sum(data['value'] == 2)
    Schao1 = len(data) + ((n1 ** 2 - 1) / (2 * n2 + 2))

    if n1 > 0 and n2 >= 0:
        varSchao1 = (n1 * (n1 - 1)) / (2 * (n2 + 1)) + (n1 * ((2 * n1 - 1) ** 2)) / (4 * ((n2 + 1) ** 2)) + ((n1 ** 2) * n2 * ((n1 - 1) ** 2)) / (4 * ((n2 + 1) ** 4))
    elif n1 > 0 and n2 == 0:
        varSchao1 = (n1 * (n1 - 1) / 2) + (n1 * ((2 * n1 - 1) ** 2)) / 4 - (n1 ** 4) / (4 * Schao1)
    else:
        varSchao1 = len(data) * np.exp(-np.sum(data['value']) / len(data)) * (1 - np.exp(-np.sum(data['value']) / len(data)))

    C = np.exp(1.96 * np.sqrt(np.log(1 + (varSchao1 / (Schao1 - len(data)) ** 2))))
    LCI = len(data) + ((Schao1 - len(data)) / C)
    UCI = len(data) + C * (Schao1 - len(data))

    return Schao1, LCI, UCI

def site_chao1(data_grouped):
    """
    Calculate Chao1 diversity estimates for grouped samples, including 95% confidence intervals.
    This function applies the `chao1` function to each group (e.g., grouped by sample age)
    to compute the Chao1 diversity estimate along with lower and upper confidence intervals.

    Parameters:
    - data_grouped (pd.DataFrameGroupBy):
        A pandas GroupBy object, where each group represents diversity data for a specific sample or time interval.
        Each group must contain:
        - 'value': Occurrence of each species or taxon in the sample.
        - 'sample_age_ma': Age of the sample.

    Returns:
    - site_chao1_df (pd.DataFrame):
        A DataFrame summarizing Chao1 diversity estimates and related information for each group, with the following columns:
        - `rawdiv`: The observed species richness (number of unique taxa) in each sample.
        - `chao1`: The Chao1 diversity estimate for each sample.
        - `lci`: The lower 95% confidence interval for the Chao1 estimate.
        - `uci`: The upper 95% confidence interval for the Chao1 estimate.
        - `Age`: The age of the sample (from 'sample_age_ma').

    Notes:
    - The input must be grouped by `sample_age_ma` or a similar grouping variable that uniquely identifies each sample (e.g., sample depth, mbsf).
    - The `chao1` function is applied to each group to calculate diversity metrics.
    """
    site_chao1_list = []
    lci_chao1_list = []
    uci_chao1_list = []
    age_list = []
    raw_list = []

    for _, group in data_grouped:
        Schao1, LCI, UCI = chao1(group)
        site_chao1_list.append(Schao1)
        lci_chao1_list.append(LCI)
        uci_chao1_list.append(UCI)
        raw_list.append(len(group))
        age_list.append(group['sample_age_ma'].iloc[0])

    site_chao1_df = pd.DataFrame({
        'rawdiv': raw_list,
        'chao1': site_chao1_list,
        'lci': lci_chao1_list,
        'uci': uci_chao1_list,
        'Age': age_list
    })

    return site_chao1_df
