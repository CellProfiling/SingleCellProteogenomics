#%% Imports
from imports import *
from Bio import SeqIO
# import requests
# import sys

#%% Import the genes names we're analyzing
ccd_transcript_regulated = np.array(pd.read_csv("output/ccd_transcript_regulated.csv")["gene"])
ccd_nontranscript_regulated = np.array(pd.read_csv("output/ccd_nontranscript_regulated.csv")["gene"])
genes_analyzed = np.array(pd.read_csv("output/gene_names.csv")["gene"])

#%% Let's take a look at protein stability
# Idea: Is there a difference in protein stability or turnover for the proteins that are
# transcriptionally regulated CCD vs non-transcriptionally regulated?
# Execution: 
# Output:

# Using the categories
def plot_protein_stabilities(filename, title, splitname):
    df = pd.read_csv(filename, delimiter="\t")
    df["ProteinName"] = df["Protein ID"].str.extract("[A-Z0-9]+_(.+)") if splitname else df["Protein ID"]
    ccd_t_stab = df[np.isin(df["ProteinName"], ccd_transcript_regulated)]
    ccd_n_stab = df[np.isin(df["ProteinName"], ccd_nontranscript_regulated)]

    df1 = df["Protein stability class"].value_counts().reindex(["non-stable", "medium", "stable", "non-melter"]).fillna(0)
    df2 = ccd_t_stab["Protein stability class"].value_counts().reindex(["non-stable", "medium", "stable", "non-melter"]).fillna(0)
    df3 = ccd_n_stab["Protein stability class"].value_counts().reindex(["non-stable", "medium", "stable", "non-melter"]).fillna(0)
    plt.plot(df1.index, df1 / np.array(df1).sum(axis=0, keepdims=True), label="All Proteins")
    plt.plot(df2.index, df2 / np.array(df2).sum(axis=0, keepdims=True), label="Transcript Reg. CCD")
    plt.plot(df3.index, df3 / np.array(df3).sum(axis=0, keepdims=True), label="Non-Transcript Reg. CCD")
    plt.title(title)
    plt.legend(loc="upper right")
    # plt.show()
    # plt.close()

plt.figure(figsize=(15,15))
plt.subplot(331)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True)
plt.subplot(332)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True)
plt.subplot(334)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True)
plt.subplot(335)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True)
plt.subplot(337)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False)
plt.subplot(338)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False)
plt.subplot(339)
plot_protein_stabilities("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False)
plt.savefig("figures/ProteinStabilities.png")
plt.show()
plt.close()

#%% [markdown]
# The stability classes don't show the greatest consistency between replicates
# or between cell lines. The melting points are a continuous category, so I might
# have more luck with it. It'll also omit the not-stable and non-melter categories,
# which seem like a bit of a crap shoot above

#%% Using melting temperatures (histograms)
# Idea: See if there is better consistency between replicates with the melting
#       points over the stability classes
# Execution: use melting points from data frames
# Output: histograms plots for the replicates and aggregate
#       statistical test for whether the transcript reg and non-transcript reg CCD
#       are different.

def weights(vals):
    '''normalizes all histogram bins to sum to 1'''
    return np.ones_like(vals)/float(len(vals))
    
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]

def add_temps(filename, title, splitname, merged):
    '''Adds melting temperature measurements from supp info files to lists'''
    df = pd.read_csv(filename, delimiter="\t")
    df["ProteinName"] = df["Protein ID"].str.extract("[A-Z0-9]+_(.+)") if splitname else df["Protein ID"]
    if len(merged) == 0: merged = df
    else: pd.merge(merged, df, on="ProteinName", suffixes=("", title))
    ccd_t_stab = df[np.isin(df["ProteinName"], ccd_transcript_regulated)]
    ccd_n_stab = df[np.isin(df["ProteinName"], ccd_nontranscript_regulated)]

    all_temps.extend(df[pd.notna(df["Melting point [°C]"])]["Melting point [°C]"])
    transcript_reg.extend(ccd_t_stab[pd.notna(ccd_t_stab["Melting point [°C]"])]["Melting point [°C]"])
    nontranscr_reg.extend(ccd_n_stab[pd.notna(ccd_n_stab["Melting point [°C]"])]["Melting point [°C]"])

def stat_tests(title):
    print(f"{title} statistical tests")
    print("Testing whether transcript reg. CCD melting points are different than non-transcript reg. CCD")
    print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
    print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
    print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
    print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
    print(scipy.stats.ttest_ind(transcript_reg, nontranscr_reg))
    print(scipy.stats.kruskal(transcript_reg, nontranscr_reg))
    print()
    print("Testing whether non-transcript reg. CCD is different than all proteins")
    print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
    print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
    print(f"{np.median(all_temps)}: All Proteins (median)")
    print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
    print(scipy.stats.ttest_ind(all_temps, nontranscr_reg))
    print(scipy.stats.kruskal(all_temps, nontranscr_reg))
    print()
    print("Testing whether transcript reg. CCD is different than all proteins")
    print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
    print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
    print(f"{np.median(all_temps)}: All Proteins (median)")
    print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
    print(scipy.stats.ttest_ind(all_temps, transcript_reg))
    print(scipy.stats.kruskal(all_temps, transcript_reg))
    print()

def temp_hist(title):
    '''Generates a histogram of melting points with bins normalized to 1'''
    bins=np.histogram(np.hstack((all_temps, transcript_reg, nontranscr_reg)), bins=40)[1] #get the bin edges
    plt.hist(all_temps, bins=bins, weights=weights(all_temps), alpha=0.5, label="All Proteins")
    plt.hist(transcript_reg, bins=bins, weights=weights(transcript_reg), alpha=0.5, label="Transcript Reg. CCD")
    plt.hist(nontranscr_reg, bins=bins, weights=weights(nontranscr_reg), alpha=0.6, label="Non-Transcript Reg. CCD")
    plt.legend(loc="upper right")
    plt.xlabel("Melting Point (°C)")
    plt.title(title)
    stat_tests(title)


# Individual histograms
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
plt.figure(figsize=(15,15))
plt.subplot(331)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True, merged)
temp_hist("A549_R1")
plt.subplot(332)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True, merged)
temp_hist("A549_R2")
plt.subplot(334)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True, merged)
temp_hist("HEK293_R1")
plt.subplot(335)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True, merged)
temp_hist("HEK293_R2")
plt.subplot(337)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False, merged)
temp_hist("HepG2_R1")
plt.subplot(338)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False, merged)
temp_hist("HepG2_R2")
plt.subplot(339)
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False, merged)
temp_hist("HepG2_R3")
plt.savefig("figures/ProteinMeltingPointsIndivid.png")
plt.show()
plt.close()

# Aggregate histogram
all_temps, transcript_reg, nontranscr_reg, merged = [],[],[],[]
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False, merged)
add_temps("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False, merged)
temp_hist("Aggregated Melting Points")
plt.savefig("figures/ProteinMeltingPoints.png")
plt.show()
plt.close()

# Statistical tests
print("Testing whether transcript reg. CCD melting points are different than non-transcript reg. CCD")
print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
print(scipy.stats.ttest_ind(transcript_reg, nontranscr_reg))
print(scipy.stats.kruskal(transcript_reg, nontranscr_reg))
print()
print("Testing whether non-transcript reg. CCD is different than all proteins")
print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
print(f"{np.mean(nontranscr_reg)} +/- {np.std(nontranscr_reg)}: Non-Transcript Reg. CCD (mean, std)")
print(f"{np.median(all_temps)}: All Proteins (median)")
print(f"{np.median(nontranscr_reg)}: Non-Transcript Reg. CCD (median)")
print(scipy.stats.ttest_ind(all_temps, nontranscr_reg))
print(scipy.stats.kruskal(all_temps, nontranscr_reg))
print()
print("Testing whether transcript reg. CCD is different than all proteins")
print(f"{np.mean(all_temps)} +/- {np.std(all_temps)}: All Proteins (mean, std)")
print(f"{np.mean(transcript_reg)} +/- {np.std(transcript_reg)}: Transcript Reg. CCD (mean, std)")
print(f"{np.median(all_temps)}: All Proteins (median)")
print(f"{np.median(transcript_reg)}: Transcript Reg. CCD (median)")
print(scipy.stats.ttest_ind(all_temps, transcript_reg))
print(scipy.stats.kruskal(all_temps, transcript_reg))

#%% [markdown]
# The replicates do look better with the melting points, and the difference is
# significant between transcript reg CCD and non-transcript reg CCD. However, the
# non-transcript regulated CCD proteins are more similar to all proteins on average,
# and that difference is not significant.

#%% Using melting temperatures (histograms, combined protein measurements)
# Idea: The above conclusions are based on tests that require independent observations,
#       but some proteins measurements are represented more than once (max of 12 times)
# Execution: use melting points from data frames; combine protein measurements
# Output: histograms for the replicates and aggregate
#       statistical test for whether the transcript reg and non-transcript reg CCD
#       are different.

atp, trp, ntp = [], [], []
all_temps, transcript_reg, nontranscr_reg = [],[],[]

def add_temps_and_names(filename, title, splitname):
    '''Adds melting temperature measurements from supp info files to lists'''
    df = pd.read_csv(filename, delimiter="\t")
    df["ProteinName"] = df["Protein ID"].str.extract("[A-Z0-9]+_(.+)") if splitname else df["Protein ID"]
    ccd_t_stab = df[np.isin(df["ProteinName"], ccd_transcript_regulated)]
    ccd_n_stab = df[np.isin(df["ProteinName"], ccd_nontranscript_regulated)]

    notna = pd.notna(df["Melting point [°C]"])
    all_temps.extend(df[notna]["Melting point [°C]"])
    transcript_reg.extend(ccd_t_stab[notna]["Melting point [°C]"])
    nontranscr_reg.extend(ccd_n_stab[notna]["Melting point [°C]"])

    atp.extend(df[notna]["ProteinName"])
    trp.extend(ccd_t_stab[notna]["ProteinName"])
    ntp.extend(ccd_n_stab[notna]["ProteinName"])

def avg_prot_temps(temps, proteinNames):
    '''Finds median of the measurements from the same protein'''
    aaa = []
    for i in range(len(temps)):
        df = pd.DataFrame({"name" : proteinNames[i], "temps": temps[i]})
        aaa.append(list(df.groupby("name")["temps"].median()))
    return aaa

# Individual histograms
plt.figure(figsize=(15,15))
plt.subplot(331)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("A549_R1")
plt.subplot(332)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("A549_R2")
plt.subplot(334)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HEK293_R1")
plt.subplot(335)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True)
temp_hist("HEK293_R2")
plt.subplot(337)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HepG2_R1")
plt.subplot(338)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HepG2_R2")
plt.subplot(339)
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("HepG2_R3")
plt.savefig("figures/ProteinMeltingPointsMedianedIndivid.png")
plt.show()
plt.close()

# Aggregate histogram
all_temps, transcript_reg, nontranscr_reg = [],[],[]
atp, trp, ntp = [], [], []
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R1.tsv", "A549_R1", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\A549_R2.tsv", "A549_R2", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R1.tsv", "HEK293_R1", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HEK293_R2.tsv", "HEK293_R2", True)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R1.tsv", "HepG2_R1", False)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R2.tsv", "HepG2_R2", False)
add_temps_and_names("C:\\Users\\antho\\Box\\ProjectData\\ProteinStability\\HepG2_R3.tsv", "HepG2_R3", False)
all_temps, transcript_reg, nontranscr_reg = avg_prot_temps([all_temps, transcript_reg, nontranscr_reg], [atp, trp, ntp])
temp_hist("Aggregated Melting Points")
plt.savefig("figures/ProteinMeltingPointsMedianed.png")
plt.show()
plt.close()

#%% [markdown]
# The replicates do look better with the melting points, and the difference is
# significant between transcript reg CCD and non-transcript reg CCD. However, the
# non-transcript regulated CCD proteins are more similar to all proteins on average,
# and that difference is not significant.

#%% Making boxplots and violin plots for the melting point data
# Idea: trying out different visualizations for the above analysis
# (boxplots work; violin not so much yet)

# data = [transcript_reg, nontranscr_reg, all_temps]
# labels=["Transcript Reg.\nCCD", "Non-Transcript Reg.\nCCD", "All Proteins"]
# plt.boxplot(data, labels=labels)
# plt.savefig("figures/ProteinMeltingPointsBoxplot.png")
# plt.show()
# plt.close()

# def adjacent_values(vals, q1, q3):
#     upper_adjacent_value = q3 + (q3 - q1) * 1.5
#     upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

#     lower_adjacent_value = q1 - (q3 - q1) * 1.5
#     lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
#     return lower_adjacent_value, upper_adjacent_value

# def set_axis_style(ax, labels):
#     ax.get_xaxis().set_tick_params(direction='out')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.set_xticks(np.arange(1, len(labels) + 1))
#     ax.set_xticklabels(labels)
#     ax.set_xlim(0.25, len(labels) + 0.75)
#     ax.set_xlabel('Sample name')

# fig, ax1 = plt.subplots(1, 1, figsize=(9,4))
# parts = ax1.violinplot(data, showmeans=True, showmedians=False, showextrema=False)
# quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75])
# whiskers = np.array([adjacent_values(sorted_array, q1, q3) for sorted_array, q1, q3 in zip(data, quartile1, quartile3)])
# whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]
# set_axis_style(ax1, labels)
# plt.show()

#%%
