import numpy as np
import pandas as pd

# Load the timeline.dat ignoring initial comments
# Typical columns are: i, j, t_first, t_last, fraction, total_time, encounters
col_names = ["res1", "res2", "t_first", "t_last", "fraction", "tot_time", "encounters"]
df = pd.read_csv("timeline1.dat", sep=r'\s+', comment='#', names=col_names)

# Filter only the rows where there was at least one contact (encounters > 0 and fraction > 0)
active_contacts = df[(df["encounters"] > 0) & (df["fraction"] > 0)].copy()

# Estimate of the effective number of observations (N_eff) based on encounters
encounters = active_contacts["encounters"]

print(f"Min encounters: {encounters.min()}")
print(f"Max encounters: {encounters.max()}")
print(f"Lower quartile (Q1) of encounters: {np.percentile(encounters, 25):.1f}")
print(f"Median encounters: {encounters.median():.1f}")
