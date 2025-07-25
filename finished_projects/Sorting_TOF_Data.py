import matplotlib.pyplot as plt
import pandas as pd
import json


json_data = {
  
  "544013": {
    "withTOF": 20.45,
    "withoutTOF": 79.55,
    "withP03": 85.57
  },
  "544028": {
    "withTOF": 25.90,
    "withoutTOF": 74.10,
    "withP03": 93.46
  },
  "544032": {
    "withTOF": 19.91,
    "withoutTOF": 80.09,
    "withP03": 86.36
  },
  "544091": {
    "withTOF": 0.00,
    "withoutTOF": 100.00,
    "withP03": 93.87
  },
  "544095": {
    "withTOF": 0.00,
    "withoutTOF": 100.00,
    "withP03": 86.84
  },
  "544098": {
    "withTOF": 20.91,
    "withoutTOF": 79.09,
    "withP03": 85.73
  },
  "544116": {
    "withTOF": 24.36,
    "withoutTOF": 75.64,
    "withP03": 93.09
  },
  "544121": {
    "withTOF": 0.00,
    "withoutTOF": 100.00,
    "withP03": 94.12
  },
  "544122": {
    "withTOF": 20.90,
    "withoutTOF": 79.10,
    "withP03": 85.60
  },
  "544123": {
    "withTOF": 31.94,
    "withoutTOF": 68.06,
    "withP03": 94.38
  },
  "544124": {
    "withTOF": 21.96,
    "withoutTOF": 78.04,
    "withP03": 85.05
  },
  "544184": {
    "withTOF": 26.51,
    "withoutTOF": 73.49,
    "withP03": 93.23
  },
  "544185": {
    "withTOF": 27.62,
    "withoutTOF": 72.38,
    "withP03": 93.35
  },
  "544389": {
    "withTOF": 26.31,
    "withoutTOF": 73.69,
    "withP03": 93.43
  },
  "544390": {
    "withTOF": 20.98,
    "withoutTOF": 79.02,
    "withP03": 85.58
  },
  "544391": {
    "withTOF": 21.21,
    "withoutTOF": 78.79,
    "withP03": 85.54
  },
  "544392": {
    "withTOF": 20.88,
    "withoutTOF": 79.12,
    "withP03": 85.39
  },
  "544451": {
    "withTOF": 0.00,
    "withoutTOF": 100.00,
    "withP03": 93.93
  },
  "544454": {
    "withTOF": 21.00,
    "withoutTOF": 79.00,
    "withP03": 85.74
  },
  "544474": {
    "withTOF": 21.76,
    "withoutTOF": 78.24,
    "withP03": 93.48
  },
  "544475": {
    "withTOF": 18.86,
    "withoutTOF": 81.14,
    "withP03": 85.82
  },
  "544476": {
    "withTOF": 19.38,
    "withoutTOF": 80.62,
    "withP03": 85.64
  },
  "544477": {
    "withTOF": 19.59,
    "withoutTOF": 80.41,
    "withP03": 85.54
  },
  "544490": {
    "withTOF": 19.86,
    "withoutTOF": 80.14,
    "withP03": 93.11
  },
  "544491": {
    "withTOF": 23.66,
    "withoutTOF": 76.34,
    "withP03": 93.62
  },
  "544492": {
    "withTOF": 19.05,
    "withoutTOF": 80.95,
    "withP03": 85.61
  },
  "544508": {
    "withTOF": 21.24,
    "withoutTOF": 78.76,
    "withP03": 93.19
  },
  "544510": {
    "withTOF": 25.64,
    "withoutTOF": 74.36,
    "withP03": 93.47
  },
  "544511": {
    "withTOF": 20.14,
    "withoutTOF": 79.86,
    "withP03": 86.10
  },
  "544512": {
    "withTOF": 20.73,
    "withoutTOF": 79.27,
    "withP03": 86.36
  },
  "544514": {
    "withTOF": 21.18,
    "withoutTOF": 78.82,
    "withP03": 86.33
  },
  "544515": {
    "withTOF": 21.18,
    "withoutTOF": 78.82,
    "withP03": 86.22
  },
  "544518": {
    "withTOF": 21.32,
    "withoutTOF": 78.68,
    "withP03": 86.12
  },
  "544548": {
    "withTOF": 23.95,
    "withoutTOF": 76.05,
    "withP03": 93.46
  },
  "544549": {
    "withTOF": 23.96,
    "withoutTOF": 76.04,
    "withP03": 93.58
  },
  "544550": {
    "withTOF": 20.22,
    "withoutTOF": 79.78,
    "withP03": 86.20
  },
  "544551": {
    "withTOF": 20.45,
    "withoutTOF": 79.55,
    "withP03": 85.88
  },
  "544564": {
    "withTOF": 23.24,
    "withoutTOF": 76.76,
    "withP03": 93.36
  },
  "544565": {
    "withTOF": 24.49,
    "withoutTOF": 75.51,
    "withP03": 93.46
  },
  "544567": {
    "withTOF": 26.47,
    "withoutTOF": 73.53,
    "withP03": 93.80
  },
  "544568": {
    "withTOF": 20.02,
    "withoutTOF": 79.98,
    "withP03": 86.32
  },
  "544580": {
    "withTOF": 24.39,
    "withoutTOF": 75.61,
    "withP03": 93.48
  },
  "544582": {
    "withTOF": 26.43,
    "withoutTOF": 73.57,
    "withP03": 93.89
  },
  "544583": {
    "withTOF": 19.58,
    "withoutTOF": 80.42,
    "withP03": 86.00
  },
  "544585": {
    "withTOF": 20.14,
    "withoutTOF": 79.86,
    "withP03": 85.99
  },
  "544614": {
    "withTOF": 19.94,
    "withoutTOF": 80.06,
    "withP03": 93.88
  },
  "544640": {
    "withTOF": 15.72,
    "withoutTOF": 84.28,
    "withP03": 92.64
  },
  "544652": {
    "withTOF": 21.10,
    "withoutTOF": 78.90,
    "withP03": 93.57
  },
  "544653": {
    "withTOF": 24.17,
    "withoutTOF": 75.83,
    "withP03": 93.58
  },
  "544672": {
    "withTOF": 23.97,
    "withoutTOF": 76.03,
    "withP03": 93.55
  },
  "544674": {
    "withTOF": 23.98,
    "withoutTOF": 76.02,
    "withP03": 93.48
  },
  "544692": {
    "withTOF": 19.01,
    "withoutTOF": 80.99,
    "withP03": 86.94
  },
  "544693": {
    "withTOF": 18.20,
    "withoutTOF": 81.80,
    "withP03": 86.74
  },
  "544694": {
    "withTOF": 19.66,
    "withoutTOF": 80.34,
    "withP03": 86.32
  },
  "544696": {
    "withTOF": 19.60,
    "withoutTOF": 80.40,
    "withP03": 86.05
  },
  "544739": {
    "withTOF": 24.59,
    "withoutTOF": 75.41,
    "withP03": 93.81
  },
  "544742": {
    "withTOF": 18.63,
    "withoutTOF": 81.37,
    "withP03": 86.34
  },
  "544754": {
    "withTOF": 18.91,
    "withoutTOF": 81.09,
    "withP03": 86.73
  },
  "544767": {
    "withTOF": 26.70,
    "withoutTOF": 73.30,
    "withP03": 93.69
  },
  "544794": {
    "withTOF": 20.32,
    "withoutTOF": 79.68,
    "withP03": 86.78
  },
  "544795": {
    "withTOF": 21.95,
    "withoutTOF": 78.05,
    "withP03": 85.86
  },
  "544797": {
    "withTOF": 22.37,
    "withoutTOF": 77.63,
    "withP03": 85.91
  },
  "544813": {
    "withTOF": 19.45,
    "withoutTOF": 80.55,
    "withP03": 86.61
  },
  "544868": {
    "withTOF": 20.45,
    "withoutTOF": 79.55,
    "withP03": 86.26
  },
  "544886": {
    "withTOF": 20.32,
    "withoutTOF": 79.68,
    "withP03": 87.16
  },
  "544887": {
    "withTOF": 21.26,
    "withoutTOF": 78.74,
    "withP03": 87.54
  },
  "544896": {
    "withTOF": 21.81,
    "withoutTOF": 78.19,
    "withP03": 86.10
  },
  "544911": {
    "withTOF": 21.82,
    "withoutTOF": 78.18,
    "withP03": 88.39
  },
  "544913": {
    "withTOF": 21.02,
    "withoutTOF": 78.98,
    "withP03": 90.92
  },
  "544914": {
    "withTOF": 20.42,
    "withoutTOF": 79.58,
    "withP03": 86.64
  },
  "544917": {
    "withTOF": 21.94,
    "withoutTOF": 78.06,
    "withP03": 86.08
  },
  "544931": {
    "withTOF": 27.96,
    "withoutTOF": 72.04,
    "withP03": 93.55
  },
  "544947": {
    "withTOF": 25.43,
    "withoutTOF": 74.57,
    "withP03": 94.57
  },
  "544961": {
    "withTOF": 20.86,
    "withoutTOF": 79.14,
    "withP03": 87.20
  },
  "544963": {
    "withTOF": 20.73,
    "withoutTOF": 79.27,
    "withP03": 86.68
  },
  "544964": {
    "withTOF": 18.42,
    "withoutTOF": 81.58,
    "withP03": 86.40
  },
  "544968": {
    "withTOF": 21.05,
    "withoutTOF": 78.95,
    "withP03": 86.16
  },
  "544991": {
    "withTOF": 19.06,
    "withoutTOF": 80.94,
    "withP03": 86.55
  },
  "544992": {
    "withTOF": 21.93,
    "withoutTOF": 78.07,
    "withP03": 86.13
  },
  "545004": {
    "withTOF": 28.65,
    "withoutTOF": 71.35,
    "withP03": 93.51
  },
  "545008": {
    "withTOF": 20.79,
    "withoutTOF": 79.21,
    "withP03": 86.89
  },
  "545009": {
    "withTOF": 21.50,
    "withoutTOF": 78.50,
    "withP03": 86.62
  },
  "545041": {
    "withTOF": 20.48,
    "withoutTOF": 79.52,
    "withP03": 87.40
  },
  "545042": {
    "withTOF": 21.26,
    "withoutTOF": 78.74,
    "withP03": 86.71
  },
  "545044": {
    "withTOF": 21.84,
    "withoutTOF": 78.16,
    "withP03": 86.31
  },
  "545047": {
    "withTOF": 21.92,
    "withoutTOF": 78.08,
    "withP03": 86.14
  },
  "545060": {
    "withTOF": 28.56,
    "withoutTOF": 71.44,
    "withP03": 93.64
  },
  "545062": {
    "withTOF": 20.88,
    "withoutTOF": 79.12,
    "withP03": 86.60
  },
  "545063": {
    "withTOF": 22.13,
    "withoutTOF": 77.87,
    "withP03": 86.21
  },
  "545064": {
    "withTOF": 21.81,
    "withoutTOF": 78.19,
    "withP03": 86.28
  },
  "545066": {
    "withTOF": 4.71,
    "withoutTOF": 95.29,
    "withP03": 78.21
  },
  "545086": {
    "withTOF": 18.42,
    "withoutTOF": 81.58,
    "withP03": 86.97
  },
  "545103": {
    "withTOF": 20.91,
    "withoutTOF": 79.09,
    "withP03": 86.75
  },
  "545117": {
    "withTOF": 21.09,
    "withoutTOF": 78.91,
    "withP03": 89.68
  },
  "545171": {
    "withTOF": 19.97,
    "withoutTOF": 80.03,
    "withP03": 86.41
  },
  "545184": {
    "withTOF": 26.80,
    "withoutTOF": 73.20,
    "withP03": 94.01
  },
  "545185": {
    "withTOF": 21.75,
    "withoutTOF": 78.25,
    "withP03": 86.05
  },
  "545210": {
    "withTOF": 19.25,
    "withoutTOF": 80.75,
    "withP03": 86.52
  },
  "545222": {
    "withTOF": 20.02,
    "withoutTOF": 79.98,
    "withP03": 86.44
  },
  "545223": {
    "withTOF": 21.19,
    "withoutTOF": 78.81,
    "withP03": 86.21
  },
  "545246": {
    "withTOF": 19.91,
    "withoutTOF": 80.09,
    "withP03": 86.95
  },
  "545249": {
    "withTOF": 21.53,
    "withoutTOF": 78.47,
    "withP03": 86.13
  },
  "545262": {
    "withTOF": 19.14,
    "withoutTOF": 80.86,
    "withP03": 86.84
  },
  "545289": {
    "withTOF": 26.86,
    "withoutTOF": 73.14,
    "withP03": 94.09
  },
  "545291": {
    "withTOF": 21.10,
    "withoutTOF": 78.90,
    "withP03": 86.50
  },
  "545294": {
    "withTOF": 21.22,
    "withoutTOF": 78.78,
    "withP03": 86.96
  },
  "545295": {
    "withTOF": 21.97,
    "withoutTOF": 78.03,
    "withP03": 88.11
  },
  "545296": {
    "withTOF": 21.62,
    "withoutTOF": 78.38,
    "withP03": 87.89
  },
  "545311": {
    "withTOF": 16.78,
    "withoutTOF": 83.22,
    "withP03": 89.76
  },
  "545312": {
    "withTOF": 21.42,
    "withoutTOF": 78.58,
    "withP03": 86.06
  },
  "545332": {
    "withTOF": 20.11,
    "withoutTOF": 79.89,
    "withP03": 87.68
  },
  "545345": {
    "withTOF": 19.05,
    "withoutTOF": 80.95,
    "withP03": 86.49
  },
  "545367": {
    "withTOF": 15.57,
    "withoutTOF": 84.43,
    "withP03": 86.49
  }
}



# Convert JSON to DataFrame
df = pd.DataFrame.from_dict(json_data, orient='index')
df.index.name = 'RunNumber'
df.reset_index(inplace=True)

# Sort by RunNumber
df['RunNumber'] = df['RunNumber'].astype(str)
df.sort_values(by='RunNumber', inplace=True)

# Estimate: % of p ≥ 0.3 tracks that have TOF
df['est_p03_with_TOF'] = (df['withTOF'] / df['withP03']) * 100

# Plot
plt.figure(figsize=(18, 7))
x = range(len(df))
width = 0.3

# Bar plots
plt.bar(x, df['withTOF'], width=width, label='% With TOF')
plt.bar([i + width for i in x], df['est_p03_with_TOF'], width=width, label='% Est_p03_with_TOF')



# Labels & legend
plt.xlabel('Run Number')
plt.ylabel('Percentage (%)')
plt.title('TOF and Momentum Distribution with Estimated Conditional TOF Rate')
plt.xticks([i + width/2 for i in x], df['RunNumber'], rotation=90, fontsize=7)
plt.ylim(0, 100)
plt.legend()
plt.tight_layout()
plt.show()