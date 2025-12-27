# Sample Data Download Instructions

## Required Data File

To run this program, you need to download the sample dataset and place it in this folder.

### Download Link

**Google Drive:** [NDT-CFPP Sample Data](https://drive.google.com/drive/folders/1YCpe9VeXKPuEkiJ6r__Pnkd2g6wBtUMZ)

### Instructions

1. Click the Google Drive link above
2. Download the file: `L16P1S10_10MHz.mat` (~2.4 GB)
3. Place the downloaded file in this `Sample Data/` folder
4. Run `Main.m` - the program will automatically detect the data file

### File Structure After Download

```
A_ Current Rev/
├── Main.m
├── Sample Data/
│   ├── README.md (this file)
│   └── L16P1S10_10MHz.mat  <-- Place downloaded file here
├── Data and Computation Scripts/
└── ...
```

### Available Datasets

| Dataset | Filename | Description |
|---------|----------|-------------|
| 4 (Default) | `L16P1S10_10MHz.mat` | 16° Defect, Sample 10 |

Additional datasets may be available upon request.

### Troubleshooting

If the program cannot find the data file:
1. Ensure the file is named exactly `L16P1S10_10MHz.mat`
2. Ensure the file is placed directly in the `Sample Data/` folder
3. Check that `DATASET = 4` is set in `Main.m` (line 22)

### Alternative Data Locations

The program also checks these locations (in order):
1. `Sample Data/` folder (recommended)
2. Environment variable `NDE_DATA_ROOT`
3. `../../Data/LaminateBeam/L0/` (relative to code)
4. `~/Documents/X_ University/Uni Research/NDT/NDE_SDSU-USD/Data/LaminateBeam/L0/`

