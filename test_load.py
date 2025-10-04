import pandas as pd

# Load Excel file
file_path = "data/data.xlsx"  # bundled data file
xls = pd.ExcelFile(file_path)

# Print sheet names
print("Sheets:", xls.sheet_names)

# Load first sheet
df = pd.read_excel(file_path, sheet_name=xls.sheet_names[0])
print(df.head())
