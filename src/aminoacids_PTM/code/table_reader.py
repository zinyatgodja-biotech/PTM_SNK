import pandas as pd
df = pd.read_excel('../dpptm_human_ptm_amino_acids.xlsx',
                 encoding='utf-8',
                 sep=';')  # или sep=','
df 