import numpy as np
import csv
import pandas as pd

if __name__ == "__main__":
    # csv_file = open("../data/STW-M-Brain-Stereo-seq-1/counts.copy.csv")
    # reader = csv.reader(csv_file)
    X = pd.read_csv("../data/STW-M-Brain-Stereo-seq-1/counts.copy.csv")
    X = X.T
    X.to_tsv("../data/STW-M-Brain-Stereo-seq-1/counts.t.txt", header=False, index=True)
