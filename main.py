import re

import bct
import matplotlib
import nibabel as nib
import numpy as np

matplotlib.use("QtAgg")

from matplotlib import pyplot as plt
from corshrink import corshrink


def hex_to_rgb(h: str):
    h = h.replace("#", "")
    h = h.replace('"', "")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))


def Main():
    path = "./Task1Data/tractor/diffusion/dti_FA.nii.gz"
    obj = nib.load(path)
    data = obj.get_fdata()
    print(data.shape)

    data[np.isnan(data)] = 0

    plt.figure(figsize=(30, 30))
    plt.suptitle("FA mapped slices at different thresholds")

    for i in range(8):
        for j in range(3):
            z = 14 * (j + 1)
            threshold = i * 0.1 + 0.1

            d = data[:, :, z].copy()
            d[d < threshold] = 0
            ax = plt.subplot(3, 8, j * 8 + i + 1)
            if i % 8 == 0:
                ax.set_ylabel(f"{z = }", rotation=90, size="large")

            if j == 0:
                ax.set_title(f"t={ threshold: .1f}")

            ax.imshow(np.flipud(d.T), cmap="grey")

    plt.show()


def GetAttributes(conn_matrix):
    den, k, n = bct.density_und(conn_matrix)
    _, path_matrix = bct.breadthdist(conn_matrix)
    msp, eff, _, _, _ = bct.charpath(path_matrix, include_infinite=False)
    clustering_coeff = np.mean(bct.clustering_coef_bu(conn_matrix))

    return np.asarray([den, msp, eff, clustering_coeff])


def BCT():

    for i in range(1, 9):
        conn_matrix = np.genfromtxt(
            f"./Task1Data/connectomes/FA.{i}_graph.csv", delimiter=","
        )
        comps, comp_sizes = bct.get_components(conn_matrix)

        final_att = np.zeros((4,))

        total_size = 0
        for idx, c in enumerate(comp_sizes):
            if c < 5:
                continue

            total_size += c
            new_matrix = conn_matrix[comps == idx + 1][:, comps == idx + 1]
            att = GetAttributes(new_matrix)
            final_att += att * c

        print(final_att / total_size)
        print(GetAttributes(conn_matrix))

        print("\n")


def parse_lut(file_path):
    data = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            # Skip comments and headers
            if not line or line.startswith("#") or line.startswith("index"):
                continue

            # Split line into fields using regex to handle multiple spaces
            parts = re.split(r"\s{2,}", line)

            if len(parts) == 7:
                index, label, nativelabel, colour, type_, hemisphere, lobe = parts
            else:
                continue  # Skip malformed lines

            # Convert index to int, handle NA values
            index = int(index)
            hemisphere = hemisphere if hemisphere != "NA" else None
            lobe = lobe if lobe != "NA" else None

            data.append([index, label, nativelabel, colour, type_, hemisphere, lobe])

    # Convert to DataFrame for easy manipulation
    return data


def WeightedFunctional():
    path = "./Task1Data/tractor/functional/data.nii.gz"
    obj = nib.load(path)
    data = obj.get_fdata()

    lut = parse_lut("./Task1Data/tractor/functional/parcellation.lut")
    parcellation_attrib = {l[0]: l[1:] for l in lut}
    parcellation = nib.load(
        "./Task1Data/tractor/functional/parcellation.nii.gz"
    ).get_fdata()

    print(data.shape)

    data[np.isnan(data)] = 0

    plt.figure(figsize=(30, 30))
    plt.suptitle("FA mapped slices at different thresholds")

    s = 17
    for i in range(4):
        d = data[:, :, s, i].copy()

        color_coded = np.zeros((d.shape[0], d.shape[1], 3))
        for x in range(d.shape[0]):
            for y in range(d.shape[1]):
                parcel = int(parcellation[x, y, s])
                if parcel in parcellation_attrib:
                    color_coded[x, y] = hex_to_rgb(parcellation_attrib[parcel][2])

        plt.subplot(2, 4, i * 2 + 1), plt.imshow(d, cmap="grey")
        plt.subplot(2, 4, i * 2 + 2), plt.imshow(color_coded)

    plt.show()


if __name__ == "__main__":
    # Main()
    # BCT()
    WeightedFunctional()
