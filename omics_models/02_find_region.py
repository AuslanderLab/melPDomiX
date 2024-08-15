import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns


def extract_first_slices(matrices):
    """Extracts the first slice from each embedding matrix

    Args:
        matrices (list): a list of numpy array matrices

    Returns:
        np.ndarray: the first slices along the third dimension of the numpy array
    """
    first_slices = []
    for matrix in matrices:
        first_slices.append(matrix[:, :, 0].flatten())
    return np.array(first_slices)


# Concatenate the first slices for PCA
def concatenate_first_slices(class1_slices, class2_slices):
    """Concatenate the GOF and LOF slice in order to perform PCA

    Args:
        class1_slices (list): class 1 (GOF) slices
        class2_slices (list): class 2 (LOF) slices

    Returns:
        np.ndarray: combined GOF and LOF slices
    """
    combined_data = []
    for i in range(len(class1_slices)):
        combined_data.append(np.concatenate((class1_slices[i], class2_slices[i]), axis=0))
    return np.array(combined_data)

# Apply PCA to the combined first slices
def perform_pca(combined_data):
    """Performs PCA on combined data to get summmed contributions

    Args:
        combined_data (np.ndarray): Combined data matrix for each of the slices

    Returns:
        summed_contributions (np.ndarray): Summed contributions from the components matrix where matrix is size (n_components, n_features)
    """
    pca = PCA()
    pca.fit(combined_data)
    loadings = pca.components_
    summed_contributions = np.sum(loadings, axis=0)
    return summed_contributions

def average_difference_in_region(matrix, x, y, region_size):
    """Slice numpy array to extract region of interest and find the mean of the matrix

    Args:
        matrix (np.ndarray): Matrix to extract
        x (int): x-coordinate of region
        y (int): y-coordinate of region
        region_size (int): the size in nxn value e.g. 5x5 = 25

    Returns:
        np.ndarray: Average of matrix array elements for the region specified
    """
    return np.mean(matrix[x:x+region_size, y:y+region_size])

def plot_original_regions(matrix_list, x, y, region_size, titles, cmap_color):
    """Plots the regions of interest

    Args:
        matrix_list (list): list of numpy arrays to plot from
        x (int): x-coordinate
        y (int): y-coordinate
        region_size (int): region size
        titles (list): titles for the individual matrices
        cmap_color (str): color palette to use
    """
    fig, axs = plt.subplots(2, 3, figsize=(24, 16))

    for ax, matrix, title in zip(axs.flatten(), matrix_list, titles):
        half_size = region_size // 2
        # using the x and y coordinates as the center, extract a region
        region = matrix[max(0, x-half_size):min(256, x+half_size), max(0, y-half_size):min(256, y+half_size)]
        
        sns.heatmap(region, cmap=cmap_color, vmin=-4000, vmax=6000, ax=ax, cbar=False)
        ax.set_xticks([])
        ax.set_yticks([])
    
    # Add a single colorbar with adjusted aspect and outline
    cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    cbar = fig.colorbar(axs[0, 0].collections[0], cbar_ax, aspect=10)
    cbar.ax.tick_params(labelsize=10)
    cbar.outline.set_linewidth(0.2)
    fig.savefig(f"./visuals/heat_gof_lof.pdf", format="pdf", transparent=True, bbox_inches="tight")
    
    # Add a single colorbar
    cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    cbar = fig.colorbar(axs[0].collections[0], cbar_ax, aspect=10)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_ticks([i*2000 for i in range(0, 4)])
    cbar.outline.set_linewidth(0.2)
    fig.savefig("./visuals/heatmap_differences_subplots.pdf", format="pdf", transparent=True, bbox_inches="tight", pad_inches = 0.25)


if __name__ == "__main__":
    # Example matrices
    rac1 = np.load("./embs/embedding_*.npy")
    braf = np.load("./embs/embedding_*.npy")
    nras = np.load("./embs/embedding_*.npy")
    alk = np.load("./embs/embedding_*.npy")
    isx = np.load("./embs/embedding_*.npy")
    pik3r2 = np.load("./embs/embedding_*.npy")

    class1_matrices = [nras, rac1, braf]
    class2_matrices = [isx, pik3r2, alk]

    str1 = ["Variant1", "Variant2", "Variant3"]
    str2 = ["Variant4", "Variant5", "Variant6"]

    class1_first_slices = extract_first_slices(class1_matrices)
    class2_first_slices = extract_first_slices(class2_matrices)

    combined_first_slices = concatenate_first_slices(class1_first_slices, class2_first_slices)

    summed_contributions = perform_pca(combined_first_slices)

    # Reshape summed contributions to 256x256 matrix
    contributions_matrix = summed_contributions[:256*256].reshape(256, 256)

    # Analyze specific square regions within the first slice
    region_size = 25

    max_diff = 0
    max_coords = (0, 0)
    for i in range(256 - region_size + 1):
        for j in range(256 - region_size + 1):
            current_diff = average_difference_in_region(contributions_matrix, i, j, region_size)
            if current_diff > max_diff:
                max_diff = current_diff
                max_coords = (i, j)

    x = 109
    y = 114
    print(f'The {region_size}x{region_size} region with the highest difference in the first slice is located at coordinates ({max_coords[0]}, {max_coords[1]}) with an average difference of {max_diff}')

    # Titles for the plots
    titles = ["Variant1", "Variant2", "Variant3",
            "Variant4", "Variant5", "Variant6"]

    # List of matrices for plotting
    matrix_list = [class1_matrices[0][:, :, 0], class1_matrices[1][:, :, 0], class1_matrices[2][:, :, 0],
                class2_matrices[0][:, :, 0], class2_matrices[1][:, :, 0], class2_matrices[2][:, :, 0]]

    # put any range of color palettes needed
    colors = ['Oranges']
    # Plot the original regions for each matrix
    for color in colors:
        plot_original_regions(matrix_list, x, y, region_size, titles, color)