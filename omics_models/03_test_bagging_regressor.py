from sklearn.ensemble import BaggingRegressor
from sklearn.tree import DecisionTreeRegressor
import pandas as pd
import numpy as np
import os, glob
from sklearn.preprocessing import RobustScaler, MinMaxScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

def load_emb(base_path: str):
    """Loads the embeddings using numpy and flattens them to one-dimension

    Args:
        base_path (str): the base path of where the embedding folder is listed.

    Returns:
        embs (pd.DataFrame): Dataframe of flattened embeddings from AMS
        gof (list): names of GOF variants in GENE_REFAAPOSALTAA format
        lof (list): names of LOF variants in GENE_REFAAPOSALTAA format
    """
    # change directory to current working directory
    os.chdir(base_path)
    fls = glob.glob('embs/*.npy')
    dct = {}
    for file in fls:
        fi = np.load(file)
        # load and flatten each embedding
        dct[file.split('embs/embedding_')[1].split('.')[0]] = np.ndarray.flatten(fi)
    df1 = pd.DataFrame(dct)
    return df1

def dimension_reduce_fit(train_df):
    """Apply dimensionality reduction using PCA and scale values

    Args:
        df1 (pd.DataFrame): Dataframe of training set embeddings
    Returns:
        scaler (RobustScaler): scaler object
        pca_obj (PCA()): PCA object
        pca_train (pd.DataFrame): Dataframe transformed by PCA
    """
    # pca of flattened embeddings
    pca_obj = PCA()
    # Apply RobustScaler
    scaler = RobustScaler()
    print("Fitting Robust Value Scaler")
    df1_scaled = scaler.fit_transform(train_df.T)
    print("Fitting PCA SVD to train set to scale test")
    pca_train = pca_obj.fit_transform(np.array(df1_scaled))
    return scaler, pca_obj, pca_train

def align_y_values(y_csv_path: str, df: pd.DataFrame):
    """Auxiliary function to map y values to appropriate indices.

    Args:
        y_csv_path (str): variant weight csv path
        df (pd.DataFrame): data columns from X to align indices

    Returns:
        ordered_y_values (list): y values for the list of genes in PCA training set
    """
    # Load y-values from the CSV file
    y_values_df = pd.read_csv(y_csv_path)
    y_values_dict = dict(zip(y_values_df['Variant'], y_values_df['Average_Matches']))
    
    # Reorder y-values based on the DataFrame index
    ordered_y_values = [y_values_dict[idx] for idx in df.columns]
    
    return ordered_y_values

if __name__ == '__main__':
    train_base_path = './train/'
    test_base_path = './test/'

    X_train_raw = load_emb(train_base_path)

    # Fit PCA and Scaler to training set
    rs, pca, X_train = dimension_reduce_fit(X_train_raw)

    # Use MOFA weights for prediction
    y_train = align_y_values("/wistar/auslander/Bryant/projects/pdx/processed/visuals/variant_weights.csv", X_train_raw)

    # Scale y_train to range [-1, 1]
    y_scaler = MinMaxScaler(feature_range=(-1, 1))
    y_train_scaled = y_scaler.fit_transform(np.array(y_train).reshape(-1, 1)).flatten()
    X_test_raw = load_emb(test_base_path)

    # Extracting variant names from the test data
    variant_test = X_test_raw.columns.tolist()

    print("Transforming test dataframe using Scaler from train set")
    X_test_scaled = rs.transform(X_test_raw.T)
    print("Applying PCA from train set to test set independently...")
    X_test = pca.transform(np.array(X_test_scaled))

    base_regressor = DecisionTreeRegressor(random_state=42)
    bagging_regressor = BaggingRegressor(estimator=base_regressor, n_estimators=10, random_state=42)

    # Train the bagging regressor
    bagging_regressor.fit(X_train, y_train_scaled)

    # Make predictions with the bagging regressor
    y_pred_regression = bagging_regressor.predict(X_test)

    # Create a DataFrame with the predictions and variant names
    predictions_df = pd.DataFrame({
        'Variant': variant_test,
        'Score': y_pred_regression
    })
    #predictions_df.to_csv('./output/test_variant_weights.csv', index=False)
    # Plot the predicted scores
    plt.figure(figsize=(10, 6))
    sns.histplot(y_pred_regression, kde=True, bins=20)
    plt.title('Distribution of Prediction Scores')
    plt.xlabel('Predicted Score')
    plt.ylabel('Frequency')

    # Save the plot as a file
    plt.savefig("./visuals/bagging_hist.pdf", format="pdf", transparent=True, bbox_inches="tight", pad_inches = 0.25)