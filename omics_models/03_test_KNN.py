import pandas as pd
import numpy as np
import os, glob
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
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
    # gof/log variants
    gof = list(pd.read_csv('gof_variants.csv')['gof_variants'])
    gof = list(set(df1.columns) & set([i.replace('_p.', '_') for i in gof]))
    lof = list(pd.read_csv('lof_variants.csv')['lof_variants'])
    lof = list(set(df1.columns) & set([i.replace('_p.', '_') for i in lof]))
    y = [1 if i in gof else 0 for i in df1.columns if i in gof or i in lof]
    return df1, y

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

if __name__ == '__main__':
    train_base_path = './train/'
    test_base_path = './test/'

    np.random.seed(42)

    X_train_raw, y_train = load_emb(train_base_path)
    # Fit PCA and Scaler to training set
    rs, pca, X_train = dimension_reduce_fit(X_train_raw)

    X_test_raw, y_test = load_emb(test_base_path)

    print("Transforming test dataframe using Scaler from train set")
    X_test_scaled = rs.transform(X_test_raw.T)
    print("Applying PCA from train set to test set independently...")
    X_test = pca.transform(np.array(X_test_scaled))

    clf = KNeighborsClassifier(n_neighbors=3)
    clf.fit(X_train, y_train)
    y_prob = clf.predict_proba(X_test)[:, 1]

    # Categorize and print predictions with their variants and binary labels
    categorized_predictions = [(X_test_raw.columns[i], y_test[i], y_prob[i]) for i in range(len(y_test))]

    print("Predicted scores categorized by binary label:")
    for variant, label, score in categorized_predictions:
        if label == score:
            print(f"Variant: {variant}, Binary Label: {label}, Predicted Score: {score}")
    
    df_categorized = pd.DataFrame(categorized_predictions, columns=['Feature', 'Binary_Label', 'Predicted_Score'])
    df_categorized.to_csv("./output/knn_test_results.csv")

    # Create histogram plot for the predicted scores across the whole dataset
    plt.figure(figsize=(10, 6))
    sns.histplot(df_categorized['Predicted_Score'], bins=30, kde=True)
    plt.title('Distribution of Prediction Scores')
    plt.xlabel('Prediction Score')
    plt.ylabel('Frequency')
    plt.savefig("./visuals/knn_hist.pdf", format="pdf", transparent=True, bbox_inches="tight", pad_inches=0.25)
