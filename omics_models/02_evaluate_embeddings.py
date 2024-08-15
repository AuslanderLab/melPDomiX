import pandas as pd
import numpy as np
import glob, os
from sklearn.preprocessing import RobustScaler   
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns

# Plot params with typesetting recommendations
plt.rcParams["font.family"] = "Arial"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.facecolor'] = 'none'


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
    fls = glob.glob('./*.npy')
    dct = {}
    for file in fls:
        fi = np.load(file)
        # load and flatten each embedding
        dct[file.split('./embs/embedding_')[1].split('.')[0]] = np.ndarray.flatten(fi)
    embs = pd.DataFrame(dct)
    # Extract GOF and LOF variant names for future use.s
    gof = list(pd.read_csv('gof_variants.csv')['gof_variants'])
    gof = list(set(embs.columns)&set([i.replace('_p.','_') for i in gof]))
    lof = list(pd.read_csv('lof_variants.csv')['lof_variants'])
    lof = list(set(embs.columns)&set([i.replace('_p.','_') for i in lof]))
    return embs, gof, lof

def dimension_reduce(df1, gof, lof):
    """Apply dimensionality reduction using PCA and scale values

    Args:
        df1 (pd.DataFrame): Dataframe of flattened embeddings from AMS
        gof (list): names of GOF variants in GENE_REFAAPOSALTAA format
        lof (list): names of LOF variants in GENE_REFAAPOSALTAA format

    Returns:
        Xt (np.ndarray): Transformed values matrix
    """
    # Applying PCA and Robust Scaler to values.
    pca = PCA()
    scaler = RobustScaler()
    print("Transforming to fit robust values scaler")
    df1_scaled = scaler.fit_transform(df1.T)
    Xt = pca.fit_transform(np.array(df1_scaled))
    y=['r' if i in gof else 'b' if i in lof else 'c' for i in df1.columns]
    return Xt, y

if __name__ == "__main__":
    df1, gof, lof = load_emb('./embs')
    Xt, y = dimension_reduce(df1, gof, lof)
    Xt2 = Xt[:,0:4]
    lb=[1 if i in gof else 0 for i in df1.columns]
    labels = np.array(lb)
    indices = [i for i in range(len(Xt2))]
    Xt2 = Xt[:,0:12]

    edge_weights = {}

    # Build 500 KNN models with 75/25 split to predict GOF/LOF based on first 12 PCs with random seed i
    accs = []
    aucs = []
    ITERS = 500
    intrain = np.zeros((Xt2.shape[0], ITERS))

    for i in range(500):
        X_train, X_test, y_train, y_test, ind_train,ind_test = train_test_split(Xt2, labels, indices, test_size = 0.25, random_state = i,stratify=y)
        intrain[ind_train,i] = 1
        clf = KNeighborsClassifier(n_neighbors=3)
        clf.fit(X_train, y_train)

        predsc=clf.predict_proba(X_test)[:,1]
        predlb = clf.predict(X_test)
        accuracy = sum(np.array(predlb) == y_test) / len(y_test)
        accs.append(accuracy)
        aucs.append(roc_auc_score(y_test, predsc))
        print(sum(aucs)/len(aucs))

        # Track which variants are in the training and test sets
        train_variants = df1.columns[ind_train]
        test_variants = df1.columns[ind_test]   
        for train_var in train_variants:
            for test_var in test_variants:
                if (train_var, test_var) not in edge_weights:
                    edge_weights[(train_var, test_var)] = []
                edge_weights[(train_var, test_var)].append(accuracy)

    # Calculate the average accuracy for each edge
    average_edge_weights = {k: np.mean(v) for k, v in edge_weights.items()}

    # Create a DataFrame for the edges and weights
    edges_df = pd.DataFrame([
        {'Source': k[0], 'Target': k[1], 'Weight': v}
        for k, v in average_edge_weights.items()
    ])

    # Save the DataFrame as a CSV file
    edges_df.to_csv('./output/variant_relationships.csv', index=False)

    print("CSV file 'variant_relationships.csv' created successfully.")

    fig, ax1 = plt.subplots()

    left, bottom, width, height = [0.15, 0.65, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])

    # Main
    sns.histplot(aucs, kde=True, ax=ax1)
    ax1.set_xlabel('ROC AUC')
    ax1.set_ylabel('Count of Cross Validation Folds')

    # Inset boxplot
    sns.boxplot(x=aucs, ax=ax2)
    ax2.set_xlabel('AUC')
    ax2.set_ylabel('')
    ax2.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_xticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'])
    fig.savefig("./visuals/roc_auc_5c.pdf", format="pdf", transparent=True, bbox_inches="tight")
