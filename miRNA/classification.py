import pandas as pd
import numpy as np
import pickle
from sklearn.ensemble import IsolationForest
import matplotlib.pyplot as plt
from sklearn import linear_model, metrics, feature_selection


def create_classification_data():
    """
    Using statistical methods, outliers and insignificant features are removed from the data.
    :return: Separated training data and labels and test data and labels
    """
    train_data = pd.read_pickle("train_data.pkl")
    train_labels = pickle.load(open('train_labels.pkl', 'rb'))
    test_data = pd.read_pickle("test_data.pkl")
    test_labels = pickle.load(open('test_labels.pkl', 'rb'))

    # Outliers detection and removal:
    clf = IsolationForest(max_samples=302, contamination=.1)
    clf.fit(train_data)
    y_pred_train = clf.predict(train_data)
    outliers_indexes = list(np.where(y_pred_train == -1)[0])
    for index in reversed(outliers_indexes):
        del train_labels[index]
    train_data_cleaned = train_data[np.where(y_pred_train == 1, True, False)]

    # Feature selection:
    best_features = feature_selection.SelectKBest(score_func=feature_selection.chi2, k='all')
    fit = best_features.fit(train_data_cleaned, train_labels)
    train_data = pd.DataFrame(fit.transform(train_data_cleaned))
    test_data = pd.DataFrame(fit.transform(test_data))

    return train_data, train_labels, test_data, test_labels


def logistic_regression(X_train, y_train, X_test, y_test):
    """
    Creates a logistic regression classifier for the miRNA data.
    :param X_train: training data - matrix that holds the amount of each miRNA for each patient
    :param y_train: training tags - 1 if the patient has Alzheimer's disease, 0 otherwise
    :param X_test: test data - matrix that holds the amount of each miRNA for each patient
    :param y_test: test tags - 1 if the patient has Alzheimer's disease, 0 otherwise
    :return: The success rate is printed, a file that holds the weights of each miRNA for the classification is
    and a confusion matrix are created
    """
    lr = linear_model.LogisticRegression(solver='saga')
    lr.fit(X_train, y_train)
    score = lr.score(X_test, y_test)
    print(score)
    lr_weights = lr.coef_
    weights_df = pd.DataFrame(lr_weights)
    weights_df.to_csv("lr_weights.csv")
    labels = [0, 1]
    metrics.plot_confusion_matrix(lr, X_test, y_test, labels=None, normalize='true', display_labels=labels,
                                  values_format='.2f')
    plt.title("Logistic Regression Confusion Matrix")
    plt.show()

