# random_forest_crc_analysis.py
# Random Forest classification analysis for CRC datasets

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.metrics import (
    roc_auc_score, accuracy_score, recall_score, precision_score,
    f1_score, confusion_matrix, roc_curve
)
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV
from sklearn.utils import resample
from scipy.stats import randint

countries = ["austria", "china", "india", "france", "italy", "usa", "japan"]
taxonomy_levels = ["genus", "species", "strain"]
base_root_prefix = "maaslin_result_"
result_root_prefix = "rf_"
fml_conditions = ["before_fml_adj", "after_fml_adj"]
n_splits = 100

param_dist = {
    'n_estimators': randint(100, 500),
    'max_depth': randint(3, 20),
    'min_samples_split': randint(2, 11),
    'min_samples_leaf': randint(1, 5),
    'max_features': ['sqrt', 'log2', None],
    'bootstrap': [True, False],
    'class_weight': [None, 'balanced']
}

def auc_ci(y_true, y_probs, n_bootstrap=1000, ci=95):
    scores = []
    for _ in range(n_bootstrap):
        idx = resample(range(len(y_true)), replace=True)
        yt, yp = y_true.iloc[idx], y_probs[idx]
        if len(np.unique(yt)) < 2:
            continue
        scores.append(roc_auc_score(yt, yp))
    if not scores:
        return np.nan, np.nan
    lower = np.percentile(scores, (100 - ci) / 2)
    upper = np.percentile(scores, ci + (100 - ci) / 2)
    return lower, upper

def save_results_level(result_root_level, country, metrics):
    cdir = os.path.join(result_root_level, country)
    os.makedirs(cdir, exist_ok=True)
    for condition in fml_conditions:
        ldir = os.path.join(cdir, condition)
        os.makedirs(ldir, exist_ok=True)
        if metrics[condition]:
            pd.DataFrame(metrics[condition]).to_csv(os.path.join(ldir, "rf_results.csv"), index=False)

def main():
    for level in taxonomy_levels:
        base_root_level = base_root_prefix + level
        result_root_level = result_root_prefix + level
        os.makedirs(result_root_level, exist_ok=True)

        for country in countries:
            metrics_res = {condition: [] for condition in fml_conditions}

            for condition in fml_conditions:
                base_path = os.path.join(base_root_level, f"{country}_rf_feature", condition)
                if not os.path.exists(base_path):
                    continue

                for split in range(1, n_splits + 1):
                    split_path = os.path.join(base_path, f"split_{split}")
                    train_file = os.path.join(split_path, "train_filtered.csv")
                    test_file = os.path.join(split_path, "test_filtered.csv")
                    if not (os.path.exists(train_file) and os.path.exists(test_file)):
                        continue

                    train_df = pd.read_csv(train_file)
                    test_df = pd.read_csv(test_file)
                    X_train = train_df.drop(columns=["Disease", "Sample_ID"])
                    y_train = train_df["Disease"].map({"Control": 0, "CRC": 1})
                    X_test = test_df.drop(columns=["Disease", "Sample_ID"])
                    y_test = test_df["Disease"].map({"Control": 0, "CRC": 1})

                    if X_train.shape[1] > 1:
                        rfecv = RFECV(
                            RandomForestClassifier(random_state=42),
                            step=1,
                            cv=StratifiedKFold(5, shuffle=True, random_state=42),
                            scoring='roc_auc',
                            n_jobs=-1
                        )
                        rfecv.fit(X_train, y_train)
                        sel_features = X_train.columns[rfecv.support_]
                        X_train_sel = X_train[sel_features]
                        X_test_sel = X_test[sel_features]
                    else:
                        sel_features = X_train.columns
                        X_train_sel = X_train
                        X_test_sel = X_test

                    search = RandomizedSearchCV(
                        RandomForestClassifier(random_state=42),
                        param_dist,
                        n_iter=50,
                        cv=5,
                        scoring='roc_auc',
                        n_jobs=-1,
                        random_state=42
                    )
                    search.fit(X_train_sel, y_train)
                    clf = search.best_estimator_

                    y_prob = clf.predict_proba(X_test_sel)[:, 1]
                    fpr, tpr, thr = roc_curve(y_test, y_prob)
                    best_thr = thr[np.argmax(tpr - fpr)]
                    y_pred = (y_prob >= best_thr).astype(int)

                    auc_val = roc_auc_score(y_test, y_prob)
                    auc_low, auc_high = auc_ci(y_test, y_prob)
                    cm = confusion_matrix(y_test, y_pred)

                    metrics = {
                        "Country": country,
                        "fml_condition": condition,
                        "Split": split,
                        "AUC": auc_val,
                        "AUC_CI_Lower": auc_low,
                        "AUC_CI_Upper": auc_high,
                        "Accuracy": accuracy_score(y_test, y_pred),
                        "Sensitivity": recall_score(y_test, y_pred),
                        "Specificity": cm[0, 0] / (cm[0, 0] + cm[0, 1]) if (cm[0, 0] + cm[0, 1]) > 0 else np.nan,
                        "PPV": precision_score(y_test, y_pred),
                        "NPV": cm[0, 0] / (cm[0, 0] + cm[1, 0]) if (cm[0, 0] + cm[1, 0]) > 0 else np.nan,
                        "F1_Score": f1_score(y_test, y_pred),
                        "Youden_Index": recall_score(y_test, y_pred) + (cm[0, 0] / (cm[0, 0] + cm[0, 1]) if (cm[0, 0] + cm[0, 1]) > 0 else 0) - 1,
                        "Best_Params": str(search.best_params_),
                        "Selected_Features_Count": len(sel_features)
                    }

                    metrics_res[condition].append(metrics)

            save_results_level(result_root_level, country, metrics_res)

if __name__ == "__main__":
    main()
