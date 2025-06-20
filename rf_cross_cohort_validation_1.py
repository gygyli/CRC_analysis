# rf_cross_cohort_validation_1.py
# Author: Guangyi Li
# Date: 2025-06-17
# Description: Perform Random Forest cross-cohort validation with feature selection and hyperparameter tuning on microbial taxonomy data; evaluate AUC on test sets and save models and features.

import os
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV
from sklearn.metrics import roc_auc_score
from scipy.stats import randint
import joblib

countries = ["austria", "china", "india", "france", "italy", "usa", "japan"]
fml_conditions = ["before_fml_adj", "after_fml_adj"]
taxonomy_levels = ["genus", "species", "strain"]

param_dist = {
    'n_estimators': randint(100, 500),
    'max_depth': randint(3, 20),
    'min_samples_split': randint(2, 11),
    'min_samples_leaf': randint(1, 5),
    'max_features': ['sqrt', 'log2', None],
    'bootstrap': [True, False],
    'class_weight': [None, 'balanced']
}

for taxonomy_level in taxonomy_levels:
    for source_country in countries:
        for fml_condition in fml_conditions:
            base_dir = f"./maaslin_result_{taxonomy_level}/{source_country}_rf_feature/{fml_condition}"
            if not os.path.exists(base_dir):
                continue

            for split in os.listdir(base_dir):
                if not split.startswith("split_"):
                    continue
                split_path = os.path.join(base_dir, split)
                split_num = split.split("_")[1]

                train_path = os.path.join(split_path, f"train_filtered.csv")
                test_path  = os.path.join(split_path, f"test_filtered.csv")
                if not os.path.exists(train_path) or not os.path.exists(test_path):
                    continue

                train = pd.read_csv(train_path)
                test = pd.read_csv(test_path)

                X_train = train.drop(columns=["Disease", "Sample_ID"])
                y_train = train["Disease"].map({'Control': 0, 'CRC': 1})
                X_test = test.drop(columns=["Disease", "Sample_ID"])
                y_test = test["Disease"].map({'Control': 0, 'CRC': 1})

                if X_train.shape[1] > 1:
                    rfecv = RFECV(RandomForestClassifier(random_state=42), step=1, cv=5, scoring='roc_auc')
                    rfecv.fit(X_train, y_train)
                    selected = X_train.columns[rfecv.support_]
                    X_train = X_train[selected]
                    X_test = X_test[selected]
                else:
                    selected = X_train.columns

                clf = RandomizedSearchCV(
                    RandomForestClassifier(random_state=42),
                    param_distributions=param_dist,
                    n_iter=50,
                    scoring='roc_auc',
                    cv=5,
                    n_jobs=-1,
                    random_state=42
                )
                clf.fit(X_train, y_train)
                auc = roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1])
                print(f"{taxonomy_level}-{source_country}-{fml_condition}-{split}: AUC = {auc:.3f}")

                save_dir = f"./cross_cohort_validation/models/{taxonomy_level}/{source_country}/{fml_condition}"
                os.makedirs(save_dir, exist_ok=True)
                joblib.dump(clf.best_estimator_, os.path.join(save_dir, f"{split}_model.pkl"))
                with open(os.path.join(save_dir, f"{split}_features.txt"), "w") as f:
                    f.writelines(f"{ftr}\n" for ftr in selected)
