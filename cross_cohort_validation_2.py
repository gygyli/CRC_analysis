import os
import pandas as pd
import numpy as np
import joblib
from sklearn.metrics import (
    roc_auc_score, accuracy_score, precision_score, recall_score,
    f1_score, confusion_matrix, roc_curve
)

countries = ["austria", "china", "india", "france", "italy", "usa", "japan"]
fml_conditions = ["before_fml_adj", "after_fml_adj"]
taxonomy_levels = ["genus", "species", "strain"]

all_results = []

for taxonomy_level in taxonomy_levels:
    for source_country in countries:
        for fml_condition in fml_conditions:
            model_dir = f"./cross_cohort_validation/models/{taxonomy_level}/{source_country}/{fml_condition}"
            if not os.path.exists(model_dir):
                continue

            target_countries = [c for c in countries if c != source_country]

            for target_country in target_countries:
                data_path = f"./cross_cohort_validation/target_country_data/{taxonomy_level}/{target_country}_data_with_metadata.tsv"
                if not os.path.exists(data_path):
                    print(f"Missing data: {data_path}")
                    continue

                df = pd.read_csv(data_path, sep="\t")
                X_ext = df.drop(columns=["Disease", "Sample_ID"], errors='ignore')
                y_ext = df["Disease"].map({'Control': 0, 'CRC': 1})

                for file in os.listdir(model_dir):
                    if not file.endswith("_model.pkl"):
                        continue

                    model_name = file.replace("_model.pkl", "")
                    model_path = os.path.join(model_dir, file)
                    feat_path = os.path.join(model_dir, f"{model_name}_features.txt")

                    try:
                        clf = joblib.load(model_path)
                        with open(feat_path) as f:
                            selected = [line.strip() for line in f]
                    except Exception as e:
                        print(f"Failed to load model or features: {model_name} - {e}")
                        continue

                    used = [f for f in selected if f in X_ext.columns]
                    missing = [f for f in selected if f not in X_ext.columns]

                    if not used:
                        print(f"All features missing: {model_name} - {target_country}")
                        continue

                    X_eval = X_ext.copy()
                    for f in missing:
                        X_eval[f] = 0
                    X_eval = X_eval[selected]

                    y_prob = clf.predict_proba(X_eval)[:, 1]
                    fpr, tpr, thresholds = roc_curve(y_ext, y_prob)
                    best_idx = np.argmax(tpr - fpr)
                    best_thresh = thresholds[best_idx]
                    y_pred = (y_prob >= best_thresh).astype(int)

                    cm = confusion_matrix(y_ext, y_pred)
                    tn, fp, fn, tp = cm.ravel()
                    result = {
                        "Taxonomy_Level": taxonomy_level,
                        "Source": source_country,
                        "fml": fml_condition,
                        "Target": target_country,
                        "Model": model_name,
                        "AUC": roc_auc_score(y_ext, y_prob),
                        "Accuracy": accuracy_score(y_ext, y_pred),
                        "Sensitivity": recall_score(y_ext, y_pred),
                        "Specificity": tn / (tn + fp) if tn + fp > 0 else 0,
                        "PPV": precision_score(y_ext, y_pred),
                        "NPV": tn / (tn + fn) if tn + fn > 0 else 0,
                        "F1": f1_score(y_ext, y_pred),
                        "Youden": tpr[best_idx] - fpr[best_idx],
                        "Threshold": best_thresh,
                        "Used_Features": len(used),
                        "Missing_Features": len(missing)
                    }
                    all_results.append(result)
                    print(f"{taxonomy_level} | {model_name} - {target_country} | AUC = {result['AUC']:.3f}")

os.makedirs("./cross_cohort_validation/results_by_taxonomy", exist_ok=True)
df = pd.DataFrame(all_results)

for level in df['Taxonomy_Level'].unique():
    df_level = df[df['Taxonomy_Level'] == level]
    save_path = f"./cross_cohort_validation/results_by_taxonomy/external_validation_results_{level}.csv"
    df_level.to_csv(save_path, index=False)
    print(f"Saved results for taxonomy level '{level}' to {save_path}")

print("All cross-cohort validations completed and results saved by taxonomy level.")
