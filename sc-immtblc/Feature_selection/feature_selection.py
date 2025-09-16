import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, classification_report, roc_curve, auc
)
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
adata = sc.read_h5ad('your_data_path/your_data.h5ad')

# Prepare features and labels
X = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
y = adata.obs["condition2"].apply(lambda x: 1 if x == "Pre-vaccination non-responder" else 0)

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

# Standard scaling
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Handle class imbalance
num_neg = (y_train == 0).sum()
num_pos = (y_train == 1).sum()
scale_pos_weight = num_neg / num_pos if num_pos > 0 else 1.0

# Hyperparameter search
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [3, 5, 7],
    'learning_rate': [0.01, 0.1, 0.2],
    'subsample': [0.8, 1.0],
    'colsample_bytree': [0.8, 1.0],
    'scale_pos_weight': [scale_pos_weight]
}

xgb = XGBClassifier(
    random_state=42,
    use_label_encoder=False,
    eval_metric='logloss'
)

random_search = RandomizedSearchCV(
    estimator=xgb,
    param_distributions=param_grid,
    n_iter=20,
    scoring='roc_auc',
    cv=5,
    random_state=42,
    n_jobs=-1
)
random_search.fit(X_train_scaled, y_train)
best_params = random_search.best_params_
print("Best Parameters:", best_params)

# Train with best parameters
best_xgb = XGBClassifier(
    **best_params,
    random_state=42,
    use_label_encoder=False,
    eval_metric='logloss'
)

# Cross-validation
cv_scores = cross_val_score(best_xgb, X_train_scaled, y_train, cv=5, scoring='roc_auc')
print("Cross-Validation AUC Scores:", cv_scores)
print("Mean AUC:", np.mean(cv_scores))

# Train final model
best_xgb.fit(X_train_scaled, y_train)

# Predict on test set
y_pred = best_xgb.predict(X_test_scaled)
y_pred_proba = best_xgb.predict_proba(X_test_scaled)[:, 1]

# Metrics
print("Accuracy:", accuracy_score(y_test, y_pred))
print("Precision:", precision_score(y_test, y_pred))
print("Recall:", recall_score(y_test, y_pred))
print("F1 Score:", f1_score(y_test, y_pred))
print("AUC-ROC:", roc_auc_score(y_test, y_pred_proba))
print("\nClassification Report:\n", classification_report(y_test, y_pred))

# Feature importance
feature_importance_df = pd.DataFrame({
    'gene': X.columns,
    'importance': best_xgb.feature_importances_
}).sort_values(by='importance', ascending=False)

# Plot top 20 important features
plt.figure(figsize=(10, 8))
sns.barplot(x='importance', y='gene', data=feature_importance_df.head(20), palette="viridis")
plt.title("Top 20 Important Genes")
plt.xlabel("Importance")
plt.ylabel("Gene")
plt.tight_layout()
plt.show()

# ROC Curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f"ROC curve (AUC = {roc_auc:.2f})")
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc="lower right")
plt.grid()
plt.show()
