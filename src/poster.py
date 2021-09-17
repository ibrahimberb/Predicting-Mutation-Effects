# poster purpose

from matplotlib import rc
import matplotlib

matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['savefig.transparent'] = True

# Evaluation metrics
plt.figure(figsize=(25, 9))  # poster purpose
sns.set(style='white', font_scale=2.5)  # poster purpose

# Less evalution metrics
predator.eval_metrics.scoring_metrics_data_melted[predator.eval_metrics.scoring_metrics_data_melted['METRIC'].isin(
    ['accuracy', 'balanced_accuracy', 'f1', 'precision', 'recall',
       'roc_auc']
)]