from .machine_learning_utils import get_default_classifier
from .mylogger import get_handler
import logging
from tqdm.notebook import tqdm

handler = get_handler()

log = logging.getLogger(__name__)
log.handlers[:] = []
log.addHandler(handler)
log.setLevel(logging.DEBUG)


class DefaultModels(list):
    def __init__(self, n_experiment):
        super().__init__()
        for _ in range(n_experiment):
            self.append(get_default_classifier(random_state=42))


class TunedModels(list):
    def __init__(self, best_estimators):
        super().__init__()
        for estimator in best_estimators:
            self.append(estimator)


class FinalizedModels(list):
    def __init__(self, tuned_models):
        super().__init__()
        for estimator in tuned_models:
            self.append(estimator)

    def fit_all(self, data_materials, determined_feature_set):
        for exp, estimator in tqdm(enumerate(self), total=len(self)):
            estimator.fit(data_materials[f"Xs_{determined_feature_set}"][exp],
                          data_materials[f"ys"][exp])
