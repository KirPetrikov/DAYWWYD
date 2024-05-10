import numpy as np
import random

from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from concurrent.futures import ProcessPoolExecutor


SEED = 111
random.seed(SEED)
np.random.seed(SEED)


class RandomForestClassifierCustom(BaseEstimator):
    """
    Custom RandomForestClassifier with multiprocessing implementation
    Based on DecisionTreeClassifier

    To use multiprocessing
    specify parameter
    n_jobs : int, default 1
    for methods
    `fit`, `predict_proba` and `predict`

    Params
    ------
    n_estimators: int, default 10
    max_depth: int, default None
    max_features: int, default None
        If not specified (None), then calculates
        depending on the number of features:
        for <3 : 1
        for 3<=, <6 : 2
        >6 : 1/3 of number of features
    random_state: int, default None
    """
    
    def __init__(self,
                 n_estimators=10,
                 max_depth=None,
                 max_features=None,
                 random_state=None
                 ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []
        self.features_idxs_by_tree = []

    def fit(self, X, y, n_jobs=1):
        self.trees = []
        self.features_idxs_by_tree = []

        features_number = X.shape[1]

        # default max_features
        if self.max_features is None:
            if features_number < 3:
                max_features = 1
            elif features_number < 6:
                max_features = 2
            else:
                max_features = features_number // 3
        else:
            if self.max_features > features_number:
                raise ValueError(
                    f'Parameter max_features ({self.max_features}) '
                    f'is greater than dataset features number ({features_number})!')
            max_features = self.max_features

        pooled_data = [_ for _ in zip([(X,
                                        y,
                                        max_features,
                                        self.max_depth,
                                        self.random_state
                                        )] * self.n_estimators,
                                      range(self.n_estimators)
                                      )
                       ]

        with ProcessPoolExecutor(n_jobs) as pool:
            results = list(pool.map(self._single_dtree_estimator,
                                    pooled_data)
                           )

        for result in results:
            self.trees.append(result[0])
            self.features_idxs_by_tree.append(result[1])

        return self

    def _single_dtree_estimator(self, dtree_params):
        X, y, features_number, depth, random_state = dtree_params[0]
        seed_modifier = dtree_params[1]

        seed = self.random_state + seed_modifier
        random.seed(seed)
        np.random.seed(seed)

        current_features_idxs = np.random.choice(X.shape[1],
                                                 features_number,
                                                 replace=False
                                                 )
        bootstrap_idxs = np.random.choice(X.shape[0], X.shape[0] // 2, replace=True)

        X_bootstrap = np.take(X, bootstrap_idxs, axis=0)
        y_bootstrap = np.take(y, bootstrap_idxs, axis=0)
        X_bootstrap = np.take(X_bootstrap, current_features_idxs, axis=1)

        dt_classifier = DecisionTreeClassifier(max_depth=depth, random_state=seed)
        dt_classifier.fit(X_bootstrap, y_bootstrap)

        return dt_classifier, current_features_idxs

    def predict_proba(self, X, n_jobs=1):
        pooled_params = [(X, _[0], _[1]) for _ in zip(self.trees,
                                                      self.features_idxs_by_tree
                                                      )
                         ]

        with ProcessPoolExecutor(n_jobs) as pool:
            pred_prob = list(pool.map(self._single_pred_proba, pooled_params))

        mean_pred = np.sum(np.array(pred_prob), axis=0) / self.n_estimators

        return mean_pred

    def _single_pred_proba(self, params):
        X, model, idxs = params
        X_feature_sample = np.take(X, idxs, axis=1)
        current_pred_proba = model.predict_proba(X_feature_sample)

        return current_pred_proba

    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions

