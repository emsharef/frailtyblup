# Contributing

Use a virtual environment and install with `python -m pip install -e '.[docs]'`. Run `python -m unittest discover -s tests -v` and the worked examples before proposing changes. Numerical changes need independent score/derivative checks and a comparison against the frozen synthetic references; update those references only for an explained estimator change.

Keep the interval API independent of clinical file locations. Preserve failure records and distinguish numerical certification from statistical validity. Document changes to score normalization, projection geometry, initialization, solver tolerances, and candidate selection because these can change the estimator delivered by a finite search.

Use synthetic data in examples and tests. Changes to the paper's theorems or statistical claims require mathematical review beyond passing tests. Code style uses Black. Development status and known statistical limits are recorded in the README and algorithm guide.
