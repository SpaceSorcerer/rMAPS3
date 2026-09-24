"""Session fixtures shared by several test modules."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))


@pytest.fixture(scope="session")
def calibrated_arm(tmp_path_factory):
    """One synthetic arm with a real v2.1 calibration and both figure layers drawn (tests/synthetic_calibrated_arm.py)."""
    import synthetic_calibrated_arm as syn
    tmp = tmp_path_factory.mktemp("calibrated_arm")
    inputs = syn.build(tmp)
    inputs["figures"] = syn.draw(inputs, tmp / "figures")
    return inputs
