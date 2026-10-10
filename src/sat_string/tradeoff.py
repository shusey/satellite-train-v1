"""Small finite-horizon comparisons; no objective weights or optimality claim."""
from __future__ import annotations


def relative_to_reference(row, reference):
    """Positive accuracy percentages mean improvement; positive cost means more."""
    result = {}
    for key in ("e_max_m", "rms_m"):
        result[key + "_reduction_pct"] = 100 * (1 - row[key] / reference[key]) if reference[key] else None
    for key in ("transitions", "extra_altitude_loss_m", "high_equivalent_total_s"):
        result[key + "_increase_pct"] = 100 * (row[key] / reference[key] - 1) if reference[key] else None
    return result


def adjacent_cost(left, right):
    """Use signed differences; a non-improving error has no cost per improvement.

    Negative cost per positive improvement means the cost also fell. It is not
    silently converted into a positive penalty. End-point error is deliberately
    excluded because its ranking depends strongly on the stopping phase.
    """
    result = dict(from_factor=left["kp_factor"], to_factor=right["kp_factor"],
                  extra_transitions=right["transitions"] - left["transitions"],
                  extra_altitude_loss_m=right["extra_altitude_loss_m"] - left["extra_altitude_loss_m"])
    for key in ("e_max_m", "rms_m"):
        improvement = left[key] - right[key]
        result[key + "_improvement_m"] = improvement
        result[key + "_extra_transitions_per_m"] = result["extra_transitions"] / improvement if improvement > 0 else None
        result[key + "_extra_altitude_loss_per_m"] = result["extra_altitude_loss_m"] / improvement if improvement > 0 else None
    return result


def dominated_by(row, candidates, accuracy_key):
    """Finite-run Pareto domination in accuracy, transitions and altitude loss."""
    keys = (accuracy_key, "transitions", "extra_altitude_loss_m")
    return [other["kp_factor"] for other in candidates
            if all(other[k] <= row[k] for k in keys) and any(other[k] < row[k] for k in keys)]
