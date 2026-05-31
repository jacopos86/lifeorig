

#   clamp method

def _clamp(value: float, min_value: float, max_value: float) -> float:
    return min(max_value, max(min_value, value))