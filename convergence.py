import numpy as np
import matplotlib.pyplot as plt


def compute_residual(old, new, norm_type='L2'):
    diff = new - old
    if norm_type == 'L2':
        return float(np.sqrt(np.mean(diff ** 2)))
    elif norm_type == 'Linf':
        return float(np.max(np.abs(diff)))
    else:
        raise ValueError(f"Unknown norm_type: {norm_type!r}")


def new_residual_history():
    return {'U1': [], 'U2': [], 'U3': []}


def update_residual_history(history, U1_old, U1_new, U2_old, U2_new, U3_old, U3_new, norm_type='L2'):
    history['U1'].append(compute_residual(U1_old, U1_new, norm_type))
    history['U2'].append(compute_residual(U2_old, U2_new, norm_type))
    history['U3'].append(compute_residual(U3_old, U3_new, norm_type))
    return history


def has_converged(history, tol):
    if any(len(v) == 0 for v in history.values()):
        return False
    return all(v[-1] < tol for v in history.values())


def plot_convergence(history, tol=None):
    plt.figure(figsize=(10, 6))

    labels = {'U1': 'Continuity (U1)', 'U2': 'Momentum (U2)', 'U3': 'Energy (U3)'}
    for key, values in history.items():
        if len(values) == 0:
            continue
        iterations = np.arange(1, len(values) + 1)
        plt.semilogy(iterations, values, label=labels.get(key, key))

    if tol is not None:
        plt.axhline(tol, color='gray', linestyle='--', linewidth=1, label=f'Tolerance ({tol:.0e})')

    plt.xlabel("Iteration")
    plt.ylabel("Residual (RMS change per step)")
    plt.title("Convergence History")
    plt.legend()
    plt.grid(True, which='both')
    plt.show()