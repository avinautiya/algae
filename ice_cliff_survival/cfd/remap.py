"""Conservative remapping of the cell temperature field between two column-structured meshes with the SAME interface s(y) (e.g. coarse -> fine nx).
Mesh: nx x ny cells; column coordinate xi in [0,1], x = xi s(y); row j spans y in [y_j, y_{j+1}] with s linear in y inside the row, so the cell (i,j) has area ∫ dy ∫ s(y) dxi = dxi_i * sbar_j * dy exactly. Overlap area of old cell i and new cell i' in row j = (xi-overlap) * sbar_j * dy,
hence the conservative remap in xi-space is EXACT in area weights: T'_{i'} = sum_i T_i overlap(i,i') / dxi'. Energy sum(rho c T V) is preserved to round-off; accuracy is first order (piecewise-constant)."""
import numpy as np
def remap_rows(T, nx_new):
    ny, nx = T.shape; xo = np.linspace(0, 1, nx + 1); xn = np.linspace(0, 1, nx_new + 1); out = np.zeros((ny, nx_new))
    for i2 in range(nx_new):
        a, b = xn[i2], xn[i2 + 1]
        for i in range(nx):
            ov = min(b, xo[i + 1]) - max(a, xo[i])
            if ov > 0: out[:, i2] += T[:, i] * ov / (b - a)
    return out
def energy(T, s_node, ny_dy_dz, nx):
    """sum_{cells} T V with V = (1/nx) sbar_j dy dz (row-wise)"""
    sbar = 0.5 * (s_node[:-1] + s_node[1:]); return float(np.sum(T * (sbar[:, None] / nx)) * ny_dy_dz)
