## Dispersion Map

In the case of [`soxs_disp_solution`](../recipes/soxs_disp_solution.md) the output dispersion map hosts the coefficients ($c_{ij}$) to two polynomials that describe the global dispersion solution for the entire detector frame:

$$
X = \sum\limits_{ij} c_{ij} \times n^i \times \lambda^j \\
$$

$$
Y = \sum\limits_{ij} c_{ij} \times n^i \times \lambda^j \\
$$

[`soxs_spatial_solution`](../recipes/soxs_spatial_solution.md), building from this dispersion solution, provides global dispersion *and* spatial solution. The dispersion map output by this recipe hosts the coefficients ($c_{ijk}$) to two polynomials (note now the inclusion of slit position):

$$
X = \sum\limits_{ijk} c_{ijk} \times n^i \times \lambda^j \times s^k \\
$$

$$
Y = \sum\limits_{ijk} c_{ijk} \times n^i \times \lambda^j \times s^k \\
$$
