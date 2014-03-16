#!/usr/bin/env python

# TODO : we will need to import the necessary functions for mathematics,
# e.g. sin, cos, tan, etc. We need a SymPy pycode function similar to the
# ccode function that handles all of the function mapping like lambdify
# does.

# external libraries
import numpy as np


def mass_forcing_matrices(constants, coordinates, speeds, specified=None,
                          mass_matrix=None, forcing_vector=None):
    """Returns the values of the mass matrix and forcing vector given the
    numerical values of the constants, coordinates, and speeds.

    M(constants, coordinates(t), speeds(t)) x'(t) = F(constants, coordinates(t), speeds(t))

    M : mass matrix
    F : forcing vector
    x : states typically the [coordinates(t), speeds(t)]

    Parameters
    ----------
    constants : array_like, shape(4,)
        The numerical values of the constants.
        [m, c, k, g]
    coordinates : array_like, shape(1,)
        The numerical values of the coordinates.
        [x(t),]
    speeds : array_like, shape(1,)
        The numerical values of the speeds.
        [v(t),]
    specified : array_like, shape(1,), optional
        The numerical values of the specified inputs.
        [f(t),]
    mass_matrix : ndarray, shape(2,2), optional
        If passed in, this array will be populated instead of creating a new
        array on each evaluation. This can increase computation speed and
        potentially memory overhead.
    forcing_vector : ndarray, shape(2,)
        If passed in, this array will be populated instead of creating a new
        array on each evaluation. This can increase computation speed and
        potentially memory overhead.

    Returns
    -------
    mass_matrix : ndarray, shape(2,2)
        The computed mass matrix.
    forcing_vector : ndarray, shape(2,)
        The computed forcing vector.

    """

    # TODO : It would be nice if these dimension checks happened outside
    # this function, i.e. only once as opposed to every evaluations in the
    # integrator loop.

    if mass_matrix is None:
        mass_matrix = np.zeros(2, 2)
        assert mass_matrix.shape == (2, 2)

    if forcing_vector is None:
        forcing_vector = np.zeros(2)
        assert forcing_matrix.shape == (2,)

    assert len(constants) == 4
    assert len(coordinates) == 1
    assert len(speeds) == 1

    if specified is not None:
        assert len(specified) == 1

    # common subexpressions
    z_0 = speeds[0]

    # mass matrix
    mass_matrix[0, 0] = 1
    mass_matrix[0, 1] = 0
    mass_matrix[1, 0] = 0
    mass_matrix[1, 1] = constants[0]

    # forcing vector
    forcing_vector[0] = z_0
    forcing_vector[1] = -constants[2]*z_0 + constants[3]*constants[0] - constants[1]*coordinates[0] + specified[0]

    return mass_matrix, forcing_vector
